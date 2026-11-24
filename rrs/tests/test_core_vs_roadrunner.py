"""The Rust core (roadrunner_nitro.RoadRunner) against libroadrunner's recorded answers on the sosbench
models and a few Antimony models: ids, values and every selection kind, setValue, reset semantics, time
courses, structure (N, Nr, L, L0, Gamma) with conserved moiety analysis off and on, steady state,
Jacobians, eigenvalues, elasticities, control coefficients, integrator settings.

Each scenario below runs through an adapter; make_reference.py runs them once on libroadrunner and saves
its answers in reference/<model>.npz (errors too: a refusal the core must mirror, or a libroadrunner
failure, where the core must refuse too or be proven right: finite differences of its own rates, its own
exact unscaled values, numpy on its own Jacobian, a twin without moiety analysis). Eigenvalues: in
libroadrunner's order (LAPACK zgeev), bit for bit on libroadrunner's own Jacobian. Also steady-state
solver settings (reference/ss_settings.npz). This test replays them on the core only (no libroadrunner).

    PYTHONPATH=<built roadrunner_nitro>:. python -m pytest rrs/tests/test_core_vs_roadrunner.py
"""
import glob, os, re
import numpy as np
import pytest
import roadrunner_nitro as rn

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'reference')
SOS = sorted(glob.glob(os.path.join(HERE, '..', '..', 'benchmarks', 'sosbench', '*', '*.xml')))
ANT = {
    'moieties': '''
        compartment C = 2
        J1: S1 -> S2; k1*S1; J2: S2 -> S3; k2*S2; J3: S3 -> S1; k3*S3
        J4: A + S1 -> B; k4*A*S1; J5: B -> A + S1; k5*B
        S1 = 5; S2 = 1; S3 = 2; A = 3; B = 0.5; k1 = 0.3; k2 = 0.2; k3 = 0.1; k4 = 0.05; k5 = 0.4
        S1 in C; S2 in C; S3 in C; A in C; B in C''',
    'chain': '''
        J0: $X -> S1; k0*X; J1: S1 -> S2; k1*S1; J2: S2 -> S3; k2*S2 - k3*S3; J3: S3 -> ; k4*S3
        X = 2; S1 = 0; S2 = 1; S3 = 0.5; k0 = 1; k1 = 0.3; k2 = 0.2; k3 = 0.05; k4 = 0.4''',
    'rules': '''
        compartment V = 1.5
        J1: S1 -> S2; k1*S1*V; J2: S2 -> S1; k2*S2*V
        S1 in V; S2 in V
        S1 = 5; S2 = k1*10; k1 = 0.3; k2 = 0.2
        p := k1*2 + S1''',
    'events': '''
        J1: S1 -> S2; k1*S1; J2: S2 -> S1; k2*S2
        S1 = 5; S2 = 0; k1 = 0.3; k2 = 0.1
        E1: at time > 2: k1 = 1, S2 = S2 + 1''',
    'raterule': '''
        J1: S1 -> S2; k1*S1*q; J2: S2 -> S1; k2*S2
        S1 = 5; S2 = 1; k1 = 0.3; k2 = 0.2; q = 1
        q' = -0.1*q''',
}
MODELS = [(os.path.basename(os.path.dirname(p)), open(p).read()) for p in SOS] + \
         [(k, rn.antimony_to_sbml(v)) for k, v in ANT.items()]
IDS = [m[0] for m in MODELS]
REFUSE = ('stoichiometry', 'No rate')  # libroadrunner refusals of getValue the core mirrors


class Core:
    """The adapter on roadrunner_nitro (make_reference.py has the libroadrunner one)."""
    def __init__(self, sbml):
        self.r = rn.RoadRunner(sbml)
    def cma(self, on): self.r.conservedMoietyAnalysis = on
    def tol(self, rtol, atol):
        self.r.setIntegratorSetting('relative_tolerance', rtol)
        self.r.setIntegratorSetting('absolute_tolerance', atol)
    def __getattr__(self, n):  # ids, values, actions with the same names on both
        return getattr(self.r, n)
    def fconc(self): return self.r.getFloatingSpeciesConcentrations()
    def finit(self): return self.r.getFloatingSpeciesInitConcentrations()
    def cms(self): return self.r.getConservedMoietyIds()
    def E(self, scaled): return self.r.getElasticityMatrix(scaled)
    def CC(self, scaled): return self.r.getConcentrationControlCoefficientMatrix(scaled)
    def FCC(self, scaled): return self.r.getFluxControlCoefficientMatrix(scaled)
    def sim(self, a, b, n): return np.asarray(self.r.simulate(a, b, n))
    def settings(self, name):
        self.r.setIntegrator(name)
        return [(k, self.r.getIntegratorSetting(k)) for k in self.r.getIntegratorSettings()]
    def ss_set(self, solver, settings):
        self.r.setSteadyStateSolver(solver)
        for k, v in settings.items():
            self.r.setSteadyStateSolverSetting(k, v)


class Rec:
    """rec(key, fn, ...): make mode (ref None) stores fn()'s value (or its error); check mode compares
    the core's fn() with the stored value. Returns the value (the adapter's own) for later inputs."""
    def __init__(self, ref=None):
        self.ref, self.out = ref, {}
        self.unverified = {}  # libroadrunner failed, the core answered: a scenario must prove it right

    def __call__(self, key, fn, rtol=1e-6, atol=1e-9, scaled=False, mirror=REFUSE, kind='close', mask=None):
        if self.ref is None:
            try:
                v = fn()
            except Exception as e:
                self.out[key] = np.array('ERR:' + ('refuse' if any(m in str(e) for m in mirror) else 'skip'))
                return None
            if isinstance(v, tuple):  # matrix (array, rows, cols)
                self.out[key], self.out[key + '.r'], self.out[key + '.c'] = np.asarray(v[0], float), np.array(v[1], str), np.array(v[2], str)
            elif isinstance(v, (list, np.ndarray)) and len(v) and isinstance(np.asarray(v).flat[0], str):
                self.out[key] = np.array(v, str)
            else:
                self.out[key] = np.asarray(v, complex if np.iscomplexobj(v) else float) if not isinstance(v, list) or v else np.array([], str)
            return v
        assert key in self.ref, 'no reference for %s (rerun make_reference.py)' % key
        r = self.ref[key]
        if r.dtype.kind == 'U' and r.ndim == 0 and str(r).startswith('ERR:'):
            if str(r) == 'ERR:refuse':
                with pytest.raises(Exception):
                    fn()
                return None
            try:  # libroadrunner failed for its own reasons: the core may refuse too, or must be proven right
                v = fn()
            except Exception:
                return None
            self.unverified[key] = v
            return v
        v = fn()
        try:
            if isinstance(v, tuple):
                rows = list(self.ref[key + '.r'])
                assert rows == list(v[1])[:len(rows)] and list(self.ref[key + '.c']) == list(v[2]), 'labels'
                cmp(r, np.asarray(v[0]).reshape(r.shape), rtol, atol, scaled, mask and mask(self.ref))
            elif r.dtype.kind == 'U':
                assert list(r) == list(v)
            else:
                cmp(r, v, rtol, atol, scaled, mask and mask(self.ref))
        except AssertionError as e:
            raise AssertionError('%s: libroadrunner %r, core %r' % (key, r, v)) from e
        return v


def eig_match(a, b, rtol, atol, extra=0.0):
    """Eigenvalues of libroadrunner's finite-difference Jacobian (a) vs the core's exact one (b).
    The order is LAPACK's on both sides (eig_exact proves it on identical input); here only
    clusters closer than the Jacobians' difference may come out permuted. `extra`: the perturbation
    bound of a (near-)defective eigenvalue, ~sqrt(|dJ| |J|)."""
    tol = atol * max(1.0, np.abs(a).max()) + rtol * np.abs(a) + extra
    assert a.shape == b.shape
    if np.all(np.abs(a - b) <= tol):
        return
    left = list(range(len(b)))
    for i, x in enumerate(a):  # each value matched, a permutation only among near-equal values
        j = min(left, key=lambda j: abs(b[j] - x))
        assert abs(b[j] - x) <= tol[i] and abs(i - left.index(j)) < len(a), (i, x, b[j])
        left.remove(j)


def cmp(a, b, rtol, atol, scaled, mask=None):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, (a.shape, b.shape)
    if mask is not None:
        a, b = a[mask], b[mask]
    if scaled and a.size:  # atol relative to the largest value (values spanning many orders of magnitude)
        atol *= max(1.0, np.nanmax(np.abs(a)))
    assert np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True), np.nanmax(np.abs(a - b))


def gps(r):
    return [p for p in r.getGlobalParameterIds() if not p.startswith('_CSUM')]


# ---------------------------------------------------------------- scenarios (r: adapter, rec: Rec)

def sc_structure(r, rec, sbml):
    for n in ('getFloatingSpeciesIds', 'getBoundarySpeciesIds', 'getGlobalParameterIds', 'getCompartmentIds',
              'getReactionIds', 'getIndependentFloatingSpeciesIds', 'getDependentFloatingSpeciesIds', 'cms',
              'getEigenValueIds', 'getRatesOfChangeIds'):
        rec(n, lambda: list(getattr(r, n)()))
    rec('tcs', lambda: list(r.timeCourseSelections))
    rec('sss', lambda: list(r.steadyStateSelections))
    rec('gpv', lambda: r.getGlobalParameterValues())
    for n in ('getFullStoichiometryMatrix', 'getNrMatrix', 'getLinkMatrix', 'getL0Matrix', 'getConservationMatrix'):
        rec(n, lambda: getattr(r, n)(), mirror=('',))


def sc_values(r, rec, sbml, cma):
    rec('sim', lambda: r.sim(0, 3, 31), atol=1e-10, scaled=True)
    rec('t', r.getCurrentTime)
    r.reset()
    ps = gps(r)
    if ps:
        r.setValue(ps[0], 0.42)
    fs, bs, rx = r.getFloatingSpeciesIds(), r.getBoundarySpeciesIds(), r.getReactionIds()
    rec('fconc', r.fconc, rtol=1e-5)
    rec('rates', r.getReactionRates, rtol=1e-5, atol=1e-8)
    if not (cma and r.getRateRuleIds()):
        rec('roc', r.getRatesOfChange, rtol=1e-5, atol=1e-8)
    twin = None
    if cma and rec.ref is not None:  # the same state without moiety analysis (where libroadrunner is no reference)
        twin = rn.RoadRunner(sbml)
        if ps:
            twin.setValue(ps[0], 0.42)
    if twin and r.getRateRuleIds():  # libroadrunner returns garbage here: the rates are the twin's
        cmp([twin.getValue(i) for i in r.getRatesOfChangeIds()], r.getRatesOfChange(), 1e-12, 1e-14, False)
    rec('finit', r.finit)
    # (libroadrunner's init() of assignment-rule species and species uec with moieties corrupt its state)
    ruled = set(re.findall(r'<assignmentRule[^>]*variable="([^"]+)"', sbml))
    fs_ode = [s for s in fs if not (cma and s in ruled)]
    sels = ['time'] + fs + ['[%s]' % s for s in fs] + bs + ['[%s]' % s for s in bs] + list(r.getGlobalParameterIds()) + \
        list(r.getCompartmentIds()) + rx + ['init(%s)' % s for s in fs_ode] + ['init([%s])' % s for s in fs_ode] + \
        ['init(%s)' % p for p in ps] + ["%s'" % s for s in fs] + ["[%s]'" % s for s in fs] + \
        ['stoich(%s, %s)' % (s, j) for s in fs[:4] for j in rx[:4]] + \
        ['uec(%s, %s)' % (j, s) for j in rx[:3] for s in fs[:3] if not cma] + ['uec(%s, %s)' % (j, p) for j in rx[:2] for p in ps[:3]]
    for s in sels:
        rec('v:' + s, lambda: r.getValue(s), rtol=1e-5, atol=1e-8)
    if twin:  # left out above (libroadrunner corrupts its own state on them): the twin's values
        for q in ['init(%s)' % s for s in fs if s in ruled] + ['init([%s])' % s for s in fs if s in ruled] + \
                ['uec(%s, %s)' % (j, s) for j in rx[:3] for s in fs[:3]]:
            cmp(twin.getValue(q), r.getValue(q), 1e-12, 1e-14, False)
    for j in rx[:3]:
        for s in fs[:3] if not cma else []:
            if abs(r.getValue(j)) > 1e-6 and abs(r.getValue(s)) > 1e-6:
                rec('ec:%s:%s' % (j, s), lambda: r.getValue('ec(%s, %s)' % (j, s)), rtol=1e-4)
    if not fs:
        return
    # libroadrunner's Jacobian is a finite difference (the core's exact): values to 1e-4, same order
    eig_exact(r, rec, 'evx', r.getFullJacobian, r.getFullEigenValues)
    ev = eig_check(r, rec, 'ev', r.getFullJacobian, r.getFullEigenValues, 'evx.J')
    for i, s in enumerate(fs[:3]):  # eigen(S): the i-th of getFullEigenValues (species index), as roadrunner
        for part in ('eigenReal', 'eigenImag'):
            v = rec('ev:%s:%s' % (part, s), lambda: r.getValue('%s(%s)' % (part, s)), 1e9, 1e9, mirror=())
            if rec.ref is not None and v is not None:
                e = r.getFullEigenValues()[i]  # checked by eig_check (or against numpy where libroadrunner fails)
                assert v == (e.real if part == 'eigenReal' else e.imag)
                rec.verified('ev:%s:%s' % (part, s))


def eig_check(r, rec, key, jac, eig, jkey):
    """getFullEigenValues / getReducedEigenValues against libroadrunner's, allowing for the
    difference of the Jacobians (libroadrunner's is a finite difference)."""
    rec(key, eig, kind='skip-compare')
    if rec.ref is None or key not in rec.ref:
        return None
    a = rec.ref[key]
    if str(a) == 'ERR:refuse':
        with pytest.raises(Exception):
            eig()
        return None
    b = np.asarray(eig())
    if a.dtype.kind == 'U':  # libroadrunner fails (its Jacobian throws): the core's against its own Jacobian
        a = np.linalg.eigvals(np.asarray(jac()[0])) if len(b) else b
        eig_match(np.sort_complex(a), np.sort_complex(b), 1e-6, 1e-10)  # non-normal J: eigenvalue condition
        return b
    extra = 0.0
    if jkey in rec.ref:
        J = rec.ref[jkey]
        extra = 2 * np.sqrt(np.linalg.norm(J - np.asarray(jac()[0]), 2) * np.linalg.norm(J, 2))
    eig_match(a, b, 1e-4, 1e-6, extra)
    return b


def eig_exact(r, rec, key, jac, eig):
    """libroadrunner's eigenvalues of its own Jacobian: the core's LAPACK port (rn.eigenvalues) on the
    same matrix must give the same values in the same order (to rounding)."""
    if rec.ref is None:
        try:
            J1, e, J2 = np.asarray(jac()[0]), np.asarray(eig()), np.asarray(jac()[0])
        except Exception:
            return
        if np.array_equal(J1, J2):  # the Jacobian behind eig() is the one recorded
            rec.rec.out[rec.prefix + key + '.J'], rec.rec.out[rec.prefix + key + '.E'] = J1, e
        return
    if key + '.J' in rec.ref:
        J, E = rec.ref[key + '.J'], rec.ref[key + '.E']
        got = rn.eigenvalues(J.tolist())
        assert got.shape == E.shape and np.allclose(got, E, rtol=1e-12, atol=1e-12 * max(1.0, np.abs(E).max())), \
            '%s: libroadrunner %r, core %r' % (key, E, got)


def sc_set(r, rec, sbml, cma):
    r.sim(0, 1, 5)
    fs, bs = r.getFloatingSpeciesIds(), r.getBoundarySpeciesIds()
    ps = gps(r)
    near = lambda s, f: f * r.getValue(s) or 0.37  # scale-relative: some models live at 1e-12
    # (libroadrunner mis-sets dependent species when constant species take part in the moieties)
    tail = fs[-1] if len(r.getFullStoichiometryMatrix()[1]) == len(fs) or not cma else r.getIndependentFloatingSpeciesIds()[-1]
    sets = [(fs[0], 1.5), ('[%s]' % tail, 0.7)] + [(p, 1.3) for p in ps[:2]] + \
        ([(bs[0], 1.3), ('[%s]' % bs[0], 0.9)] if bs else []) + [(c, 1.2) for c in r.cms()[:1]]
    for k, (s, f) in enumerate(sets):
        x = near(s, f)
        if rec('set%d' % k, lambda: r.setValue(s, x) or 0.0, mirror=('',)) is None:
            continue  # refused (libroadrunner, and so the core)
        rec('set%d:get' % k, lambda: r.getValue(s))
        rec('set%d:fconc' % k, r.fconc)
        rec('set%d:gpv' % k, r.getGlobalParameterValues)
    rec('sim2', lambda: r.sim(1, 2, 5), atol=1e-10, scaled=True, mirror=('',))
    if rec.ref is not None and tail != fs[-1]:  # the species libroadrunner mis-sets: the core sets it
        x = 0.7 * r.getValue('[%s]' % fs[-1]) or 0.37
        r.setValue('[%s]' % fs[-1], x)
        assert abs(r.getValue('[%s]' % fs[-1]) - x) <= 1e-12 * max(1, abs(x)) and abs(r.fconc()[-1] - x) <= 1e-12 * max(1, abs(x))


def sc_reset(r, rec, sbml):
    fs, ps = r.getFloatingSpeciesIds(), r.getGlobalParameterIds()
    if not ps:
        return
    n = [0]
    def snap():
        n[0] += 1
        rec('s%d:fconc' % n[0], r.fconc)
        rec('s%d:gpv' % n[0], r.getGlobalParameterValues)
        rec('s%d:t' % n[0], r.getCurrentTime)
    r.sim(0, 4, 9); r.setValue(ps[0], 0.123); r.reset(); snap()  # species back, parameter stays
    r.sim(0, 4, 9); r.setValue(ps[0], 0.321); r.resetAll(); snap()  # parameters back too
    r.sim(0, 4, 9); r.setValue('init([%s])' % fs[0], 4.2); snap()  # reloads at the initial values
    r.sim(0, 4, 9); r.reset(); snap()
    r.setValue('init(%s)' % ps[0], 0.77); snap()
    r.resetToOrigin(); snap()
    # init(...) set after a parameter change: libroadrunner reloads at the initial values (parameters too)
    r.sim(0, 4, 9); r.setValue(ps[0], 0.555); r.setValue('init([%s])' % fs[0], 2.5); snap()
    rec('sim', lambda: r.sim(0, 2, 5), atol=1e-10, scaled=True)


def sc_volume(r, rec, sbml):  # 'rules' model: initial assignment re-evaluated, compartment change
    r.sim(0, 3, 4); r.setValue('k1', 0.5); r.reset()
    rec('a', r.fconc)
    rec('b', lambda: r.getValue('init([S2])'))
    r.sim(0, 3, 4); r.setValue('V', 4.0)
    rec('c', r.fconc)
    rec('d', lambda: r.sim(3, 6, 7), atol=1e-10, scaled=True)
    r.reset(); rec('e', r.fconc)
    r.resetAll(); rec('f', r.fconc)
    rec('g', lambda: r.getValue('V'))


def fd_check(r, sbml, what):
    """The core's Jacobian ('J', d[S]'/d[S]) or unscaled elasticities ('E', dv/d[S]) against central
    differences of its own rates, at the current state (where libroadrunner gives no reference)."""
    c = rn.RoadRunner(sbml)
    ids = c.getFloatingSpeciesIds()
    for s, x in zip(r.getFloatingSpeciesIds(), r.fconc()):
        c.setValue('[%s]' % s, x)
    for p in gps(r):
        try:
            c.setValue(p, r.getValue(p))
        except ValueError:  # set by an assignment rule: follows from the rest
            pass
    m, rows, cols = r.getFullJacobian() if what == 'J' else r.E(False)
    f = (lambda: np.array([c.getValue("[%s]'" % s) for s in ids])) if what == 'J' else (lambda: np.asarray(c.getReactionRates()))
    out = list(ids) if what == 'J' else list(c.getReactionIds())
    for j, s in enumerate(cols):
        x = c.getValue('[%s]' % s)
        h = 1e-6 * max(abs(x), 1e-3)
        c.setValue('[%s]' % s, x + h); a = f()
        c.setValue('[%s]' % s, x - h); b = f()
        c.setValue('[%s]' % s, x)
        fd = (a - b) / (2 * h)
        got = np.array([m[rows.index(k)][j] if k in rows else 0.0 for k in out])
        cmp(fd, got, 1e-5, 1e-8, True)


def sc_steady(r, rec, sbml):
    if rec('ss', lambda: r.steadyState() * 0, mirror=()) is None:
        return  # both find none
    if rec.verified('ss') is not None:  # libroadrunner finds none, the core one: it must be one
        S = np.asarray(r.fconc())
        assert np.abs(r.getRatesOfChange()).max() <= 1e-9 * max(1, np.abs(S).max()) and S.min() >= -1e-10
        return
    rec('fconc', r.fconc, rtol=1e-6, atol=1e-10)
    rec('ssv', r.getSteadyStateValues, rtol=1e-6, atol=1e-10)
    # libroadrunner's Jacobian, elasticities and MCA are finite differences (the core's exact): 1e-4
    J = rec('J', r.getFullJacobian, 1e-4, 1e-6, scaled=True, mirror=())
    if rec.verified('J') is not None:  # libroadrunner's Jacobian fails here
        fd_check(r, sbml, 'J')
    elif J is not None:
        rec('Jr', r.getReducedJacobian, 1e-4, 1e-6, scaled=True)
        eig_exact(r, rec, 'evx', r.getFullJacobian, r.getFullEigenValues)
        eig_exact(r, rec, 'evxr', r.getReducedJacobian, r.getReducedEigenValues)
        eig_check(r, rec, 'getFullEigenValues', r.getFullJacobian, r.getFullEigenValues, 'evx.J')
        eig_check(r, rec, 'getReducedEigenValues', r.getReducedJacobian, r.getReducedEigenValues, 'evxr.J')
    if len(r.getFloatingSpeciesIds()) > 20:
        return  # libroadrunner's finite-difference MCA takes minutes on the big models
    E = rec('E', lambda: r.E(False), 1e-4, 1e-6, scaled=True, mirror=())
    if rec.verified('E') is not None:  # libroadrunner fails on some rule models
        fd_check(r, sbml, 'E')
        return
    if E is None:
        return  # both refuse
    rec('v', lambda: np.abs(r.getReactionRates()))
    fine = lambda ref: (ref['v'][:, None] > max(1e-12, 1e-9 * ref['v'].max()))  # near-zero fluxes: ill-conditioned
    rec('Es', lambda: r.E(True), 1e-4, 1e-6, scaled=True, mask=lambda ref: fine(ref) & np.isfinite(ref['Es']))
    rec('uCC', lambda: r.CC(False), 1e-4, 1e-6, scaled=True)
    rec('S', lambda: np.abs(r.fconc()))
    rec('CC', lambda: r.CC(True), 1e-4, 1e-6, scaled=True,
        mask=lambda ref: (ref['S'][:, None] > 1e-12) & (ref['v'][None, :] > 1e-12) & np.isfinite(ref['CC']))
    rec('uFCC', lambda: r.FCC(False), 1e-4, 1e-6, scaled=True)
    rec('FCC', lambda: r.FCC(True), 1e-4, 1e-6, scaled=True)
    if rec.ref is not None:
        # entries masked above (libroadrunner's finite differences are meaningless at ~0 flux or
        # concentration): the core's scaled values must be its exact unscaled ones times x / v,
        # undefined only where the divisor is exactly 0
        v, S = np.asarray(r.getReactionRates()), np.asarray(r.fconc())
        uE, sE = np.asarray(r.E(False)[0]), np.asarray(r.E(True)[0])
        ok = (v[:, None] != 0) & (uE != 0)
        cmp((uE * S[None, :] / np.where(v[:, None] != 0, v[:, None], 1))[ok], sE[ok], 1e-12, 1e-300, False)
        uC, sC = np.asarray(r.CC(False)[0]), np.asarray(r.CC(True)[0])
        ok = (S[:, None] != 0) & np.isfinite(sC)
        cmp((uC * v[None, :] / np.where(S[:, None] != 0, S[:, None], 1))[ok], sC[ok], 1e-12, 1e-300, False)
    fs, rx, ps = r.getFloatingSpeciesIds(), r.getReactionIds(), gps(r)
    for x in fs[:2] + rx[:2]:
        for p in ps[:2]:
            rec('ucc:%s:%s' % (x, p), lambda: r.getuCC(x, p), 1e-4, 1e-8)
            xv = r.getValue(x if x in rx else '[%s]' % x)
            if abs(xv) > 1e-12:  # scaled by ~0: libroadrunner's finite difference is ill-conditioned
                rec('cc:%s:%s' % (x, p), lambda: r.getCC(x, p), 1e-4, 1e-8)
                rec('ccs:%s:%s' % (x, p), lambda: r.getValue('cc(%s, %s)' % (x, p)), 1e-4, 1e-8)
            elif rec.ref is not None and xv != 0:  # the core's own: ucc * p / x, exactly
                cmp(r.getuCC(x, p) * r.getValue(p) / xv, r.getCC(x, p), 1e-10, 1e-300, False)
    for j in rx[:2]:
        for x in fs[:2] + ps[:2]:
            rec('uee:%s:%s' % (j, x), lambda: r.getuEE(j, x), 1e-4, 1e-8)
            if abs(r.getValue(j)) > 1e-8:
                rec('ee:%s:%s' % (j, x), lambda: r.getEE(j, x), 1e-4, 1e-8)
            elif rec.ref is not None and r.getValue(j) != 0:  # ~0 flux: the core's own uee * x / v, exactly
                xv = r.getValue('[%s]' % x if x in fs else x)
                cmp(r.getuEE(j, x) * xv / r.getValue(j), r.getEE(j, x), 1e-10, 1e-300, False)


def run_all(make, sbml, name, rec):
    """Every scenario on a fresh adapter (make(sbml)); keys prefixed by scenario."""
    def sub(prefix, fn, cma=None, *args):
        r = make(sbml)
        r.tol(1e-10, 1e-12)
        if cma is not None:
            sr = Sub(rec, prefix)
            if sr('cma', lambda: r.cma(cma) or 0.0, mirror=('',)) is None:
                return  # libroadrunner refuses moiety analysis here; so must the core
            return fn(r, sr, sbml, *args)
        return fn(r, Sub(rec, prefix), sbml, *args)
    for cma in (False, True):
        sub('structure%d/' % cma, sc_structure, cma)
        sub('values%d/' % cma, sc_values, cma, cma)
        sub('set%d/' % cma, sc_set, cma, cma)
        if name != 'raterule':
            sub('steady%d/' % cma, sc_steady, cma)
    sub('reset/', sc_reset)
    if name == 'rules':
        sub('volume/', sc_volume)
    if name == 'chain':
        r = make(sbml)
        for n in ('cvode', 'gillespie', 'rk45', 'euler', 'rk4'):
            kv = r.settings(n)
            rec('settings/%s' % n, lambda: [k for k, _ in kv])
            for k, v in kv:
                if k != 'seed':
                    rec('settings/%s/%s' % (n, k), lambda: v)


# ---------------------------------------------------------------- steady-state solver settings

SS_ANT = {
    # Newton from S2 = 0 meets a singular Jacobian: nleq1 / newton need presimulation, nleq2 not
    'slow': 'J1: S1 -> S2; k1*S1 - k2*S2; J2: S2 -> ; k3*S2^2; J0: -> S1; k0; '
            'k0 = 1; k1 = 1; k2 = 0.5; k3 = 1e-3; S1 = 1e4; S2 = 0',
    # no steady state, dS/dt = 1e-8: only the approximation (|dS/dt| <= approx_tolerance) succeeds
    'const': 'J0: -> S1; k; k = 1e-8; S1 = 1',
    # no steady state at all
    'grow': 'J0: -> S1; k; k = 0.5; S1 = 1',
    # Newton from S1 = 1e3 on a flat Hill curve fails; presimulation reaches S1 = 1
    'hill': 'J0: -> S1; v0; J1: S1 -> ; Vm*S1^4/(K^4 + S1^4); v0 = 1; Vm = 2; K = 1; S1 = 1e3',
}
SS_SOLVERS = ('nleq1', 'nleq2', 'newton', 'newton_linesearch')
SS_SETTINGS = [{'allow_presimulation': p, 'allow_approx': a} for p in (False, True) for a in (False, True)] + [
    {'allow_presimulation': True, 'presimulation_times': [1000.0]},
    {'allow_presimulation': False, 'allow_approx': True, 'approx_tolerance': 1e-12},
    {'allow_presimulation': False, 'allow_approx': True, 'approx_time': 100.0, 'approx_maximum_steps': 100},
]
SS_EXTRA = {'nleq': [{'maximum_iterations': 1}, {'relative_tolerance': 1e-6}, {'minimum_damping': 0.5}],
            'newton': [{'num_max_iters': 1}, {'allow_negative': True}, {'func_norm_tol': 1e-4}]}


def ss_cases():
    for name, ant in SS_ANT.items():
        for solver in SS_SOLVERS:
            for k, st in enumerate(SS_SETTINGS + SS_EXTRA[solver[:-1] if solver.startswith('nleq') else 'newton']):
                yield '%s/%s/%d' % (name, solver, k), ant, solver, st


def ss_run(make, sbml, solver, settings):
    """(steadyState value, time, concentrations) or the error message."""
    r = make(sbml)
    try:
        r.ss_set(solver, settings)
        v = r.steadyState()
        return float(v), float(r.getCurrentTime()), np.asarray(r.fconc(), float)
    except Exception as e:
        return str(e)


def ss_check(key, ref, got, sbml, settings, solver):
    """The core's answer against libroadrunner's (ref) for one solver x settings case. Where they
    disagree on success, the side that succeeded must be right: a true steady state, or (approx) a
    point whose |d/dt| sum is within approx_tolerance."""
    def rates(conc):
        r = rn.RoadRunner(sbml)
        for s, c in zip(r.getFloatingSpeciesIds(), conc):
            r.setValue('[%s]' % s, c)
        return np.abs(r.getRatesOfChange())
    approx = settings.get('allow_approx', solver.startswith('newton')) and settings.get('approx_time', 10000.0)
    if isinstance(ref, str) and isinstance(got, str):
        return 'both fail'
    if isinstance(got, str):  # libroadrunner "succeeded": only acceptable at a non-steady state
        v, t, conc = ref
        assert rates(conc).max() > 1e-12 * max(1, conc.max()), '%s: core failed (%s), libroadrunner found %r' % (key, got, conc)
        return 'libroadrunner accepts a non-steady state (KINSOL tolerance)'
    v, t, conc = got
    if approx and t == approx:  # the approximation: libroadrunner's criterion, its end point
        assert v <= settings.get('approx_tolerance', 1e-6)
        if not isinstance(ref, str):
            if ref[1] != t and rates(ref[2]).max() > 1e-12 * max(1, ref[2].max()):
                return 'libroadrunner accepts a non-steady state (KINSOL tolerance); core approximates'
            assert ref[1] == t and np.allclose(conc, ref[2], rtol=1e-5, atol=1e-10), (key, ref, got)
            return 'match (approximation)'
        return 'approximation (libroadrunner: NLEQ errors skip its approximation decorator)'
    ftol = max(1e-9 * max(1, np.abs(conc).max()), settings.get('func_norm_tol', 0))  # a loose tolerance asked for
    assert rates(conc).max() <= ftol and conc.min() >= -1e-10, (key, got)
    if isinstance(ref, str):
        return 'core finds the steady state libroadrunner misses (%s)' % ref[:40]
    assert np.allclose(conc, ref[2], rtol=1e-5 if not settings.get('func_norm_tol') else 1e-3, atol=1e-9), (key, ref, got)
    assert t == ref[1], (key, 'time', ref[1], t)
    return 'match'


def test_steady_state_settings():
    with np.load(os.path.join(REF, 'ss_settings.npz'), allow_pickle=True) as z:
        ref = z['cases'].item()
    seen = {}
    for key, ant, solver, st in ss_cases():
        sbml = rn.antimony_to_sbml(ant)
        out = ss_check(key, ref[key], ss_run(Core, sbml, solver, st), sbml, st, solver)
        seen[key] = out
    # the behaviour the settings must produce (as libroadrunner)
    # settings index: 0 neither, 1 approx only, 2 presimulation only, 3 both (SS_SETTINGS)
    # no presimulation: libroadrunner's Newton meets a singular Jacobian; Nitro's engine solver finds it (verified in ss_check)
    assert seen['slow/nleq1/0'].startswith('core finds')
    assert seen['slow/nleq1/2'] == 'match' and seen['slow/newton/2'] == 'match'  # presimulation rescues
    assert seen['slow/nleq2/0'] == 'match'  # nleq2's rank-deficient steps need none
    assert seen['const/nleq2/3'] == 'match (approximation)'  # approximation rescues a model with no steady state
    assert seen['const/newton/1'].startswith('libroadrunner accepts a non-steady state')
    assert seen['const/nleq2/5'] == 'both fail'  # approx_tolerance 1e-12 < 1e-8
    assert all(seen['grow/%s/%d' % (s, k)] == 'both fail' for s in SS_SOLVERS for k in range(7))
    assert seen['hill/nleq2/2'].startswith('core finds')  # libroadrunner's presimulated NLEQ2 fails at S1 = 1


class Sub:
    """A Rec view with a key prefix."""
    def __init__(self, rec, prefix):
        self.rec, self.prefix = rec, prefix
    def verified(self, key):
        """The core's answer where libroadrunner failed (if any), now proven right by the caller."""
        return self.rec.unverified.pop(self.prefix + key, None)
    @property
    def ref(self):
        return None if self.rec.ref is None else _Prefixed(self.rec.ref, self.prefix)
    def __call__(self, key, *a, **k):
        if k.get('kind') == 'skip-compare':  # recorded for other checks, not compared itself
            if self.rec.ref is None:
                return self.rec(self.prefix + key, a[0])
            return None
        return self.rec(self.prefix + key, *a, **k) if 'mask' not in k or k['mask'] is None else \
            self.rec(self.prefix + key, *a, **dict(k, mask=lambda ref: k['mask'](_Prefixed(ref, self.prefix))))


class _Prefixed(dict):
    def __init__(self, d, p):
        super().__init__()
        self.d, self.p = d, p
    def __getitem__(self, k): return self.d[self.p + k]
    def __setitem__(self, k, v): self.d[self.p + k] = v
    def __contains__(self, k): return self.p + k in self.d


@pytest.mark.parametrize('name,sbml', MODELS, ids=IDS)
def test_against_libroadrunner(name, sbml):
    print(rn.__file__)
    with np.load(os.path.join(REF, name + '.npz')) as z:
        ref = dict(z)
    rec = Rec(ref)
    run_all(Core, sbml, name, rec)
    assert not rec.unverified, 'core answers where libroadrunner failed, unchecked: %s' % sorted(rec.unverified)


def test_integrator_registry():
    b = rn.RoadRunner(rn.antimony_to_sbml(ANT['chain']))
    for name in ('rk4', 'euler'):  # crate::fixed (values vs libroadrunner: test_fixed_vs_roadrunner.py)
        b.setIntegrator(name)
        assert np.asarray(b.simulate(0, 1, 5)).shape[0] == 5
    b.setIntegrator('rk45')
    y = np.asarray(b.simulate(0, 1, 5))  # variable_step_size=True (default): a row per accepted step
    assert y.shape[0] > 5 and y[-1, 0] == 1
    b.setIntegratorSetting('variable_step_size', False)
    assert np.asarray(b.simulate(0, 1, 5)).shape[0] == 5
    b.setIntegrator('cvode')
    b.reset()
    b.setIntegratorSetting('variable_step_size', True)
    y = np.asarray(b.simulate(0, 1, 5))
    assert y.shape[0] > 5 and y[0, 0] == 0 and y[-1, 0] == 1
    b.setIntegratorSetting('maximum_num_steps', 5000)
    with pytest.raises(ValueError):
        b.setIntegrator('lsoda')
    b.setSteadyStateSolver('newton')
    assert b.getSteadyStateSolverSetting('allow_presimulation') is True
    b.setSteadyStateSolverSetting('presimulation_time', 10)
    assert b.getSteadyStateSolverSetting('presimulation_time') == 10.0


def test_cvode_step_settings():
    sbml = rn.antimony_to_sbml("S1 -> S2; 0.4*S1; S1 = 10; S2 = 0")
    b = rn.RoadRunner(sbml)
    b.setIntegratorSetting('variable_step_size', True)
    b.setIntegratorSetting('initial_time_step', 1e-4)
    b.setIntegratorSetting('minimum_time_step', 1e-4)
    b.setIntegratorSetting('maximum_time_step', 0.04)
    y = np.asarray(b.simulate(0, 1, 5))
    assert y[0, 0] == 0 and y[-1, 0] == 1
    dt = np.diff(y[:, 0])
    assert dt[0] == pytest.approx(1e-4)
    assert np.all(dt[:-1] >= 1e-4 - 1e-15)
    assert np.max(dt) == pytest.approx(0.04)
    ref = rn.RoadRunner(sbml)
    ref.setIntegratorSetting('relative_tolerance', 1e-12)
    ref.setIntegratorSetting('absolute_tolerance', 1e-14)
    tight = np.asarray(ref.simulateTimes(y[:, 0].tolist()))
    np.testing.assert_allclose(y[:, 1:], tight[:, 1:], rtol=2e-5, atol=1e-6)

    b.reset()
    b.setIntegratorSetting('max_output_rows', 3)
    truncated = np.asarray(b.simulate(0, 1, 5))
    np.testing.assert_array_equal(truncated, y[:3])
    assert b.getCurrentTime() == truncated[-1, 0]

    b.reset()
    b.setIntegratorSetting('variable_step_size', False)
    b.setIntegratorSetting('maximum_num_steps', 1)
    with pytest.raises(ValueError, match='integration failed'):
        b.simulate(0, 1, 2)
    b.setIntegratorSetting('multiple_steps', True)
    assert np.asarray(b.simulate(0, 1, 2))[-1, 0] == 1


def test_gillespie_step_settings():
    b = rn.RoadRunner(rn.antimony_to_sbml("S1 -> ; 100; S1 = 1"))
    b.setIntegrator('gillespie')
    b.setSeed(17)
    b.setIntegratorSetting('nonnegative', True)
    b.setIntegratorSetting('minimum_time_step', 0.25)
    y = np.asarray(b.simulate(0, 1, 2))
    np.testing.assert_array_equal(y[:, 0], [0, 0.25, 0.5, 0.75, 1])
    np.testing.assert_array_equal(y[:, 1], [1, 0, 0, 0, 0])
    b.reset()
    b.setIntegratorSetting('maximum_num_steps', 1)
    with pytest.raises(ValueError, match='maximum_num_steps'):
        b.simulate(0, 1, 2)

    b = rn.RoadRunner(rn.antimony_to_sbml("S1 -> ; 0; S1 = 1"))
    b.setIntegrator('gillespie')
    b.setIntegratorSetting('maximum_time_step', 0.25)
    y = np.asarray(b.simulate(0, 1, 2))
    np.testing.assert_array_equal(y[:, 0], [0, 0.25, 0.5, 0.75, 1])
    b.reset()
    b.setIntegratorSetting('variable_step_size', False)
    np.testing.assert_array_equal(np.asarray(b.simulate(0, 1, 2))[:, 0], [0, 1])


def test_gillespie_and_state(tmp_path):
    b = rn.RoadRunner(rn.antimony_to_sbml('J1: S1 -> S2; k1*S1; S1 = 100; S2 = 0; k1 = 0.5'))
    b.setIntegrator('gillespie')
    b.setSeed(7)
    y1 = b.gillespie(0, 5, 11)
    b.reset()
    b.setSeed(7)
    assert np.array_equal(y1, b.gillespie(0, 5, 11))  # same seed: same run
    assert np.all(y1[:, 1] + y1[:, 2] == 100)
    b.reset()
    v = b.gillespie(0, 5)  # a row per firing
    assert v.shape[0] > 50 and v[-1, 0] == 5
    b.reset()
    mean, sd = b.gillespieMeanSD(2000, 0, 5, 6)
    cmp(mean[:, 1], 100 * np.exp(-0.5 * mean[:, 0]), 0.02, 1e-9, False)
    b.setIntegrator('cvode')
    b.reset()
    b.simulate(0, 2, 3)
    b.setValue('k1', 0.7)
    p = str(tmp_path / 'state.txt')
    b.saveState(p)
    c = rn.RoadRunner()
    c.loadState(p)
    assert c.getCurrentTime() == 2 and c.getValue('k1') == 0.7
    cmp(c.simulate(2, 4, 5), b.simulate(2, 4, 5), 1e-6, 1e-9, False)


def test_reset_rate_rule_compartment():
    # libroadrunner 2.11: reset() brings the rate-rule compartment and [S1] back together (1, 1, 1)
    r = rn.RoadRunner(rn.antimony_to_sbml("compartment c = 1; c' = 0.1; S1 in c; S1 = 1; J1: S1 -> ; k*S1; k = 0.2"))
    r.simulate(0, 5, 6)
    r.reset()
    assert r.getValue('c') == 1 and abs(r.getValue('[S1]') - 1) < 1e-12 and abs(r.getValue('S1') - 1) < 1e-12


def test_big_model_load_defers_structural_analysis():
    import time
    src = open(os.path.join(HERE, '..', '..', 'benchmarks', 'bigmodel', 'synthetic_1000_species.ant')).read()
    sbml = rn.antimony_to_sbml(src)
    t = time.perf_counter()
    r = rn.RoadRunner(sbml)
    assert time.perf_counter() - t < 0.74  # libstruct's QR runs on the first structural query
    assert len(r.getIndependentFloatingSpeciesIds()) + len(r.getDependentFloatingSpeciesIds()) == len(r.getFloatingSpeciesIds())


@pytest.mark.parametrize('volume', ["c = 1; c' = 0.1", "c := 1 + 2*p; p = 0; p' = 0.05"])
def test_amount_rate_variable_compartment(volume):
    # S1' is d(amount)/dt = reaction rate k*[S1] (libroadrunner agrees); V d[S1]/dt alone misses C dV/dt
    r = rn.RoadRunner(rn.antimony_to_sbml(f"compartment c; {volume}; S1 in c; S1 = 1; J1: S1 -> ; k*S1; k = 0.2"))
    r.simulate(0, 2, 3)
    assert abs(r.getValue("S1'") + 0.2 * r.getValue('[S1]')) < 1e-12


def test_amount_rate_rule_in_assignment_compartment():
    r = rn.RoadRunner(rn.antimony_to_sbml(
        "compartment c; c := 1 + 0.1*time; substanceOnly species S1 in c; S1 = 1; S1' = 2"))
    r.setIntegratorSetting('relative_tolerance', 1e-11)
    r.setIntegratorSetting('absolute_tolerance', 1e-13)
    r.simulate(0, 2, 3)
    assert r.getValue("S1'") == pytest.approx(2.0, abs=1e-12)
    assert r.getValue('S1') == pytest.approx(5.0, abs=1e-8)


def test_steady_state_found_where_newton_alone_fails():
    # sosbench 00009 (MAPK cascade) after simulate(0, 10): Newton fails from there with default settings
    # (libroadrunner: "Maximum iterations exceeded"), but the model has a steady state. Nitro finds it.
    f = next(p for p in SOS if os.sep + '00009' + os.sep in p)
    r = rn.RoadRunner(open(f).read())
    r.simulate(0, 10, 6)
    r.steadyState()
    x = np.array([r.getValue(i.rstrip("'")) for i in r.getRatesOfChangeIds()])
    assert np.max(np.abs(r.getRatesOfChange())) <= 1e-9 * np.max(np.abs(x))  # as ss_check: a true steady state


def test_reset_reevaluates_initial_assignment_with_function():
    # libroadrunner 2.11: reset() re-evaluates S1 = f(k1) with the new k1 (10)
    r = rn.RoadRunner(rn.antimony_to_sbml('function f(x)\n 2*x\nend\nJ1: S1 -> S2; k2*S1; S1 = f(k1); S2 = 0; k1 = 3; k2 = 0.1'))
    r.simulate(0, 2, 3)
    r.setValue('k1', 5)
    r.reset()
    assert r.getValue('S1') == 10


def test_save_load_state_keeps_delay_history(tmp_path):
    # a continued delay model reads its history: a loaded state must carry it
    b = rn.RoadRunner(rn.antimony_to_sbml('J1: -> S1; k*delay(S1, 1.5); J2: S1 -> ; 0.5*S1; S1 = 1; k = 0.4'))
    b.simulate(0, 3, 7)
    p = str(tmp_path / 'delay.txt')
    b.saveState(p)
    c = rn.RoadRunner()
    c.loadState(p)
    cmp(c.simulate(3, 5, 9), b.simulate(3, 5, 9), 1e-12, 1e-14, False)


def test_get_sbml_explicit_stoichiometry():
    # libroadrunner 2.11's getSBML writes stoichiometry="1" on references that leave it implicit
    sbml = ('<?xml version="1.0" encoding="UTF-8"?><sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">'
            '<model id="m"><listOfCompartments><compartment id="c" size="1"/></listOfCompartments><listOfSpecies>'
            '<species id="A" compartment="c" initialConcentration="1"/><species id="B" compartment="c" initialConcentration="0"/>'
            '</listOfSpecies><listOfParameters><parameter id="k" value="0.1"/></listOfParameters><listOfReactions>'
            '<reaction id="J"><listOfReactants><speciesReference species="A"/></listOfReactants><listOfProducts>'
            '<speciesReference species="B"/></listOfProducts><kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">'
            '<apply><times/><ci>k</ci><ci>A</ci></apply></math></kineticLaw></reaction></listOfReactions></model></sbml>')
    out = rn.RoadRunner(sbml).getSBML()
    assert re.findall(r'<speciesReference [^>]*>', out) == ['<speciesReference species="A" stoichiometry="1"/>',
                                                           '<speciesReference species="B" stoichiometry="1"/>']
    cmp(rn.RoadRunner(out).simulate(0, 5, 6), rn.RoadRunner(sbml).simulate(0, 5, 6), 1e-14, 1e-16, False)

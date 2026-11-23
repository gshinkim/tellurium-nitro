"""The Rust core (roadrunner_nitro.RoadRunner) against libroadrunner's recorded answers on the sosbench
models and a few Antimony models: ids, values and every selection kind, setValue, reset semantics, time
courses, structure (N, Nr, L, L0, Gamma) with conserved moiety analysis off and on, steady state,
Jacobians, eigenvalues, elasticities, control coefficients, integrator settings.

Each scenario below runs through an adapter; make_reference.py runs them once on libroadrunner and saves
its answers in reference/<model>.npz (errors too: a refusal the core must mirror, or a libroadrunner bug
to skip). This test replays them on the core only (no libroadrunner needed).

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


class Rec:
    """rec(key, fn, ...): make mode (ref None) stores fn()'s value (or its error); check mode compares
    the core's fn() with the stored value. Returns the value (the adapter's own) for later inputs."""
    def __init__(self, ref=None):
        self.ref, self.out = ref, {}

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
        v = fn()
        try:
            if isinstance(v, tuple):
                rows = list(self.ref[key + '.r'])
                assert rows == list(v[1])[:len(rows)] and list(self.ref[key + '.c']) == list(v[2]), 'labels'
                cmp(r, np.asarray(v[0]).reshape(r.shape), rtol, atol, scaled, mask and mask(self.ref))
            elif r.dtype.kind == 'U':
                assert list(r) == list(v)
            elif kind == 'member':
                assert np.min(np.abs(r.real - v)) < rtol * max(1, np.abs(r).max())
            else:
                cmp(r, v, rtol, atol, scaled, mask and mask(self.ref))
        except AssertionError as e:
            raise AssertionError('%s: libroadrunner %r, core %r' % (key, r, v)) from e
        return v


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
    if not (cma and r.getRateRuleIds()):  # libroadrunner returns garbage there
        rec('roc', r.getRatesOfChange, rtol=1e-5, atol=1e-8)
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
    for j in rx[:3]:
        for s in fs[:3] if not cma else []:
            if abs(r.getValue(j)) > 1e-6 and abs(r.getValue(s)) > 1e-6:
                rec('ec:%s:%s' % (j, s), lambda: r.getValue('ec(%s, %s)' % (j, s)), rtol=1e-4)
    ev = rec('ev', r.getFullEigenValues, kind='skip-compare') if fs else None  # eigen(): own order; members
    if rec.ref is not None and fs and 'ev' in rec.ref and rec.ref['ev'].dtype.kind != 'U':
        for s in fs[:3]:
            rec.ref['ev:' + s] = rec.ref['ev']
            rec('ev:' + s, lambda: r.getValue('eigenReal(%s)' % s), rtol=1e-6, kind='member')
    return ev


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


def sc_steady(r, rec, sbml):
    if rec('ss', lambda: r.steadyState() * 0, mirror=()) is None:  # libroadrunner finds none: skip
        return
    rec('fconc', r.fconc, rtol=1e-6, atol=1e-10)
    rec('ssv', r.getSteadyStateValues, rtol=1e-6, atol=1e-10)
    # libroadrunner's Jacobian, elasticities and MCA are finite differences (the core's exact): 1e-4
    if rec('J', r.getFullJacobian, 1e-4, 1e-6, scaled=True, mirror=()) is not None:
        rec('Jr', r.getReducedJacobian, 1e-4, 1e-6, scaled=True)
        for f in ('getFullEigenValues', 'getReducedEigenValues'):
            rec(f + '.re', lambda: np.sort_complex(getattr(r, f)()).real, 1e-4, 1e-6, scaled=True)
            rec(f + '.im', lambda: np.sort(np.abs(getattr(r, f)().imag)), 1e-4, 1e-6, scaled=True)
    if len(r.getFloatingSpeciesIds()) > 20:
        return  # libroadrunner's finite-difference MCA takes minutes on the big models
    if rec('E', lambda: r.E(False), 1e-4, 1e-6, scaled=True, mirror=()) is None:
        return  # libroadrunner fails on some rule models
    rec('v', lambda: np.abs(r.getReactionRates()))
    fine = lambda ref: (ref['v'][:, None] > max(1e-12, 1e-9 * ref['v'].max()))  # near-zero fluxes: ill-conditioned
    rec('Es', lambda: r.E(True), 1e-4, 1e-6, scaled=True, mask=lambda ref: fine(ref) & np.isfinite(ref['Es']))
    rec('uCC', lambda: r.CC(False), 1e-4, 1e-6, scaled=True)
    rec('S', lambda: np.abs(r.fconc()))
    rec('CC', lambda: r.CC(True), 1e-4, 1e-6, scaled=True,
        mask=lambda ref: (ref['S'][:, None] > 1e-12) & (ref['v'][None, :] > 1e-12) & np.isfinite(ref['CC']))
    rec('uFCC', lambda: r.FCC(False), 1e-4, 1e-6, scaled=True)
    rec('FCC', lambda: r.FCC(True), 1e-4, 1e-6, scaled=True)
    fs, rx, ps = r.getFloatingSpeciesIds(), r.getReactionIds(), gps(r)
    for x in fs[:2] + rx[:2]:
        for p in ps[:2]:
            rec('ucc:%s:%s' % (x, p), lambda: r.getuCC(x, p), 1e-4, 1e-8)
            if abs(r.getValue(x if x in rx else '[%s]' % x)) > 1e-12:  # scaled by ~0: ill-conditioned
                rec('cc:%s:%s' % (x, p), lambda: r.getCC(x, p), 1e-4, 1e-8)
                rec('ccs:%s:%s' % (x, p), lambda: r.getValue('cc(%s, %s)' % (x, p)), 1e-4, 1e-8)
    for j in rx[:2]:
        for x in fs[:2] + ps[:2]:
            rec('uee:%s:%s' % (j, x), lambda: r.getuEE(j, x), 1e-4, 1e-8)
            if abs(r.getValue(j)) > 1e-8:
                rec('ee:%s:%s' % (j, x), lambda: r.getEE(j, x), 1e-4, 1e-8)


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
        if name not in ('events', 'raterule'):
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


class Sub:
    """A Rec view with a key prefix."""
    def __init__(self, rec, prefix):
        self.rec, self.prefix = rec, prefix
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
    run_all(Core, sbml, name, Rec(ref))


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
    with pytest.raises(ValueError):
        b.setIntegratorSetting('variable_step_size', True)  # not honoured: refused
    b.setIntegratorSetting('maximum_num_steps', 5000)
    with pytest.raises(ValueError):
        b.setIntegrator('lsoda')
    b.setSteadyStateSolver('newton')
    assert b.getSteadyStateSolverSetting('allow_presimulation') is True
    b.setSteadyStateSolverSetting('presimulation_time', 10)
    assert b.getSteadyStateSolverSetting('presimulation_time') == 10.0


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

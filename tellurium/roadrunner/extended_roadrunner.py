"""
Extending the RoadRunner object with additional fields.
"""

from __future__ import print_function, division, absolute_import
import os
import roadrunner
import warnings
import copy
import numpy as np

# ---------------------------------------------------------------------
# Turbo: route roadrunner's expensive calls to the Rust engine roadrunner_rs
# ---------------------------------------------------------------------
try:
    import roadrunner_rs  # noqa: F401
    _turbo = os.environ.get('TELLURIUM_ENGINE', '').strip().lower() != 'roadrunner'
except ImportError:
    _turbo = False

_NamedArray = roadrunner._roadrunner.NamedArray


class _Opts:  # the SimulateOptions fields _rs_simulate reads, for the simulate fast path
    times, output_file, reset_model = (), '', False

    def __init__(self, start, end, steps):
        self.start, self.end, self.steps = start, end, steps
# integrator settings the Rust engine honours only at these values (others fall back to roadrunner).
# cvode: tolerances are passed through; stiff/BDF/Adams orders/maximum_num_steps don't change the
# solution within tolerance (Rodas4 is stiff-stable and has its own step cap).
_RS_INTEGRATORS = {
    'cvode': {'maximum_time_step': 0.0, 'minimum_time_step': 0.0, 'initial_time_step': 0.0,
              'multiple_steps': False, 'variable_step_size': False},
    'gillespie': {'minimum_time_step': 0.0, 'maximum_time_step': 0.0, 'nonnegative': False, 'maximum_num_steps': 0},
}
_RS_INTEGRATORS = {k: (tuple(v), tuple(v.values())) for k, v in _RS_INTEGRATORS.items()}  # (keys, values)


def turbo(enable=None):
    """Use the Rust engine (roadrunner_rs) for simulate/steadyState/MCA when it supports the model.

    :param enable: True/False to switch globally; None only queries.
    :returns: whether turbo is on. On by default when roadrunner_rs is installed;
        env ``TELLURIUM_ENGINE=roadrunner`` starts it off. Unsupported models,
        integrators or settings always fall back to roadrunner.
    """
    global _turbo
    if enable is not None:
        if enable:
            try:
                import roadrunner_rs  # noqa: F401
            except ImportError:
                warnings.warn("te.turbo: roadrunner_rs is not installed; using roadrunner.")
                enable = False
        _turbo = bool(enable)
    return _turbo


def _named(a, rows=None, cols=None):
    a = np.ascontiguousarray(a, dtype=float).view(_NamedArray)
    if rows is not None:
        a.rownames = list(rows)
    if cols is not None:
        a.colnames = list(cols)
    return a


# ---------------------------------------------------------------------
# Extended RoadRunner class
# ---------------------------------------------------------------------
class ExtendedRoadRunner(roadrunner.RoadRunner):

    def __init__(self, *args, **kwargs):
        super(ExtendedRoadRunner, self).__init__(*args, **kwargs)

    # ---------------------------------------------------------------------
    # Turbo (Rust engine) overrides; each falls back to roadrunner via super()
    # ---------------------------------------------------------------------
    def _rs(self):
        """Compiled Rust model + index maps for the current model, or None (off/unsupported)."""
        mdl = self.model
        if not _turbo or mdl is None:
            return None
        c = self.__dict__.get('_rs_cache')
        # rebuild on a CMA switch or when a compartment volume baked into the Rust model changes
        cma = self.conservedMoietyAnalysis
        if c is None or c[0] != (cma, mdl.getCompartmentVolumes()[c[2]].tobytes()):
            st = self._rs_build()
            cc = st['cconst'] if st is not None else list(range(self.model.getNumCompartments()))
            c = self.__dict__['_rs_cache'] = ((cma, mdl.getCompartmentVolumes()[cc].tobytes()), st, cc)
        return c[1]

    def _rs_build(self):
        from ..analysis.rrs_bridge import from_roadrunner
        info = {}
        try:
            m = from_roadrunner(self, info)
        except (NotImplementedError, ValueError):
            return None
        mdl = self.model
        ids, rids = list(mdl.getFloatingSpeciesIds()), list(mdl.getReactionIds())
        gids, cids, bids = list(mdl.getGlobalParameterIds()), list(mdl.getCompartmentIds()), list(mdl.getBoundarySpeciesIds())
        states, outs = list(m.state_ids), set(m.output_ids)
        if not set(ids) <= outs or sorted(rids) != sorted(m.reactions):
            return None

        def known(k):  # Rust parameter or state (not a rule output / moiety total)
            try:
                m.steady_state_batch([k], np.zeros((0, 1)))
                return True
            except ValueError:
                return False
        pidx, cidx, bidx, sidx = ([i for i, k in enumerate(l) if known(k)] for l in (gids, cids, bids, ids))
        keys = [gids[i] for i in pidx] + [cids[i] for i in cidx] + [bids[i] for i in bidx] + [ids[i] for i in sidx]
        comp_of, general = info['comp_of'], info['general']
        rx = [list(m.reactions).index(j) for j in rids]
        NR, L, _ = m.conservation()
        Ngen = np.asarray(L) @ np.asarray(NR)  # amount stoichiometry, states x (reactions + rate rules)
        vols = mdl.getCompartmentVolumes()
        sbml = self.getCurrentSBML()
        st = dict(m=m, kres=m.resolve_keys(keys), ids=ids, rids=rids, gids=gids, cids=cids, bids=bids, pidx=pidx, cidx=cidx, bidx=bidx,
                  sidx=sidx, keys=keys, kpos={k: i for i, k in enumerate(keys)}, off=len(keys) - len(sidx),
                  bamt=np.array([bids[i] in info['amount_params'] for i in bidx], dtype=bool),
                  amount_params=info['amount_params'], comp_of=comp_of, general=general, outs=outs, states=states,
                  rx=rx, ev=(), cols={}, sp=[states.index(s) if s in states else None for s in ids], asc=np.ones(len(states)),
                  V=np.array([vols[cids.index(comp_of[s])] for s in ids]),
                  # species volumes the Rust model bakes in (rebuild when they change)
                  cconst=sorted({cids.index(comp_of[s]) for s in m.species if s not in general}),
                  # steady state: roadrunner refuses events, and time-dependent rates have none
                  nss='symbols/time' in sbml or mdl.getNumEvents() > 0)
        st['spi'] = [i for i, j in enumerate(st['sp']) if j is not None]
        st['spj'] = [st['sp'][i] for i in st['spi']]
        # what _rs_sync writes back: species, and parameters/compartments/boundary species that can change
        chg = info['changing']
        ba = [i for i, a in zip(bidx, st['bamt']) if bids[i] in chg and a]
        bc = [i for i, a in zip(bidx, st['bamt']) if bids[i] in chg and not a]
        oid = {n: i for i, n in enumerate(m.output_ids)}
        nb = len(bids) if bidx else 0  # _rs_x gathers [globals | volumes | boundary conc | boundary amounts | floating]
        o1, o2 = len(gids), len(gids) + len(cids)
        st['xidx'] = np.array(pidx + [o1 + i for i in cidx] + [o2 + i + nb * a for i, a in zip(bidx, st['bamt'])] +
                              [o2 + 2 * nb + i for i in sidx], dtype=int)
        st['sync'] = [(f, idx, [names[i] for i in idx], np.array([oid[names[i]] for i in idx])) for f, idx, names in (
            ('setCompartmentVolumes', [i for i in cidx if cids[i] in chg], cids),  # volumes first: species keep conc
            ('setGlobalParameterValues', [i for i in pidx if gids[i] in chg], gids),
            ('amount', ba, bids), ('conc', bc, bids),  # boundary species: via setValue (no array setter)
            ('setFloatingSpeciesConcentrations', sidx, ids)) if idx]
        # reaction-network views (MCA, Jacobian, elasticities) over the floating species; rate-rule
        # states are not floating species in roadrunner and are held fixed, as roadrunner does
        st['plain'] = all(s in m.species and s not in general for s in ids)
        st['mca'] = st['plain'] and not st['nss']
        if st['plain']:
            st['N'] = Ngen[st['sp']][:, rx]
        if st['mca']:
            Lm = self.getLinkMatrix()
            if set(ids) <= set(Lm.rownames):
                st.update(L=np.array(Lm)[[list(Lm.rownames).index(s) for s in ids]], ind=list(Lm.colnames))
            else:
                st['mca'] = False
        return st

    def _rs_x(self, st, mdl=None):
        """Current values of st['keys']: parameters, compartments, boundary species (SBML units), species conc."""
        mdl = self.model if mdl is None else mdl
        b = [mdl.getBoundarySpeciesConcentrations(), mdl.getBoundarySpeciesAmounts()] if st['bidx'] else []
        return np.concatenate([mdl.getGlobalParameterValues(), mdl.getCompartmentVolumes(), *b,
                               mdl.getFloatingSpeciesConcentrations()])[st['xidx']]

    def _rs_sync(self, st, d, mdl=None):
        """Write Rust values (dict id -> value; species as concentrations) back into the roadrunner model."""
        mdl = self.model if mdl is None else mdl
        arr = isinstance(d, np.ndarray)  # all outputs, in output_ids order (else a dict id -> value)
        for f, idx, names, oidx in st['sync']:
            if arr and f not in ('amount', 'conc'):
                getattr(mdl, f)(idx, d[oidx])
                continue
            if arr:
                have = list(zip(idx, names, d[oidx]))
            else:
                have = [(i, n, d[n]) for i, n in zip(idx, names) if n in d]
            if not have:
                continue
            if f in ('amount', 'conc'):
                for _, n, x in have:
                    roadrunner.RoadRunner.setValue(self, n if f == 'amount' else '[%s]' % n, float(x))
            else:
                i, _, x = zip(*have)
                getattr(mdl, f)(list(i), np.array(x, dtype=float))

    def _rs_mcast(self):
        st = self._rs()
        return st if st is not None and st['mca'] else None

    def simulate(self, start=None, end=None, points=None, selections=None, output_file=None, steps=None, times=None):
        # Turbo fast path for simulate(start, end, points or steps[, selections]) with CVODE: what
        # roadrunner's wrapper does (its simulate options are private and reset by every call) without
        # its per-call overhead. Anything else, or a model/setting Rust can't run, takes the wrapper
        # (which also turns gillespie's variable_step_size off for such calls).
        if _turbo and output_file is None and times is None and start is not None and end is not None \
                and (points is None) != (steps is None) and type(points) is not list:
            n = points - 1 if steps is None else steps
            integ = self.integrator
            st = self._rs() if n >= 1 and integ.getName() == 'cvode' else None
            sel = None if st is None else list(self.timeCourseSelections)
            # new selections go through roadrunner's wrapper once (it validates and sets them)
            if st is not None and (selections is None or (isinstance(selections, (list, tuple)) and sel == list(selections))):
                o = _Opts(start=start, end=end, steps=n)
                res = self.__dict__['_rs_data'] = self._rs_simulate(st, o, sel, integ)
                if res is not None:
                    return res
        return super(ExtendedRoadRunner, self).simulate(start, end, points, selections, output_file, steps, times)
    simulate.__doc__ = roadrunner.RoadRunner.simulate.__doc__

    def _simulate(self, o):
        st = self._rs()
        res = None if st is None else self._rs_simulate(st, o)
        self.__dict__['_rs_data'] = res
        if res is None and st is not None:
            st['ev'] = ()  # roadrunner moves the event triggers itself; Rust restarts them next time
        return super(ExtendedRoadRunner, self)._simulate(o) if res is None else res

    def _rs_cols(self, st, sel):
        """(Rust outputs, plan turning their columns into the selections), or None if Rust can't."""
        key = tuple(sel)
        if key in st['cols']:
            return st['cols'][key]
        outs = []

        def o(n):
            if n not in outs:
                outs.append(n)
            return outs.index(n)
        sp = set(st['ids']) | set(st['bids'])
        I, mul, div, lin, init = [], [], [], [], []
        for k, s in enumerate(sel):
            br = s[:1] == '[' and s[-1:] == ']'
            v = s[1:-1] if br else s
            if v in sp:
                if v not in st['outs']:
                    return None
                I.append(o(v))
                amt = v in st['amount_params']  # boundary species whose Rust value is an amount
                if br and amt:
                    div.append((k, o(st['comp_of'][v])))
                elif not br and not amt:
                    mul.append((k, o(st['comp_of'][v])))
            elif s in st['outs']:  # time, parameters, compartments, reaction rates, rule variables
                I.append(o(s))
            elif s.startswith('init(') and s.endswith(')'):
                I.append(0)
                init.append((k, s))
            elif s.endswith("'") and st['plain'] and s[:-1].strip('[]') in st['ids']:
                i = st['ids'].index(s[:-1].strip('[]'))
                I.append(0)  # amount (or conc: constant compartment, plain) rate from the reaction rates
                lin.append((k, [o(j) for j in st['rids']], st['N'][i] / (st['V'][i] if s[0] == '[' else 1.0)))
            else:
                return None  # ponytail: eigen(), cc(), ee() ... selections -> roadrunner
        ident = I == list(range(len(outs))) and not (mul or div or lin or init)
        st['cols'][key] = r = (outs, (ident, I, mul, div, lin, init), st['m'].resolve_outputs(outs))
        return r

    def _rs_apply(self, plan, Y, pre=None):
        ident, I, mul, div, lin, init = plan
        if ident:
            return Y
        out = Y[:, I]
        for ops, f in ((mul, np.multiply), (div, np.divide)):
            if ops:
                k, j = zip(*ops)
                out[:, list(k)] = f(out[:, list(k)], Y[:, list(j)])
        for k, j, c in lin:
            out[:, k] = Y[:, j] @ c
        for (k, _), x in zip(init, pre or ()):
            out[:, k] = x
        return out

    def _rs_simulate(self, st, o, sel=None, integ=None):
        integ = self.integrator if integ is None else integ
        name = integ.getName()
        if o.output_file or name not in _RS_INTEGRATORS:
            return None
        cfg = integ.getSettingsMap()
        keys, vals = _RS_INTEGRATORS[name]
        if tuple(map(cfg.get, keys)) != vals:
            return None
        if len(o.times):
            times = np.asarray(o.times, dtype=float)
        else:
            tk = (o.start, o.end, o.steps)
            if st.get('tk') != tk:
                st['tk'], st['times'] = tk, np.linspace(o.start, o.end, o.steps + 1)
            times = st['times']
        sel = list(self.timeCourseSelections) if sel is None else sel
        cols = self._rs_cols(st, sel)
        if cols is None:
            return None
        if o.reset_model:
            self.reset()
        outs, plan, ores = cols
        pre = [roadrunner.RoadRunner.getValue(self, s) for _, s in plan[5]]  # init(...) values before the run
        if name == 'gillespie':
            r = self._rs_ssa(st, o, times, outs)
        else:
            try:
                rtol, atol = float(cfg['relative_tolerance']), float(cfg['absolute_tolerance'])
            except TypeError:  # per-species tolerance vector
                return None
            mdl = self.model
            x, asc = self._rs_x(st, mdl).tolist(), self._rs_asc(st, mdl)
            for a in (atol, min(atol, 1e-12)):  # a loose atol can let tiny species go astray: retry tighter
                try:
                    y, final, ev = st['m'].simulate_one(st['kres'], x, times, outputs=ores, rtol=rtol, atol=a,
                                                        event_state=st['ev'], final_array=True, atol_scale=asc)
                    break
                except RuntimeError:
                    pass
            else:
                return None
            r = (np.asarray(y), final, times[-1], ev)
        if r is None or not np.isfinite(r[0]).all():
            return None
        Y, final, tend, st['ev'] = r
        mdl = self.model
        self._rs_sync(st, final, mdl)
        mdl.setTime(tend)
        return _named(self._rs_apply(plan, Y, pre), cols=sel)

    def _rs_asc(self, st, mdl):
        """Per-state atol multipliers: roadrunner's CVODE multiplies a scalar atol by each species'
        initial amount (by its volume if that is zero), i.e. atol * |initial concentration| (atol if
        zero). Rate-rule states keep atol."""
        ic = mdl.getFloatingSpeciesInitConcentrations()
        k = ic.tobytes()
        if st.get('ick') != k:
            a = np.abs(ic[st['spi']])
            st['ick'], st['asc'][st['spj']] = k, np.where(a < 1e-30, 1.0, a)
        return st['asc']

    def _rs_ssa(self, st, o, times, outs):
        """Gillespie in Rust: fixed output points, or every reaction firing when variable_step_size.
        Returns (Y, final values, end time, event state) or None."""
        integ = self.integrator
        # roadrunner draws successive runs from one seeded stream: same seed -> same sequence of runs
        seed = int(integ.getValue('seed')) % 2 ** 64
        c = self.__dict__.get('_rs_seed')
        if c is None or c[0] != seed:
            c = self.__dict__['_rs_seed'] = [seed, 0]
        c[1] += 1
        run_seed = (seed + c[1] * 0x9E3779B97F4A7C15) % 2 ** 64
        m, x, k = st['m'], self._rs_x(st), st['keys']
        try:
            if integ.getValue('variable_step_size'):
                R = roadrunner_rs.gillespie_steps(m, k, x.tolist(), o.start, o.end, run_seed, outputs=outs + k,
                                                  max_rows=int(integ.getValue('max_output_rows')))
                Y, tend = R[:, 1:], R[-1, 0]
            else:
                Y, tend = m.gillespie_batch(k, x[None], times, seed=run_seed, threads=1, outputs=outs + k)[0], times[-1]
        except ValueError:  # rate rules (roadrunner refuses them too)
            return None
        return Y[:, :len(outs)], dict(zip(k, Y[-1, len(outs):])), tend, ()

    def timeSeriesSensitivities(self, start, stop, num, params=(), species=(), k=0):
        st = self._rs()
        r = None if st is None or k not in (0, 1) else self._rs_tss(st, start, stop, num, list(params), list(species), k)
        return super(ExtendedRoadRunner, self).timeSeriesSensitivities(start, stop, num, list(params), list(species), k) \
            if r is None else r
    timeSeriesSensitivities.__doc__ = roadrunner.RoadRunner.timeSeriesSensitivities.__doc__

    def _rs_tss(self, st, start, stop, num, params, species, k):
        """Forward sensitivities in Rust: exact symbolic derivatives, parameters integrated in parallel."""
        cma = self.conservedMoietyAnalysis
        cs = list(self.getConservedMoietyIds()) if cma else []
        # roadrunner: rows = global parameters (+ moiety totals, + rate-rule targets as zero rows);
        # columns = rate-rule states, then (independent) floating species, in model order
        rr_ids = [v for v in self.getRateRuleIds() if v in st['states']]
        base = rr_ids + (list(self.getIndependentFloatingSpeciesIds()) if cma else st['ids'])
        pkeys = set(st['keys'][:st['off']]) - set(st['bids'])
        if any(p not in pkeys and p not in cs for p in params) or \
                any(s in st['general'] or s not in st['states'] for s in base):
            return None
        params, species = params or st['gids'], [s for s in base if s in species or not species]
        solver = self.getSensitivitySolver()
        try:
            rtol, atol = float(solver.getValue('relative_tolerance')), float(solver.getValue('absolute_tolerance'))
        except TypeError:
            return None
        zero = set(rr_ids)
        rp = [p for p in params if p not in cs and p not in zero]
        dep = list(self.getDependentFloatingSpeciesIds()) if any(p in cs for p in params) else []
        times = np.linspace(start, stop, int(num))
        try:
            X, S = roadrunner_rs.time_sensitivities(st['m'], st['keys'], self._rs_x(st).tolist(), times, rp + dep,
                                                    rtol=rtol, atol=atol, order=k)
        except ValueError:  # integration failure
            return None
        rows = []
        if dep:  # moiety total T_k moves the dependent species: d conc_dep / d T = inv(Gamma[:, dep]) / V_dep
            G = self.getConservationMatrix()
            gcols = list(G.colnames)
            C = np.linalg.inv(np.asarray(G)[:, [gcols.index(d) for d in dep]]) / st['V'][[st['ids'].index(d) for d in dep]][:, None]
        for p in params:
            if p in zero:
                rows.append(np.zeros((len(times), S.shape[2])))
            elif p in cs:
                rows.append(np.einsum('tdx,d->tx', S[:, len(rp):], C[:, cs.index(p)]))
            else:
                rows.append(S[:, rp.index(p)])
        S = np.stack(rows, axis=1) if rows else S[:, :0]
        if k:
            S[0] = 0.0  # roadrunner reports a zero first row for k = 1 too
        self._rs_sync(st, dict(zip(st['states'], X[-1])))
        self.model.setTime(times[-1])
        idx = [st['states'].index(s) for s in species]
        V = np.array([st['V'][st['ids'].index(s)] if s in st['ids'] else 1.0 for s in species])
        return times, S[:, :, idx] * V, params, species  # roadrunner reports amount sensitivities

    def oneStep(self, currentTime, stepSize, reset=True):
        st = self._rs()
        integ = self.integrator
        if st is None or integ.getName() != 'cvode' or \
                any(integ.getValue(k) != v for k, v in zip(*_RS_INTEGRATORS['cvode'])):
            return super(ExtendedRoadRunner, self).oneStep(currentTime, stepSize, reset)
        try:
            rtol, atol = float(integ.getValue('relative_tolerance')), float(integ.getValue('absolute_tolerance'))
            _, final, st['ev'] = st['m'].simulate_one(st['keys'], self._rs_x(st).tolist(),
                                                      np.array([currentTime, currentTime + stepSize]),
                                                      rtol=rtol, atol=atol, event_state=st['ev'],
                                                      atol_scale=self._rs_asc(st, self.model))
        except (TypeError, RuntimeError):
            return super(ExtendedRoadRunner, self).oneStep(currentTime, stepSize, reset)
        self._rs_sync(st, final)
        self.model.setTime(currentTime + stepSize)
        return currentTime + stepSize

    def getSimulationData(self):
        d = self.__dict__.get('_rs_data')
        return super(ExtendedRoadRunner, self).getSimulationData() if d is None else d

    def steadyState(self, dict=None):
        st = self._rs()
        if st is None or dict is not None or st['nss'] or self.steadyStateSolver.getName() != 'nleq2':
            return super(ExtendedRoadRunner, self).steadyState(dict)
        S, v, _, ok = st['m'].steady_state_batch(st['keys'], self._rs_x(st)[None])
        if not ok[0]:
            return super(ExtendedRoadRunner, self).steadyState()
        self._rs_sync(st, {k: val for k, val in zip(st['states'], S[0])})
        dx = (st['N'] @ v[0][st['rx']] / st['V']) if st['plain'] else 0.0
        return float(np.linalg.norm(dx))

    def _rs_mca(self):
        """(st, unscaled CCC, unscaled FCC, v, S) at steady state, model left there (as roadrunner does)."""
        st = self._rs_mcast()
        if st is None:
            return None
        res = st['m'].mca_batch(st['keys'], self._rs_x(st)[None])
        if not res['ok'][0]:
            return None
        sp, rx = st['sp'], st['rx']
        self._rs_sync(st, dict(zip(st['states'], res['S'][0])))
        return st, res['uccc'][0][sp][:, rx], res['ufcc'][0][rx][:, rx], res['v'][0][rx], res['S'][0][sp]

    def getUnscaledConcentrationControlCoefficientMatrix(self):
        r = self._rs_mca()
        if r is None:
            return super(ExtendedRoadRunner, self).getUnscaledConcentrationControlCoefficientMatrix()
        return _named(r[1], r[0]['ids'], r[0]['rids'])

    def getScaledConcentrationControlCoefficientMatrix(self):
        r = self._rs_mca()
        if r is None:
            return super(ExtendedRoadRunner, self).getScaledConcentrationControlCoefficientMatrix()
        st, C, _, v, S = r
        return _named(C * v[None, :] / S[:, None], st['ids'], st['rids'])

    def getUnscaledFluxControlCoefficientMatrix(self):
        r = self._rs_mca()
        if r is None:
            return super(ExtendedRoadRunner, self).getUnscaledFluxControlCoefficientMatrix()
        return _named(r[2], r[0]['rids'], r[0]['rids'])

    def getScaledFluxControlCoefficientMatrix(self):
        r = self._rs_mca()
        if r is None:
            return super(ExtendedRoadRunner, self).getScaledFluxControlCoefficientMatrix()
        st, _, F, v, _ = r
        vi = np.where(np.abs(v) < self.getFluxThreshold(), 0.0, v)  # roadrunner zeroes near-zero flux rows
        with np.errstate(divide='ignore', invalid='ignore'):
            out = np.where(vi[:, None] != 0, F * v[None, :] / vi[:, None], 0.0)
        return _named(out, st['rids'], st['rids'])

    def _rs_rd(self):
        """(st, v, dv/dS, S) for the reactions at the current state (exact derivatives)."""
        st = self._rs()
        if st is None or not st['plain']:
            return None
        x = self._rs_x(st)
        if 'rd' not in st:  # compiled once per Rust model
            st['rd'] = roadrunner_rs.RateDerivatives(st['m'])
        v, dx, _, _ = st['rd'](st['keys'], x.tolist(), self.model.getTime())
        return st, v[st['rx']], dx[st['rx']][:, st['sp']], x[st['off']:]

    def getUnscaledElasticityMatrix(self):
        r = self._rs_rd()
        if r is None:
            return super(ExtendedRoadRunner, self).getUnscaledElasticityMatrix()
        return _named(r[2], r[0]['rids'], r[0]['ids'])

    def getScaledElasticityMatrix(self):
        r = self._rs_rd()
        if r is None:
            return super(ExtendedRoadRunner, self).getScaledElasticityMatrix()
        st, v, E, S = r
        with np.errstate(divide='ignore', invalid='ignore'):
            return _named(E * S[None, :] / v[:, None], st['rids'], st['ids'])

    def getSteadyStateValues(self):
        if self._rs() is None:
            return super(ExtendedRoadRunner, self).getSteadyStateValues()
        self.steadyState()  # Rust when it can, roadrunner otherwise
        return np.array([self.getValue(s) for s in self.steadyStateSelections], dtype=float)

    def getSteadyStateValuesNamedArray(self):
        if self._rs() is None:
            return super(ExtendedRoadRunner, self).getSteadyStateValuesNamedArray()
        return _named(self.getSteadyStateValues()[None], cols=self.steadyStateSelections)

    def steadyStateNamedArray(self, dict=None):
        if self._rs() is None or dict is not None:
            return super(ExtendedRoadRunner, self).steadyStateNamedArray(dict)
        return self.getSteadyStateValuesNamedArray()

    def mcaSteadyState(self, *args):
        if self._rs_mcast() is None:
            return super(ExtendedRoadRunner, self).mcaSteadyState(*args)
        return self.steadyState()

    def _rs_jac(self):
        """(st, Jacobian of the floating species in concentrations) at the current state and time."""
        st = self._rs_mcast()
        if st is None:
            return None
        sp = st['sp']
        J = np.asarray(st['m'].jacobian(st['keys'], self._rs_x(st).tolist(), self.model.getTime()))
        return st, J[sp][:, sp]

    def getFullJacobian(self):
        r = self._rs_jac()
        if r is None:
            return super(ExtendedRoadRunner, self).getFullJacobian()
        return _named(r[1], r[0]['ids'], r[0]['ids'])

    def _rs_rjac(self):
        r = self._rs_jac()
        if r is None:
            return None
        st, J = r
        if not self.conservedMoietyAnalysis:
            return st['ids'], J
        # d(dependent conc)/d(independent conc): the link matrix is in amounts
        V, ind = st['V'], [st['ids'].index(s) for s in st['ind']]
        Lc = st['L'] / V[:, None] * V[ind][None, :]
        return st['ind'], (J @ Lc)[ind]

    def getReducedJacobian(self, h=-1.0):
        r = self._rs_rjac()
        if r is None:
            return super(ExtendedRoadRunner, self).getReducedJacobian(h)
        return _named(r[1], r[0], r[0])

    def getFullEigenValues(self):
        r = self._rs_jac()
        if r is None:
            return super(ExtendedRoadRunner, self).getFullEigenValues()
        return np.linalg.eigvals(r[1]).astype(complex)

    def getReducedEigenValues(self):
        r = self._rs_rjac()
        if r is None:
            return super(ExtendedRoadRunner, self).getReducedEigenValues()
        return np.linalg.eigvals(r[1]).astype(complex)

    def getFullEigenValuesNamedArray(self):
        r = self._rs_jac()
        if r is None:
            return super(ExtendedRoadRunner, self).getFullEigenValuesNamedArray()
        e = np.linalg.eigvals(r[1])
        return _named(np.column_stack([e.real, e.imag]), r[0]['ids'], ['real', 'imaginary'])

    def getReducedEigenValuesNamedArray(self):
        r = self._rs_rjac()
        if r is None:
            return super(ExtendedRoadRunner, self).getReducedEigenValuesNamedArray()
        e = np.linalg.eigvals(r[1])
        return _named(np.column_stack([e.real, e.imag]), r[0], ['real', 'imaginary'])

    def _rs_ucc(self, variableName, parameterName):
        """(unscaled CC, parameter value, variable value now, is flux) with roadrunner's 5-point stencil."""
        st = self._rs_mcast()
        var = variableName.replace('[', '').replace(']', '')
        if st is None or st['kpos'].get(parameterName, st['off']) >= st['off'] or \
                (var not in st['ids'] and var not in st['rids']):
            return None  # species / moiety totals as parameters -> roadrunner
        x = self._rs_x(st)
        j = st['kpos'][parameterName]
        p0 = x[j]
        h = self.getDiffStepSize() * p0
        if abs(h) < 1e-12:
            h = self.getDiffStepSize()
        X = np.repeat(x[None], 4, axis=0)
        X[:, j] += [h, 2 * h, -h, -2 * h]
        S, v, _, ok = st['m'].steady_state_batch(st['keys'], X)
        if not ok.all():
            return None
        if var in st['ids']:
            i = st['ids'].index(var)
            f, now = S[:, st['sp'][i]], self.model.getFloatingSpeciesConcentrations()[i]
        else:
            i = st['rids'].index(var)
            f, now = v[:, st['rx'][i]], self.model.getReactionRates()[i]
        return (f[3] + 8 * f[0] - 8 * f[2] - f[1]) / (12 * h), p0, now, var in st['rids']

    def getuCC(self, variableName, parameterName):
        r = self._rs_ucc(variableName, parameterName)
        return super(ExtendedRoadRunner, self).getuCC(variableName, parameterName) if r is None else r[0]

    def getCC(self, variableName, parameterName):
        r = self._rs_ucc(variableName, parameterName)
        if r is None:
            return super(ExtendedRoadRunner, self).getCC(variableName, parameterName)
        u, p, now, flux = r
        if flux and abs(now) < self.getFluxThreshold():
            return 0.0
        with np.errstate(divide='ignore', invalid='ignore'):
            return float(np.float64(u) * p / now)

    # ---------------------------------------------------------------------
    # Model access
    # ---------------------------------------------------------------------
    # def getBoundarySpeciesConcentrations(self):
    #     return self.model.getBoundarySpeciesConcentrations()
    # getBoundarySpeciesConcentrations.__doc__ = roadrunner.ExecutableModel.getBoundarySpeciesConcentrations.__doc__

    # These model functions are attached after class creation
    _model_functions = [
        'getBoundarySpeciesIds',
        'getBoundarySpeciesConcentrationIds',
        'getBoundarySpeciesConcentrations',
        'getBoundarySpeciesAmounts',
        'getNumBoundarySpecies',

        'getFloatingSpeciesIds',
        'getFloatingSpeciesConcentrationIds',
        'getFloatingSpeciesConcentrations',
        'getFloatingSpeciesAmounts',
        'getNumFloatingSpecies',

        'getGlobalParameterIds',
        'getGlobalParameterValues',
        'getNumGlobalParameters',

        'getCompartmentIds',
        'getCompartmentVolumes',
        'getNumCompartments',

        'getConservedMoietyValues',
        'getNumConservedMoieties',

        'getReactionIds',
        'getReactionRates',
        'getNumReactions',
        'getNumEvents',
        'getNumRateRules'
     ]

    # ---------------------------------------------------------------------
    # Jarnac compatibility layer
    # ---------------------------------------------------------------------
    def fjac(self):
        return self.getFullJacobian()
    fjac.__doc__ = roadrunner.RoadRunner.getFullJacobian.__doc__

    def sm(self):
        return self.getFullStoichiometryMatrix()
    sm.__doc__ = roadrunner.RoadRunner.getFullStoichiometryMatrix.__doc__

    def rs(self):
        return self.model.getReactionIds()
    rs.__doc__ = roadrunner.ExecutableModel.getReactionIds.__doc__

    def fs(self):
        return self.model.getFloatingSpeciesIds()
    fs.__doc__ = roadrunner.ExecutableModel.getFloatingSpeciesIds.__doc__

    def bs(self):
        return self.model.getBoundarySpeciesIds()
    bs.__doc__ = roadrunner.ExecutableModel.getBoundarySpeciesIds.__doc__

    def ps(self):
        return self.model.getGlobalParameterIds()
    ps.__doc__ = roadrunner.ExecutableModel.getGlobalParameterIds.__doc__

    def vs(self):
        return self.model.getCompartmentIds()
    vs.__doc__ = roadrunner.ExecutableModel.getCompartmentIds.__doc__

    def dv(self):
        return self.getRatesOfChange()
    dv.__doc__ = roadrunner.RoadRunner.getRatesOfChange.__doc__

    def rv(self):
        return self.model.getReactionRates()
    rv.__doc__ = roadrunner.ExecutableModel.getReactionRates.__doc__

    def sv(self):
        return self.model.getFloatingSpeciesConcentrations()
    sv.__doc__ = roadrunner.ExecutableModel.getFloatingSpeciesConcentrations.__doc__

    # ---------------------------------------------------------------------
    # Export Utilities
    # ---------------------------------------------------------------------
    def __getSBML(self, current):
        if current is True:
            return self.getCurrentSBML()
        else:
            return self.getSBML()

    def getAntimony(self, current=False, removeFunctionDefinitions=None):
        """ Antimony string of the original model loaded into roadrunner.

        :param current: return current model state
        :param removeFunctionDefinitions: remove function definitions
        :type current: bool
        :return: Antimony
        :rtype: str
        """
        sbml = self.__getSBML(current)
        from .. import sbmlToAntimony
        return sbmlToAntimony(sbml, removeFunctionDefinitions)

    def getCurrentAntimony(self, removeFunctionDefinitions=None):
        """ Antimony string of the current model state.

        See also: :func:`getAntimony`
        :returns: Antimony string
        :rtype: str
        """
        return self.getAntimony(current=True)

    def getCellML(self, current=False):
        """ CellML string of the original model loaded into roadrunner.

        :param current: return current model state
        :type current: bool
        :returns: CellML string
        :rtype: str
        """
        sbml = self.__getSBML(current)
        from .. import sbmlToCellML
        return sbmlToCellML(sbml)

    def getCurrentCellML(self):
        """ CellML string of current model state.

        See also: :func:`getCellML`
        :returns: CellML string
        :rtype: str
        """
        return self.getCellML(current=True)

    def getMatlab(self, current=False):
        """ Matlab string of the original model loaded into roadrunner.

        See also: :func:`getCurrentMatlab`
        :returns: Matlab string
        :rtype: str
        """
        try:
            from sbml2matlab import sbml2matlab
            sbml = self.__getSBML(current)
            return sbml2matlab.sbml2matlab(sbml)
        except ImportError:
            warnings.warn("'sbml2matlab' could not be imported, no support for Matlab code generation",
                          RuntimeWarning, stacklevel=2)
            return ""

    def getCurrentMatlab(self):
        """ Matlab string of current model state.

        :param current: return current model state
        :type current: bool
        :returns: Matlab string
        :rtype: str
        """
        return self.getMatlab(current=True)

    def exportToSBML(self, filePath, current=True):
        """ Save current model as SBML file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of SBML file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.__getSBML(current))

    def exportToAntimony(self, filePath, current=True):
        """ Save current model as Antimony file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of Antimony file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getAntimony(current))

    def exportToCellML(self, filePath, current=True):
        """ Save current model as CellML file.

        :param current: export current model state
        :type current: bool
        :param filePath: file path of CellML file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getCellML(current))

    def exportToMatlab(self, filePath, current=True):
        """ Save current model as Matlab file.
        To save the original model loaded into roadrunner use
        current=False.

        :param self: RoadRunner instance
        :type self: RoadRunner.roadrunner
        :param filePath: file path of Matlab file
        :type filePath: str
        """
        with open(filePath, 'w', encoding="utf-8") as f:
            f.write(self.getMatlab(current))

    # ---------------------------------------------------------------------
    # Plotting Utilities
    # ---------------------------------------------------------------------
    def draw(self, **kwargs):
        """ Draws an SBMLDiagram of the current model.

        To set the width of the output plot provide the 'width' argument.
        Species are drawn as white circles (boundary species
        shaded in blue), reactions as grey squares.
        Currently only the drawing of medium-size networks is supported.
        """
        import shutil
        # PATH variable may not be present - cannot use os.environ['PATH']
        # but can use shutil.which for Python 3.3+
        if hasattr(shutil, 'which'):
            if shutil.which('dot') is None:
                warnings.warn("Graphviz is not installed in your machine or could not be found. 'draw' command cannot produce a diagram.",
                    Warning, stacklevel=2)
                return
        elif not 'PATH' in os.environ or any([ os.access( os.path.join( p, 'dot' ), os.X_OK ) for p in os.environ['PATH'].split( os.pathsep )]):
            warnings.warn("Graphviz is not installed in your machine or could not be found. 'draw' command cannot produce a diagram.",
                Warning, stacklevel=2)
            return

        from tellurium import SBMLDiagram
        diagram = SBMLDiagram(self.getSBML())
        diagram.draw(**kwargs)


    def plot(self, result=None, show=True,
             xlabel=None, ylabel=None, title=None, linewidth=2, xlim=None, ylim=None, logx=False, logy=False,
             xscale='linear', yscale='linear', grid=False, ordinates=None, tag=None, labels=None, 
             figsize=(6,4), savefig=None, dpi=80, alpha=1.0, **kwargs):
        """ Plot roadrunner simulation data.

        Plot is called with simulation data to plot as the first argument. If no data is provided the data currently
        held by roadrunner generated in the last simulation is used. The first column is considered the x axis and
        all remaining columns the y axis.
        If the result array has no names, then the current r.selections are used for naming. In this case the
        dimension of the r.selections has to be the same like the number of columns of the result array.

        Curves are plotted in order of selection (columns in result).

        In addition to the listed keywords plot supports all matplotlib.pyplot.plot keyword arguments,
        like color, alpha, linewidth, linestyle, marker, ...
        ::

            sbml = te.getTestModel('feedback.xml')
            r = te.loadSBMLModel(sbml)
            s = r.simulate(0, 100, 201)
            r.plot(s, loc="upper right", linewidth=2.0, lineStyle='-', marker='o', markersize=2.0,
                   alpha=0.8, title="Feedback Oscillation", xlabel="time", ylabel="concentration",
                   xlim=[0,100], ylim=[-1, 4])

        :param result: results data to plot (numpy array)
        :param show: show=True (default) shows the plot, use show=False to plot multiple simulations in one plot
        :type show: bool
        :param xlabel: x-axis label
        :type xlabel: str
        :param ylabel: y-axis label
        :type ylabel: str
        :param title: plot title
        :type title: str
        :param linewidth: linewidth of the plot
        :type linewidth: float
        :param xlim: limits on x-axis (tuple [start, end])
        :param ylim: limits on y-axis
        :param logx: use log scale for x-axis
        :type logx: bool
        :param logy: use log scale for y-axis
        :type logy: bool
        :param xscale: 'linear' or 'log' scale for x-axis
        :param yscale: 'linear' or 'log' scale for y-axis
        :param grid: show grid
        :type grid: bool
        :param ordinates: If supplied, only these selections will be plotted (see RoadRunner selections)
        :param tag: If supplied, all traces with the same tag will be plotted with the same color/style
        :param labels: 'id' to use species IDs
        :param figsize: If supplied, customize the size of the figure (width,height)
        :param savefig: If supplied, saves the figure to specified location
        :type savefig: str
        :param dpi: Change the dpi of the saved figure
        :type dpi: int
        :param alpha: Change the alpha value of the figure
        :type alpha: float
        :param kwargs: additional matplotlib keywords like marker, lineStyle, ...
        """
        if result is None:
            simData = self.getSimulationData()
            result = copy.copy(simData)
            result.colnames = simData.colnames

        from .. import getPlottingEngine

        if xlabel:
            kwargs['xlabel'] = xlabel
        if ylabel:
            kwargs['ylabel'] = ylabel
        if title:
            kwargs['title'] = title
        if linewidth:
            kwargs['linewidth'] = linewidth
        if xlim:
            kwargs['xlim'] = xlim
        if ylim:
            kwargs['ylim'] = ylim
        if logx:
            kwargs['logx'] = logx
        if logy:
            kwargs['logy'] = logy
        if xscale:
            kwargs['xscale'] = xscale
        if yscale:
            kwargs['yscale'] = yscale
        if grid:
            kwargs['grid'] = grid
        if ordinates:
            kwargs['ordinates'] = ordinates
        if tag:
            kwargs['tag'] = tag
        if labels:
            kwargs['labels'] = labels
        if figsize:
            kwargs['figsize'] = figsize
        if savefig:
            kwargs['savefig'] = savefig
        if dpi:
            kwargs['dpi'] = dpi
        if alpha:
            kwargs['alpha'] = alpha            

        # FIXME: provide the additional parameters to the plotting engine

        engine = getPlottingEngine()
        if show:
            # show the plot immediately
            engine.plotTimecourse(result, **kwargs)
        else:
            # otherwise, accumulate the traces
            engine.accumulateTimecourse(result, **kwargs)

    def show(self, reset=True):
        """ Show plot. """
        from .. import getPlottingEngine
        getPlottingEngine().show(reset=reset)

    # ---------------------------------------------------------------------
    # Stochastic Simulation Methods
    # ---------------------------------------------------------------------

    @property
    def seed(self):
        """ Getter for global Configuration seed. """
        return self.getSeed()

    @seed.setter
    def seed(self, value):
        """ Setter for model and global Configuration option seed. """
        return self.setSeed(value)

    def gillespie(self, *args, **kwargs):
        """ Run a Gillespie stochastic simulation.

        Sets the integrator to gillespie and performs simulation.
        ::

            rr = te.loada ('S1 -> S2; k1*S1; k1 = 0.1; S1 = 40')
            # Simulate from time zero to 40 time units using variable step sizes (classic Gillespie)
            result = rr.gillespie (0, 40)
            # Simulate on a grid with 10 points from start 0 to end time 40
            rr.reset()
            result = rr.gillespie (0, 40, 10)
            # Simulate from time zero to 40 time units using variable step sizes with given selection list
            # This means that the first column will be time and the second column species S1
            rr.reset()
            result = rr.gillespie (0, 40, selections=['time', 'S1'])
            # Simulate on a grid with 20 points from time zero to 40 time units
            # using the given selection list
            rr.reset()
            result = rr.gillespie (0, 40, 20, ['time', 'S1'])
            rr.plot(result)

        :param seed: seed for gillespie
        :type seed: int
        :param args: parameters for simulate
        :param kwargs: parameters for simulate
        :returns: simulation results
        """

        integratorName = self.integrator.getName()
        vss = self.integrator.variable_step_size
        if (integratorName != 'gillespie'):
            self.setIntegrator('gillespie')
            self.integrator.variable_step_size = True        
        if (len(args) > 3):
            self.integrator.variable_step_size = False
        if (len(args) == 3):
            if (type(args[2]) != list):
                self.integrator.variable_step_size = False
        elif ('points' in kwargs or 'steps' in kwargs):
            self.integrator.variable_step_size = False
        if "seed" in kwargs:
            self.integrator.setValue("seed", kwargs["seed"])
            del kwargs["seed"]

        s = self.simulate(*args, **kwargs)
        self.setIntegrator(integratorName)
        self.integrator.variable_step_size = vss
        return s


def _resets(f):  # reset() / resetAll() restart event triggers: so does the Rust engine
    def g(self, *args, **kwargs):
        c = self.__dict__.get('_rs_cache')
        if c is not None and c[1] is not None:
            c[1]['ev'] = ()
        return f(self, *args, **kwargs)
    g.__name__, g.__doc__ = f.__name__, f.__doc__
    return g


for _n in ('reset', 'resetAll'):
    setattr(ExtendedRoadRunner, _n, _resets(getattr(roadrunner.RoadRunner, _n)))


# Structure-changing roadrunner methods drop the cached Rust model (rebuilt on next use)
def _invalidates(f):
    def g(self, *args, **kwargs):
        self.__dict__.pop('_rs_cache', None)
        return f(self, *args, **kwargs)
    g.__name__, g.__doc__ = f.__name__, f.__doc__
    return g


for _n in dir(roadrunner.RoadRunner):
    if _n.startswith(('add', 'remove', 'load', 'clear', 'regenerate', 'resetToOrigin', 'setBoundary', 'setConstant',
                      'setHasOnlySubstanceUnits', 'setInitStoichiometryValue', 'setKineticLaw', 'setReversible')):
        setattr(ExtendedRoadRunner, _n, _invalidates(getattr(roadrunner.RoadRunner, _n)))

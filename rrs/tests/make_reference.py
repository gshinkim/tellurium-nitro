"""Record libroadrunner's answers for test_core_vs_roadrunner.py into reference/<model>.npz.

    python rrs/tests/make_reference.py [model ...]     (needs libroadrunner; only missing entries are computed)
    python rrs/tests/make_reference.py --fixed          (fixed-step reference courses)
"""
import os, sys
import numpy as np
import roadrunner

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_core_vs_roadrunner as T  # noqa: E402

MODEL_FNS = {'getFloatingSpeciesIds', 'getBoundarySpeciesIds', 'getGlobalParameterIds', 'getCompartmentIds',
             'getReactionIds', 'getGlobalParameterValues', 'getReactionRates', 'getFloatingSpeciesConcentrations'}


class RR:
    """The adapter on libroadrunner (same names as test_core_vs_roadrunner.Core)."""
    def __init__(self, sbml):
        self.r = roadrunner.RoadRunner(sbml)
    def cma(self, on): self.r.conservedMoietyAnalysis = on
    def tol(self, rtol, atol):
        self.r.integrator.relative_tolerance, self.r.integrator.absolute_tolerance = rtol, atol
    def __getattr__(self, n):
        f = getattr(self.r.model, n) if n in MODEL_FNS else getattr(self.r, n)
        if not callable(f):
            return f
        def call(*a):
            v = f(*a)
            return (np.asarray(v), list(v.rownames), list(v.colnames)) if hasattr(v, 'rownames') else v
        return call
    timeCourseSelections = property(lambda s: s.r.timeCourseSelections)
    steadyStateSelections = property(lambda s: s.r.steadyStateSelections)
    def getConservedMoietyIds(self): return list(self.r.model.getConservedMoietyIds())
    def getRatesOfChangeIds(self): return list(self.r.getRatesOfChangeNamedArray().colnames)
    cms = getConservedMoietyIds
    def fconc(self): return self.r.model.getFloatingSpeciesConcentrations()
    def finit(self): return self.r.model.getFloatingSpeciesInitConcentrations()
    def _m(self, v): return np.asarray(v), list(v.rownames), list(v.colnames)
    def E(self, scaled): return self._m(self.r.getScaledElasticityMatrix() if scaled else self.r.getUnscaledElasticityMatrix())
    def CC(self, scaled):
        return self._m(self.r.getScaledConcentrationControlCoefficientMatrix() if scaled else self.r.getUnscaledConcentrationControlCoefficientMatrix())
    def FCC(self, scaled):
        return self._m(self.r.getScaledFluxControlCoefficientMatrix() if scaled else self.r.getUnscaledFluxControlCoefficientMatrix())
    def sim(self, a, b, n): return np.asarray(self.r.simulate(a, b, n))
    def sids(self): return list(self.r.model.getStoichiometryIds())
    def stoich(self, i): return self.r.model.getStoichiometry(i)
    def setStoichiometry(self, i, v): self.r.model.setStoichiometry(i, v)
    def ids(self, t, model): return list((self.r.model if model else self.r).getIds(t))
    def npw(self): return self.r.model.getNumPiecewiseTriggers()
    def events(self): return list(self.r.model.getEventIds())
    def settings(self, name):
        self.r.setIntegrator(name)
        return [(k, self.r.integrator.getValue(k)) for k in self.r.integrator.getSettings()]
    def ss_set(self, solver, settings):
        self.r.setSteadyStateSolver(solver)
        for k, v in settings.items():
            setattr(self.r.steadyStateSolver, k, v)


def load(path):
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=True) as z:
        return dict(z)


def fixed_references():
    import test_fixed_vs_roadrunner as F

    def simulate(sbml, sel, method, settings, *args, **kwargs):
        r = roadrunner.RoadRunner(sbml)
        r.setIntegrator(method)
        for k, v in settings.items():
            r.integrator.setValue(k, v)
        r.timeCourseSelections = sel
        return np.asarray(r.simulate(*args, **kwargs))

    ref = load(F.REF)
    cv = dict(relative_tolerance=1e-12, absolute_tolerance=1e-24, maximum_num_steps=10**7)
    for name, sbml, end, n in F.models():
        key = f'{name}/ids'
        if key not in ref:
            rr = roadrunner.RoadRunner(sbml)
            ref[key] = np.array(rr.model.getFloatingSpeciesIds())
        ids = ref[key]
        sel = ['time'] + [f'[{i}]' for i in ids]
        k = f'{name}/cvode'
        if k not in ref:
            ref[k] = simulate(sbml, sel, 'cvode', cv, 0, end, n)
        for label, integ, settings, method, _ in F.CASES:
            key = f'{name}/{label}'
            if key not in ref:
                if method == 'rk45' and name in F.TIGHT:
                    ref[key] = ref[k]
                else:
                    stop, points = F.fixed_grid(name, method, end, n)
                    ref[key] = simulate(sbml, sel, integ, settings, 0, stop, points)
        if name in F.ANT:
            for i, (eps, hmax, rows) in enumerate(F.VARIABLE):
                key = f'{name}/variable{i}'
                if key not in ref:
                    ref[key] = simulate(sbml, sel, 'rk45', dict(epsilon=eps, maximum_time_step=hmax, max_output_rows=rows), 0, 10)
                if key + '/cvode' not in ref:
                    ref[key + '/cvode'] = simulate(sbml, sel, 'cvode', cv, times=ref[key][:, 0].tolist())
        print(name, flush=True)
    np.savez_compressed(F.REF, **ref)


class ExtraRecorder(T.Rec):
    """Use the saved values to drive each scenario and record only missing answers."""
    def __init__(self, existing):
        super().__init__()
        self.existing = existing

    def __call__(self, key, fn, *args, **kwargs):
        if key not in self.existing:
            return super().__call__(key, fn, *args, **kwargs)
        value = self.existing[key]
        if value.ndim == 0 and str(value).startswith('ERR:'):
            return None
        return fn() if key.endswith('/cma') else value


def extra_metadata(ref, name, sbml):
    """Record compatibility details once, including labels preserved by NamedArray."""
    def add(key, get):
        if key not in ref:
            ref[key] = np.asarray(get())

    for cma in (False, True):
        prefix = 'extra/%d/' % cma
        if str(ref.get(prefix + 'cma', '')) == 'ERR:refuse':
            continue
        r = RR(sbml)
        r.cma(cma)
        for matrix in ('getKMatrix', 'getExtendedStoichiometryMatrix',
                       'getReactantsStoichiometryMatrix', 'getProductsStoichiometryMatrix'):
            rows, cols = (prefix + matrix + '/' + label for label in ('rownames', 'colnames'))
            if rows not in ref or cols not in ref:
                result = getattr(r.r, matrix)()
                add(rows, lambda: result.rownames)
                add(cols, lambda: result.colnames)
        if name == 'chain' and not cma:
            add(prefix + 'supported', r.r.getSupportedIdTypes)
            add(prefix + 'model_supported', r.r.model.getSupportedIdTypes)
            add(prefix + 'ids_all', lambda: list(r.r.getIds(-1)))
            add(prefix + 'threshold', r.r.getSteadyStateThreshold)
            add(prefix + 'selection_flags', lambda: [
                (key, str(getattr(roadrunner.SelectionRecord, key)))
                for key in dir(roadrunner.SelectionRecord) if key.isupper()
            ])
            model = r.r.model
            add(prefix + 'stoich_index', lambda: [
                model.getStoichiometryIndex(s) for s in [*r.sids(), 'missing']
            ])
            add(prefix + 'initial_stoich', lambda: [
                model.getInitStoichiometry(i) for i in range(len(r.sids()))
            ])
            add(prefix + 'stoich_pair', lambda: [
                [model.getStoichiometry(i, j) for j in range(len(r.getReactionIds()))]
                for i in range(len(r.getFloatingSpeciesIds()))
            ])
            add(prefix + 'initial_stoich_pair', lambda: [
                [model.getInitStoichiometry(i, j) for j in range(len(r.getReactionIds()))]
                for i in range(len(r.getFloatingSpeciesIds()))
            ])
            add(prefix + 'stoich_pair_index', lambda: [
                [model.getStoichiometryIndex(sp, rx) for rx in r.getReactionIds()]
                for sp in r.getFloatingSpeciesIds()
            ])
            def set_stoich():
                result = model.setStoichiometry(0, 2.5)
                current = model.getStoichiometry(0)
                initial = model.getInitStoichiometry(0)
                r.r.reset()
                after_reset = model.getStoichiometry(0)
                r.r.resetAll()
                after_reset_all = model.getStoichiometry(0)
                r.r.resetToOrigin()
                return [result, current, initial, after_reset, after_reset_all,
                        model.getStoichiometry(0)]
            add(prefix + 'set_stoich_reset', set_stoich)

            def frequency_labels():
                result = r.r.getFrequencyResponse(0.1, 2, 5, 'k0', 'S1', False, False)
                return list(getattr(result, 'colnames', ()))
            add(prefix + 'frequency_columns', frequency_labels)

            def reset_seed_state():
                r.r.setSeed(42)
                r.r.simulate(0, 1, 2)
                r.r.resetSeed()
                return [r.r.getCurrentTime(), r.r['S1']]
            add(prefix + 'reset_seed_state', reset_seed_state)

            def internal_step_support():
                out = []
                for name in ('cvode', 'rk45', 'gillespie', 'rk4', 'euler'):
                    r.r.setIntegrator(name)
                    r.r.resetToOrigin()
                    try:
                        r.r.internalOneStep(0, 1)
                    except (ValueError, RuntimeError):
                        out.append(False)
                    else:
                        out.append(True)
                return out
            add(prefix + 'internal_step_support', internal_step_support)


if __name__ == '__main__':
    if sys.argv[1:] == ['--fixed']:
        fixed_references()
        sys.exit()
    # libroadrunner is slow: only what is missing is computed; existing entries are kept as they are.
    # ss_settings: per case. Models: only those without a reference file, or named on the command line
    # (their scenarios rerun in order, since each call depends on the state the previous ones left;
    # new keys are added, existing ones kept).
    os.makedirs(T.REF, exist_ok=True)
    path = os.path.join(T.REF, 'ss_settings.npz')
    cases = load(path).get('cases', np.array({}, dtype=object)).item()
    todo = [c for c in T.ss_cases() if c[0] not in cases]
    for key, ant, solver, st in todo:
        cases[key] = T.ss_run(RR, T.rn.antimony_to_sbml(ant), solver, st)
    if todo:
        np.savez_compressed(path, cases=np.array(cases, dtype=object))
    print('ss_settings', len(todo), 'new cases', flush=True)
    metadata = '--extra-metadata' in sys.argv
    extra = '--extra' in sys.argv or metadata  # add test_extra_vs_roadrunner.py's keys to the existing files
    names = [a for a in sys.argv[1:] if a not in ('--extra', '--extra-metadata')]
    for name, sbml in T.MODELS:
        if names and name not in names:
            continue
        path = os.path.join(T.REF, name + '.npz')
        old = load(path)
        if extra:
            import test_extra_vs_roadrunner as X
            if not metadata:
                rec = ExtraRecorder(old)
                X.run_extra(RR, sbml, name, rec)
                old.update({k: v for k, v in rec.out.items() if k not in old})
            extra_metadata(old, name, sbml)
            np.savez_compressed(path, **old)
            print(name, 'extra entries recorded', flush=True)
            continue
        if old and not names:
            continue
        rec = T.Rec()
        T.run_all(RR, sbml, name, rec)
        new = {k: v for k, v in rec.out.items() if k not in old}
        np.savez_compressed(path, **{**old, **new})
        print(name, len(new), 'new entries', flush=True)

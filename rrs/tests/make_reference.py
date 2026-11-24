"""Record libroadrunner's answers for test_core_vs_roadrunner.py into reference/<model>.npz.

    python rrs/tests/make_reference.py [model ...]     (needs libroadrunner; slow on the big models)
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
    os.makedirs(T.REF, exist_ok=True)
    metadata = '--extra-metadata' in sys.argv
    extra = '--extra' in sys.argv or metadata  # add test_extra_vs_roadrunner.py's keys to the existing files
    names = [a for a in sys.argv[1:] if a not in ('--extra', '--extra-metadata')]
    for name, sbml in T.MODELS:
        if names and name not in names:
            continue
        if extra:
            import test_extra_vs_roadrunner as X
            p = os.path.join(T.REF, name + '.npz')
            with np.load(p) as z:
                old = dict(z)
            if not metadata:
                rec = ExtraRecorder(old)
                X.run_extra(RR, sbml, name, rec)
                old.update({k: v for k, v in rec.out.items() if k not in old})
            extra_metadata(old, name, sbml)
            np.savez_compressed(p, **old)
            print(name, 'extra entries recorded', flush=True)
            continue
        rec = T.Rec()
        T.run_all(RR, sbml, name, rec)
        np.savez_compressed(os.path.join(T.REF, name + '.npz'), **rec.out)
        print(name, len(rec.out), 'entries', flush=True)

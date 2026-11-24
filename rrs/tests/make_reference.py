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
    for name, sbml in T.MODELS:
        path = os.path.join(T.REF, name + '.npz')
        old = load(path)
        if old and name not in sys.argv[1:]:
            continue
        rec = T.Rec()
        T.run_all(RR, sbml, name, rec)
        new = {k: v for k, v in rec.out.items() if k not in old}
        np.savez_compressed(path, **{**old, **new})
        print(name, len(new), 'new entries', flush=True)

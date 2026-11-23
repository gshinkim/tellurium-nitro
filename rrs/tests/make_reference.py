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
    def settings(self, name):
        self.r.setIntegrator(name)
        return [(k, self.r.integrator.getValue(k)) for k in self.r.integrator.getSettings()]


if __name__ == '__main__':
    os.makedirs(T.REF, exist_ok=True)
    for name, sbml in T.MODELS:
        if len(sys.argv) > 1 and name not in sys.argv[1:]:
            continue
        rec = T.Rec()
        T.run_all(RR, sbml, name, rec)
        np.savez_compressed(os.path.join(T.REF, name + '.npz'), **rec.out)
        print(name, len(rec.out), 'entries', flush=True)

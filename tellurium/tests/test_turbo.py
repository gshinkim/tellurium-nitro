"""te.turbo(): the Rust engine behind simulate/steadyState/MCA must match plain roadrunner."""
import unittest
from unittest import mock

import numpy as np
import tellurium as te
from tellurium.roadrunner import ExtendedRoadRunner

try:
    import roadrunner_rs  # noqa: F401
except ImportError:
    roadrunner_rs = None

# conserved moiety (S4+S5), non-unit compartment, an amount-only species
MOIETY = '''
compartment c = 2; S1 in c; S2 in c; S3 in c; S4 in c; S5 in c
J0: $X0 -> S1; k0*X0
J1: S1 -> S2; k1*S1
J2: S2 -> S3; k2*S2*S4
J3: S3 -> $X1; k3*S3
J4: S4 -> S5; k4*S4
J5: S5 -> S4; k5*S5*S2
X0 = 1; k0 = 1; k1 = 0.5; k2 = 0.3; k3 = 0.2; k4 = 0.4; k5 = 0.1
S1 = 1; S2 = 1; S3 = 1; S4 = 3; S5 = 1
'''

# assignment rule in a rate law, Michaelis-Menten, two compartments
RULES = '''
compartment a = 1.5, b = 0.5; S1 in a; S2 in b; S3 in b
J0: $X0 -> S1; k0*X0
J1: S1 -> S2; Vm*S1/(Km + S1)
J2: S2 -> S3; k2*S2
J3: S3 -> $X1; k3*S3
Vm := 2*k1
X0 = 3; k0 = 0.4; k1 = 1; Km = 0.5; k2 = 0.7; k3 = 0.9; S1 = 0; S2 = 0; S3 = 0
'''

# events Rust carries across simulate calls, a rate rule, rule/rate/parameter selections
EVENTS_RULES = '''
compartment c = 2; S1 in c; S2 in c; substanceOnly species S3 in c
J1: S1 -> S2; k1*S1
J2: S2 -> S1 + S3; k2*S2
x := 2*k1
p' = 0.1*S1
k1 = 0.5; k2 = 0.2; S1 = 10; S2 = 0; S3 = 0; p = 0
at time > 3: S1 = S1 + 5, k1 = 0.8
at S1 < 4: k2 = 0.05
'''

EVENT = '''
J1: S1 -> S2; k1*S1
J2: S2 -> S1; k2*S2
k1 = 0.5; k2 = 0.2; S1 = 10; S2 = 0
at time > 3, t0 = false: k1 = 2  # t0=false: Rust can't carry trigger state -> roadrunner
'''


def chain(n):
    lines = ["J0: $X0 -> S1; k0*X0"]
    for i in range(1, n):
        lines.append(f"J{i}: S{i} -> S{i+1}; k{i}*S{i} - k{i}r*S{i+1}")
    lines.append(f"J{n}: S{n} -> $X1; k{n}*S{n}")
    lines += [f"B{i}: S{i} -> $W; kb{i}*S{i}/(1+S{i})" for i in range(5, n, 5)]
    lines.append("X0=10; k0=1")
    lines += [f"k{i}={1+0.1*(i%7)}; k{i}r=0.2; S{i}=1" for i in range(1, n + 1)]
    lines += [f"kb{i}=0.3" for i in range(5, n, 5)]
    return "\n".join(lines)


GETTERS = ['getUnscaledElasticityMatrix', 'getScaledElasticityMatrix', 'getFullJacobian', 'getReducedJacobian',
           'getFullEigenValues', 'getReducedEigenValues',
           'getUnscaledConcentrationControlCoefficientMatrix', 'getScaledConcentrationControlCoefficientMatrix',
           'getUnscaledFluxControlCoefficientMatrix', 'getScaledFluxControlCoefficientMatrix']
SUPER = ['_simulate', 'steadyState', 'getuCC', 'getCC'] + GETTERS


def _raise(*args, **kwargs):
    raise AssertionError("turbo fell back to roadrunner")


@unittest.skipIf(roadrunner_rs is None, "roadrunner_rs not built")
class TestTurbo(unittest.TestCase):

    def setUp(self):
        self.was = te.turbo()

    def tearDown(self):
        te.turbo(self.was)

    def both(self, f, strict=True):
        """f() with turbo off (roadrunner) and on; with strict, roadrunner's own methods are made to raise."""
        te.turbo(False)
        ref = f()
        te.turbo(True)
        if not strict:
            return ref, f()
        with mock.patch.multiple(ExtendedRoadRunner.__bases__[0],  # the SWIG class
                                  **{n: _raise for n in SUPER}):
            return ref, f()

    def same(self, a, b, rtol=1e-5, atol=1e-8, msg=None):
        if hasattr(b, 'colnames'):
            self.assertEqual(list(a.colnames), list(b.colnames), msg)
            self.assertEqual(list(getattr(a, 'rownames', [])), list(getattr(b, 'rownames', [])), msg)
            self.assertIs(type(a), type(b), msg)
        np.testing.assert_allclose(np.asarray(a), np.asarray(b), rtol=rtol, atol=atol, err_msg=str(msg))

    def test_simulate(self):
        for ant in (MOIETY, RULES, chain(12)):
            for cma in (False, True):
                def f():
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    sel = ['time'] + r.getFloatingSpeciesIds()[:2] + ['[%s]' % r.getFloatingSpeciesIds()[-1]]
                    out = [r.simulate(0, 10, 51), r.simulate(10, 20, 11, sel), r.simulate(times=[20, 21, 25]),
                           r.simulate(25, 30, steps=5, selections=sel)]
                    r.timeCourseSelections = sel[::-1]
                    out += [r.simulate(30, 35), np.array(r.getSimulationData()),  # copy: roadrunner's is a view into r r.model.getFloatingSpeciesConcentrations(),
                            [r.getValue('time')]]
                    return out
                for i, (a, b) in enumerate(zip(*self.both(f))):
                    self.same(a, b, rtol=1e-4, atol=1e-7, msg=(cma, i))

    def test_events_rules_outputs(self):
        def f():
            r = te.loada(EVENTS_RULES)
            sel = ['time', 'S1', '[S3]', 'S3', 'x', 'p', 'J1', 'k1', 'c']
            out = [r.simulate(0, 2, 11, sel), r.simulate(2, 6, 21, sel), r.simulate(6, 30, 25, sel)]
            return out + [[r.k1, r.k2, r.p, r.S1, r['[S3]'], r.getValue('time')]]
        ref, got = self.both(f)
        for i, (a, b) in enumerate(zip(ref, got)):
            self.same(a, b, rtol=1e-4, atol=1e-7, msg=i)

    def test_state_sync(self):
        def f():
            r = te.loada(MOIETY)
            a = r.simulate(0, 5, 11)
            vals = [r.getValue('S1'), r['[S2]'], r.S5, r.getValue('time')]
            r.k1 = 2.0  # changed parameter is honoured by the next simulate
            b = r.simulate(5, 10, 11)
            r.reset()
            return a, vals, b, r.simulate(0, 1, 3)
        ref, got = self.both(f)
        for a, b in zip(ref, got):
            self.same(a, b, rtol=1e-4, atol=1e-7)

    def test_steady_state_and_getters(self):
        for ant in (MOIETY, RULES, chain(12)):
            for cma in (False, True):
                def f():
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    r.setDiffStepSize(1e-4)  # roadrunner's FD is the inexact side
                    out = [r.steadyState() < 1e-6, r.model.getFloatingSpeciesConcentrations()]
                    out += [getattr(r, g)() for g in GETTERS]
                    r.reset()
                    out += [r.getFullJacobian(), r.getReducedJacobian()]  # off steady state too
                    out.append(r.getScaledFluxControlCoefficientMatrix())  # runs steadyState first
                    out.append(r.model.getFloatingSpeciesConcentrations())
                    r.reset()
                    p, s, j = r.getGlobalParameterIds()[0], r.getFloatingSpeciesIds()[0], r.getReactionIds()[1]
                    out += [[r.getCC(s, p), r.getuCC(s, p), r.getCC(j, p), r.getuCC('[%s]' % s, p)],
                            r.model.getFloatingSpeciesConcentrations()]
                    return out
                ref, got = self.both(f)
                for i, (a, b) in enumerate(zip(ref, got)):
                    if np.iscomplexobj(a):
                        a, b = (sorted(np.asarray(x), key=lambda z: (round(z.real, 6), round(z.imag, 6))) for x in (a, b))
                    self.same(a, b, rtol=1e-5, atol=1e-7, msg=(ant[:20], cma, i))

    def test_event_state_carries_over(self):  # t0=false trigger and a delayed event across simulate calls
        def f():
            r = te.loada(EVENT + '\nat 3 after time > 4: k2 = 1')
            out = [r.simulate(0, 5, 11), r.simulate(5, 10, 11)]
            r.reset()
            return out + [r.simulate(0, 10, 21)]
        ref, got = self.both(f)
        for a, b in zip(ref, got):
            self.same(a, b, rtol=1e-4, atol=1e-7)

    def test_fallback_settings(self):
        def f():
            r = te.loada(MOIETY)
            r.integrator.variable_step_size = True
            a = r.simulate(0, 5)
            r.integrator.variable_step_size = False
            r.setIntegrator('rk4')
            return a, r.simulate(0, 5, 11)
        ref, got = self.both(f, strict=False)
        for a, b in zip(ref, got):
            self.same(a, b, rtol=1e-12)

    def test_structure_change_invalidates(self):
        def f():
            r = te.loada(MOIETY)
            r.simulate(0, 1, 2)
            r.addParameter('k6', 0.3)
            r.addReaction('J6', ['S3'], ['S1'], 'k6*S3')
            r.X0 = 2.0
            r.reset()
            return r.simulate(0, 10, 11), r.getFullJacobian() if r.steadyState() < 1e-6 else None
        ref, got = self.both(f, strict=False)
        for a, b in zip(ref, got):
            self.same(a, b, rtol=1e-4, atol=1e-7)

    def test_switch(self):
        te.turbo(False)
        self.assertFalse(te.turbo())
        self.assertIsNone(te.loada(MOIETY)._rs())
        te.turbo(True)
        self.assertTrue(te.turbo())
        self.assertIsNotNone(te.loada(MOIETY)._rs())
        self.assertTrue(issubclass(type(te.loada(MOIETY)), ExtendedRoadRunner))


if __name__ == '__main__':
    unittest.main()


# ---------------------------------------------------------------------------------------------
# Sensitivities, elasticities, rates of change, steady-state values, eigenvalues, Gillespie, oneStep
# ---------------------------------------------------------------------------------------------
@unittest.skipIf(roadrunner_rs is None, "roadrunner_rs not built")
class TestTurboFacade(unittest.TestCase):
    setUp, tearDown, same = TestTurbo.setUp, TestTurbo.tearDown, TestTurbo.same

    def strict(self, f, names):
        """(roadrunner result, turbo result) with roadrunner's `names` made to raise under turbo."""
        te.turbo(False)
        ref = f()
        te.turbo(True)
        with mock.patch.multiple(ExtendedRoadRunner.__bases__[0], **{n: _raise for n in names}):
            return ref, f()

    def compare(self, ref, got, rtol=1e-5, atol=1e-8, msg=None):
        for i, (a, b) in enumerate(zip(ref, got)):
            if isinstance(a, list) and a and isinstance(a[0], str):
                self.assertEqual(a, b, (msg, i))
            else:
                self.same(a, b, rtol=rtol, atol=atol, msg=(msg, i))

    def test_time_series_sensitivities(self):
        for ant in (MOIETY, chain(8)):
            for cma in (False, True):
                def f():
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    t, s, p, x = r.timeSeriesSensitivities(0, 10, 11)
                    # independent species only: roadrunner leaves CMA-dependent ones slightly stale here
                    out = [t, s, p, x, r.getIndependentFloatingSpeciesConcentrationsV(), [r.getCurrentTime()]]
                    # fresh model: roadrunner carries its sensitivity state over between calls, even across reset()
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    sp = r.getFloatingSpeciesIds()
                    t, s, p, x = r.timeSeriesSensitivities(0, 5, 6, [], [sp[-1], sp[0]])  # columns: model order
                    return out + [t, s, p, x]
                ref, got = self.strict(f, ['timeSeriesSensitivities'])
                self.compare(ref, got, rtol=2e-4, atol=1e-6 * max(1, np.abs(ref[1]).max()), msg=(cma, ant[:10]))

    def test_elasticity_matrices(self):
        # scalar getters (getuEE, getRatesOfChange, ...) stay on roadrunner: ~1 us there
        names = ['getUnscaledElasticityMatrix', 'getScaledElasticityMatrix']
        for ant in (MOIETY, RULES, EVENTS_RULES.split('x :=')[0] + 'k1 = 0.5; k2 = 0.2; S1 = 10'):
            for cma in (False, True):
                def f():
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    for i, s in enumerate(r.getFloatingSpeciesIds()):  # off steady state, same state for both engines
                        r.setValue('[%s]' % s, 0.3 + 0.7 * i)
                    return [r.getUnscaledElasticityMatrix(), r.getScaledElasticityMatrix()]
                ref, got = self.strict(f, names)
                self.compare(ref, got, rtol=1e-5, atol=1e-7, msg=(ant[:10], cma))

    def test_steady_state_values_and_eigenvalues(self):
        names = ['getSteadyStateValues', 'getSteadyStateValuesNamedArray', 'steadyStateNamedArray', 'mcaSteadyState',
                 'getFullEigenValuesNamedArray', 'getReducedEigenValuesNamedArray', 'steadyState']
        for ant in (MOIETY, RULES):
            for cma in (False, True):
                def f():
                    r = te.loada(ant)
                    r.conservedMoietyAnalysis = cma
                    r.steadyStateSelections = ['[S2]', 'J1'] + r.getFloatingSpeciesIds()[:2]
                    out = [r.getSteadyStateValues()]
                    r.reset()
                    out += [r.getSteadyStateValuesNamedArray()]
                    r.reset()
                    out += [r.steadyStateNamedArray(), [r.mcaSteadyState() < 1e-6]]
                    for e in (r.getFullEigenValuesNamedArray(), r.getReducedEigenValuesNamedArray()):
                        a = np.array(sorted(map(tuple, np.round(np.asarray(e), 9))))
                        out += [list(e.colnames), len(e.rownames), a]
                    return out
                ref, got = self.strict(f, names)
                self.compare(ref, got, rtol=1e-6, atol=1e-8, msg=(ant[:10], cma))

    def test_gillespie(self):
        ant = 'compartment c = 2; S1 in c; S2 in c; J1: S1 -> S2; k1*S1; J2: S2 -> ; k2*S2; S1 = 50; k1 = 0.5; k2 = 0.2'

        def runs(n, seed=None, **kw):
            r = te.loada(ant)
            r.setIntegrator('gillespie')
            if seed is not None:
                r.integrator.seed = seed
            out = []
            for _ in range(n):
                r.reset()
                out.append(r.simulate(**kw))
            return r, out

        te.turbo(True)
        with mock.patch.multiple(ExtendedRoadRunner.__bases__[0], _simulate=_raise):
            r, fixed = runs(400, seed=7, start=0, end=4, points=5, selections=['time', 'S1', '[S2]'])
            _, again = runs(400, seed=7, start=0, end=4, points=5, selections=['time', 'S1', '[S2]'])
            r2, var = runs(50, seed=3, start=0, end=4)
        self.assertEqual(list(fixed[0].colnames), ['time', 'S1', '[S2]'])
        np.testing.assert_array_equal(np.array(fixed), np.array(again))  # same seed, same runs
        self.assertGreater(len({tuple(np.asarray(a)[:, 1]) for a in fixed}), 300)  # successive runs differ
        a = np.array(fixed)
        np.testing.assert_array_equal(a[..., 1], np.round(a[..., 1]))  # amounts are whole molecules
        self.assertEqual(r.getCurrentTime(), 4.0)
        self.assertEqual(r.getValue('S1'), fixed[-1][-1, 1])  # state synced to the last row
        te.turbo(False)
        _, ref = runs(400, start=0, end=4, points=5, selections=['time', 'S1', '[S2]'])
        ref = np.array(ref)
        se = np.sqrt(ref.var(0) / len(ref) + a.var(0) / len(a)) + 1e-12
        self.assertLess(np.max(np.abs(ref.mean(0) - a.mean(0)) / se), 5)
        for v in var:  # variable step: a row at start, after every firing, at end
            v = np.asarray(v)
            self.assertEqual((v[0, 0], v[-1, 0]), (0.0, 4.0))
            self.assertTrue((np.diff(v[:, 0]) >= 0).all() and len(v) > 10)
            self.assertTrue((np.abs(np.diff(v[:, 1:] * 2, axis=0)).sum(1)[:-1] > 0.5).all())  # one firing per row
        self.assertAlmostEqual(np.mean([np.asarray(v)[-1, 1] for v in var]), 50 * np.exp(-1), delta=2.0)  # rate law in concentrations: k_eff = k1 / 2

    def test_one_step_and_integrator_settings(self):
        def f():
            r = te.loada(MOIETY)
            t = 0.0
            for _ in range(5):
                t = r.oneStep(t, 0.5)
            out = [[t], r.model.getFloatingSpeciesConcentrations()]
            r.reset()
            r.integrator.stiff = False
            r.integrator.maximum_num_steps = 5000
            return out + [r.simulate(0, 5, 11)]
        ref, got = self.strict(f, ['oneStep', '_simulate'])
        self.compare(ref, got, rtol=1e-4, atol=1e-7)


def _edit(ant, f):
    """SBML of `ant` after f(libsbml model) (features antimony cannot write)."""
    import libsbml
    d = libsbml.readSBMLFromString(te.antimonyToSBML(ant))
    f(d.getModel())
    return libsbml.writeSBMLToString(d)


def _cf(m, model_level, constant=True):
    p = m.createParameter()
    p.setId('cf'), p.setValue(2.0), p.setConstant(constant)
    (m if model_level else m.getSpecies('S2')).setConversionFactor('cf')


def _vstoich(m, how):
    import libsbml
    sr = m.getReaction('J1').getProduct(0)
    sr.setId('n2'), sr.setConstant(False)
    if how == 'event':
        e = m.createEvent()
        e.setUseValuesFromTriggerTime(True)
        t = e.createTrigger()
        t.setMath(libsbml.parseL3Formula('time > 3')), t.setPersistent(True), t.setInitialValue(True)
        ea = e.createEventAssignment()
        ea.setVariable('n2'), ea.setMath(libsbml.parseL3Formula('3'))
    else:
        rule = m.createAssignmentRule() if how == 'rule' else m.createRateRule()
        rule.setVariable('n2'), rule.setMath(libsbml.parseL3Formula('1 + time/5' if how == 'rule' else '0.2'))


def _stoich_math():
    import libsbml
    d = libsbml.readSBMLFromString(te.antimonyToSBML("J1: S1 -> S2; k*S1; k=1; S1=10; S2=0; n = 2"))
    d.setLevelAndVersion(2, 4, False)
    d.getModel().getReaction('J1').getProduct(0).createStoichiometryMath().setMath(libsbml.parseL3Formula('n + time/10'))
    return libsbml.writeSBMLToString(d)


BASE = "J1: S1 -> S2; k*S1; k=1; S1=10; S2=0"
VC = "compartment c=2; species S1 in c, S2 in c; J1: S1 -> S2; k*S1; k=1; S1=10; S2=0; "
# name -> (antimony or SBML, extra outputs as (roadrunner selection, roadrunner_rs output))
FEATURES = {
    'compartment rate rule': (VC + "c' = 1", ['c']),
    'compartment assignment rule': (VC + "c := 2 + sin(time)*0.5 + k", ['c']),
    'compartment event': (VC + "at time > 1: c = 4", ['c']),
    'compartment and species event': (VC + "at time > 1: S1 = 7, c = 4", ['c']),
    'compartment rule, event on its parameter': (VC + "c := 1 + kc; kc = 1; at time > 1: kc = 3", ['c']),
    'compartment rate rule and event': (VC + "c' = 0.2; at time > 1: c = c*2", ['c']),
    'amount species, variable compartment': (
        "compartment c=2; species S1 in c; substanceOnly species S2 in c; J1: S1 -> S2; k*S1*S2; k=0.1; "
        "S1=10; S2=1; c' = 0.3", ['c']),
    'boundary species, variable compartment': (
        "compartment c=2; species S1 in c, $B in c; J1: B -> S1; k*B; k=0.1; S1=10; B=1; c' = 0.3", ['c', 'B']),
    'rate-rule species, variable compartment': (
        "compartment c=2; species S1 in c, S2 in c; S1' = -0.5*S1; J1: -> S2; 1; S1=10; c' = 0.3", ['c']),
    'conversion factor (species)': (_edit(BASE, lambda m: _cf(m, False)), []),
    'conversion factor (model)': (_edit(BASE, lambda m: _cf(m, True)), []),
    'conversion factor (event-changed)': (_edit(BASE + "; at time>2: k=2", lambda m: (
        _cf(m, True, False), m.getEvent(0).createEventAssignment().setVariable('cf'),
        m.getEvent(0).getEventAssignment(1).setMath(__import__('libsbml').parseL3Formula('3')))), []),
    'stoichiometry: assignment rule': (_edit(BASE, lambda m: _vstoich(m, 'rule')), []),
    'stoichiometry: rate rule': (_edit(BASE, lambda m: _vstoich(m, 'rate')), []),
    'stoichiometry: event': (_edit(BASE, lambda m: _vstoich(m, 'event')), []),
    'stoichiometryMath (L2)': (_stoich_math(), []),
    'assignment-rule species': ("species S3, S4; " + BASE + "; S3 := 2*S1; S4 := time", ['S3']),
    'assignment-rule amount species': (
        "compartment c=3; substanceOnly species S3 in c; species S1 in c, S2 in c; " + BASE + "; S3 := 2*S1*time",
        [('[S3]', 'S3')]),
    'rules in reverse order': ("J1: S1 -> S2; a*S1; S1=10; a := b/2; b := c + 1; c := 0.5*time", ['a', 'b']),
    'function definitions (nested, swapped args)': (
        "function f(x, y)\n x*y\nend\nfunction g(y, x)\n f(x+1, y)\nend\nJ1: S1 -> S2; g(k, S1)*0.1; k=1; S1=10", []),
    'avogadro': ("J1: S1 -> S2; k*S1*avogadro/6.02214179e23; k=1; S1=10", []),
    'rateOf': (BASE + "; p := rateOf(S1); x = 0; x' = time; q := rateOf(x)*2", ['p', 'q']),
    'rateOf in kinetic law': ("J1: S1 -> S2; k*S1; J2: -> S3; -rateOf(S1); k=1; S1=10", []),
    'reaction without kinetic law': ("J1: S1 -> S2; ; J2: S1 -> S3; k*S1; k=1; S1=10", ['J1']),
    'reaction id in math': (BASE + "; p := J1*2", ['p']),
}

EVENTS = '''
J1: S1 -> S2; k1*S1
J2: S2 -> S1; k2*S2
k1 = 0.5; k2 = 0.2; S1 = 10; S2 = 0; p = 0; q = 0
E1: at time > 1.5, t0 = false: k1 = 2
E2: at S1 < 5: k2 = k2*1.5
E3: at 1.0 after time > 2.9, fromTrigger = false: p = p + S2
E4: at 0.7 after S2 > 6, persistent = false: q = q + 1
E5: at time > 3.5, priority = 2: S1 = S1 + 3
E6: at time > 3.5, priority = 1: S1 = S1 * 2
E7: at 2 after time > 4.5: S2 = S2 + 1, p = 100
'''


@unittest.skipIf(roadrunner_rs is None, "roadrunner_rs not built")
class TestBridgeFeatures(unittest.TestCase):
    """from_roadrunner + roadrunner_rs vs plain libroadrunner, one SBML feature at a time."""

    def rr(self, src):
        import roadrunner
        return roadrunner.roadrunner.RoadRunner(src if src.lstrip().startswith('<') else te.antimonyToSBML(src))

    def test_features(self):
        from tellurium.analysis.rrs_bridge import from_roadrunner
        for name, (src, extra) in FEATURES.items():
            with self.subTest(name):
                r = self.rr(src)
                m = from_roadrunner(r)
                extra = [e if isinstance(e, tuple) else (e, e) for e in extra]
                ids = [s for s in r.model.getFloatingSpeciesIds()]
                ref = r.simulate(0, 10, 101, ['time'] + ['[%s]' % s for s in ids] + [a for a, _ in extra])
                got, _ = m.simulate_one([], [], np.linspace(0, 10, 101), outputs=['time'] + ids + [b for _, b in extra])
                np.testing.assert_allclose(got, np.asarray(ref), rtol=1e-4, atol=1e-6, err_msg=name)

    def test_events_across_simulate_calls(self):
        from tellurium.analysis.rrs_bridge import from_roadrunner
        r = self.rr(EVENTS)
        m = from_roadrunner(r)
        sel, outs, keys = ['time', '[S1]', '[S2]', 'p', 'q', 'k1', 'k2'], ['time', 'S1', 'S2', 'p', 'q', 'k1', 'k2'], \
            ['S1', 'S2', 'k1', 'k2', 'p', 'q']
        state, x = (), None
        for a, b, n in [(0, 2, 21), (2, 3.5, 16), (3.5, 5, 7), (5, 6, 3), (6, 12, 61)]:  # events pending across calls
            ref = r.simulate(a, b, n, sel)
            y, fin, state = m.simulate_one(keys if x else [], x or [], np.linspace(a, b, n), outputs=outs,
                                           event_state=state)
            x = [fin[k] for k in keys]
            np.testing.assert_allclose(y, np.asarray(ref), rtol=1e-4, atol=1e-6, err_msg=str((a, b)))

    def test_mca_with_rate_rules(self):
        from tellurium.analysis.rrs_bridge import from_roadrunner
        r = self.rr('''J0: $X -> S1; k0*X; J1: S1 -> S2; k1*S1*p; J2: S2 -> S3; k2*S2; J3: S3 ->; k3*S3*S1
            X = 2; k0=1; k1=1; k2=0.7; k3=0.5; S1=1; S2=1; S3=1; p=1; p' = 0.1*(2 + S2 - p); q = 0; q' = 0.3*(1 - q*S1)''')
        m = from_roadrunner(r)
        res = m.mca_batch([], np.zeros((1, 0)))
        self.assertTrue(res['ok'][0])
        r.setDiffStepSize(1e-5)
        r.steadyState()
        np.testing.assert_allclose(res['S'][0][m.state_ids.index('p')], r.p, rtol=1e-8)
        for key, get in (('ccc', r.getScaledConcentrationControlCoefficientMatrix),
                         ('fcc', r.getScaledFluxControlCoefficientMatrix),
                         ('uccc', r.getUnscaledConcentrationControlCoefficientMatrix)):
            M = get()
            rows = [(m.reactions if key == 'fcc' else m.species).index(x) for x in M.rownames]
            cols = [m.reactions.index(x) for x in M.colnames]
            np.testing.assert_allclose(res[key][0][np.ix_(rows, cols)], np.asarray(M), rtol=1e-7, atol=1e-9, err_msg=key)


# ---------------------------------------------------------------------------------------------
# Facade coverage: rule species, variable compartments, boundary/init()/rate-of-change selections,
# rate-rule MCA and sensitivities, event sensitivities, variable-step Gillespie, sosbench cases
# ---------------------------------------------------------------------------------------------
VARS = '''
compartment c = 1; c' = 0.1
species S1 in c, S2 in c, S3 in c, $B in c
J1: $B -> S1; k1*B
J2: S1 -> S2; k2*S1
S3 := 2*S1
B = 2; k1 = 0.5; k2 = 0.3; S1 = 1; S2 = 0
E1: at time > 2, t0 = false: B = 5
E2: at 1 after time > 3: k2 = 0.6
'''

RATE = '''
compartment c = 2; S1 in c; S2 in c; S3 in c
J0: $X0 -> S1; k0*X0
J1: S1 -> S2; k1*S1
J2: S2 -> ; k2*S2*p
p' = 0.1 - 0.2*p*S1
S3' = 0.5 - S3
X0 = 1; k0 = 1; k1 = 2; k2 = 0.3; S1 = 1; S2 = 1; p = 1; S3 = 0
'''


@unittest.skipIf(roadrunner_rs is None, "roadrunner_rs not built")
class TestTurboCoverage(unittest.TestCase):
    setUp, tearDown, same = TestTurbo.setUp, TestTurbo.tearDown, TestTurbo.same
    strict, compare = TestTurboFacade.strict, TestTurboFacade.compare

    def test_rules_compartments_boundary_init(self):
        def f():
            r = te.loada(VARS)
            sel = ['time', 'S1', '[S1]', 'S3', '[S3]', 'B', '[B]', 'c', 'J1', 'k2', 'init(S1)', 'init([S2])']
            out = [r.simulate(0, 2.5, 6, sel), r.simulate(2.5, 6, 8, sel)]
            return out + [[r.S1, r.c, r.B, r.k2, r.getValue('S2')]]
        ref, got = self.strict(f, ['_simulate'])
        self.compare(ref, got, rtol=1e-4, atol=1e-7)

    def test_rate_of_change_selections(self):
        def f():
            return [te.loada(RATE).simulate(0, 1, 3, ['time', "S1'", "[S1]'", 'X0', '[X0]', 'init(p)', 'J1'])]
        ref, got = self.strict(f, ['_simulate'])
        self.compare(ref, got, rtol=1e-4, atol=1e-7)

    def test_rate_rule_mca_and_sensitivities(self):
        def f():
            r = te.loada(RATE)
            t, s, p, x = r.timeSeriesSensitivities(0, 2, 5)
            # fresh model for k=1: roadrunner carries its sensitivity state over between calls
            out = [s, p, x, te.loada(RATE).timeSeriesSensitivities(0, 2, 5, [], [], 1)[1]]
            r.reset()
            r.steadyState()
            return out + [r.getFullJacobian(), r.getReducedJacobian(), r.getUnscaledElasticityMatrix(),
                          r.getScaledConcentrationControlCoefficientMatrix(), r.getScaledFluxControlCoefficientMatrix(),
                          r.model.getFloatingSpeciesConcentrations(), [r.p, r.S3], [r.getuCC('S1', 'k1'), r.getCC('J1', 'k2')]]
        names = ['timeSeriesSensitivities', 'steadyState', 'getFullJacobian', 'getReducedJacobian', 'getUnscaledElasticityMatrix',
                 'getScaledConcentrationControlCoefficientMatrix', 'getScaledFluxControlCoefficientMatrix', 'getuCC', 'getCC']
        ref, got = self.strict(f, names)
        self.compare(ref, got, rtol=1e-3, atol=1e-5)

    def test_sensitivities_across_events(self):
        # S1' = -k1 S1, S1 += 3 at t = 1: dS1/dk1 = dS1(1)/dk1 e^{-k1 (t-1)} - (t-1) S1(t) after the event
        # (roadrunner's CVODES run drifts ~1% from this after the event, so compare to the exact value)
        te.turbo(True)
        r = te.loada('J1: S1 -> ; k1*S1; S1 = 10; k1 = 0.5\nat time > 1: S1 = S1 + 3')
        with mock.patch.multiple(ExtendedRoadRunner.__bases__[0], timeSeriesSensitivities=_raise):
            t, s, p, x = r.timeSeriesSensitivities(0, 3, 7, ['k1'])
        S = np.where(t <= 1, 10 * np.exp(-0.5 * t), (10 * np.exp(-0.5) + 3) * np.exp(-0.5 * (t - 1)))
        ds = np.where(t <= 1, -t * S, -10 * np.exp(-0.5) * np.exp(-0.5 * (t - 1)) - (t - 1) * S)
        np.testing.assert_allclose(s[:, 0, 0], ds, rtol=1e-4, atol=1e-6)

    def test_gillespie_variable_step_events_rules(self):
        ant = 'J1: S1 -> S2; k1*S1; S3 := 2*S1; S1 = 40; k1 = 0.5\nat time > 1: S1 = S1 + 20'
        te.turbo(True)
        r = te.loada(ant)
        r.setIntegrator('gillespie')
        r.integrator.seed = 5
        with mock.patch.multiple(ExtendedRoadRunner.__bases__[0], _simulate=_raise):
            a = np.asarray(r.simulate(0, 3, selections=['time', 'S1', 'S3', 'k1']))
        self.assertEqual((a[0, 0], a[-1, 0]), (0.0, 3.0))
        np.testing.assert_array_equal(a[:, 2], 2 * a[:, 1])  # rule column
        jump = np.flatnonzero(np.diff(a[:, 1]) > 0)  # the event adds 20 molecules once
        self.assertEqual(len(jump), 1)
        self.assertGreater(a[jump[0] + 1, 0], 1.0)

    def test_sosbench_cases_run_in_rust(self):
        import os
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'benchmarks'))
        import bench
        for case in ('00009', '00014', '00022', '00033', '10001', 'repressilator'):
            cfg = {}
            for line in open(os.path.join(bench.HERE, 'sosbench', case, '%s-settings.txt' % case)):
                k, _, v = line.partition(':')
                cfg[k.strip()] = v.strip()
            atol = float(cfg['absolute'])
            out = {}
            for turbo in (False, True):
                te.turbo(turbo)
                run, _ = bench.sos(case)
                r = [c.cell_contents for c in run.__closure__ if isinstance(c.cell_contents, ExtendedRoadRunner)][0]
                try:
                    out[turbo] = run()
                except RuntimeError:  # roadrunner's CVODE gives up on the repressilator (CV_TOO_MUCH_WORK)
                    out[turbo] = None
                if turbo:
                    self.assertIsNotNone(r.__dict__.get('_rs_data'), case)
            self.assertTrue(np.isfinite(out[True]).all(), case)
            if out[False] is not None:  # both within their (loose) tolerances of each other
                np.testing.assert_allclose(out[True], out[False], rtol=1e-4, atol=300 * atol, err_msg=case)


# ---------------------------------------------------------------------------------------------
# Lazy loading: te.loada/loadSBMLModel build only the Rust model; roadrunner compiles on demand
# ---------------------------------------------------------------------------------------------
IA = '''
compartment c = 2; S1 in c; S2 in c; S3 in c; substanceOnly species S4 in c; $B in c
J1: S1 -> S2; k1*S1*B; J2: S2 -> S3; k1*S2; J3: S3 -> S4; k2*S3; J4: S4 -> S1; 0.1*S4
k1 = 1; k2 = 3; S1 = 1; S2 = S1*2; S3 = k2 + S1; S4 = k2*S1; B = k2/2; p = k2*2; p' = -0.1*p; q = S2*2
'''


def _lazy_models():
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'benchmarks'))
    import bench
    out = {k: te.antimonyToSBML(v) for k, v in dict(MOIETY=MOIETY, RULES=RULES, EVENTS_RULES=EVENTS_RULES, IA=IA,
                                                    VARS=VARS, RATE=RATE, chain=chain(12)).items()}
    for case in ('00009', '00014', '00022', '00033', '10000', '10001', 'repressilator'):
        d = os.path.join(bench.HERE, 'sosbench', case)
        out[case] = open(os.path.join(d, [f for f in os.listdir(d) if f.endswith('.xml')][0])).read()
    return out


@unittest.skipIf(roadrunner_rs is None, "roadrunner_rs not built")
class TestLazy(unittest.TestCase):
    setUp, tearDown, same = TestTurbo.setUp, TestTurbo.tearDown, TestTurbo.same

    def values(self, r):
        m = r.model
        return [list(m.getFloatingSpeciesIds()), list(m.getBoundarySpeciesIds()), list(m.getGlobalParameterIds()),
                list(m.getCompartmentIds()), list(m.getReactionIds()), list(r.timeCourseSelections),
                list(r.steadyStateSelections), m.getFloatingSpeciesConcentrations(), m.getFloatingSpeciesAmounts(),
                m.getBoundarySpeciesConcentrations(), m.getBoundarySpeciesAmounts(), m.getGlobalParameterValues(),
                m.getCompartmentVolumes(), m.getReactionRates(), [r.getValue('time')]]

    def workflow(self, r, sbml):
        """What the lazy model serves without roadrunner; r is lazy or plain roadrunner."""
        import warnings
        out = self.values(r) + [r.getCurrentSBML(), r.getSBML()]
        m = r.model
        p, s = (list(m.getGlobalParameterIds()) or ['time'])[0], m.getFloatingSpeciesIds()[0]  # rule species last
        r.setValue(p, r.getValue(p) * 1.5)
        r.setValue('[%s]' % s, 0.7)
        out += self.values(r) + [r.getCurrentSBML()]
        T = 0.5 if 'lorenz' in sbml else 5
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out.append(r.simulate(0, T, 11))
            sel = ['time', s, '[%s]' % s, p] + list(m.getReactionIds())[:1]
            out += [r.simulate(T, 2 * T, 6, sel)] + self.values(r)
        r.reset()
        out += self.values(r)
        r.resetAll()
        return out + self.values(r)

    def compare(self, ref, got, msg):
        for i, (a, b) in enumerate(zip(ref, got)):
            if isinstance(a, str) or (isinstance(a, list) and a and isinstance(a[0], str)):
                self.assertEqual(a, b, (msg, i))
            else:
                self.same(a, b, rtol=2e-3, atol=1e-6, msg=(msg, i))

    def test_lazy_matches_roadrunner(self):
        import roadrunner
        te.turbo(True)
        base = roadrunner.roadrunner.RoadRunner
        for name, sbml in _lazy_models().items():
            with self.subTest(name):
                # VARS: reset() with a variable compartment needs roadrunner (loaded from the lazy state).
                # Compared with turbo on an eagerly loaded model: after Rust ran the events, roadrunner's
                # reset() doesn't see their t0=false trigger state, unlike plain roadrunner.
                ref = self.workflow(ExtendedRoadRunner(sbml) if name == 'VARS' else base(sbml), sbml)
                if name == 'VARS':
                    r = te.loadSBMLModel(sbml)
                    got = self.workflow(r, sbml)
                    self.assertNotIn('_lz', r.__dict__)
                else:
                    with mock.patch.object(base, 'load', _raise), mock.patch.object(base, '_load', _raise):
                        r = te.loadSBMLModel(sbml)
                        got = self.workflow(r, sbml)
                        self.assertIn('_lz', r.__dict__)  # roadrunner never compiled the model
                self.compare(ref, got, name)

    def test_attributes_mca_and_load_on_demand(self):
        import roadrunner
        te.turbo(True)
        for ant in (MOIETY, RULES, chain(12)):  # (roadrunner's Jacobian fails with initial assignments)
            ref = roadrunner.roadrunner.RoadRunner(te.antimonyToSBML(ant))
            with mock.patch.object(roadrunner.roadrunner.RoadRunner, '_load', _raise):
                r = te.loada(ant)
                s, p = r.getFloatingSpeciesIds()[0], r.getGlobalParameterIds()[0]
                self.assertEqual((getattr(r, s), getattr(r, s + '_amt'), getattr(r, p), r['[%s]' % s]),
                                 (ref['[%s]' % s], ref[s], ref[p], ref['[%s]' % s]))
                setattr(r, p, 0.3)
                ref.setValue(p, 0.3)
                r.simulate(0, 20, 3)  # near the steady state: Rust's Newton solver has no presimulation
                ref.simulate(0, 20, 3)
                r.steadyState()
                ref.steadyState()
                got = [r.model.getFloatingSpeciesConcentrations(), r.getScaledConcentrationControlCoefficientMatrix(),
                       r.getScaledFluxControlCoefficientMatrix(), r.getFullJacobian(), r.getReducedJacobian(),
                       r.getScaledElasticityMatrix(), r.getFullStoichiometryMatrix(), [r.getCC(s, p)],
                       r.getRatesOfChange()]
                sens = te.loada(ant).timeSeriesSensitivities(0, 2, 5, [p])[1]
                self.assertIn('_lz', r.__dict__)
            want = [ref.model.getFloatingSpeciesConcentrations(), ref.getScaledConcentrationControlCoefficientMatrix(),
                    ref.getScaledFluxControlCoefficientMatrix(), ref.getFullJacobian(), ref.getReducedJacobian(),
                    ref.getScaledElasticityMatrix(), ref.getFullStoichiometryMatrix(), [ref.getCC(s, p)],
                    ref.getRatesOfChange()]
            for i, (a, b) in enumerate(zip(want, got)):
                self.same(a, b, rtol=1e-4, atol=1e-7, msg=(ant[:20], i))
            fresh = roadrunner.roadrunner.RoadRunner(te.antimonyToSBML(ant))
            np.testing.assert_allclose(sens, fresh.timeSeriesSensitivities(0, 2, 5, [p])[1], rtol=1e-3, atol=1e-6)
            # a call Nitro doesn't serve compiles the model in roadrunner, from the current state
            ref.simulate(20, 23, 3)
            r.simulate(20, 23, 3)
            self.same(r.getLinkMatrix(), ref.getLinkMatrix())
            self.assertNotIn('_lz', r.__dict__)
            self.same(r.model.getFloatingSpeciesConcentrations(), ref.model.getFloatingSpeciesConcentrations(), rtol=1e-4)
            self.assertEqual((r.getValue(p), r.getValue('time')), (ref.getValue(p), ref.getValue('time')))
            self.same(r.simulate(23, 26, 7), ref.simulate(23, 26, 7), rtol=1e-4, atol=1e-7)

    def test_current_sbml_round_trip_and_switch(self):
        te.turbo(True)
        r = te.loada(IA)
        r.simulate(0, 3, 4)
        r2 = te.loadSBMLModel(r.getCurrentSBML())
        for a, b in zip(self.values(r)[:-1], self.values(r2)):  # not the time
            if not (isinstance(a, list) and a and isinstance(a[0], str)):
                self.same(a, b, rtol=1e-12)
        self.assertEqual(te.loada(IA).getAntimony(), te.loada(IA)._load_rr().getAntimony())
        te.turbo(False)
        self.assertNotIn('_lz', te.loada(IA).__dict__)  # turbo off: eager roadrunner, as before
        te.turbo(True)
        r = te.loada(MOIETY)
        r.conservedMoietyAnalysis = True  # libstruct's moieties: roadrunner
        self.assertNotIn('_lz', r.__dict__)
        r = te.loada(MOIETY)
        r.addReaction('J9', ['S1'], ['S3'], 'k1*S1')  # model edits: roadrunner, Rust model rebuilt
        self.assertNotIn('_lz', r.__dict__)
        self.assertIn('J9', r.getReactionIds())

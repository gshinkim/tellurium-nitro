"""tellurium_nitro's RoadRunner (Rust engine) against libroadrunner's recorded answers.

Each workflow below takes ``load`` (SBML -> runner) and returns a list of results. make_reference.py
runs them once on libroadrunner and saves its answers in reference/turbo.npz; this test replays them on
tellurium_nitro only (no libroadrunner needed). Where libroadrunner is approximate (finite-difference
elasticities and control coefficients) the workflow tightens its step so its answer is the accurate one.
"""
import os
import unittest

import numpy as np
import tellurium_nitro as te
from tellurium_nitro.roadrunner import ExtendedRoadRunner

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'reference', 'turbo.npz')
SOS = os.path.join(HERE, '..', '..', 'benchmarks', 'sosbench')

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

# events carried across simulate calls, a rate rule, rule/rate/parameter selections
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
at time > 3, t0 = false: k1 = 2
at 3 after time > 4: k2 = 1
'''

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

IA = '''
compartment c = 2; S1 in c; S2 in c; S3 in c; substanceOnly species S4 in c; $B in c
J1: S1 -> S2; k1*S1*B; J2: S2 -> S3; k1*S2; J3: S3 -> S4; k2*S3; J4: S4 -> S1; 0.1*S4
k1 = 1; k2 = 3; S1 = 1; S2 = S1*2; S3 = k2 + S1; S4 = k2*S1; B = k2/2; p = k2*2; p' = -0.1*p; q = S2*2
'''

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

ELAST = EVENTS_RULES.split('x :=')[0] + 'k1 = 0.5; k2 = 0.2; S1 = 10'


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


ANT = dict(MOIETY=MOIETY, RULES=RULES, EVENTS_RULES=EVENTS_RULES, EVENT=EVENT, VARS=VARS, RATE=RATE, IA=IA,
           EVENTS=EVENTS, ELAST=ELAST, chain=chain(12), chain8=chain(8))


def sbml_of(name):
    """SBML of a model: an Antimony model above or a sosbench case."""
    if name in ANT:
        return te.antimonyToSBML(ANT[name])
    d = os.path.join(SOS, name)
    return open(os.path.join(d, [f for f in os.listdir(d) if f.endswith('.xml')][0])).read()


GETTERS = ['getUnscaledElasticityMatrix', 'getScaledElasticityMatrix', 'getFullJacobian', 'getReducedJacobian',
           'getFullEigenValues', 'getReducedEigenValues',
           'getUnscaledConcentrationControlCoefficientMatrix', 'getScaledConcentrationControlCoefficientMatrix',
           'getUnscaledFluxControlCoefficientMatrix', 'getScaledFluxControlCoefficientMatrix']


def eig(e):  # eigenvalues in a fixed order
    e = np.asarray(e)
    return np.array(sorted(e, key=lambda z: (round(z.real, 6), round(z.imag, 6))))


# ---------------------------------------------------------------- workflows: load(sbml) -> runner

def w_simulate(load, name, cma):
    r = load(sbml_of(name))
    r.conservedMoietyAnalysis = cma
    sel = ['time'] + list(r.model.getFloatingSpeciesIds())[:2] + ['[%s]' % list(r.model.getFloatingSpeciesIds())[-1]]
    out = [r.simulate(0, 10, 51), r.simulate(10, 20, 11, sel), r.simulate(times=[20, 21, 25]),
           r.simulate(25, 30, steps=5, selections=sel)]
    r.timeCourseSelections = sel[::-1]
    return out + [r.simulate(30, 35), np.array(r.getSimulationData()), r.model.getFloatingSpeciesConcentrations(),
                  [r.getValue('time')]]


def w_events_rules(load):
    r = load(sbml_of('EVENTS_RULES'))
    sel = ['time', 'S1', '[S3]', 'S3', 'x', 'p', 'J1', 'k1', 'c']
    out = [r.simulate(0, 2, 11, sel), r.simulate(2, 6, 21, sel), r.simulate(6, 30, 25, sel)]
    return out + [[r.k1, r.k2, r.p, r.S1, r['[S3]'], r.getValue('time')]]


def w_state_sync(load):
    r = load(sbml_of('MOIETY'))
    a = r.simulate(0, 5, 11)
    vals = [r.getValue('S1'), r['[S2]'], r.S5, r.getValue('time')]
    r.k1 = 2.0
    b = r.simulate(5, 10, 11)
    r.reset()
    return [a, vals, b, r.simulate(0, 1, 3)]


def w_steady_state(load, name, cma):
    r = load(sbml_of(name))
    r.conservedMoietyAnalysis = cma
    r.setDiffStepSize(1e-5)  # libroadrunner's finite differences: small step, so its answer is the accurate one
    out = [[float(r.steadyState() < 1e-6)], r.model.getFloatingSpeciesConcentrations()]
    for g in GETTERS:
        v = getattr(r, g)()
        out.append(eig(v) if 'Eigen' in g else v)
    r.reset()
    out += [r.getFullJacobian(), r.getReducedJacobian()]  # off steady state too
    out.append(r.getScaledFluxControlCoefficientMatrix())  # runs steadyState first
    out.append(r.model.getFloatingSpeciesConcentrations())
    r.reset()
    p, s, j = r.model.getGlobalParameterIds()[0], r.model.getFloatingSpeciesIds()[0], r.model.getReactionIds()[1]
    return out + [[r.getCC(s, p), r.getuCC(s, p), r.getCC(j, p), r.getuCC('[%s]' % s, p)],
                  r.model.getFloatingSpeciesConcentrations()]


def w_event_state(load):
    r = load(sbml_of('EVENT'))
    out = [r.simulate(0, 5, 11), r.simulate(5, 10, 11)]
    r.reset()
    return out + [r.simulate(0, 10, 21)]


def w_rk4(load):
    r = load(sbml_of('MOIETY'))
    r.setIntegrator('rk4')
    return [r.simulate(0, 5, 11)]


def w_edit(load):
    r = load(sbml_of('MOIETY'))
    r.simulate(0, 1, 2)
    r.addParameter('k6', 0.3)
    r.addReaction('J6', ['S3'], ['S1'], 'k6*S3')
    r.X0 = 2.0
    r.reset()
    return [r.simulate(0, 10, 11), r.getFullJacobian() if r.steadyState() < 1e-6 else [0.0],
            list(r.model.getReactionIds()), list(r.model.getGlobalParameterIds())]


def w_sensitivities(load, name, cma):
    r = load(sbml_of(name))
    r.conservedMoietyAnalysis = cma
    t, s, p, x = r.timeSeriesSensitivities(0, 10, 11)
    out = [t, s, list(p), list(x), [r['[%s]' % i] for i in r.getIndependentFloatingSpeciesIds()], [r.getCurrentTime()]]
    r = load(sbml_of(name))  # fresh model: libroadrunner carries its sensitivity state over between calls
    r.conservedMoietyAnalysis = cma
    sp = r.model.getFloatingSpeciesIds()
    t, s, p, x = r.timeSeriesSensitivities(0, 5, 6, [], [sp[-1], sp[0]])
    return out + [t, s, list(p), list(x)]


def w_elasticities(load, name, cma):
    r = load(sbml_of(name))
    r.conservedMoietyAnalysis = cma
    r.setDiffStepSize(1e-5)
    for i, s in enumerate(r.model.getFloatingSpeciesIds()):  # off steady state
        r.setValue('[%s]' % s, 0.3 + 0.7 * i)
    return [r.getUnscaledElasticityMatrix(), r.getScaledElasticityMatrix()]


def w_ss_values(load, name, cma):
    r = load(sbml_of(name))
    r.conservedMoietyAnalysis = cma
    r.steadyStateSelections = ['[S2]', 'J1'] + list(r.model.getFloatingSpeciesIds())[:2]
    out = [r.getSteadyStateValues()]
    r.reset()
    out += [r.getSteadyStateValuesNamedArray()]
    r.reset()
    out += [r.steadyStateNamedArray(), [float(r.mcaSteadyState() < 1e-6)]]
    for e in (r.getFullEigenValuesNamedArray(), r.getReducedEigenValuesNamedArray()):
        a = np.array(sorted(map(tuple, np.round(np.asarray(e), 9))))
        out += [list(e.colnames), [len(e.rownames)], a]
    return out


def w_one_step(load):
    r = load(sbml_of('MOIETY'))
    t = 0.0
    for _ in range(5):
        t = r.oneStep(t, 0.5)
    out = [[t], r.model.getFloatingSpeciesConcentrations()]
    r.reset()
    r.integrator.stiff = False
    r.integrator.maximum_num_steps = 5000
    return out + [r.simulate(0, 5, 11)]


def w_gillespie(load):
    r = load(te.antimonyToSBML('compartment c = 2; S1 in c; S2 in c; J1: S1 -> S2; k1*S1; J2: S2 -> ; k2*S2; '
                               'S1 = 50; k1 = 0.5; k2 = 0.2'))
    r.setIntegrator('gillespie')
    r.integrator.seed = 11
    a = []
    for _ in range(400):
        r.reset()
        a.append(np.asarray(r.simulate(0, 4, 5, ['time', 'S1', '[S2]'])))
    a = np.array(a)
    return [a.mean(0), a.var(0)]


def w_rules_init(load):
    r = load(sbml_of('VARS'))
    sel = ['time', 'S1', '[S1]', 'S3', '[S3]', 'B', '[B]', 'c', 'J1', 'k2', 'init(S1)', 'init([S2])']
    out = [r.simulate(0, 2.5, 6, sel), r.simulate(2.5, 6, 8, sel)]
    return out + [[r.S1, r.c, r.B, r.k2, r.getValue('S2')]]


def w_rate_selections(load):
    return [load(sbml_of('RATE')).simulate(0, 1, 3, ['time', "S1'", "[S1]'", 'X0', '[X0]', 'init(p)', 'J1'])]


def w_rate_rule_mca(load):
    r = load(sbml_of('RATE'))
    r.setDiffStepSize(1e-5)
    t, s, p, x = r.timeSeriesSensitivities(0, 2, 5)
    out = [s, list(p), list(x), load(sbml_of('RATE')).timeSeriesSensitivities(0, 2, 5, [], [], 1)[1]]
    r.reset()
    r.steadyState()
    return out + [r.getFullJacobian(), r.getReducedJacobian(), r.getUnscaledElasticityMatrix(),
                  r.getScaledConcentrationControlCoefficientMatrix(), r.getScaledFluxControlCoefficientMatrix(),
                  r.model.getFloatingSpeciesConcentrations(), [r.p, r.S3], [r.getuCC('S1', 'k1'), r.getCC('J1', 'k2')]]


def values(r):
    m = r.model
    return [list(m.getFloatingSpeciesIds()), list(m.getBoundarySpeciesIds()), list(m.getGlobalParameterIds()),
            list(m.getCompartmentIds()), list(m.getReactionIds()), list(r.timeCourseSelections),
            list(r.steadyStateSelections), m.getFloatingSpeciesConcentrations(), m.getFloatingSpeciesAmounts(),
            m.getBoundarySpeciesConcentrations(), m.getBoundarySpeciesAmounts(), m.getGlobalParameterValues(),
            m.getCompartmentVolumes(), m.getReactionRates(), [r.getValue('time')]]


def w_workflow(load, name):
    """Ids, values, setValue, simulate, reset, resetAll on one model."""
    sbml = sbml_of(name)
    r = load(sbml)
    out = values(r)
    m = r.model
    p, s = (list(m.getGlobalParameterIds()) or ['time'])[0], m.getFloatingSpeciesIds()[0]
    r.setValue(p, r.getValue(p) * 1.5)
    r.setValue('[%s]' % s, 0.7)
    out += values(r)
    T = 0.5 if 'lorenz' in sbml else 5
    out.append(r.simulate(0, T, 11))
    sel = ['time', s, '[%s]' % s, p] + list(m.getReactionIds())[:1]
    out += [r.simulate(T, 2 * T, 6, sel)] + values(r)
    r.reset()
    out += values(r)
    r.resetAll()
    return out + values(r)


def w_attributes_mca(load, name):
    r = load(sbml_of(name))
    r.setDiffStepSize(1e-5)
    s, p = r.model.getFloatingSpeciesIds()[0], r.model.getGlobalParameterIds()[0]
    out = [[getattr(r, s), getattr(r, s + '_amt'), getattr(r, p), r['[%s]' % s]]]
    setattr(r, p, 0.3)
    r.simulate(0, 20, 3)
    r.steadyState()
    out += [r.model.getFloatingSpeciesConcentrations(), r.getScaledConcentrationControlCoefficientMatrix(),
            r.getScaledFluxControlCoefficientMatrix(), r.getFullJacobian(), r.getReducedJacobian(),
            r.getScaledElasticityMatrix(), r.getFullStoichiometryMatrix(), [r.getCC(s, p)], r.getRatesOfChange(),
            load(sbml_of(name)).timeSeriesSensitivities(0, 2, 5, [p])[1]]
    r.simulate(20, 23, 3)
    return out + [r.getLinkMatrix(), r.model.getFloatingSpeciesConcentrations(), [r.getValue(p), r.getValue('time')],
                  r.simulate(23, 26, 7)]


def w_sos(load, case):
    d = os.path.join(SOS, case)
    cfg = {}
    for line in open(os.path.join(d, '%s-settings.txt' % case)):
        k, _, v = line.partition(':')
        cfg[k.strip()] = [x.strip() for x in v.split(',') if x.strip()]
    r = load(sbml_of(case))
    sp = set(r.model.getFloatingSpeciesIds()) | set(r.model.getBoundarySpeciesIds())
    sel = ['time'] + ['[%s]' % v if v in sp and v in cfg.get('concentration', []) else v for v in cfg['variables']]
    r.integrator.absolute_tolerance = 1e-14  # tight: the reference is then within the comparison tolerance
    r.integrator.relative_tolerance = 1e-10
    r.integrator.maximum_num_steps = 10 ** 6
    t0, dur, n = float(cfg['start'][0]), float(cfg['duration'][0]), int(cfg['steps'][0])
    return [r.simulate(t0, t0 + dur, n + 1, sel)]


CMA = (False, True)
WORKFLOWS = {}
for _n in ('MOIETY', 'RULES', 'chain'):
    for _c in CMA:
        WORKFLOWS['simulate/%s/%d' % (_n, _c)] = (lambda l, n=_n, c=_c: w_simulate(l, n, c), 1e-4, 1e-7)
        WORKFLOWS['steady_state/%s/%d' % (_n, _c)] = (lambda l, n=_n, c=_c: w_steady_state(l, n, c), 1e-5, 1e-7)
for _n in ('MOIETY', 'chain8'):
    for _c in CMA:
        WORKFLOWS['sensitivities/%s/%d' % (_n, _c)] = (lambda l, n=_n, c=_c: w_sensitivities(l, n, c), 2e-4, 1e-5)
for _n in ('MOIETY', 'RULES', 'ELAST'):
    for _c in CMA:
        WORKFLOWS['elasticities/%s/%d' % (_n, _c)] = (lambda l, n=_n, c=_c: w_elasticities(l, n, c), 1e-5, 1e-7)
for _n in ('MOIETY', 'RULES'):
    for _c in CMA:
        WORKFLOWS['ss_values/%s/%d' % (_n, _c)] = (lambda l, n=_n, c=_c: w_ss_values(l, n, c), 1e-6, 1e-8)
for _n in ('MOIETY', 'RULES', 'chain'):
    WORKFLOWS['attributes_mca/' + _n] = (lambda l, n=_n: w_attributes_mca(l, n), 1e-4, 1e-7)
for _n in ('MOIETY', 'RULES', 'EVENTS_RULES', 'IA', 'VARS', 'RATE', 'chain', '00009', '00014', '00022', '00033',
           '10000', '10001', 'repressilator'):
    WORKFLOWS['workflow/' + _n] = (lambda l, n=_n: w_workflow(l, n), 2e-3, 1e-6)
for _n in ('00009', '00014', '00022', '00033', '10001', 'repressilator'):
    WORKFLOWS['sos/' + _n] = (lambda l, n=_n: w_sos(l, n), 1e-5, 1e-8)
WORKFLOWS.update({
    'events_rules': (w_events_rules, 1e-4, 1e-7), 'state_sync': (w_state_sync, 1e-4, 1e-7),
    'event_state': (w_event_state, 1e-4, 1e-7), 'rk4': (w_rk4, 1e-12, 1e-14), 'edit': (w_edit, 1e-4, 1e-7),
    'one_step': (w_one_step, 1e-4, 1e-7), 'gillespie': (w_gillespie, None, None),
    'rules_init': (w_rules_init, 1e-4, 1e-7), 'rate_selections': (w_rate_selections, 1e-4, 1e-7),
    'rate_rule_mca': (w_rate_rule_mca, 1e-3, 1e-5),
})


# ---------------------------------------------------------------- the Rust batch model (rrs_bridge) per SBML feature

def _vc(extra):
    return "compartment c=2; species S1 in c, S2 in c; J1: S1 -> S2; k*S1; k=1; S1=10; S2=0; " + extra


BASE = "J1: S1 -> S2; k*S1; k=1; S1=10; S2=0"
# name -> (Antimony, extra outputs as (roadrunner selection, roadrunner_nitro output)); the SBML-edited
# features (conversion factors, variable stoichiometry) are in make_reference.py, their SBML in turbo.npz
FEATURES = {
    'compartment rate rule': (_vc("c' = 1"), ['c']),
    'compartment assignment rule': (_vc("c := 2 + sin(time)*0.5 + k"), ['c']),
    'compartment event': (_vc("at time > 1: c = 4"), ['c']),
    'compartment and species event': (_vc("at time > 1: S1 = 7, c = 4"), ['c']),
    'compartment rule, event on its parameter': (_vc("c := 1 + kc; kc = 1; at time > 1: kc = 3"), ['c']),
    'compartment rate rule and event': (_vc("c' = 0.2; at time > 1: c = c*2"), ['c']),
    'amount species, variable compartment': (
        "compartment c=2; species S1 in c; substanceOnly species S2 in c; J1: S1 -> S2; k*S1*S2; k=0.1; "
        "S1=10; S2=1; c' = 0.3", ['c']),
    'boundary species, variable compartment': (
        "compartment c=2; species S1 in c, $B in c; J1: B -> S1; k*B; k=0.1; S1=10; B=1; c' = 0.3", ['c', 'B']),
    'rate-rule species, variable compartment': (
        "compartment c=2; species S1 in c, S2 in c; S1' = -0.5*S1; J1: -> S2; 1; S1=10; c' = 0.3", ['c']),
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
SBML_FEATURES = {  # SBML from make_reference.py (libsbml edits Antimony cannot write)
    'conversion factor (species)': [], 'conversion factor (model)': [], 'conversion factor (event-changed)': [],
    'stoichiometry: assignment rule': [], 'stoichiometry: rate rule': [], 'stoichiometry: event': [],
    'stoichiometryMath (L2)': [],
}


def feature_sbml(name, ref=None):
    return te.antimonyToSBML(FEATURES[name][0]) if name in FEATURES else str(ref['feature/%s/sbml' % name])


def f_feature(load, sbml, extra):
    r = load(sbml)
    ids = list(r.model.getFloatingSpeciesIds())
    return [list(ids), r.simulate(0, 10, 101, ['time'] + ['[%s]' % s for s in ids] + [a for a, _ in extra])]


def f_events(load):
    r = load(sbml_of('EVENTS'))
    return [r.simulate(a, b, n, ['time', '[S1]', '[S2]', 'p', 'q', 'k1', 'k2'])
            for a, b, n in [(0, 2, 21), (2, 3.5, 16), (3.5, 5, 7), (5, 6, 3), (6, 12, 61)]]


RATE_MCA = '''J0: $X -> S1; k0*X; J1: S1 -> S2; k1*S1*p; J2: S2 -> S3; k2*S2; J3: S3 ->; k3*S3*S1
    X = 2; k0=1; k1=1; k2=0.7; k3=0.5; S1=1; S2=1; S3=1; p=1; p' = 0.1*(2 + S2 - p); q = 0; q' = 0.3*(1 - q*S1)'''


def f_rate_mca(load):
    r = load(te.antimonyToSBML(RATE_MCA))
    r.setDiffStepSize(1e-5)
    r.steadyState()
    return [[r.p], r.getScaledConcentrationControlCoefficientMatrix(), r.getScaledFluxControlCoefficientMatrix(),
            r.getUnscaledConcentrationControlCoefficientMatrix()]


# ---------------------------------------------------------------- comparison

def flatten(name, out):
    """Workflow results -> {key: array} (strings as str arrays, NamedArray names as .r / .c)."""
    d = {}
    for i, v in enumerate(out):
        k = '%s/%d' % (name, i)
        if isinstance(v, list) and v and isinstance(v[0], str):
            d[k] = np.array(v, str)
            continue
        if hasattr(v, 'colnames'):
            d[k + '.c'] = np.array(list(v.colnames), str)
            d[k + '.r'] = np.array(list(getattr(v, 'rownames', [])), str)
        a = np.asarray(v)
        d[k] = a.astype(complex if np.iscomplexobj(a) else float)
    return d


class TestAgainstLibroadrunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = dict(np.load(REF))

    def check(self, name, out, rtol, atol):
        got = flatten(name, out)
        want = {k: v for k, v in self.ref.items() if k.startswith(name + '/') and not k.endswith('/sbml')}
        self.assertTrue(want, 'no reference for %s (run make_reference.py)' % name)
        self.assertEqual(sorted(want), sorted(got), name)
        for k, a in want.items():
            b = got[k]
            if a.dtype.kind == 'U':
                self.assertEqual(list(a), list(b), k)
            else:
                np.testing.assert_allclose(b.reshape(a.shape), a, rtol=rtol, atol=atol * max(1.0, np.abs(a).max(initial=0)),
                                           err_msg=k)

    def test_workflows(self):
        for name, (f, rtol, atol) in WORKFLOWS.items():
            if rtol is None:
                continue
            with self.subTest(name):
                self.check(name, f(te.loadSBMLModel), rtol, atol)

    def test_gillespie_statistics(self):
        mean, var = f_gillespie_ref = WORKFLOWS['gillespie'][0](te.loadSBMLModel)
        rm, rv = self.ref['gillespie/0'], self.ref['gillespie/1']
        se = np.sqrt(rv / 400 + var / 400) + 1e-12
        self.assertLess(np.max(np.abs(rm - mean) / se), 5)
        np.testing.assert_allclose(var, rv, rtol=0.4, atol=1e-9)

    def test_batch_model_features(self):
        """rrs_bridge.from_sbml (the Rust batch model) per SBML feature against libroadrunner."""
        from tellurium_nitro.analysis.rrs_bridge import from_sbml
        for name in list(FEATURES) + list(SBML_FEATURES):
            with self.subTest(name):
                extra = [e if isinstance(e, tuple) else (e, e) for e in FEATURES.get(name, (None, []))[1]]
                m = from_sbml(feature_sbml(name, self.ref))
                ids = list(self.ref['feature/%s/0' % name])
                got, _ = m.simulate_one([], [], np.linspace(0, 10, 101), outputs=['time'] + ids + [b for _, b in extra])
                np.testing.assert_allclose(got, self.ref['feature/%s/1' % name], rtol=1e-4, atol=1e-6, err_msg=name)
                # and the RoadRunner
                self.check('feature/' + name, f_feature(te.loadSBMLModel, feature_sbml(name, self.ref), extra), 1e-4, 1e-6)

    def test_batch_model_events_across_calls(self):
        from tellurium_nitro.analysis.rrs_bridge import from_sbml
        m = from_sbml(sbml_of('EVENTS'))
        outs, keys = ['time', 'S1', 'S2', 'p', 'q', 'k1', 'k2'], ['S1', 'S2', 'k1', 'k2', 'p', 'q']
        state, x = (), None
        for i, (a, b, n) in enumerate([(0, 2, 21), (2, 3.5, 16), (3.5, 5, 7), (5, 6, 3), (6, 12, 61)]):
            y, fin, state = m.simulate_one(keys if x else [], x or [], np.linspace(a, b, n), outputs=outs, event_state=state)
            x = [fin[k] for k in keys]
            np.testing.assert_allclose(y, self.ref['events/%d' % i], rtol=1e-4, atol=1e-6, err_msg=str((a, b)))
        self.check('events', f_events(te.loadSBMLModel), 1e-4, 1e-6)

    def test_batch_model_mca_with_rate_rules(self):
        from tellurium_nitro.analysis.rrs_bridge import from_sbml
        m = from_sbml(te.antimonyToSBML(RATE_MCA))
        res = m.mca_batch([], np.zeros((1, 0)))
        self.assertTrue(res['ok'][0])
        np.testing.assert_allclose(res['S'][0][m.state_ids.index('p')], self.ref['rate_mca/0'][0], rtol=1e-8)
        for i, key in ((1, 'ccc'), (2, 'fcc'), (3, 'uccc')):
            rows, cols = list(self.ref['rate_mca/%d.r' % i]), list(self.ref['rate_mca/%d.c' % i])
            r = [(m.reactions if key == 'fcc' else m.species).index(x) for x in rows]
            c = [m.reactions.index(x) for x in cols]
            np.testing.assert_allclose(res[key][0][np.ix_(r, c)], self.ref['rate_mca/%d' % i], rtol=1e-6, atol=1e-9, err_msg=key)
        self.check('rate_mca', f_rate_mca(te.loadSBMLModel), 1e-5, 1e-8)


# ---------------------------------------------------------------- properties needing no reference

class TestRoadRunner(unittest.TestCase):

    def test_turbo_is_always_on(self):
        self.assertTrue(te.turbo())
        self.assertTrue(te.turbo('strict'))
        with self.assertWarns(UserWarning):
            self.assertTrue(te.turbo(False))
        self.assertTrue(te.turbo())

    def test_gillespie(self):
        ant = 'compartment c = 2; S1 in c; S2 in c; J1: S1 -> S2; k1*S1; J2: S2 -> ; k2*S2; S1 = 50; k1 = 0.5; k2 = 0.2'

        def runs(n, seed, **kw):
            r = te.loada(ant)
            r.setIntegrator('gillespie')
            r.integrator.seed = seed
            out = []
            for _ in range(n):
                r.reset()
                out.append(r.simulate(**kw))
            return r, out

        r, fixed = runs(200, 7, start=0, end=4, points=5, selections=['time', 'S1', '[S2]'])
        _, again = runs(200, 7, start=0, end=4, points=5, selections=['time', 'S1', '[S2]'])
        _, var = runs(50, 3, start=0, end=4)
        self.assertEqual(list(fixed[0].colnames), ['time', 'S1', '[S2]'])
        np.testing.assert_array_equal(np.array(fixed), np.array(again))  # same seed, same runs
        self.assertGreater(len({tuple(np.asarray(a)[:, 1]) for a in fixed}), 150)  # successive runs differ
        a = np.array(fixed)
        np.testing.assert_array_equal(a[..., 1], np.round(a[..., 1]))  # amounts are whole molecules
        self.assertEqual(r.getCurrentTime(), 4.0)
        self.assertEqual(r.getValue('S1'), fixed[-1][-1, 1])  # state at the last row
        for v in var:  # variable step: a row at start, after every firing, at end
            v = np.asarray(v)
            self.assertEqual((v[0, 0], v[-1, 0]), (0.0, 4.0))
            self.assertTrue((np.diff(v[:, 0]) >= 0).all() and len(v) > 10)

    def test_gillespie_variable_step_events_rules(self):
        r = te.loada('J1: S1 -> S2; k1*S1; S3 := 2*S1; S1 = 40; k1 = 0.5\nat time > 1: S1 = S1 + 20')
        r.setIntegrator('gillespie')
        r.integrator.seed = 5
        a = np.asarray(r.simulate(0, 3, selections=['time', 'S1', 'S3', 'k1']))
        self.assertEqual((a[0, 0], a[-1, 0]), (0.0, 3.0))
        np.testing.assert_array_equal(a[:, 2], 2 * a[:, 1])  # rule column
        jump = np.flatnonzero(np.diff(a[:, 1]) > 0)  # the event adds 20 molecules once
        self.assertEqual(len(jump), 1)
        self.assertGreater(a[jump[0] + 1, 0], 1.0)

    def test_sensitivities_across_events(self):
        # S1' = -k1 S1, S1 += 3 at t = 1: exact dS1/dk1 (libroadrunner's CVODES drifts ~1% after the event)
        r = te.loada('J1: S1 -> ; k1*S1; S1 = 10; k1 = 0.5\nat time > 1: S1 = S1 + 3')
        t, s, p, x = r.timeSeriesSensitivities(0, 3, 7, ['k1'])
        S = np.where(t <= 1, 10 * np.exp(-0.5 * t), (10 * np.exp(-0.5) + 3) * np.exp(-0.5 * (t - 1)))
        ds = np.where(t <= 1, -t * S, -10 * np.exp(-0.5) * np.exp(-0.5 * (t - 1)) - (t - 1) * S)
        np.testing.assert_allclose(s[:, 0, 0], ds, rtol=1e-4, atol=1e-6)

    def test_last_run(self):
        r = te.loada(chain(12))
        self.assertIsNone(r.lastRun)
        r.simulate(0, 10, 11)
        self.assertEqual((r.lastRun['method'], r.lastRun['engine']), ('simulate', 'rust'))
        r.steadyState()
        self.assertEqual((r.lastRun['method'], r.lastRun['engine']), ('steadyState', 'rust'))

    def test_simulate_batch(self):
        ant = chain(12) + "\nat time > 3: k1 = 2"
        r = te.loada(ant)
        ks, s3 = [0.5, 1.0, 2.0], [0.1, 1.0, 3.0]
        got = r.simulateBatch({'k0': ks, 'S3': s3}, 0, 10, 21)
        self.assertEqual(got.shape, (3, 21, len(r.timeCourseSelections)))
        for i in range(3):  # same as one simulate per set from this state
            q = te.loada(ant)
            q.k0, q['[S3]'] = ks[i], s3[i]
            np.testing.assert_allclose(got[i], np.asarray(q.simulate(0, 10, 21)), rtol=1e-5, atol=1e-9)
        self.assertEqual(r.getCurrentTime(), 0.0)  # the model itself is untouched
        # amounts, rates of change: the per-set path
        sel = ['time', 'S1', "S2'", 'J1']
        got = r.simulateBatch({'k0': ks}, 0, 5, 6, sel)
        for i in range(3):
            q = te.loada(ant)
            q.k0 = ks[i]
            np.testing.assert_allclose(got[i], np.asarray(q.simulate(0, 5, 6, sel)), rtol=1e-5, atol=1e-9)

    def test_current_sbml_round_trip(self):
        r = te.loada(IA)
        r.simulate(0, 3, 4)
        r2 = te.loadSBMLModel(r.getCurrentSBML())
        for a, b in zip(values(r)[:-1], values(r2)):  # not the time
            if not (isinstance(a, list) and a and isinstance(a[0], str)):
                np.testing.assert_allclose(a, b, rtol=1e-12)
        r = te.loada(MOIETY)
        r.addReaction('J9', ['S1'], ['S3'], 'k1*S1')
        self.assertIn('J9', r.model.getReactionIds())

    def test_dropin_module(self):
        import tellurium_nitro.nitro as rr
        r = rr.RoadRunner(te.antimonyToSBML(chain(12)))
        self.assertIsInstance(r, ExtendedRoadRunner)
        self.assertIsInstance(r.simulate(0, 10, 11), rr.NamedArray)
        rr.Config.setValue('x', 1)
        self.assertEqual(rr.Config.getValue('x'), 1)
        rr.Logger.setLevel(rr.Logger.LOG_WARNING)

    def test_named_array_and_pickle(self):
        import copy
        import pickle
        r = te.loada('S1 -> S2; k1*S1; k1 = 0.5; S1 = 10')
        a = r.simulate(0, 3, 4)
        self.assertEqual(str(a).splitlines()[0].split(','), ['    time', '    [S1]', '    [S2]'])
        np.testing.assert_array_equal(a['[S1]'], np.asarray(a)[:, 1])
        b = pickle.loads(pickle.dumps(a))
        self.assertEqual(b.colnames, a.colnames)
        r.integrator.relative_tolerance = 1e-8
        r.k1 = 0.7
        c = copy.deepcopy(r)
        self.assertEqual((c.k1, c.getCurrentTime(), c.integrator.relative_tolerance), (0.7, 3.0, 1e-8))
        np.testing.assert_allclose(c.simulate(3, 5, 3), r.simulate(3, 5, 3), rtol=1e-12)

    def test_model_editing(self):
        r = te.loada('compartment c = 1; S1 in c; J1: S1 -> S2; k1*S1; k1 = 0.5; S1 = 10; S2 = 0; p = k1*2')
        r.addSpeciesConcentration('S3', 'c', 2.0)
        r.addParameter('k2', 0.25)
        r.addReaction('J2', ['S2'], ['S3'], 'k2*S2')
        r.setKineticLaw('J1', 'k1*S1/2')
        self.assertEqual(r.getKineticLaw('J1').replace(' ', ''), 'k1*S1/2')
        r.removeInitialAssignment('p')
        self.assertEqual(r.getInitialAssignmentIds(), [])
        self.assertEqual(r.p, 1.0)
        r.setInitValue('k1', 0.1)
        self.assertEqual(r.getValue('init(k1)'), 0.1)
        r.setModelName('edited')
        self.assertEqual(r.getModelName(), 'edited')
        r.resetAll()
        y = np.asarray(r.simulate(0, 4, 5, ['time', '[S1]', '[S3]']))
        np.testing.assert_allclose(y[:, 1], 10 * np.exp(-0.05 * y[:, 0]), rtol=1e-5)
        self.assertGreater(y[-1, 2], 2.0)

    def test_value_vectors_follow_their_ids(self):
        # libroadrunner (cma off) returns getIndependentFloatingSpeciesConcentrationsV in model order, not in
        # getIndependentFloatingSpeciesIds order; here values and ids always correspond
        r = te.loada(MOIETY)
        for cma in (False, True):
            r.conservedMoietyAnalysis = cma
            ind = r.getIndependentFloatingSpeciesIds()
            np.testing.assert_array_equal(r.getIndependentFloatingSpeciesConcentrationsV(), [r['[%s]' % i] for i in ind])
            dep = r.getDependentFloatingSpeciesIds()
            np.testing.assert_array_equal(r.getDependentFloatingSpeciesAmountsV(), [r[i] for i in dep])

    def test_sensitivity_solver_tolerances(self):
        r = te.loada('J1: S1 -> ; k1*S1; S1 = 10; k1 = 0.5')
        r.getSensitivitySolver().relative_tolerance = 1e-10
        r.getSensitivitySolver().absolute_tolerance = 1e-14
        t, s, p, x = r.timeSeriesSensitivities(0, 2, 5, ['k1'])
        np.testing.assert_allclose(s[:, 0, 0], -t * 10 * np.exp(-0.5 * t), rtol=1e-8, atol=1e-12)

    def test_model_counts_and_ids(self):
        r = te.loada(EVENTS_RULES)
        m = r.model
        self.assertEqual([r.getNumFloatingSpecies(), r.getNumBoundarySpecies(), r.getNumGlobalParameters(),
                          r.getNumCompartments(), r.getNumReactions(), r.getNumEvents(), r.getNumRateRules()],
                         [3, 0, 4, 1, 2, 2, 1])
        self.assertEqual([m.getFloatingSpeciesId(1), m.getReactionId(0), m.getCompartmentId(0)], ['S2', 'J1', 'c'])
        self.assertEqual(m.getAssignmentRuleIds(), ['x'])
        self.assertEqual(m.getStoichiometry(2, 1), 1.0)
        self.assertTrue(r.getHasOnlySubstanceUnits('S3'))
        self.assertAlmostEqual(m.getFloatingSpeciesAmountRate(1), r.getValue("S2'"))

    def test_test_models(self):
        self.assertIn('feedback.xml', te.listTestModels())
        for m in te.listTestModels():
            self.assertTrue(te.getTestModel(m).lstrip().startswith('<'))
        r = te.loadTestModel('feedback.xml')
        self.assertEqual(np.asarray(r.simulate(0, 10, 11)).shape[0], 11)


if __name__ == '__main__':
    unittest.main()

"""Fixed-step integrators against cached libroadrunner courses and tight CVODE."""
import glob
import os

import antimony
import numpy as np
import pytest
import roadrunner_nitro as rn

HERE = os.path.dirname(os.path.abspath(__file__))
SOS = sorted(glob.glob(os.path.join(HERE, '..', '..', 'benchmarks', 'sosbench', '*', '*.xml')))

ANT = {
    'rate_rule_time': '''
        compartment C = 2.5
        species S1 in C, S2 in C
        J1: S1 -> S2; k1*S1*(1 + sin(time))
        J2: S2 -> S1; k2*S2
        p' = -0.3*p + S1
        k1 = 0.4; k2 = 0.15; S1 = 3; S2 = 1; p = 0.5
    ''',
    'mass_action_volumes': '''
        compartment A = 0.2, B = 4
        species X in A, Y in B, Z in B
        J0: $X0 -> X; v0
        J1: X -> Y; k1*X
        J2: Y + Y -> Z; k2*Y^2
        J3: Z -> ; k3*Z
        X0 = 1; v0 = 0.5; k1 = 1.2; k2 = 0.05; k3 = 0.3; X = 0.1; Y = 0; Z = 0
    ''',
    'events': '''
        S1 -> S2; k*S1
        k = 0.5; S1 = 10; S2 = 0
        E1: at time > 2.05: S1 = 10
        E0: at time >= 0, t0 = false: S2 = S2 + 100
        E2: at 0.3 after S1 < 5: k = k*1.1
        E3: at time > 2.65, priority = 10, fromTrigger = false: S2 = 2*S2
        E4: at time > 2.65, priority = 1, fromTrigger = false: S2 = S2 + 3
        E5: at 0.4 after time > 3.1 && time < 3.3, persistent = false: S1 = 99
    ''',
}

CASES = [  # (label, roadrunner integrator, settings, nitro method, steps per output interval)
    ('euler', 'euler', {}, 'euler', 1),
    ('euler/4', 'euler', {'subdivision_steps': 4}, 'euler', 4),
    ('rk4', 'rk4', {}, 'rk4', 1),
    ('rk45', 'rk45', {'variable_step_size': False}, 'rk45', None),
    ('rk45 eps1e-8 hmax0.25', 'rk45', {'variable_step_size': False, 'epsilon': 1e-8, 'maximum_time_step': 0.25}, 'rk45', None),
]


def sbml_of(src):
    antimony.clearPreviousLoads()
    assert antimony.loadAntimonyString(src) >= 0, antimony.getLastError()
    return antimony.getSBMLString(antimony.getMainModuleName())


def models():
    for f in SOS:
        yield os.path.basename(f).split("-")[0], open(f).read(), 10.0, 1001
    for k in ('rate_rule_time', 'mass_action_volumes'):
        yield k, sbml_of(ANT[k]), 4.0, 41


# rk45's epsilon is absolute (on amounts). 00009 (compartment 4e-12, amounts ~1e-17): there it is
# meaningless (libroadrunner's rk45 moves by 63% when one initial value changes by 1e-15); 00033
# (values ~1e5): rounding never meets it, libroadrunner's steps fall to minimum_time_step and grind
# for hours. Nitro clamps the bound to the state's scale (fixed.rs rk45_tol); those two are checked
# against a tight CVODE run instead.
TIGHT = {'00009', '00033'}


REF = os.path.join(HERE, 'reference', 'fixed.npz')
VARIABLE = [(1e-12, 1.0, 100000), (1e-8, 0.25, 100000), (1e-10, 1.0, 7)]


def fixed_grid(name, method, end, points):
    if method != 'rk45':
        # Resolve the fast transient within the stability range of the explicit method.
        end = {'00009': 0.001, '10001': 0.001, '00033': 1.0}.get(name, end)
    return end, points


def error(got, want):
    return (np.abs(got - want) / (np.abs(want).max(axis=0) + 1e-300)).max()


def test_fixed_vs_roadrunner():
    with np.load(REF) as ref:
        for name, sbml, end, n in models():
            m, _ = rn.from_sbml(sbml)
            ids = ref[f'{name}/ids'].tolist()
            for label, _, settings, method, sub in CASES:
                key = f'{name}/{label}'
                y = ref[key]
                times = np.ascontiguousarray(y[:, 0])
                if method == 'rk45':
                    step, atol = settings.get('maximum_time_step', 1.0), settings.get('epsilon', 1e-12)
                else:
                    step, atol = (times[1] - times[0]) / sub, 1e-12
                got = rn.simulate_fixed(m, method, times, step, 1e-6, atol, ids)
                if method == 'rk45' and name in TIGHT:
                    assert error(got, y[:, 1:]) <= 1e4 * atol, key
                elif method == 'rk45' and error(got, y[:, 1:]) > 1e-9:
                    # The scale-aware epsilon can choose a different step sequence.
                    exact = ref[f'{name}/cvode'][:, 1:]
                    assert error(got, exact) <= error(y[:, 1:], exact), key
                else:
                    assert error(got, y[:, 1:]) <= 1e-9, key


def tight(sbml, sel, times):
    r = rn.RoadRunner(sbml)
    r.setIntegratorSetting('relative_tolerance', 1e-12)
    r.setIntegratorSetting('absolute_tolerance', 1e-14)
    r.timeCourseSelections = sel
    return np.asarray(r.simulateTimes(list(times))).reshape(len(times), len(sel))


@pytest.mark.parametrize('name', ['rate_rule_time', 'mass_action_volumes'])
@pytest.mark.parametrize('index', range(len(VARIABLE)))
def test_rk45_variable_step(name, index):
    sbml = sbml_of(ANT[name])
    eps, hmax, rows = VARIABLE[index]
    key = f'{name}/variable{index}'
    with np.load(REF) as ref:
        sel = ['time'] + [f'[{i}]' for i in ref[f'{name}/ids']]
        ya = ref[key]
        exact = ref[key + '/cvode'][:, 1:]
        scale = np.abs(exact).max(axis=0) + 1e-300
        ea = (np.abs(ya[:, 1:] - exact) / scale).max()
    b = rn.RoadRunner(sbml)
    b.setIntegrator('rk45')
    for k, v in dict(epsilon=eps, maximum_time_step=hmax, max_output_rows=rows).items():
        b.setIntegratorSetting(k, v)
    b.timeCourseSelections = sel
    yb = np.asarray(b.simulate(0, 10)).reshape(-1, len(sel))
    assert yb[0, 0] == 0.0
    assert np.all((0 < np.diff(yb[:, 0])) & (np.diff(yb[:, 0]) <= hmax + 1e-14))
    if len(ya) == rows:
        assert len(yb) == rows and yb[-1, 0] < 10
    else:
        assert yb[-1, 0] == 10
    if name == 'mass_action_volumes':
        # The initial amount norm is below one, so Nitro's epsilon clamp tightens the bound.
        assert yb[1, 0] < ya[1, 0] and len(yb) >= len(ya)
    else:
        assert len(yb) == len(ya)
        assert yb[1, 0] == pytest.approx(ya[1, 0], rel=1e-9)
    eb = (np.abs(yb[:, 1:] - tight(sbml, sel, yb[:, 0])[:, 1:]) / scale).max()
    assert eb <= max(2 * ea, 1e-9)


@pytest.mark.parametrize('method,step,tol', [('euler', 0.0005, 2e-4), ('rk4', 0.01, 1e-8), ('rk45', 0.08, 1e-8)])
def test_fixed_events(method, step, tol):
    sbml = sbml_of(ANT['events'])
    m, _ = rn.from_sbml(sbml)
    times = np.linspace(0, 4, 41)
    ids = ['S1', 'S2', 'k']
    got = rn.simulate_fixed(m, method, times, step, 1e-6, 1e-12, ids)
    exact = tight(sbml, ['[S1]', '[S2]', 'k'], times)
    assert error(got, exact) <= tol


def test_rk45_variable_events_across_one_step_calls():
    sbml = sbml_of(ANT['events'])
    r = rn.RoadRunner(sbml)
    r.setIntegrator('rk45')
    r.timeCourseSelections = ['time', '[S1]', '[S2]', 'k']
    y = np.asarray(r.simulate(0, 4)).reshape(-1, 4)
    # Event rows are post-event. Sample the final state on the same grid as CVODE.
    sel = r.timeCourseSelections
    expected = tight(sbml, sel, np.array([0.0, 4.0]))[-1]
    np.testing.assert_allclose(y[-1], expected, rtol=1e-8, atol=1e-9)
    r.resetAll()
    t = 0.0
    for _ in range(4 * len(y)):
        if t == 4.0:
            break
        tn = r.oneStep(t, 4.0 - t)
        assert t < tn <= 4.0
        t = tn
    assert t == 4.0
    np.testing.assert_allclose([r.getValue(s) for s in sel], expected, rtol=1e-8, atol=1e-9)

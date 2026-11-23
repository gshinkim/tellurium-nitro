"""roadrunner_nitro.simulate_fixed vs libroadrunner's own euler / rk4 / rk45 (same method, same steps).
Run: <venv>/bin/python rrs/tests/test_fixed_vs_roadrunner.py   (or pytest)
"""
import glob
import os

import antimony
import numpy as np
import roadrunner
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
    # events: libroadrunner fires them at the start time for all three; then euler checks triggers at
    # every step end and rk4 / rk45 ignore them
    'events': '''
        S1 -> S2; k*S1
        k = 0.5; S1 = 10; S2 = 0
        E1: at time > 2.05: S1 = 10
        E0: at time >= 0, t0 = false: S2 = S2 + 100
        E2: at 0.3 after S1 < 5: k = k*1.1
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
    for k, src in ANT.items():
        yield k, sbml_of(src), 4.0, 41


def compare(sbml, t_end, n, integ, settings, method, sub):
    r = roadrunner.RoadRunner(sbml)
    r.setIntegrator(integ)
    for k, v in settings.items():
        r.integrator.setValue(k, v)
    ids = list(r.model.getFloatingSpeciesIds())
    r.timeCourseSelections = ['time'] + [f'[{i}]' for i in ids]  # nitro states are concentrations
    ref = r.simulate(0, t_end, n)
    times = np.ascontiguousarray(ref[:, 0])
    ref = np.asarray(ref[:, 1:])
    dt = times[1] - times[0]
    if method == 'rk45':
        step, atol = settings.get('maximum_time_step', 1.0), settings.get('epsilon', 1e-12)
    else:
        step, atol = dt / sub, 1e-12
    m, _ = rn.from_sbml(sbml)
    got = rn.simulate_fixed(m, method, times, step, 1e-6, atol, ids)
    # explicit methods blow up on stiff models (both engines); compare the rows before that
    y0 = np.abs(ref[0]).max()
    bad = ~np.isfinite(ref).all(1) | (np.abs(ref).max(1) > 1e3 * max(y0, 1e-300))
    rows = int(np.argmax(bad)) if bad.any() else len(ref)
    ref, got = ref[:rows], got[:rows]
    scale = np.abs(ref).max(0)  # per species
    d = np.abs(got - ref)
    rel = d / np.maximum(np.abs(ref), 1e-6 * scale + 1e-300)
    return d.max(), rel.max(), (d / (scale + 1e-300)).max(), rows


def test_fixed_vs_roadrunner(verbose=False):
    worst = []
    for name, sbml, t_end, n in models():
        for label, integ, settings, method, sub in CASES:
            # rk45's epsilon is absolute (on amounts). 00009 (compartment 4e-12, amounts ~1e-17): it is
            # meaningless and libroadrunner's own rk45 result moves by 63% when one initial value changes
            # by 1e-15. 00033 (values ~1e5): steps fall to minimum_time_step, both engines grind for hours.
            if name in ('00009', '00033') and method == 'rk45':
                continue
            da, dr, dn, rows = compare(sbml, t_end, n, integ, settings, method, sub)
            if verbose:
                print(f'{name:20s} {label:22s} rows {rows:5d}  max abs {da:.1e}  max rel {dr:.1e}  max abs/colmax {dn:.1e}', flush=True)
            worst.append((dn, name, label))
            if not verbose: assert dn <= 1e-9, (name, label, dn)
    return max(worst)


if __name__ == '__main__':
    print('roadrunner_nitro:', rn.__file__, ' roadrunner:', roadrunner.__version__)
    print('worst (abs diff / max |y|):', test_fixed_vs_roadrunner(verbose=True))

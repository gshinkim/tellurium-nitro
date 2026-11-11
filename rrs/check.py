"""roadrunner_rs vs roadrunner: MCA to 1e-8, simulation (rules, events) to 1e-4, SSA statistics,
plus timings.
Run: ~/venvs/tellurium/bin/python rrs/check.py
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))  # this repo's tellurium
import roadrunner
import tellurium as te
from tellurium.analysis.rrs_bridge import from_roadrunner


def load(ant):
    # plain libroadrunner as the reference: tellurium replaces roadrunner.RoadRunner with a facade
    # that may itself run on roadrunner_rs
    return roadrunner.roadrunner.RoadRunner(te.antimonyToSBML(ant))

SMALL = '''
J0: $X0 -> S1; k0*X0
J1: S1 -> S2; k1*S1 - k1r*S2
J2: S2 + S3 -> S4; k2*S2*S3
J3: S4 -> S3 + S2; k3*S4
J4: S2 -> $X1; k4*S2
X0=10; k0=1; k1=2; k1r=0.5; k2=0.3; k3=1.1; k4=0.7; S1=1; S2=1; S3=2; S4=1
'''


def model(n):
    lines = ["J0: $X0 -> S1; k0*X0"]
    for i in range(1, n):
        lines.append(f"J{i}: S{i} -> S{i+1}; k{i}*S{i} - k{i}r*S{i+1}")
    lines.append(f"J{n}: S{n} -> $X1; k{n}*S{n}")
    lines += [f"B{i}: S{i} -> $W; kb{i}*S{i}/(1+S{i})" for i in range(5, n, 5)]
    lines.append("X0=10; k0=1")
    lines += [f"k{i}={1+0.1*(i%7)}; k{i}r=0.2; S{i}=1" for i in range(1, n + 1)]
    lines += [f"kb{i}=0.3" for i in range(5, n, 5)]
    return "\n".join(lines)


def check(name, ant, k1s):
    m = from_roadrunner(load(ant))
    res = m.mca_batch(['k1'], np.array(k1s)[:, None])
    assert res['ok'].all(), name
    times = np.linspace(0, 10, 51)
    sim = m.simulate_batch(['k1'], np.array(k1s)[:, None], times)
    for i, k in enumerate(k1s):
        r = load(ant)
        r.k1 = k
        ref = r.simulate(0, 10, 51)
        cols = [c.strip('[]') for c in ref.colnames[1:]]
        idx = [m.species.index(c) for c in cols]
        np.testing.assert_allclose(sim[i][:, idx], np.array(ref)[:, 1:], rtol=1e-4, atol=1e-8)
        r.resetAll()
        r.k1 = k
        # roadrunner's elasticities are finite differences (default step 0.02 gives ~5e-8 error in
        # fcc on the bench model); roadrunner_rs differentiates exactly, so tighten the reference.
        r.setDiffStepSize(1e-4)
        r.steadyState()
        for key, get in (('ccc', r.getScaledConcentrationControlCoefficientMatrix),
                         ('fcc', r.getScaledFluxControlCoefficientMatrix)):
            M = get()
            rows = [(m.species if key == 'ccc' else m.reactions).index(x) for x in M.rownames]
            colsr = [m.reactions.index(x) for x in M.colnames]
            np.testing.assert_allclose(res[key][i][np.ix_(rows, colsr)], np.array(M), rtol=1e-8, atol=1e-8)
    print(f"PASS {name}: mca (ccc, fcc) rtol=atol=1e-8 and simulate rtol=1e-4 for {len(k1s)} samples")


def timing(ant, n=1000):
    r = load(ant)
    nref = 20
    t0 = time.perf_counter()
    for k in np.linspace(0.5, 3, nref):
        r.resetAll()
        r.k1 = k
        # roadrunner's elasticities are finite differences (default step 0.02 gives ~5e-8 error in
        # fcc on the bench model); roadrunner_rs differentiates exactly, so tighten the reference.
        r.setDiffStepSize(1e-4)
        r.steadyState()
        r.getUnscaledElasticityMatrix()
    t_rr = (time.perf_counter() - t0) / nref
    t0 = time.perf_counter()
    m = from_roadrunner(r)
    t_build = time.perf_counter() - t0
    vals = np.linspace(0.5, 3, n)[:, None]
    out = {}
    for label, fn, th in (('steady_state_batch', m.steady_state_batch, 1), ('steady_state_batch', m.steady_state_batch, None),
                          ('mca_batch', m.mca_batch, 1), ('mca_batch', m.mca_batch, None)):
        fn(['k1'], vals[:10], threads=th)  # warm-up
        t0 = time.perf_counter()
        res = fn(['k1'], vals, threads=th)
        out[(label, th)] = (time.perf_counter() - t0) / n
        ok = res[3] if isinstance(res, tuple) else res['ok']
        assert ok.all()
    print(f"\nTiming, 40-species model, {n} samples (from_roadrunner build: {t_build*1e3:.1f} ms)")
    print(f"{'engine':42s} {'per sample':>12s} {'speedup':>9s}")
    print(f"{'roadrunner steadyState+elasticities (1 thr)':42s} {t_rr*1e6:10.1f} us {1:8.0f}x")
    for (label, th), t in out.items():
        lab = f"roadrunner_rs {label} ({'1 thr' if th else 'all thr'})"
        print(f"{lab:42s} {t*1e6:10.2f} us {t_rr/t:8.0f}x")
    t_ss1 = out[('steady_state_batch', 1)]
    assert t_rr / t_ss1 >= 100, "less than 100x faster single-threaded"


RULES = '''
compartment C = 2
species S1 in C, S2 in C, S3 in C
substanceOnly species A in C
J0: S1 -> S2; C*k1*S1*(1 + 0.1*p)
J1: S2 -> ; C*k2*S2
S1 = 10; S2 = 0; S3 = 1; A = 3; k1 = 0.5; k2 = 0.2
S3' = -0.1*S3 + 0.01*S2
A' = -0.3*A
p' = 0.2 - 0.05*p
p = 1
S4 := 2*S1 + A
q := k1*p
'''

EVENTS = {
    'threshold': '''
J0: S1 -> S2; k1*S1
S1 = 10; S2 = 0; k1 = 1; n = 0
E1: at (S1 < 4): S1 = 10, n = n + 1
''',
    'delayed (fromTrigger true/false)': '''
J0: S1 -> S2; k1*S1
S1 = 10; S2 = 0; k1 = 1; a = 0; b = 0
E1: at 0.7 after (S1 < 5): a = S1
E2: at 0.7 after (S1 < 5), fromTrigger=false: b = S1, S1 = S1 + 8
''',
    'persistent=false cancel': '''
J0: -> S1; 1
J1: S1 -> ; k*S1
S1 = 0; k = 1; a = 0; b = 0
E0: at (time > 1.5): k = 20
E1: at 1 after (S1 > 0.7), persistent=false: a = 1
E2: at 1 after (S1 > 0.7), persistent=true: b = 1
''',
    'simultaneous priority': '''
J0: S1 -> S2; k1*S1
S1 = 10; S2 = 0; k1 = 0.3
E1: at (time > 2), priority=1, fromTrigger=false: S1 = S1 + 1
E2: at (time > 2), priority=5, fromTrigger=false: S1 = 2*S1
E3: at (time > 4), priority=10, fromTrigger=false: S2 = S2 + 1
E4: at (time > 4), priority=-3, fromTrigger=false: S2 = 3*S2
''',
    'cascade': '''
J0: S1 -> S2; k1*S1
S1 = 10; S2 = 0; k1 = 0.3; p = 0
E1: at (time > 1.3): p = 1
E2: at (p > 0.5): S1 = 20, k1 = 0.6
E3: at (S1 > 15): S2 = S2 + 100
''',
    't0 initialValue=false': '''
J0: S1 -> S2; k1*S1
S1 = 10; S2 = 0; k1 = 0.3
E1: at (S1 > 1), t0=false: S2 = 5
E2: at (S1 > 1), t0=true: S1 = 50
''',
}


def compare(name, ant, T=10, N=101):
    """roadrunner_rs simulate_one/simulate_batch vs roadrunner simulate, all outputs, tight tolerances."""
    r = load(ant)
    m = from_roadrunner(r)
    ids = [i for i in m.output_ids if i in r.model.getFloatingSpeciesIds() + r.model.getGlobalParameterIds()
           + r.model.getReactionIds() + r.model.getBoundarySpeciesIds() + r.model.getCompartmentIds()]
    sel = ['time'] + ['[%s]' % i if i in r.model.getFloatingSpeciesIds() + r.model.getBoundarySpeciesIds() else i for i in ids]
    r.integrator.relative_tolerance = 1e-10
    r.integrator.absolute_tolerance = 1e-12
    ref = np.array(r.simulate(0, T, N, sel))
    times = np.linspace(0, T, N)
    traj, fin = m.simulate_one([], [], times, ['time'] + ids, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(traj, ref, rtol=1e-4, atol=1e-8, err_msg=name)
    batch = m.simulate_batch([], np.zeros((3, 0)), times, rtol=1e-10, outputs=['time'] + ids)
    np.testing.assert_allclose(batch, np.broadcast_to(traj, batch.shape), rtol=1e-12, atol=1e-14)
    for k, v in fin.items():  # final values incl. event-changed parameters
        if k in ids:
            np.testing.assert_allclose(v, ref[-1, 1 + ids.index(k)], rtol=1e-4, atol=1e-8, err_msg=name)
    print(f"PASS {name}: {len(ids)} outputs x {N} points match roadrunner (rtol=1e-4, atol=1e-8)")
    return r, m


def check_rules():
    r, m = compare('rate + assignment rules (C=2, substanceOnly)', RULES)
    assert m.state_ids == ['S1', 'S2', 'S3', 'A', 'p'], m.state_ids
    # Jacobian at the initial state vs roadrunner full Jacobian (floating species incl. rate-rule species)
    r.reset()
    J = r.getFullJacobian()
    ours = m.jacobian([], [])
    idx = [m.state_ids.index(c) for c in J.colnames]
    # roadrunner Jacobian is in amounts for substanceOnly A; identical here since A has no coupling
    np.testing.assert_allclose(ours[np.ix_(idx, idx)], np.array(J), rtol=1e-6, atol=1e-10)
    # steady state solves rate-rule states too (p* = 4, A* = 0)
    S, v, E, ok = m.steady_state_batch([], np.zeros((1, 0)))
    assert ok[0] and abs(S[0, 4] - 4) < 1e-9 and abs(S[0, 3]) < 1e-9, S
    print("PASS jacobian vs getFullJacobian, steady state with rate-rule states")


def check_ssa():
    ant = '''
J0: -> S1; k0
J1: S1 -> ; k1*S1
J2: S1 -> S2; k2*S1
S1 = 5; S2 = 0; k0 = 10; k1 = 0.1; k2 = 0.05
E1: at (S2 > 3): k0 = 30
'''
    m = from_roadrunner(load(ant))
    times = np.linspace(0, 10, 11)
    n_rs, n_rr = 4000, 1000
    ours = m.gillespie_batch([], np.zeros((n_rs, 0)), times, seed=1)
    a = ours.copy()
    again = m.gillespie_batch([], np.zeros((n_rs, 0)), times, seed=1, threads=3)
    assert np.array_equal(a, again), "SSA not reproducible across thread counts"
    r = load(ant)
    r.integrator = 'gillespie'
    ref = []
    for _ in range(n_rr):
        r.resetAll()  # the event changes k0; reset() would keep it
        ref.append(np.array(r.simulate(0, 10, 11, ['[S1]', '[S2]'])))
    ref = np.array(ref)
    for k in range(2):
        mu1, mu2 = a[:, :, k].mean(0), ref[:, :, k].mean(0)
        se = np.sqrt(a[:, :, k].var(0) / n_rs + ref[:, :, k].var(0) / n_rr) + 1e-12
        z = np.abs(mu1 - mu2) / se
        assert z[1:].max() < 4.5, (k, z)
        v1, v2 = a[:, -1, k].var(), ref[:, -1, k].var()
        assert abs(v1 / v2 - 1) < 0.2, (k, v1, v2)
    print(f"PASS gillespie (event S2 > 3): means at 10 times within 4.5 SE of roadrunner ({n_rs} vs {n_rr} runs), "
          f"final variances within 20%; max z = {z[1:].max():.2f}")


def sim_timing():
    rows = []
    times = np.linspace(0, 10, 51)
    for label, ant in (('40-species bench', model(40)), ('event model (threshold)', EVENTS['threshold'])):
        r = load(ant)
        m = from_roadrunner(r)
        nref = 50
        t0 = time.perf_counter()
        for _ in range(nref):
            r.resetAll()
            r.simulate(0, 10, 51)
        t_rr = (time.perf_counter() - t0) / nref
        for n, th in ((1, 1), (1000, 1), (1000, None)):
            vals = np.linspace(0.5, 3, n)[:, None]
            m.simulate_batch(['k1'], vals[:4], times, threads=th)
            t0 = time.perf_counter()
            reps = 200 if n == 1 else 1
            for _ in range(reps):
                out = m.simulate_batch(['k1'], vals, times, threads=th)
            t = (time.perf_counter() - t0) / reps
            assert np.isfinite(out).all()
            rows.append((label, n, th, t_rr, t / n))
    print(f"\nTiming simulate(0, 10, 51), default tolerances (roadrunner CVODE vs roadrunner_rs Rodas4)")
    print(f"{'model':26s} {'samples':>7s} {'threads':>7s} {'roadrunner/sample':>18s} {'rs/sample':>11s} {'speedup':>8s}")
    for label, n, th, t_rr, t in rows:
        print(f"{label:26s} {n:7d} {th or 'all'!s:>7s} {t_rr*1e6:15.1f} us {t*1e6:8.2f} us {t_rr/t:7.1f}x")


if __name__ == '__main__':
    check('test_gpu model', SMALL, [0.5, 1.0, 2.0, 4.0])
    check('bench model n=40', model(40), list(np.linspace(0.5, 3, 50)))
    check_rules()
    for name, ant in EVENTS.items():
        compare('event: ' + name, ant)
    check_ssa()
    timing(model(40))
    sim_timing()

"""One benchmark task in the current Python environment; prints one JSON line.

    python benchmarks/bench.py <task> <out.npz>

Tasks: sos:<case> (roadrunner's sosbench, see sosbench/), ss, fcc, sim, sens, cc200 (40-species model).
Reports median wall time per call, CPU time / wall time (cores busy) and peak RSS (this process + children).
Run all three engines with benchmarks/run.py.
"""
import json
import os
import resource
import sys
import time

import numpy as np
import tellurium as te

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.environ.get('BENCH_REF') == '1'  # high-precision reference run (plain roadrunner, tight tolerances)


def tight(r):
    r.integrator.relative_tolerance, r.integrator.absolute_tolerance = 1e-12, 1e-14
    r.integrator.maximum_num_steps = 10 ** 7
    r.setDiffStepSize(1e-5)  # finite-difference elasticities: smaller step, ~1e-10 error
    r.steadyStateSolver.relative_tolerance = 1e-14
    ss = r.getSensitivitySolver()
    ss.relative_tolerance, ss.absolute_tolerance, ss.maximum_num_steps = 1e-10, 1e-12, 10 ** 7  # tighter takes >15 min


def chain(n=40):  # the Tellurium Turbo benchmark model (same as rrs/check.py)
    lines = ["J0: $X0 -> S1; k0*X0"]
    for i in range(1, n):
        lines.append(f"J{i}: S{i} -> S{i+1}; k{i}*S{i} - k{i}r*S{i+1}")
    lines.append(f"J{n}: S{n} -> $X1; k{n}*S{n}")
    lines += [f"B{i}: S{i} -> $W; kb{i}*S{i}/(1+S{i})" for i in range(5, n, 5)]
    lines.append("X0=10; k0=1")
    lines += [f"k{i}={1+0.1*(i%7)}; k{i}r=0.2; S{i}=1" for i in range(1, n + 1)]
    lines += [f"kb{i}=0.3" for i in range(5, n, 5)]
    return "\n".join(lines)


def cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def timed(f, budget=2.0, minrep=3, maxrep=200):
    """Median seconds per call after one warm-up call; CPU/wall over the timed calls."""
    out = f()  # first call's output is the one checked: roadrunner's sensitivities leak state into later calls
    ts, c0, w0 = [], cpu(), time.perf_counter()
    while len(ts) < minrep or (sum(ts) < budget and len(ts) < maxrep):
        t = time.perf_counter()
        f()
        ts.append(time.perf_counter() - t)
        if ts[-1] > 10:
            break
    return float(np.median(ts)), (cpu() - c0) / (time.perf_counter() - w0), len(ts), out


def sos(case):
    d = os.path.join(HERE, 'sosbench', case)
    xml = [f for f in os.listdir(d) if f.endswith('.xml')][0]
    cfg = {}
    for line in open(os.path.join(d, f'{case}-settings.txt')):
        k, _, v = line.partition(':')
        cfg[k.strip()] = [s.strip() for s in v.split(',') if s.strip()]
    t = time.perf_counter()
    r = te.loadSBMLModel(os.path.join(d, xml))
    load = time.perf_counter() - t
    sp = set(r.getFloatingSpeciesIds()) | set(r.getBoundarySpeciesIds())
    conc = set(cfg.get('concentration', []))
    sel = ['time'] + [f'[{v}]' if v in sp and v in conc else v for v in cfg['variables']]
    r.integrator.absolute_tolerance = float(cfg['absolute'][0])
    r.integrator.relative_tolerance = float(cfg['relative'][0])
    if REF:
        tight(r)
    t0, dur, n = float(cfg['start'][0]), float(cfg['duration'][0]), int(cfg['steps'][0])

    def run():
        r.resetAll()
        return np.array(r.simulate(t0, t0 + dur, n + 1, sel))
    return run, dict(load_s=load)


def task(name):
    if name.startswith('sos:'):
        return sos(name[4:])
    r = te.loada(chain())
    if REF:
        tight(r)
    if name == 'ss':
        def run():
            r.resetAll()
            r.steadyState()
            return np.array(r.getFloatingSpeciesConcentrations())
    elif name == 'fcc':
        r.steadyState()
        run = lambda: np.array(r.getScaledFluxControlCoefficientMatrix())
    elif name == 'sim':
        def run():
            r.resetAll()
            return np.array(r.simulate(0, 100, 1001))
    elif name == 'sens':
        def run():
            r.resetAll()
            return np.array(r.timeSeriesSensitivities(0, 50, 101)[1])
    elif name == 'cc200':
        ks = np.linspace(0.5, 3, 200)
        if hasattr(te, 'controlCoefficients') and not REF:
            run = lambda: te.controlCoefficients(r, {'k1': ks})['fcc'][0]
        else:  # regular tellurium: per-set loop, timed on 20 sets and scaled to 200
            def run():
                for k in ks[:20]:
                    r.resetAll()
                    r.k1 = k
                    r.steadyState()
                    f = np.array(r.getScaledFluxControlCoefficientMatrix())
                    r.getScaledConcentrationControlCoefficientMatrix()
                    if k == ks[0]:
                        first = f
                return first
            return run, dict(scale=10)  # 20 sets timed, x10 for 200
        return run, {}
    return run, {}


if __name__ == '__main__':
    name, out = sys.argv[1], sys.argv[2]
    run, extra = task(name)
    slow = REF or name in ('cc200', 'sens')
    t, ratio, reps, y = timed(run, budget=0 if slow else 2.0, minrep=1 if slow else 3)
    t *= extra.get('scale', 1)
    np.savez(out, y=np.asarray(y, dtype=float))
    rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    print(json.dumps(dict(task=name, s=t, cores=ratio, reps=reps, rss_mb=rss / 2 ** 20, **extra)))

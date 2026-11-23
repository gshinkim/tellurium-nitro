"""Big-model benchmark: regular tellurium (roadrunner) vs Tellurium Nitro on a 1000-species, 3400-reaction model.

python benchmarks/bigmodel/run.py      (regular tellurium is timed too when the `tellurium` package is installed)
"""
import os
import time
import numpy as np
import tellurium_nitro

MODEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'synthetic_1000_species.ant')


def bench(te, name):
    t0 = time.perf_counter()
    r = te.loada(MODEL)
    t1 = time.perf_counter()
    res = r.simulate(0, 100, 1000)
    t2 = time.perf_counter()
    r.reset()
    r.simulate(0, 100, 1000)  # second run: model already built/compiled
    t3 = time.perf_counter()
    print(f"{name:>10}: load {t1 - t0:7.3f} s | simulate {t2 - t1:7.3f} s | re-simulate {t3 - t2:7.3f} s")
    return np.asarray(res)


if __name__ == '__main__':
    nx = bench(tellurium_nitro, 'nitro')
    try:
        import tellurium
    except ImportError:
        tellurium = None
    if tellurium is not None:
        print("max abs diff:", np.nanmax(np.abs(bench(tellurium, 'roadrunner') - nx)))

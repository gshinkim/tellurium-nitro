"""Run bench.py for regular tellurium, Tellurium Turbo and Tellurium Nitro; write results.json + results.md.

    python benchmarks/run.py   # env vars TE_REGULAR / TE_TURBO / TE_NITRO = python executables
    ONLY=nitro python benchmarks/run.py [tasks]   # rerun one engine; others kept from results.json
References are cached in benchmarks/ref/ (computed with the regular engine when missing).
Each task runs in a fresh process. Accuracy = max |x - x_ref| / max |x_ref|, x_ref = plain roadrunner at rtol 1e-12, atol 1e-14, FD step 1e-5.
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HOME = os.path.expanduser('~')
ENGINES = {
    'regular': os.environ.get('TE_REGULAR', f'{HOME}/venvs/te-regular/bin/python'),
    'turbo': os.environ.get('TE_TURBO', f'{HOME}/venvs/te-turbo/bin/python'),
    'nitro': os.environ.get('TE_NITRO', f'{HOME}/venvs/tellurium/bin/python'),
}
if __name__ == '__main__':
    CASES = sorted(os.listdir(os.path.join(HERE, 'sosbench')))
    TASKS = [f'sos:{c}' for c in CASES] + ['ss', 'fcc', 'sim', 'sens', 'cc200', 'load']

    RUN = {e: p for e, p in ENGINES.items() if os.environ.get('ONLY', e) == e}
    NOREF = {'sens'}  # tight-tolerance roadrunner sensitivities take >15 min
    TASKS = [] if sys.argv[1:] == ['--table'] else (sys.argv[1:] or TASKS)  # rerun selected tasks; others are kept from results.json
    OUT = os.path.join(HERE, 'results.json')
    res, tmp = (json.load(open(OUT)) if sys.argv[1:] or len(RUN) < 3 else {}), tempfile.mkdtemp()
    os.makedirs(os.path.join(HERE, 'ref'), exist_ok=True)
    for task in TASKS:
        ref = os.path.join(HERE, 'ref', f"ref-{task.replace(':', '_')}.npz")
        p = None if os.path.exists(ref) or task in NOREF or 'regular' not in RUN else subprocess.run([ENGINES['regular'], os.path.join(HERE, 'bench.py'), task, ref], capture_output=True, text=True,
                           cwd=tmp, env={**os.environ, 'PYTHONWARNINGS': 'ignore', 'BENCH_REF': '1'})
        if p is not None and p.returncode:
            print(task, 'reference FAILED', p.stderr.strip().splitlines()[-1:], file=sys.stderr)
        for eng, py in RUN.items():
            out = os.path.join(tmp, f"{eng}-{task.replace(':', '_')}.npz")
            p = subprocess.run([py, os.path.join(HERE, 'bench.py'), task, out], capture_output=True, text=True,
                               cwd=tmp, env={**os.environ, 'PYTHONWARNINGS': 'ignore'})
            if p.returncode:
                print(task, eng, 'FAILED', p.stderr.strip().splitlines()[-1:], file=sys.stderr)
                continue
            r = json.loads(p.stdout.strip().splitlines()[-1])
            if task in NOREF and eng == 'regular':  # no tight reference (too slow): regular's output, marked *
                import shutil; shutil.copy(out, ref)
            rf = ref
            if task in NOREF and eng != 'regular':
                r['vs_regular'] = True
            if os.path.exists(rf):
                a, b = np.load(out)['y'], np.load(rf)['y']
                r['err'] = float(np.abs(a - b).max() / np.abs(b).max()) if a.shape == b.shape else None
            res.setdefault(task, {})[eng] = r
            print(task, eng, r, file=sys.stderr)

    json.dump(res, open(OUT, 'w'), indent=1)


    def fmt(s):
        return f'{s * 1e3:.3g} ms' if s < 1 else f'{s:.3g} s'


    rows = ['| Task | Regular tellurium | Tellurium Turbo | Tellurium Nitro | Nitro vs regular | Nitro vs Turbo | Error (regular / Turbo / Nitro) |',
            '|---|---|---|---|---|---|---|']
    for task, r in res.items():
        g, t, n = (r.get(e) for e in ENGINES)
        cell = lambda x: '—' if x is None else f"{fmt(x['s'])} · {x['cores']:.1f} cores · {x['rss_mb']:.0f} MB"
        sp = lambda a, b: '—' if a is None or b is None else f"**{a['s'] / b['s']:,.0f}x**" if a['s'] / b['s'] >= 10 else f"**{a['s'] / b['s']:.1f}x**"
        err = ' / '.join('—' if x is None or x.get('err') is None else f"{x['err']:.0e}" + ('*' if x.get('vs_regular') else '') for x in (g, t, n))
        rows.append(f'| {task} | {cell(g)} | {cell(t)} | {cell(n)} | {sp(g, n)} | {sp(t, n)} | {err} |')
    open(os.path.join(HERE, 'results.md'), 'w').write('\n'.join(rows) + '\n')
    print('\n'.join(rows))

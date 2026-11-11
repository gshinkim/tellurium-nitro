"""Run the same upstream tellurium pytest suite against regular tellurium, Turbo and Nitro.

    python benchmarks/tests.py <upstream tellurium checkout> > benchmarks/tests.md
Tests import whichever tellurium the interpreter has installed (--import-mode=importlib), so all three
engines run identical test files. Reports wall time, passed/failed/skipped, and tests whose outcome differs
from regular tellurium.
"""
import os
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

from run import ENGINES

src = os.path.join(os.path.abspath(sys.argv[1]), 'tellurium', 'tests')
tmp = tempfile.mkdtemp()
outcome, rows = {}, []
for eng, py in ENGINES.items():
    xml = os.path.join(tmp, f'{eng}.xml')
    t = time.perf_counter()
    subprocess.run([py, '-m', 'pytest', src, '-q', '-p', 'no:cacheprovider', '--import-mode=importlib',
                    '--timeout=600', '--continue-on-collection-errors', f'--junitxml={xml}'], cwd=tmp, capture_output=True,
                   env={**os.environ, 'MPLBACKEND': 'Agg', 'PYTHONWARNINGS': 'ignore'})
    wall = time.perf_counter() - t
    res = {}
    for tc in ET.parse(xml).iter('testcase'):
        kind = next((c.tag for c in tc if c.tag in ('failure', 'error', 'skipped')), 'passed')
        res[f"{tc.get('classname')}::{tc.get('name')}"] = {'error': 'failed', 'failure': 'failed'}.get(kind, kind)
    outcome[eng] = res
    n = lambda k: sum(v == k for v in res.values())
    rows.append(f"| {eng} | {wall:.1f} s | {n('passed')} | {n('failed')} | {n('skipped')} |")

print('| Engine | Wall time | Passed | Failed | Skipped |\n|---|---|---|---|---|')
print('\n'.join(rows))
for eng in ('turbo', 'nitro'):
    diff = sorted(k for k, v in outcome[eng].items() if outcome['regular'].get(k) != v)
    print(f'\n{eng}: {len(diff)} tests differ from regular' + ''.join(
        f"\n- `{k}`: regular {outcome['regular'].get(k)}, {eng} {outcome[eng][k]}" for k in diff))

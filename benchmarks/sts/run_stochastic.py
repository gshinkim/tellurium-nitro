"""SBML Test Suite (stochastic, DSMTS) on Tellurium Nitro's Rust SSA.

usage: python sts_stoch.py <cases/stochastic> <out.jsonl> [n_replicates]
Per case: n Gillespie runs (Rust, batched), then the DSMTS tests at every time point:
  Z = sqrt(n) (mean - mu) / sigma in meanRange, Y = sqrt(n / 2) (s^2 / sigma^2 - 1) in sdRange.
"""
import glob, json, os, sys, time, warnings
import numpy as np

warnings.filterwarnings('ignore')
import tellurium  # noqa: F401  (import order: libsedml before libsbml)
from tellurium.analysis.rrs_bridge import from_sbml

LEVELS = ['l3v2', 'l3v1', 'l2v5', 'l2v4', 'l2v3', 'l2v2', 'l2v1']


def settings(path):
    s = {}
    for line in open(path):
        if ':' in line:
            k, v = line.split(':', 1)
            s[k.strip()] = v.strip()
    rng = lambda k: tuple(float(x) for x in s[k].strip('()').split(','))
    return dict(start=float(s['start']), duration=float(s['duration']), steps=int(s['steps']),
                variables=[x.strip() for x in s['variables'].split(',') if x.strip()],
                mean=rng('meanRange'), sd=rng('sdRange'))


def run(root, cid, n):
    d = os.path.join(root, cid)
    tags = {}
    for line in open(os.path.join(d, cid + '-model.m')):
        for k in ('componentTags', 'testTags'):
            if line.startswith(k + ':'):
                tags[k] = [x.strip() for x in line.split(':', 1)[1].split(',') if x.strip()]
    rec = dict(id=cid, **tags)
    lv = next((l for l in LEVELS if os.path.exists(os.path.join(d, '%s-sbml-%s.xml' % (cid, l)))), None)
    st = settings(os.path.join(d, cid + '-settings.txt'))
    exp = np.genfromtxt(os.path.join(d, cid + '-results.csv'), delimiter=',', names=True, deletechars='')
    times = np.linspace(st['start'], st['start'] + st['duration'], st['steps'] + 1)
    t0 = time.perf_counter()
    try:
        info = {}
        m = from_sbml(open(os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)), encoding='utf-8').read(), info, d)
        comp_of = info['comp_of']
        outs = st['variables'] + sorted({comp_of[v] for v in st['variables'] if v in comp_of})
        Y = np.asarray(m.gillespie_batch([], np.zeros((n, 0)), times, seed=12345, outputs=outs))
    except Exception as e:
        return dict(rec, status='error', reason='%s: %s' % (type(e).__name__, str(e)[:200]))
    rec['seconds'] = round(time.perf_counter() - t0, 3)
    col = {o: k for k, o in enumerate(outs)}
    zbad = ybad = checks = 0
    worst = (0.0, None)
    for k, v in enumerate(st['variables']):
        x = Y[:, :, k] * (Y[:, :, col[comp_of[v]]] if v in comp_of else 1.0)  # amounts
        mu, sig = exp[v + '-mean'], exp[v + '-sd']
        mean, var = x.mean(0), x.var(0, ddof=1)
        for q in range(len(times)):
            checks += 2
            if sig[q] == 0:  # deterministic point: must match exactly
                ok = abs(mean[q] - mu[q]) < 1e-9 and var[q] < 1e-18
                zbad += not ok
                continue
            z = np.sqrt(n) * (mean[q] - mu[q]) / sig[q]
            y = np.sqrt(n / 2) * (var[q] / sig[q] ** 2 - 1)
            zbad += not (st['mean'][0] <= z <= st['mean'][1])
            ybad += not (st['sd'][0] <= y <= st['sd'][1])
            worst = max(worst, (abs(z), '%s Z=%.2f at t=%g' % (v, z, times[q])), (abs(y) * 3 / 5, '%s Y=%.2f at t=%g' % (v, y, times[q])),
                        key=lambda w: w[0])
    bad = zbad + ybad
    return dict(rec, status='pass' if bad == 0 else 'fail', z_out=int(zbad), y_out=int(ybad), checks=checks,
                worst=worst[1], n=n)


if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
    with open(out, 'w') as f:
        for cid in sorted(x for x in os.listdir(root) if x.isdigit()):
            r = run(root, cid, n)
            f.write(json.dumps(r) + '\n')
            f.flush()
            print(cid, r['status'], r.get('z_out', ''), r.get('y_out', ''), r.get('seconds', ''), r.get('reason', r.get('worst', '')), flush=True)

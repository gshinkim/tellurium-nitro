"""SBML Test Suite (stochastic, DSMTS) on Tellurium Nitro's Rust SSA.

usage: python sts_stoch.py <cases/stochastic> <out.jsonl> [n_replicates]
Per case: n replicates (Rust, batched: Gillespie, or one ODE run each when nothing reacts), then the
suite's tests at every time point for the outputs settings.txt requests (X or ln(X)):
  Z = sqrt(n) (mean - mu) / sigma in meanRange, Y = sqrt(n / 2) (s^2 / sigma^2 - 1) in sdRange;
  a mean without an expected sd is the median for an untruncated Cauchy (distribution-free 3-sigma
  interval for the sample median) and else a t-test with the sample sd (truncated Cauchy).
  Y assumes normal kurtosis (3); a Y outside sdRange is re-tested with the kurtosis-aware standard
  error of s^2, sqrt((m4 - s^4 (n-3)/(n-1)) / n); such points are counted in 'y_kurtosis'.
Pass rule (the suite's README): Z and Y each fail at most once per test. With 2 or more failures
the case is re-run with two other seeds: 'pass-recheck' if the failures do not recur there.
"""
import glob, json, os, sys, time, warnings
import numpy as np

warnings.filterwarnings('ignore')
import tellurium_nitro  # noqa: F401  (import order: libsedml before libsbml)
from tellurium_nitro.analysis.rrs_bridge import from_sbml

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
                output=[x.strip() for x in s['output'].split(',') if x.strip()],
                mean=rng('meanRange'), sd=rng('sdRange'))


def run(root, cid, n, seed=12345):
    d = os.path.join(root, cid)
    tags = {}
    for line in open(os.path.join(d, cid + '-model.m')):
        for k in ('componentTags', 'testTags'):
            if line.startswith(k + ':'):
                tags[k] = [x.strip() for x in line.split(':', 1)[1].split(',') if x.strip()]
    rec = dict(id=cid, **tags)
    synopsis = open(os.path.join(d, cid + '-model.m')).read().split('synopsis:')[1].split('\n')[0].lower()
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
        if m.reactions:  # Gillespie
            Y = np.asarray(m.gillespie_batch([], np.zeros((n, 0)), times, seed=seed, outputs=outs))
        else:  # nothing fires at random: one ODE run per replicate (distrib draws in events, own stream each)
            Y = np.asarray(m.simulate_batch([], np.zeros((n, 0)), times, None, 1e-10, 1e-12, outs, seed=seed))
    except Exception as e:
        return dict(rec, status='error', reason='%s: %s' % (type(e).__name__, str(e)[:200]))
    rec['seconds'] = round(time.perf_counter() - t0, 3)
    col = {o: k for k, o in enumerate(outs)}
    zbad, ybad, ykurt, checks, where = 0, 0, 0, 0, []
    names = list(exp.dtype.names)
    for item in st['output']:
        name, stat = item.rsplit('-', 1)
        v = name[3:-1] if name.startswith('ln(') else name
        x = Y[:, :, col[v]] * (Y[:, :, col[comp_of[v]]] if v in comp_of else 1.0)  # amounts
        if name != v:
            with np.errstate(all='ignore'):
                x = np.log(x)
        mu = exp[name + '-mean']
        sig = exp[name + '-sd'] if name + '-sd' in names else None
        for q in range(len(times)):
            checks += 1
            if stat == 'mean' and sig is None and 'truncation' in synopsis:  # a defined mean: t-test
                z = np.sqrt(n) * (x[:, q].mean() - mu[q]) / x[:, q].std(ddof=1)
                ok = st['mean'][0] <= z <= st['mean'][1]
            elif stat == 'mean' and sig is None:  # median: order-statistic 3-sigma interval
                srt = np.sort(x[:, q])
                k = int(np.floor(n / 2 - 1.5 * np.sqrt(n)))
                ok = srt[max(k, 0)] <= mu[q] <= srt[min(n - 1 - k, n - 1)]
            elif sig[q] == 0:  # deterministic point: exact
                ok = abs(x[:, q].mean() - mu[q]) < 1e-9 and x[:, q].var() < 1e-18
            elif stat == 'mean':
                z = np.sqrt(n) * (x[:, q].mean() - mu[q]) / sig[q]
                ok = st['mean'][0] <= z <= st['mean'][1]
            else:
                s2 = x[:, q].var(ddof=1)
                y = np.sqrt(n / 2) * (s2 / sig[q] ** 2 - 1)
                ok = st['sd'][0] <= y <= st['sd'][1]
                if not ok:  # the same test with the sample's own kurtosis
                    m4 = ((x[:, q] - x[:, q].mean()) ** 4).mean()
                    yk = (s2 - sig[q] ** 2) / np.sqrt(max(m4 - s2 ** 2 * (n - 3) / (n - 1), 1e-300) / n)
                    ok = st['sd'][0] <= yk <= st['sd'][1]
                    ykurt += ok
            if not ok:
                where.append('%s-%s@%g' % (name, stat, times[q]))
                if stat == 'mean':
                    zbad += 1
                else:
                    ybad += 1
    return dict(rec, status='pass' if zbad <= 1 and ybad <= 1 else 'fail', z_out=zbad, y_out=ybad, checks=checks,
                failed_at=where[:6], y_kurtosis=int(ykurt), n=n, seed=seed)


def judge(root, cid, n):
    r = run(root, cid, n)
    if r['status'] != 'fail':
        return r
    again = [run(root, cid, n, seed) for seed in (424242, 777)]
    first = set(r['failed_at'])
    recurs = any(a['status'] != 'pass' or first & set(a['failed_at']) for a in again)
    return dict(r, status='fail' if recurs else 'pass-recheck', recheck=[(a['z_out'], a['y_out']) for a in again])

if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
    with open(out, 'w') as f:
        for cid in sorted(x for x in os.listdir(root) if x.isdigit()):
            r = judge(root, cid, n)
            f.write(json.dumps(r) + '\n')
            f.flush()
            print(cid, r['status'], r.get('z_out', ''), r.get('y_out', ''), r.get('reason', r.get('failed_at', '')), r.get('recheck', ''), flush=True)

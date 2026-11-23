"""SBML Test Suite (stochastic, DSMTS) on the Rust core alone: roadrunner_nitro.RoadRunner's
gillespieRuns (no tellurium_nitro). Same tests and pass rule as run_stochastic.py.

usage: python run_core_stochastic.py <cases/stochastic> <out.jsonl> [n_replicates]
"""
import json, os, sys, time, warnings
import numpy as np
import roadrunner_nitro

warnings.filterwarnings('ignore')
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
        r = roadrunner_nitro.RoadRunner(os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)))
        r.setIntegratorSetting('relative_tolerance', 1e-10)  # reaction-free models: one ODE run per seed
        r.setIntegratorSetting('absolute_tolerance', 1e-12)
        r.timeCourseSelections = st['variables']  # species ids select amounts
        Y = np.asarray(r.gillespieRuns(n, times[0], times[-1], len(times), seed=seed))
    except Exception as e:
        return dict(rec, status='error', reason='%s: %s' % (type(e).__name__, str(e)[:200]))
    rec['seconds'] = round(time.perf_counter() - t0, 3)
    col = {o: k for k, o in enumerate(st['variables'])}
    zbad, ybad, ykurt, checks, where = 0, 0, 0, 0, []
    names = list(exp.dtype.names)
    for item in st['output']:
        name, stat = item.rsplit('-', 1)
        v = name[3:-1] if name.startswith('ln(') else name
        x = Y[:, :, col[v]]
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

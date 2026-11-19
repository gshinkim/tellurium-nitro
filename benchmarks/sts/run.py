"""SBML Test Suite (semantic) on Tellurium Nitro.

worker: python sts_run.py <cases_dir> <out.jsonl> <id> [<id> ...]   -> one JSON line per case
driver: python sts_run.py <cases_dir> <out.jsonl> --all [nproc]      -> shards, restarts after crashes/hangs
"""
import json, os, subprocess, sys, time

LEVELS = ['l3v2', 'l3v1', 'l2v5', 'l2v4', 'l2v3', 'l2v2', 'l2v1', 'l1v2']
RTOL, ATOL = 1e-8, 1e-10  # integrator tolerances (the suite's own pass tolerances come from settings)


def settings(path):
    s = {}
    for line in open(path):
        if ':' in line:
            k, v = line.split(':', 1)
            s[k.strip()] = v.strip()
    lst = lambda k: [x.strip() for x in s.get(k, '').split(',') if x.strip()]
    return dict(start=float(s['start']), duration=float(s['duration']), steps=int(s['steps']),
                variables=lst('variables'), amount=set(lst('amount')), concentration=set(lst('concentration')),
                atol=float(s['absolute']), rtol=float(s['relative']))


def tags(path):
    out = {}
    for line in open(path):
        for k in ('componentTags', 'testTags', 'testType'):
            if line.startswith(k + ':'):
                out[k] = [x.strip() for x in line.split(':', 1)[1].split(',') if x.strip()]
    return out


def fba_case(d, cid, rec):
    """FluxBalanceSteadyState: te.fba on the newest level, compared like a time course row."""
    import numpy as np
    import tellurium as te
    st = {k.strip(): v.strip() for k, v in (l.split(':', 1) for l in open(os.path.join(d, cid + '-settings.txt')) if ':' in l)}
    names = [x.strip() for x in st['variables'].split(',') if x.strip()]
    atol, rtol = float(st['absolute']), float(st['relative'])
    lv = next((l for l in LEVELS if os.path.exists(os.path.join(d, '%s-sbml-%s.xml' % (cid, l)))), None)
    rec['level'] = lv
    exp = np.genfromtxt(os.path.join(d, cid + '-results.csv'), delimiter=',', names=True, deletechars='')
    try:
        got = te.fba(os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)))
    except Exception as e:
        return dict(rec, status='error', stage='fba', reason='%s: %s' % (type(e).__name__, str(e)[:200]))
    rec['engine'] = 'fba'
    bad = [(n, got.get(n), float(exp[n])) for n in names
           if not (np.isnan(float(exp[n])) and np.isnan(got.get(n, 0.0))
                   or abs(got.get(n, np.nan) - float(exp[n])) <= atol + rtol * abs(float(exp[n])))]
    if bad:
        n, g, e = bad[0]
        return dict(rec, status='fail', first_bad=dict(column=n, t=0.0, got=g, expected=e), n_bad=len(bad))
    return dict(rec, status='pass')


def run_case(root, cid):
    import numpy as np
    import tellurium as te
    d = os.path.join(root, cid)
    rec = dict(id=cid, **tags(os.path.join(d, cid + '-model.m')))
    if 'FluxBalanceSteadyState' in rec.get('testType', []):
        return fba_case(d, cid, rec)
    if 'TimeCourse' not in rec.get('testType', []):
        return dict(rec, status='skip', reason='not a time course: %s' % rec.get('testType'))
    lv = next((l for l in LEVELS if os.path.exists(os.path.join(d, '%s-sbml-%s.xml' % (cid, l)))), None)
    if lv is None:
        return dict(rec, status='skip', reason='no SBML file')
    rec['level'] = lv
    st = settings(os.path.join(d, cid + '-settings.txt'))
    exp = np.genfromtxt(os.path.join(d, cid + '-results.csv'), delimiter=',', names=True, deletechars='')
    t0 = time.perf_counter()
    try:
        r = te.loadSBMLModel(os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)))
    except Exception as e:
        return dict(rec, status='error', stage='load', reason='%s: %s' % (type(e).__name__, str(e)[:300]))
    # the time column is the requested grid (a model id may itself be named 'time')
    sel = ['[%s]' % v if v in st['concentration'] else v for v in st['variables']]
    try:
        r.integrator.relative_tolerance = RTOL
        r.integrator.absolute_tolerance = ATOL
        res = np.asarray(r.simulate(st['start'], st['start'] + st['duration'], st['steps'] + 1, selections=sel))
        res = np.column_stack([np.linspace(st['start'], st['start'] + st['duration'], st['steps'] + 1), res.reshape(st['steps'] + 1, -1)])
    except Exception as e:
        return dict(rec, status='error', stage='simulate', reason='%s: %s' % (type(e).__name__, str(e)[:300]),
                    engine=(r.lastRun or {}).get('engine'), fallback=(r.lastRun or {}).get('reason'))
    run = r.lastRun or {}
    rec.update(engine=run.get('engine'), fallback=run.get('reason'), solver=run.get('solver'),
               seconds=round(time.perf_counter() - t0, 4))
    want = np.column_stack([exp[n] for n in exp.dtype.names])
    if want.shape != res.shape:
        return dict(rec, status='fail', reason='shape %s vs expected %s' % (res.shape, want.shape))
    ok = np.isclose(res, want, rtol=0, atol=0, equal_nan=True) | \
        (np.abs(res - want) <= st['atol'] + st['rtol'] * np.abs(want))
    if ok.all():
        return dict(rec, status='pass')
    bad = np.argwhere(~ok)
    i, j = bad[0]
    return dict(rec, status='fail', n_bad=int(len(bad)), first_bad=dict(
        column=(['time'] + st['variables'])[j], t=float(want[i, 0]),
        got=float(res[i, j]), expected=float(want[i, j])),
        max_rel=float(np.nanmax(np.abs(res - want)[~ok] / (np.abs(want[~ok]) + st['atol']))))


def worker(root, out, ids):
    import warnings
    warnings.filterwarnings('ignore')
    with open(out, 'a') as f:
        for cid in ids:
            with open(out + '.cur', 'w') as c:
                c.write(cid)
            try:
                rec = run_case(root, cid)
            except Exception as e:  # harness bug, not a model result
                rec = dict(id=cid, status='error', stage='harness', reason='%s: %s' % (type(e).__name__, e))
            f.write(json.dumps(rec) + '\n')
            f.flush()


def driver(root, out, nproc, timeout=120):
    ids = sorted(x for x in os.listdir(root) if x.isdigit())
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)['id'] for l in open(out) if l.strip()}
    todo = [c for c in ids if c not in done]
    shards = [todo[k::nproc] for k in range(nproc)]
    procs = {}

    def launch(k):
        rest = [c for c in shards[k] if c not in done]
        if not rest:
            return
        p = subprocess.Popen([sys.executable, '-I', __file__, root, out + '.%d' % k] + rest,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        procs[k] = dict(p=p, cur=None, since=time.time())

    for k in range(nproc):
        launch(k)
    while procs:
        time.sleep(0.5)
        for k in list(procs):
            s = procs[k]
            try:
                cur = open(out + '.%d.cur' % k).read().strip() or None
            except OSError:
                cur = None
            if cur != s['cur']:
                s['cur'], s['since'] = cur, time.time()
            dead = s['p'].poll() is not None
            hung = not dead and s['cur'] and time.time() - s['since'] > timeout
            if dead or hung:
                if hung:
                    s['p'].kill()
                    s['p'].wait()
                part = out + '.%d' % k
                got = {json.loads(l)['id'] for l in open(part) if l.strip()} if os.path.exists(part) else set()
                done |= got
                if s['cur'] and s['cur'] not in got:  # crashed or hung on this case
                    with open(part, 'a') as f:
                        f.write(json.dumps(dict(id=s['cur'], status='error', stage='timeout' if hung else 'crash',
                                                reason='no result after %ds' % timeout if hung else
                                                'process died (exit %s)' % s['p'].returncode)) + '\n')
                    done.add(s['cur'])
                del procs[k]
                launch(k)
    with open(out, 'a') as f:
        for k in range(nproc):
            part = out + '.%d' % k
            if os.path.exists(part):
                f.write(open(part).read())
                os.remove(part)
            if os.path.exists(part + '.cur'):
                os.remove(part + '.cur')


if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    if sys.argv[3] == '--all':
        driver(root, out, int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count())
    else:
        worker(root, out, sys.argv[3:])

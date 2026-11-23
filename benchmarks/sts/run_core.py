"""SBML Test Suite (semantic time courses) on the Rust core alone: roadrunner_nitro.RoadRunner, no
tellurium_nitro, no libroadrunner. Same judging as run.py; flux-balance cases are skipped (te.fba).

worker: python run_core.py <cases_dir> <out.jsonl> <id> [<id> ...]
driver: python run_core.py <cases_dir> <out.jsonl> --all [nproc]
Workers inherit PYTHONPATH (no -I), so a locally built roadrunner_nitro can be put first.
"""
import json, os, subprocess, sys, time

import run as base  # settings(), tags(), LEVELS, RTOL / ATOL and the worker loop

LEVELS, RTOL, ATOL, settings, tags = base.LEVELS, base.RTOL, base.ATOL, base.settings, base.tags


def run_case(root, cid):
    import numpy as np
    import roadrunner_nitro
    d = os.path.join(root, cid)
    rec = dict(id=cid, **tags(os.path.join(d, cid + '-model.m')))
    if 'TimeCourse' not in rec.get('testType', []) or 'FluxBalanceSteadyState' in rec.get('testType', []):
        return dict(rec, status='skip', reason='not a time course: %s' % rec.get('testType'))
    lv = next((l for l in LEVELS if os.path.exists(os.path.join(d, '%s-sbml-%s.xml' % (cid, l)))), None)
    if lv is None:
        return dict(rec, status='skip', reason='no SBML file')
    rec['level'] = lv
    st = settings(os.path.join(d, cid + '-settings.txt'))
    exp = np.genfromtxt(os.path.join(d, cid + '-results.csv'), delimiter=',', names=True, deletechars='')
    t0 = time.perf_counter()
    try:
        r = roadrunner_nitro.RoadRunner(os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)))
    except Exception as e:
        return dict(rec, status='error', stage='load', reason='%s: %s' % (type(e).__name__, str(e)[:300]))
    sel = ['[%s]' % v if v in st['concentration'] else v for v in st['variables']]
    try:
        r.setIntegratorSetting('relative_tolerance', RTOL)
        r.setIntegratorSetting('absolute_tolerance', ATOL)
        res = np.asarray(r.simulate(st['start'], st['start'] + st['duration'], st['steps'] + 1, selections=sel))
        res = np.column_stack([np.linspace(st['start'], st['start'] + st['duration'], st['steps'] + 1), res.reshape(st['steps'] + 1, -1)])
    except Exception as e:
        return dict(rec, status='error', stage='simulate', reason='%s: %s' % (type(e).__name__, str(e)[:300]))
    rec.update(engine='core', seconds=round(time.perf_counter() - t0, 4))
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
        got=float(res[i, j]), expected=float(want[i, j])))


base.run_case = run_case  # base.worker runs these cases


def driver(root, out, nproc, timeout=120):
    """run.py's driver, launching this file's workers without -I."""
    ids = sorted(x for x in os.listdir(root) if x.isdigit())
    done = {json.loads(l)['id'] for l in open(out) if l.strip()} if os.path.exists(out) else set()
    shards = [[c for c in ids if c not in done][k::nproc] for k in range(nproc)]
    procs = {}

    def launch(k):
        rest = [c for c in shards[k] if c not in done]
        if rest:
            procs[k] = dict(p=subprocess.Popen([sys.executable, __file__, root, out + '.%d' % k] + rest,
                                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL), cur=None, since=time.time())

    for k in range(nproc):
        launch(k)
    while procs:
        time.sleep(0.5)
        for k in list(procs):
            s, part = procs[k], out + '.%d' % k
            try:
                cur = open(part + '.cur').read().strip() or None
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
                got = {json.loads(l)['id'] for l in open(part) if l.strip()} if os.path.exists(part) else set()
                done |= got
                if s['cur'] and s['cur'] not in got:
                    with open(part, 'a') as f:
                        f.write(json.dumps(dict(id=s['cur'], status='error', stage='timeout' if hung else 'crash')) + '\n')
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
        base.worker(root, out, sys.argv[3:])

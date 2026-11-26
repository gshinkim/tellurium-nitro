"""SBML Test Suite (semantic time courses) through the roadrunner_nitro C API alone (the libroadrunner
rrc_api.h drop-in, libroadrunner_nitro/capi) via ctypes: no Python engine. Same judging as run.py / run_core.py;
flux-balance cases are skipped (te.fba).

worker: python run_capi.py <cases_dir> <out.jsonl> <id> [<id> ...]
driver: python run_capi.py <cases_dir> <out.jsonl> --all [nproc]
Library: $RRC_LIB, else ../libroadrunner_nitro/target/release/libroadrunner_nitro_c_api.{dylib,so}.
"""
import ctypes as C, os, sys, time

import run as base
import run_core

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.environ.get('RRC_LIB') or os.path.join(HERE, '..', '..', '..', 'libroadrunner_nitro', 'target', 'release',
                                                  {'darwin': 'libroadrunner_nitro_c_api.dylib', 'win32': 'roadrunner_nitro_c_api.dll'}
                                                  .get(sys.platform, 'libroadrunner_nitro_c_api.so'))


class RRCData(C.Structure):
    _fields_ = [('RSize', C.c_int), ('CSize', C.c_int), ('Data', C.POINTER(C.c_double)),
                ('Weights', C.POINTER(C.c_double)), ('ColumnHeaders', C.POINTER(C.c_char_p))]


def lib():
    L = C.CDLL(LIB)
    for name, res, args in [
        ('createRRInstance', C.c_void_p, []), ('freeRRInstance', C.c_bool, [C.c_void_p]),
        ('loadSBMLFromFileE', C.c_bool, [C.c_void_p, C.c_char_p, C.c_bool]),
        ('setCurrentIntegratorParameterDouble', C.c_int, [C.c_void_p, C.c_char_p, C.c_double]),
        ('setTimeCourseSelectionList', C.c_bool, [C.c_void_p, C.c_char_p]),
        ('simulateEx', C.POINTER(RRCData), [C.c_void_p, C.c_double, C.c_double, C.c_int]),
        ('getRRDataNumRows', C.c_int, [C.POINTER(RRCData)]), ('getRRDataNumCols', C.c_int, [C.POINTER(RRCData)]),
        ('getRRCDataElement', C.c_bool, [C.POINTER(RRCData), C.c_int, C.c_int, C.POINTER(C.c_double)]),
        ('freeRRCData', C.c_int, [C.POINTER(RRCData)]), ('getLastError', C.c_char_p, []),
    ]:
        f = getattr(L, name)
        f.restype, f.argtypes = res, args
    return L


_L = None


def run_case(root, cid):
    import numpy as np
    global _L
    _L = _L or lib()
    L = _L
    d = os.path.join(root, cid)
    rec = dict(id=cid, **base.tags(os.path.join(d, cid + '-model.m')))
    if 'TimeCourse' not in rec.get('testType', []) or 'FluxBalanceSteadyState' in rec.get('testType', []):
        return dict(rec, status='skip', reason='not a time course: %s' % rec.get('testType'))
    lv = next((l for l in base.LEVELS if os.path.exists(os.path.join(d, '%s-sbml-%s.xml' % (cid, l)))), None)
    if lv is None:
        return dict(rec, status='skip', reason='no SBML file')
    rec['level'] = lv
    st = base.settings(os.path.join(d, cid + '-settings.txt'))
    exp = np.genfromtxt(os.path.join(d, cid + '-results.csv'), delimiter=',', names=True, deletechars='')
    err = lambda: L.getLastError().decode(errors='replace').strip()[:300]
    t0 = time.perf_counter()
    h = L.createRRInstance()
    try:
        if not L.loadSBMLFromFileE(h, os.path.join(d, '%s-sbml-%s.xml' % (cid, lv)).encode(), True):
            return dict(rec, status='error', stage='load', reason=err())
        sel = ['[%s]' % v if v in st['concentration'] else v for v in st['variables']]
        ok = L.setCurrentIntegratorParameterDouble(h, b'relative_tolerance', base.RTOL) == 1 and \
            L.setCurrentIntegratorParameterDouble(h, b'absolute_tolerance', base.ATOL) == 1 and \
            L.setTimeCourseSelectionList(h, ', '.join(sel).encode())
        n = st['steps'] + 1
        res = L.simulateEx(h, st['start'], st['start'] + st['duration'], n) if ok else None
        if not res:
            return dict(rec, status='error', stage='simulate', reason=err())
        nr, nc = L.getRRDataNumRows(res), L.getRRDataNumCols(res)
        v = C.c_double()
        got = np.array([[L.getRRCDataElement(res, i, j, C.byref(v)) and v.value for j in range(nc)] for i in range(nr)], float).reshape(nr, nc)
        L.freeRRCData(res)
    finally:
        L.freeRRInstance(h)
    # the time column is the requested grid (a model id may itself be named 'time')
    got = np.column_stack([np.linspace(st['start'], st['start'] + st['duration'], nr), got])
    rec.update(engine='capi', seconds=round(time.perf_counter() - t0, 4))
    want = np.column_stack([exp[k] for k in exp.dtype.names])
    if want.shape != got.shape:
        return dict(rec, status='fail', reason='shape %s vs expected %s' % (got.shape, want.shape))
    okm = np.isclose(got, want, rtol=0, atol=0, equal_nan=True) | (np.abs(got - want) <= st['atol'] + st['rtol'] * np.abs(want))
    if okm.all():
        return dict(rec, status='pass')
    bad = np.argwhere(~okm)
    i, j = bad[0]
    return dict(rec, status='fail', n_bad=int(len(bad)), first_bad=dict(
        column=(['time'] + st['variables'])[j], t=float(want[i, 0]), got=float(got[i, j]), expected=float(want[i, j])))


base.run_case = run_case
run_core.__file__ = os.path.abspath(__file__)  # the driver launches this file's workers

if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    if sys.argv[3] == '--all':
        run_core.driver(root, out, int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count())
    else:
        base.worker(root, out, sys.argv[3:])

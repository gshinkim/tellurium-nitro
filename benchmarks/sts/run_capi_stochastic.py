"""SBML Test Suite (stochastic, DSMTS) through the roadrunner_nitro C API alone (ctypes): each replicate
is reset + gillespieOnGridEx on one instance seeded once with setSeed. Same tests and pass rule as
run_core_stochastic.py (whose judging this reuses).

usage: python run_capi_stochastic.py <cases/stochastic> <out.jsonl> [n_replicates]
"""
import ctypes as C, json, sys

import numpy as np

import run_capi
import run_core_stochastic as base


class _RR:
    """The slice of roadrunner_nitro.RoadRunner run_core_stochastic.run uses, on the C API."""
    L = None

    def __init__(self, path):
        L = _RR.L = _RR.L or self._lib()
        self.h = L.createRRInstance()
        if not L.loadSBMLFromFileE(self.h, path.encode(), True):
            raise RuntimeError(L.getLastError().decode())

    @staticmethod
    def _lib():
        L = run_capi.lib()
        for name, res, args in [('setSeed', C.c_bool, [C.c_void_p, C.c_int64]), ('reset', C.c_bool, [C.c_void_p]),
                                ('gillespieOnGridEx', C.POINTER(run_capi.RRCData), [C.c_void_p, C.c_double, C.c_double, C.c_int])]:
            f = getattr(L, name)
            f.restype, f.argtypes = res, args
        return L

    def _check(self, ok):
        if not ok:
            raise RuntimeError(self.L.getLastError().decode())

    def setIntegratorSetting(self, key, v):  # the current integrator (cvode): reaction-free models
        self._check(self.L.setCurrentIntegratorParameterDouble(self.h, key.encode(), v) == 1)

    @property
    def timeCourseSelections(self):
        raise AttributeError

    @timeCourseSelections.setter
    def timeCourseSelections(self, sel):
        self._check(self.L.setTimeCourseSelectionList(self.h, ', '.join(sel).encode()))

    def gillespieRuns(self, n, t0, t1, points, seed):
        L = self.L
        self._check(L.setSeed(self.h, seed))
        out = None
        for i in range(n):
            self._check(L.reset(self.h))
            p = L.gillespieOnGridEx(self.h, t0, t1, points)
            self._check(p)
            d = p.contents
            y = np.ctypeslib.as_array(d.Data, shape=(d.RSize, d.CSize))
            if out is None:
                out = np.empty((n, d.RSize, d.CSize))
            out[i] = y
            L.freeRRCData(p)
        return out

    def __del__(self):
        if self.L and getattr(self, 'h', None):
            self.L.freeRRInstance(self.h)


class _Module:
    RoadRunner = _RR


base.roadrunner_nitro = _Module

if __name__ == '__main__':
    root, out = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
    import os
    with open(out, 'w') as f:
        for cid in sorted(x for x in os.listdir(root) if x.isdigit()):
            r = base.judge(root, cid, n)
            f.write(json.dumps(r) + '\n')
            f.flush()
            print(cid, r['status'], r.get('z_out', ''), r.get('y_out', ''), r.get('reason', r.get('failed_at', '')), r.get('recheck', ''), flush=True)

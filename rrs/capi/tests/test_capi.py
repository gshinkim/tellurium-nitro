"""The C API (libroadrunner_nitro_c_api, rrc_api.h) through ctypes, against libroadrunner's recorded
answers (rrs/tests/reference/*.npz, replayed with test_core_vs_roadrunner's scenarios on a C API adapter)
and against the Rust core (roadrunner_nitro.RoadRunner, same process) for the C-only surface: by-index
access, initial conditions, symbol lists, integrator / solver settings, gillespie variants, reset
family, save / load state, editing + regenerate, SBML out, ownership and error paths.

    cargo build --release -p roadrunner_nitro_c_api   (in rrs/)
    PYTHONPATH=<built roadrunner_nitro> python -m pytest rrs/capi/tests/test_capi.py
Library: $RRC_LIB, else rrs/target/release/libroadrunner_nitro_c_api.{dylib,so}.
"""
import ctypes as C
import os
import sys

import numpy as np
import pytest
import roadrunner_nitro as rn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'tests'))
import test_core_vs_roadrunner as T  # noqa: E402

LIB = os.environ.get('RRC_LIB') or os.path.join(HERE, '..', '..', 'target', 'release',
                                                  'libroadrunner_nitro_c_api.' + ('dylib' if sys.platform == 'darwin' else 'so'))
L = C.CDLL(LIB)
H, S, D, I, B, V, P = C.c_void_p, C.c_char_p, C.c_double, C.c_int, C.c_bool, None, C.c_void_p
pD, pI = C.POINTER(C.c_double), C.POINTER(C.c_int)


class RRVector(C.Structure):
    _fields_ = [('Count', C.c_int), ('Data', pD)]


class RRStringArray(C.Structure):
    _fields_ = [('Count', C.c_int), ('String', C.POINTER(C.c_char_p))]


class RRDoubleMatrix(C.Structure):
    _fields_ = [('RSize', C.c_int), ('CSize', C.c_int), ('Data', pD)]


class RRCData(C.Structure):
    _fields_ = [('RSize', C.c_int), ('CSize', C.c_int), ('Data', pD), ('Weights', pD), ('ColumnHeaders', C.POINTER(C.c_char_p))]


class RRList(C.Structure):
    pass


class RRListItem(C.Structure):
    _fields_ = [('ItemType', C.c_int), ('data', C.c_void_p)]  # union: 8 bytes (pointer or double)


RRList._fields_ = [('Count', C.c_int), ('Items', C.POINTER(C.POINTER(RRListItem)))]
pV, pSA, pM, pCD, pL = (C.POINTER(t) for t in (RRVector, RRStringArray, RRDoubleMatrix, RRCData, RRList))

PROTOS = {  # name: (restype, argtypes); handle-first functions list only the arguments after it
    'createRRInstance': (H, []), 'freeRRInstance': (B, [H]), 'hasError': (B, []), 'getLastError': (S, []),
    'freeText': (I, [P]), 'freeVector': (I, [pV]), 'freeStringArray': (I, [pSA]), 'freeMatrix': (I, [pM]),
    'freeRRCData': (I, [pCD]), 'freeRRList': (V, [pL]),
    'getNumRegisteredIntegrators': (I, []), 'getRegisteredIntegratorName': (P, [I]), 'getRegisteredIntegratorHint': (P, [I]),
    'getNumRegisteredSteadyStateSolvers': (I, []), 'getRegisteredSteadyStateSolverName': (P, [I]),
}
HANDLE = {
    'loadSBML': (B, [S]), 'loadSBMLFromFile': (B, [S]), 'isModelLoaded': (B, []), 'clearModel': (B, []),
    'getSBML': (P, []), 'getCurrentSBML': (P, []), 'saveState': (B, [S]), 'loadState': (B, [S]),
    'setComputeAndAssignConservationLaws': (B, [B]), 'getComputeAndAssignConservationLaws': (I, [pI]),
    'setTimeCourseSelectionList': (B, [S]), 'getTimeCourseSelectionList': (pSA, []),
    'setSteadyStateSelectionList': (B, [S]), 'getSteadyStateSelectionList': (pSA, []),
    'simulate': (pCD, []), 'simulateEx': (pCD, [D, D, I]), 'simulateTimes': (pCD, [pD, I]), 'getSimulationResult': (pCD, []),
    'getSimulationResultAsDoubleMatrix': (pM, []), 'oneStep': (B, [D, D, pD]),
    'setTimeStart': (B, [D]), 'setTimeEnd': (B, [D]), 'setNumPoints': (B, [I]),
    'getTimeStart': (B, [pD]), 'getTimeEnd': (B, [pD]), 'getNumPoints': (B, [pI]),
    'steadyState': (B, [pD]), 'computeSteadyStateValues': (pV, []), 'getValue': (B, [S, pD]), 'setValue': (B, [S, D]),
    'getFloatingSpeciesConcentrations': (pV, []), 'getFloatingSpeciesAmounts': (pV, []), 'getRateRuleValues': (pV, []),
    'getIndependentFloatingSpeciesConcentrations': (pV, []), 'getDependentFloatingSpeciesConcentrations': (pV, []),
    'getIndependentFloatingSpeciesAmounts': (pV, []), 'getDependentFloatingSpeciesAmounts': (pV, []),
    'getBoundarySpeciesConcentrations': (pV, []), 'getBoundarySpeciesAmounts': (pV, []), 'getGlobalParameterValues': (pV, []),
    'getFloatingSpeciesInitialConcentrations': (pV, []), 'getReactionRates': (pV, []), 'getRatesOfChange': (pV, []),
    'getReactionRatesEx': (pV, [pV]), 'getRatesOfChangeEx': (pV, [pV]),
    'setFloatingSpeciesConcentrations': (B, [pV]), 'setBoundarySpeciesConcentrations': (B, [pV]),
    'setFloatingSpeciesInitialConcentrations': (B, [pV]),
    'setBoundarySpeciesByIndex': (B, [C.c_uint, D]), 'getBoundarySpeciesByIndex': (B, [I, pD]),
    'setFloatingSpeciesByIndex': (B, [I, D]), 'getFloatingSpeciesByIndex': (B, [I, pD]),
    'setFloatingSpeciesInitialConcentrationByIndex': (B, [I, D]), 'getFloatingSpeciesInitialConcentrationByIndex': (B, [I, pD]),
    'setGlobalParameterByIndex': (B, [I, D]), 'getGlobalParameterByIndex': (B, [I, pD]),
    'setCompartmentByIndex': (B, [I, D]), 'getCompartmentByIndex': (B, [I, pD]),
    'getReactionRate': (B, [I, pD]), 'getRateOfChange': (B, [I, pD]), 'evalModel': (B, []),
    'getFullJacobian': (pM, []), 'getReducedJacobian': (pM, []), 'getEigenvalues': (pM, []),
    'getStoichiometryMatrix': (pM, []), 'getLinkMatrix': (pM, []), 'getNrMatrix': (pM, []), 'getConservationMatrix': (pM, []),
    'getUnscaledElasticityMatrix': (pM, []), 'getScaledElasticityMatrix': (pM, []),
    'getUnscaledConcentrationControlCoefficientMatrix': (pM, []), 'getScaledConcentrationControlCoefficientMatrix': (pM, []),
    'getUnscaledFluxControlCoefficientMatrix': (pM, []), 'getScaledFluxControlCoefficientMatrix': (pM, []),
    'getuCC': (B, [S, S, pD]), 'getCC': (B, [S, S, pD]), 'getuEE': (B, [S, S, pD]), 'getEE': (B, [S, S, pD]),
    'getScaledFloatingSpeciesElasticity': (B, [S, S, pD]), 'getDiffStepSize': (B, [pD]), 'setDiffStepSize': (B, [D]),
    'getNumberOfReactions': (I, []), 'getNumberOfCompartments': (I, []), 'getNumberOfBoundarySpecies': (I, []),
    'getNumberOfFloatingSpecies': (I, []), 'getNumberOfGlobalParameters': (I, []), 'getNumberOfDependentSpecies': (I, []),
    'getNumberOfIndependentSpecies': (I, []),
    'getReactionIds': (pSA, []), 'getBoundarySpeciesIds': (pSA, []), 'getBoundarySpeciesConcentrationIds': (pSA, []),
    'getFloatingSpeciesIds': (pSA, []), 'getDependentFloatingSpeciesIds': (pSA, []), 'getIndependentFloatingSpeciesIds': (pSA, []),
    'getFloatingSpeciesConcentrationIds': (pSA, []), 'getGlobalParameterIds': (pSA, []), 'getCompartmentIds': (pSA, []),
    'getEigenvalueIds': (pSA, []), 'getRatesOfChangeIds': (pSA, []), 'getFloatingSpeciesInitialConditionIds': (pSA, []),
    'getFloatingSpeciesInitialConcentrationIds': (pSA, []),
    'getAvailableTimeCourseSymbols': (pL, []), 'getAvailableSteadyStateSymbols': (pL, []),
    'getElasticityCoefficientIds': (pL, []), 'getUnscaledElasticityCoefficientIds': (pL, []),
    'getFluxControlCoefficientIds': (pL, []), 'getUnscaledFluxControlCoefficientIds': (pL, []),
    'getConcentrationControlCoefficientIds': (pL, []), 'getUnscaledConcentrationControlCoefficientIds': (pL, []),
    'getSeed': (B, [C.POINTER(C.c_int64)]), 'setSeed': (B, [C.c_int64]),
    'gillespie': (pCD, []), 'gillespieEx': (pCD, [D, D]), 'gillespieOnGrid': (pCD, []), 'gillespieOnGridEx': (pCD, [D, D, I]),
    'gillespieMeanOnGridEx': (pCD, [D, D, I, I]), 'gillespieMeanSDOnGridEx': (pCD, [D, D, I, I]),
    'reset': (B, []), 'resetAll': (B, []), 'resetToOrigin': (B, []), 'resetParameter': (B, []),
    'setCurrentIntegrator': (I, [S]), 'getCurrentIntegratorName': (P, []), 'getNumInstantiatedIntegrators': (I, []),
    'getNumberOfCurrentIntegratorParameters': (I, []), 'getCurrentIntegratorNthParameterName': (P, [I]),
    'getCurrentIntegratorNthParameterType': (I, [I]), 'getListOfCurrentIntegratorParameterNames': (pSA, []),
    'getCurrentIntegratorParameterType': (I, [S]), 'getCurrentIntegratorParameterHint': (P, [S]),
    'getCurrentIntegratorParameterInt': (I, [S]), 'setCurrentIntegratorParameterInt': (I, [S, I]),
    'getCurrentIntegratorParameterUInt': (C.c_uint, [S]), 'setCurrentIntegratorParameterUInt': (I, [S, C.c_uint]),
    'getCurrentIntegratorParameterDouble': (D, [S]), 'setCurrentIntegratorParameterDouble': (I, [S, D]),
    'getCurrentIntegratorParameterString': (P, [S]), 'setCurrentIntegratorParameterString': (I, [S, S]),
    'getCurrentIntegratorParameterBoolean': (I, [S]), 'setCurrentIntegratorParameterBoolean': (I, [S, I]),
    'resetCurrentIntegratorParameters': (I, []), 'setCurrentIntegratorIndividualTolerance': (I, [S, D]),
    'getCurrentIntegratorParameterDoubleArray': (I, [S, C.POINTER(pD), pI]),
    'setCurrentSteadyStateSolver': (I, [S]), 'getCurrentSteadyStateSolverName': (P, []),
    'getCurrentSteadyStateSolverParameterDouble': (D, [S]), 'setCurrentSteadyStateSolverParameterDouble': (I, [S, D]),
    'getCurrentSteadyStateSolverParameterBoolean': (I, [S]), 'getListOfCurrentSteadyStateSolverParameterNames': (pSA, []),
    'getCurrentSteadyStateSolverNthParameterType': (I, [I]), 'resetCurrentSteadyStateSolverParameters': (I, []),
    'getConfigurationXML': (P, []), 'setConfigurationXML': (B, [S]), 'loadSimulationSettings': (B, [S]),
    'addParameter': (B, [S, D]), 'addParameterNoRegen': (B, [S, D]), 'addSpeciesConcentration': (B, [S, S, D, B, B]),
    'addReaction': (B, [S, C.POINTER(C.c_char_p), I, C.POINTER(C.c_char_p), I, S]),
    'addReactionNoRegen': (B, [S, C.POINTER(C.c_char_p), I, C.POINTER(C.c_char_p), I, S]),
    'setKineticLaw': (B, [S, S]), 'getKineticLaw': (P, [S]), 'removeReaction': (B, [S]), 'regenerate': (B, []),
    'addRateRuleNoRegen': (B, [S, S]), 'addEvent': (B, [S, B, S]), 'addEventAssignment': (B, [S, S, S]),
    'validateCurrentSBML': (B, []), 'getParamPromotedSBML': (P, [S]),
    '_getTime': (D, []), '_getNumIndFloatingSpecies': (I, []), '_getNumRateRules': (I, []), '_getStateVector': (I, []),
    '_getStateVectorRate': (V, [D, pD]), 'writeRRData': (I, [S]),
    'setTimeCourseSelectionListEx': (B, [I, C.POINTER(C.c_char_p)]),
}
for _n, (_r, _a) in list(PROTOS.items()) + [(k, (r, [H] + a)) for k, (r, a) in HANDLE.items()]:
    _f = getattr(L, _n)
    _f.restype, _f.argtypes = _r, _a


def err():
    return L.getLastError().decode()


def take_text(p):
    assert p, err()
    s = C.string_at(p).decode()
    L.freeText(p)
    return s


def take_vector(p):
    if not p:
        return None
    v = p.contents
    out = [v.Data[i] for i in range(v.Count)]
    L.freeVector(p)
    return np.array(out)


def take_strings(p):
    if not p:
        return []
    a = p.contents
    out = [a.String[i].decode() for i in range(a.Count)]
    L.freeStringArray(p)
    return out


def take_matrix(p):
    if not p:
        return None
    m = p.contents
    out = np.array([m.Data[i] for i in range(m.RSize * m.CSize)]).reshape(m.RSize, m.CSize)
    L.freeMatrix(p)
    return out


def take_cdata(p, weights=False):
    assert p, err()
    d = p.contents
    y = np.array([d.Data[i] for i in range(d.RSize * d.CSize)]).reshape(d.RSize, d.CSize)
    cols = [d.ColumnHeaders[i].decode() for i in range(d.CSize)]
    w = np.array([d.Weights[i] for i in range(d.RSize * d.CSize)]).reshape(d.RSize, d.CSize) if weights else None
    L.freeRRCData(p)
    return (y, cols, w) if weights else (y, cols)


def list_value(p):
    """RRList -> nested Python lists (strings / numbers)."""
    if not p:
        return None
    out = []
    for i in range(p.contents.Count):
        it = p.contents.Items[i].contents
        if it.ItemType == 0:
            out.append(C.string_at(it.data).decode())
        elif it.ItemType == 3:
            out.append(list_value(C.cast(it.data, pL)) if it.data else None)
        else:
            out.append(C.cast(C.pointer(C.c_void_p(it.data)), pD).contents.value if it.ItemType == 2 else it.data)
    return out


def take_list(p):
    v = list_value(p)
    L.freeRRList(p)
    return v


def dval(fn, *a):
    x = C.c_double()
    assert fn(*a, C.byref(x)), err()
    return x.value


def vec(xs):
    a = (C.c_double * len(xs))(*xs)
    return C.pointer(RRVector(len(xs), C.cast(a, pD))), a


def cstrs(xs):
    return (C.c_char_p * len(xs))(*[x.encode() for x in xs])


class Capi:
    """test_core_vs_roadrunner's adapter interface, on the C API only."""
    def __init__(self, sbml):
        self.h, self.sbml = L.createRRInstance(), sbml
        assert L.loadSBML(self.h, sbml.encode()), err()

    def __del__(self):
        L.freeRRInstance(self.h)

    def _ok(self, ok):
        if not ok:
            raise RuntimeError(err())

    def cma(self, on): self._ok(L.setComputeAndAssignConservationLaws(self.h, on))

    def tol(self, rtol, atol):
        assert L.setCurrentIntegratorParameterDouble(self.h, b'relative_tolerance', rtol) == 1
        assert L.setCurrentIntegratorParameterDouble(self.h, b'absolute_tolerance', atol) == 1

    def _v(self, name, *a):
        p = getattr(L, name)(self.h, *a)
        if not p and L.hasError() and 'RoadRunner exception' in err():
            raise RuntimeError(err())
        v = take_vector(p)
        return np.array([]) if v is None else v

    def _ids(self, name): return take_strings(getattr(L, name)(self.h))
    def getFloatingSpeciesIds(self): return self._ids('getFloatingSpeciesIds')
    def getBoundarySpeciesIds(self): return self._ids('getBoundarySpeciesIds')
    def getGlobalParameterIds(self): return self._ids('getGlobalParameterIds')
    def getCompartmentIds(self): return self._ids('getCompartmentIds')
    def getReactionIds(self): return self._ids('getReactionIds')
    def getIndependentFloatingSpeciesIds(self): return self._ids('getIndependentFloatingSpeciesIds')
    def getDependentFloatingSpeciesIds(self): return self._ids('getDependentFloatingSpeciesIds')
    def getEigenValueIds(self): return self._ids('getEigenvalueIds')
    def getRatesOfChangeIds(self): return self._ids('getRatesOfChangeIds')
    def cms(self): return [p for p in self.getGlobalParameterIds() if p.startswith('_CSUM')]
    def getRateRuleIds(self): return [None] * L._getNumRateRules(self.h)
    timeCourseSelections = property(lambda s: take_strings(L.getTimeCourseSelectionList(s.h)))
    steadyStateSelections = property(lambda s: take_strings(L.getSteadyStateSelectionList(s.h)))
    def getGlobalParameterValues(self): return self._v('getGlobalParameterValues')
    def fconc(self): return self._v('getFloatingSpeciesConcentrations')
    def finit(self): return self._v('getFloatingSpeciesInitialConcentrations')
    def getReactionRates(self): return self._v('getReactionRates')
    def getRatesOfChange(self): return self._v('getRatesOfChange')
    def getSteadyStateValues(self): return self._v('computeSteadyStateValues')
    def getValue(self, s): return dval(L.getValue, self.h, s.encode()) if L.getValue(self.h, s.encode(), C.byref(C.c_double())) else self._raise()
    def _raise(self): raise RuntimeError(err())
    def setValue(self, s, v): self._ok(L.setValue(self.h, s.encode(), v))
    def reset(self): self._ok(L.reset(self.h))
    def resetAll(self): self._ok(L.resetAll(self.h))
    def resetToOrigin(self): self._ok(L.resetToOrigin(self.h))
    def getCurrentTime(self): return L._getTime(self.h)
    def steadyState(self): return dval(L.steadyState, self.h) if True else None

    def sim(self, a, b, n):
        p = L.simulateEx(self.h, a, b, n)
        if not p:
            raise RuntimeError(err())
        return take_cdata(p)[0]

    # Matrix labels are not part of the C API: they come from a shadow core instance (same model and
    # moiety setting); the values are the C API's.
    def _m(self, name, labels):
        p = getattr(L, name)(self.h)
        if not p and L.hasError() and 'RoadRunner exception' in err():
            raise RuntimeError(err())
        m = take_matrix(p)
        _, rows, cols = labels()
        return (np.zeros((len(rows), len(cols))) if m is None else m), list(rows), list(cols)

    def _cma(self):
        x = C.c_int()
        L.getComputeAndAssignConservationLaws(self.h, C.byref(x))
        return x.value

    def _shadow(self):
        if getattr(self, 'sh', None) is None:
            self.sh = rn.RoadRunner(self.sbml)
        self.sh.conservedMoietyAnalysis = bool(self._cma())
        return self.sh

    def getFullStoichiometryMatrix(self): return self._m('getStoichiometryMatrix', lambda: self._shadow().getFullStoichiometryMatrix())
    def getNrMatrix(self): return self._m('getNrMatrix', lambda: self._shadow().getNrMatrix())
    def getLinkMatrix(self): return self._m('getLinkMatrix', lambda: self._shadow().getLinkMatrix())
    def getL0Matrix(self):  # no C call: the dependent rows of the link matrix
        m, rows, cols = self.getLinkMatrix()
        ni = len(cols)
        return m[ni:], rows[ni:], cols
    def getConservationMatrix(self): return self._m('getConservationMatrix', lambda: self._shadow().getConservationMatrix())
    def getFullJacobian(self): return self._m('getFullJacobian', lambda: (None, self.getFloatingSpeciesIds(), self.getFloatingSpeciesIds()))
    def getReducedJacobian(self):
        f = self.getIndependentFloatingSpeciesIds() if self._cma() else self.getFloatingSpeciesIds()
        return self._m('getReducedJacobian', lambda: (None, f, f))
    def getFullEigenValues(self):
        m = take_matrix(L.getEigenvalues(self.h))
        return np.array([]) if m is None else m[:, 0] + 1j * m[:, 1]
    def getReducedEigenValues(self):  # no C call: eigenvalues of getReducedJacobian
        m = self.getReducedJacobian()[0]
        return np.linalg.eigvals(m) if m.size else np.array([])
    def E(self, scaled):
        return self._m('getScaledElasticityMatrix' if scaled else 'getUnscaledElasticityMatrix',
                       lambda: (None, self.getReactionIds(), self.getFloatingSpeciesIds()))
    def CC(self, scaled):
        return self._m('getScaledConcentrationControlCoefficientMatrix' if scaled else 'getUnscaledConcentrationControlCoefficientMatrix',
                       lambda: (None, self.getFloatingSpeciesIds(), self.getReactionIds()))
    def FCC(self, scaled):
        r = self.getReactionIds()
        return self._m('getScaledFluxControlCoefficientMatrix' if scaled else 'getUnscaledFluxControlCoefficientMatrix', lambda: (None, r, r))
    def _2(self, name, a, b):
        x = C.c_double()
        self._ok(getattr(L, name)(self.h, a.encode(), b.encode(), C.byref(x)))
        return x.value
    def getuCC(self, a, b): return self._2('getuCC', a, b)
    def getCC(self, a, b): return self._2('getCC', a, b)
    def getuEE(self, a, b): return self._2('getuEE', a, b)
    def getEE(self, a, b): return self._2('getEE', a, b)

    def settings(self, name):
        assert L.setCurrentIntegrator(self.h, name.encode()) == 1, err()
        out = []
        for k in take_strings(L.getListOfCurrentIntegratorParameterNames(self.h)):
            t, kb = L.getCurrentIntegratorParameterType(self.h, k.encode()), k.encode()
            v = {2: lambda: bool(L.getCurrentIntegratorParameterBoolean(self.h, kb)),
                 3: lambda: L.getCurrentIntegratorParameterInt(self.h, kb),
                 6: lambda: L.getCurrentIntegratorParameterUInt(self.h, kb),
                 8: lambda: L.getCurrentIntegratorParameterDouble(self.h, kb)}[t]()
            out.append((k, v))
        return out


@pytest.mark.parametrize('name,sbml', T.MODELS, ids=T.IDS)
def test_against_libroadrunner(name, sbml):
    """Every scenario of test_core_vs_roadrunner (ids, values, selections, setValue, reset, time courses,
    N / Nr / L / L0 / Gamma with moiety analysis off and on, steady state, Jacobians, eigenvalues,
    elasticities, control coefficients, getCC / getEE, integrator settings) through the C API."""
    with np.load(os.path.join(T.REF, name + '.npz')) as z:
        ref = dict(z)
    T.run_all(Capi, sbml, name, T.Rec(ref))


# ---------------------------------------------------------------- C API surface vs the Rust core

MOI = rn.antimony_to_sbml(T.ANT['moieties'])
CHAIN = rn.antimony_to_sbml(T.ANT['chain'])


def pair(sbml):
    c = Capi(sbml)
    r = rn.RoadRunner(sbml)
    return c, r


def test_ids_counts_values_and_indexing():
    c, r = pair(CHAIN)
    h = c.h
    assert L.getNumberOfFloatingSpecies(h) == 3 and L.getNumberOfReactions(h) == 4 and L.getNumberOfBoundarySpecies(h) == 1
    assert L.getNumberOfGlobalParameters(h) == len(r.getGlobalParameterIds()) and L.getNumberOfCompartments(h) == 1
    assert c._ids('getFloatingSpeciesConcentrationIds') == ['[S1]', '[S2]', '[S3]']
    assert c._ids('getBoundarySpeciesConcentrationIds') == ['[X]']
    assert c._ids('getFloatingSpeciesInitialConditionIds') == ['init([S1])', 'init([S2])', 'init([S3])']
    assert c.getRatesOfChangeIds() == r.getRatesOfChangeIds()
    c.sim(0, 2, 5)
    r.simulate(0, 2, 5)
    for i, s in enumerate(r.getFloatingSpeciesIds()):
        assert dval(L.getFloatingSpeciesByIndex, h, i) == pytest.approx(r.getValue('[%s]' % s), rel=1e-9)
        assert dval(L.getRateOfChange, h, i) == pytest.approx(r.getValue("%s'" % s), rel=1e-9)
        assert dval(L.getFloatingSpeciesInitialConcentrationByIndex, h, i) == r.getValue('init([%s])' % s)
    for i, p in enumerate(r.getGlobalParameterIds()):
        assert dval(L.getGlobalParameterByIndex, h, i) == r.getValue(p)
    for i, j in enumerate(r.getReactionIds()):
        assert dval(L.getReactionRate, h, i) == pytest.approx(r.getValue(j), rel=1e-9)
    assert dval(L.getBoundarySpeciesByIndex, h, 0) == 2.0 and dval(L.getCompartmentByIndex, h, 0) == 1.0
    np.testing.assert_allclose(c.fconc(), r.getFloatingSpeciesConcentrations(), rtol=1e-9)
    np.testing.assert_allclose(c._v('getFloatingSpeciesAmounts'), r.getFloatingSpeciesAmounts(), rtol=1e-9)
    np.testing.assert_allclose(c._v('getBoundarySpeciesConcentrations'), [2.0])
    np.testing.assert_allclose(c._v('getBoundarySpeciesAmounts'), [2.0])
    np.testing.assert_allclose(c.getRatesOfChange(), r.getRatesOfChange(), rtol=1e-9)
    np.testing.assert_allclose(c._v('getIndependentFloatingSpeciesConcentrations'), c.fconc()[:L.getNumberOfIndependentSpecies(h)])
    assert len(c._v('getDependentFloatingSpeciesAmounts')) == L.getNumberOfDependentSpecies(h) == 0
    # setters by index
    assert L.setFloatingSpeciesByIndex(h, 1, 0.75) and dval(L.getValue, h, b'[S2]') == 0.75
    assert L.setBoundarySpeciesByIndex(h, 0, 3.0) and dval(L.getValue, h, b'[X]') == 3.0
    assert L.setGlobalParameterByIndex(h, 0, 1.25) and dval(L.getValue, h, r.getGlobalParameterIds()[0].encode()) == 1.25
    comp = c.getCompartmentIds()[0].encode()
    assert L.setCompartmentByIndex(h, 0, 2.0) and dval(L.getValue, h, comp) == 2.0
    # index errors (upstream message format)
    x = C.c_double(-7)
    assert not L.getFloatingSpeciesByIndex(h, 3, C.byref(x)) and 'Index in getFloatingSpeciesByIndex out of range: [3]' in err()
    assert not L.getFloatingSpeciesInitialConcentrationByIndex(h, 9, C.byref(x)) and x.value == -7
    assert not L.getReactionRate(h, -1, C.byref(x)) and 'out of range' in err()
    # *Ex: set concentrations, then rates
    v, keep = vec([1.0, 2.0, 3.0])
    rates = take_vector(L.getReactionRatesEx(h, v))
    r.setValue('[S1]', 1.0); r.setValue('[S2]', 2.0); r.setValue('[S3]', 3.0)
    r.setValue('X', 3.0); r.setValue(r.getGlobalParameterIds()[0], 1.25); r.setValue(comp.decode(), 2.0)
    r.setValue('[S1]', 1.0); r.setValue('[S2]', 2.0); r.setValue('[S3]', 3.0)
    np.testing.assert_allclose(rates, r.getReactionRates(), rtol=1e-12)
    np.testing.assert_allclose(take_vector(L.getRatesOfChangeEx(h, v)), r.getRatesOfChange(), rtol=1e-12)
    assert L.setFloatingSpeciesConcentrations(h, v) and L.setBoundarySpeciesConcentrations(h, vec([4.0])[0])
    assert dval(L.getValue, h, b'[X]') == 4.0
    assert L.evalModel(h)


def test_initial_conditions():
    c, r = pair(CHAIN)
    h = c.h
    v, keep = vec([0.5, 0.25])
    assert L.setFloatingSpeciesInitialConcentrations(h, v)
    np.testing.assert_allclose(c.finit()[:2], [0.5, 0.25])
    np.testing.assert_allclose(c.fconc()[:2], [0.5, 0.25])
    assert L.setFloatingSpeciesInitialConcentrationByIndex(h, 2, 0.125)
    assert dval(L.getFloatingSpeciesInitialConcentrationByIndex, h, 2) == 0.125
    assert L.reset(h)
    np.testing.assert_allclose(c.fconc(), [0.5, 0.25, 0.125])
    r.setValue('init([S1])', 0.5); r.setValue('init([S2])', 0.25); r.setValue('init([S3])', 0.125); r.reset()
    np.testing.assert_allclose(c.sim(0, 3, 7), np.asarray(r.simulate(0, 3, 7)), rtol=1e-9)
    assert not L.setFloatingSpeciesInitialConcentrations(h, None) and 'NULL vector' in err()


def test_initial_conditions_keep_parameters():
    """Upstream: model write (+ reset() for ByIndex); parameters and boundary species keep their values,
    unlike setValue('init(...)'), which reloads everything at the initial values (as libroadrunner)."""
    c = Capi(CHAIN)
    h = c.h
    c.sim(0, 2, 5)
    c.setValue('k1', 0.9); c.setValue('X', 5)
    assert L.setFloatingSpeciesInitialConcentrationByIndex(h, 0, 0.25)
    assert c.getValue('k1') == 0.9 and c.getValue('X') == 5 and L._getTime(h) == 0  # reset(): species, time
    np.testing.assert_allclose(c.fconc(), [0.25, 1, 0.5])
    c.sim(0, 2, 5)
    s3 = c.getValue('[S3]')
    v, keep = vec([0.1, 0.2])
    assert L.setFloatingSpeciesInitialConcentrations(h, v)  # changeInitialConditions: no reset at all
    assert L._getTime(h) == 2 and c.getValue('k1') == 0.9 and c.getValue('X') == 5
    np.testing.assert_allclose(c.fconc(), [0.1, 0.2, s3])
    np.testing.assert_allclose(c.finit(), [0.1, 0.2, 0.5])
    c.setValue('init([S1])', 0.3)  # setValue: libroadrunner reloads at the initial values, parameters too
    assert c.getValue('k1') == 0.3 and c.getValue('X') == 2


@pytest.mark.parametrize('rule', ["V' = 0.1", 'V := 1 + 0.1*time', 'at time > 1: V = 1.25'])
def test_initial_amounts_use_initial_volume(rule):
    """init(S) of a species in a changing compartment: the initial amount, in the initial volume
    (libroadrunner 2.11: 2, 2 at t = 2.5 with V = 1.25; set init(S) = 1 -> 1, 1)."""
    c = Capi(rn.antimony_to_sbml('compartment V = 1\n%s\nspecies S in V\nJ1: S -> ; k*S\nS = 2; k = 0.1' % rule))
    h = c.h
    c.sim(0, 2.5, 6)
    assert c.getValue('V') == pytest.approx(1.25)
    assert c.getValue('init(S)') == 2 and c.getValue('init([S])') == 2
    assert L.setFloatingSpeciesInitialConcentrationByIndex(h, 0, 1.0)
    assert c.getValue('init([S])') == 1 and c.getValue('init(S)') == 1
    c.sim(0, 2.5, 6)
    c.setValue('init(S)', 1.0)
    assert c.getValue('init(S)') == 1 and c.getValue('init([S])') == 1 and c.getValue('S') == 1


TOL = """compartment C = 2
species A in C, B in C
J1: A -> ; k*A
J2: A -> B; eps*A
J3: B -> ; kb*B
A = 1000; B = 0; k = 0.5; eps = 1e-12; kb = 2
p' = -0.1*p; p = 3"""


def test_individual_tolerance():
    """CVODEIntegrator::setIndividualTolerance: amount tolerances in libroadrunner's state order (rate
    rules, then floating species; values recorded from libroadrunner 2.11), and a tight one on the tiny
    species B makes it accurate where the scaled scalar tolerance leaves it unresolved."""
    sbml = rn.antimony_to_sbml(TOL)
    c, r = Capi(sbml), rn.RoadRunner(sbml)
    h = c.h
    c.tol(1e-6, 1e-6)
    sel = b'time, [A], [B]'
    assert L.setTimeCourseSelectionList(h, sel)
    t = np.linspace(0, 4, 9)
    exact = 1e-12 * 1000 / (2 - 0.5) * (np.exp(-0.5 * t / 2) - np.exp(-2 * t / 2))  # rates are amount / time, V = 2
    loose = take_cdata(L.simulateEx(h, 0, 4, 9))[0][:, 2]
    e_loose = np.abs(loose - exact).max() / exact.max()  # B (~1e-10) far below its 2e-6 amount tolerance
    assert L.setCurrentIntegratorIndividualTolerance(h, b'B', 1e-20) == 1, err()
    assert L.getCurrentIntegratorParameterType(h, b'absolute_tolerance') == 11  # double vector
    p, n = pD(), C.c_int()
    assert L.getCurrentIntegratorParameterDoubleArray(h, b'absolute_tolerance', C.byref(p), C.byref(n)) == 1
    vals = [p[i] for i in range(n.value)]
    C.CDLL(None).free(p)
    np.testing.assert_allclose(vals, [3e-6, 2e-3, 1e-20], rtol=1e-15)  # libroadrunner: [3e-06, 0.002, 1e-18 as set]
    L.reset(h)
    tight = take_cdata(L.simulateEx(h, 0, 4, 9))[0][:, 2]
    e_tight = np.abs(tight - exact).max() / exact.max()
    assert e_tight < 1e-5 and e_tight < e_loose / 100, (e_tight, e_loose)
    # the same through the core's PyO3 class
    r.setIntegratorSetting('absolute_tolerance', 1e-6)
    np.testing.assert_allclose(r.getAbsoluteToleranceVector(), [3e-6, 2e-3, 2e-6], rtol=1e-15)
    r.setIndividualTolerance('p', 1e-15)
    np.testing.assert_allclose(r.getIntegratorSetting('absolute_tolerance'), [1e-15, 2e-3, 2e-6], rtol=1e-15)
    with pytest.raises(ValueError):
        r.setIndividualTolerance('k', 1e-9)
    assert L.setCurrentIntegrator(h, b'rk4') == 1 and L.setCurrentIntegratorIndividualTolerance(h, b'B', 1e-9) == -1


def test_time_course_controls_and_results(tmp_path):
    c, r = pair(CHAIN)
    h = c.h
    assert L.setTimeStart(h, 1.0) and L.setTimeEnd(h, 4.0) and L.setNumPoints(h, 7)
    assert dval(L.getTimeStart, h) == 1.0 and dval(L.getTimeEnd, h) == 4.0
    n = C.c_int()
    assert L.getNumPoints(h, C.byref(n)) and n.value == 7
    assert L.setTimeStart(h, 2.0) and dval(L.getTimeEnd, h) == 5.0  # upstream keeps the duration
    assert L.setTimeCourseSelectionList(h, b'time, [S1] S2')
    assert c.timeCourseSelections == ['time', '[S1]', 'S2']
    assert L.setSteadyStateSelectionList(h, b'[S1], cc(S1, k1) ,J1')  # commas inside () stay
    assert c.steadyStateSelections == ['[S1]', 'cc(S1, k1)', 'J1']
    assert not L.setTimeCourseSelectionList(h, b'time, nope') and 'nope' in err()
    assert L.setTimeCourseSelectionListEx(h, 3, cstrs(['time', '[S1]', "S1'"]))
    y, cols = take_cdata(L.simulate(h))
    assert cols == ['time', '[S1]', "S1'"] and y.shape == (7, 3) and y[0, 0] == 2.0 and y[-1, 0] == 5.0
    r.timeCourseSelections = cols
    np.testing.assert_allclose(y, np.asarray(r.simulate(2, 5, 7)).reshape(7, 3), rtol=1e-9)
    y2, _ = take_cdata(L.getSimulationResult(h))
    assert np.array_equal(y, y2)
    m = L.getSimulationResultAsDoubleMatrix(h)  # owned by the instance: not freed here
    assert m.contents.RSize == 7 and m.contents.Data[3 * 3 + 1] == y[3, 1]
    p = str(tmp_path / 'out.csv')
    assert L.writeRRData(h, p.encode()) == 1
    assert open(p).readline().strip() == "time,[S1],S1'" and np.allclose(np.loadtxt(p, delimiter=',', skiprows=1), y)
    ts = (C.c_double * 4)(0, 0.5, 2, 3)
    y3, _ = take_cdata(L.simulateTimes(h, ts, 4))
    assert list(y3[:, 0]) == [0, 0.5, 2, 3]
    out = C.c_double()
    assert L.oneStep(h, 3.0, 0.5, C.byref(out)) and out.value == 3.5 and L._getTime(h) == 3.5
    # steady state selections
    assert L.setSteadyStateSelectionList(h, b'[S1], J1')
    assert c.steadyStateSelections == ['[S1]', 'J1']
    ss = c.getSteadyStateValues()
    r.steadyStateSelections = ['[S1]', 'J1']
    np.testing.assert_allclose(ss, r.getSteadyStateValues(), rtol=1e-8)
    # settings file
    f = tmp_path / 'settings.txt'
    f.write_text('start: 0\nduration: 10\nsteps: 20\nvariables: S1\nabsolute: 1e-9\nrelative: 1e-7\n')
    assert L.loadSimulationSettings(h, str(f).encode())
    assert dval(L.getTimeEnd, h) == 10 and L.getNumPoints(h, C.byref(n)) and n.value == 21
    assert L.getCurrentIntegratorParameterDouble(h, b'relative_tolerance') == 1e-7
    assert not L.loadSimulationSettings(h, b'/nonexistent/x.txt') and err() == 'The file /nonexistent/x.txt was not found'


def test_rk45_variable_steps_and_one_step():
    c, r = pair(CHAIN)
    h = c.h
    assert L.setCurrentIntegrator(h, b'rk45') == 1
    assert L.getCurrentIntegratorParameterBoolean(h, b'variable_step_size') == 1
    y, _ = take_cdata(L.simulateEx(h, 0, 5, 11))  # libroadrunner's C API: variable steps despite 11 points
    assert y.shape[0] > 11 and y[-1, 0] == 5 and np.all(np.diff(y[:, 0]) > 0)
    r.setIntegrator('rk45')
    np.testing.assert_allclose(y, np.asarray(r.simulate(0, 5, 11)).reshape(y.shape), rtol=1e-12)
    for name in ('rk4', 'euler', 'gillespie'):
        assert L.setCurrentIntegrator(h, name.encode()) == 1
        out = C.c_double()
        assert L.oneStep(h, 5.0, 0.25, C.byref(out)), (name, err())
        assert out.value == 5.25 if name != 'gillespie' else 5.0 < out.value <= 5.25  # variable steps: one firing
        assert L._getTime(h) == out.value
    assert L.setCurrentIntegratorParameterBoolean(h, b'variable_step_size', 0) == 1
    assert L.oneStep(h, 6.0, 0.25, C.byref(out)) and out.value == 6.25
    assert L.getNumInstantiatedIntegrators(h) == 5


def test_matrices_and_mca_with_moieties():
    c, r = pair(MOI)
    h = c.h
    for x in (c, r):
        x.cma(True) if x is c else setattr(r, 'conservedMoietyAnalysis', True)
    assert c._cma() == 1
    assert L.getNumberOfDependentSpecies(h) == 2 and L.getNumberOfIndependentSpecies(h) == 3
    for name, core in [('getStoichiometryMatrix', r.getFullStoichiometryMatrix), ('getLinkMatrix', r.getLinkMatrix),
                       ('getNrMatrix', r.getNrMatrix), ('getConservationMatrix', r.getConservationMatrix)]:
        np.testing.assert_array_equal(take_matrix(getattr(L, name)(h)), np.asarray(core()[0]), err_msg=name)
    assert L.steadyState(h, C.byref(C.c_double()))
    r.steadyState()
    for name, core in [('getFullJacobian', r.getFullJacobian), ('getReducedJacobian', r.getReducedJacobian),
                       ('getUnscaledElasticityMatrix', lambda: r.getElasticityMatrix(False)),
                       ('getScaledElasticityMatrix', lambda: r.getElasticityMatrix(True)),
                       ('getUnscaledConcentrationControlCoefficientMatrix', lambda: r.getConcentrationControlCoefficientMatrix(False)),
                       ('getScaledConcentrationControlCoefficientMatrix', lambda: r.getConcentrationControlCoefficientMatrix(True)),
                       ('getUnscaledFluxControlCoefficientMatrix', lambda: r.getFluxControlCoefficientMatrix(False)),
                       ('getScaledFluxControlCoefficientMatrix', lambda: r.getFluxControlCoefficientMatrix(True))]:
        np.testing.assert_allclose(take_matrix(getattr(L, name)(h)), np.asarray(core()[0]), rtol=1e-9, atol=1e-12, err_msg=name)
    ev = take_matrix(L.getEigenvalues(h))
    core_ev = np.asarray(r.getFullEigenValues())
    assert ev.shape == (5, 2)
    np.testing.assert_allclose(np.sort_complex(ev[:, 0] + 1j * ev[:, 1]), np.sort_complex(core_ev), atol=1e-10)
    assert c._ids('getEigenvalueIds') == r.getEigenValueIds()
    assert dval(L.getCC, h, b'S1', b'k1') == pytest.approx(r.getCC('S1', 'k1'), rel=1e-9)
    assert dval(L.getuCC, h, b'J1', b'k2') == pytest.approx(r.getuCC('J1', 'k2'), rel=1e-9)
    assert dval(L.getEE, h, b'J1', b'S1') == pytest.approx(r.getEE('J1', 'S1'), rel=1e-12)
    assert dval(L.getuEE, h, b'J4', b'k4') == pytest.approx(r.getuEE('J4', 'k4'), rel=1e-12)
    assert dval(L.getScaledFloatingSpeciesElasticity, h, b'J1', b'S1') == pytest.approx(r.getEE('J1', 'S1', False), rel=1e-12)
    assert not L.getScaledFloatingSpeciesElasticity(h, b'J1', b'k1', C.byref(C.c_double())) and 'invalid species name: k1' in err()
    assert L.setDiffStepSize(h, 0.01) and dval(L.getDiffStepSize, h) == 0.01
    # symbol lists: {label, {items}} groups; nested coefficient groups
    tc = take_list(L.getAvailableTimeCourseSymbols(h))
    assert [g[0] for g in tc] == ['Time', 'Floating Species', 'Boundary Species', 'Floating Species (amount)', 'Boundary Species (amount)',
                                  'Global Parameters', 'Fluxes', 'Rates of Change', 'Volumes', 'Elasticity Coefficients',
                                  'Unscaled Elasticity Coefficients', 'Eigenvalues']
    assert tc[0] == ['Time', ['time']] and tc[2] == ['Boundary Species', None]  # empty group: NULL list (as upstream)
    assert tc[6][1] == r.getReactionIds() and tc[9][1][0][0] == 'J1' and tc[9][1][0][1][0] == 'ec(J1, %s)' % r.getFloatingSpeciesIds()[0]
    st = take_list(L.getAvailableSteadyStateSymbols(h))
    assert [g[0] for g in st][5:10] == ['Volumes', 'Fluxes', 'Flux Control Coefficients', 'Concentration Control Coefficients',
                                        'Unscaled Concentration Control Coefficients']
    fcc = take_list(L.getFluxControlCoefficientIds(h))
    assert fcc[0][0] == 'J1' and fcc[0][1] == ['cc(J1, %s)' % p for p in ['k1', 'k2', 'k3', 'k4', 'k5', '_CSUM0', '_CSUM1']]
    ucc = take_list(L.getUnscaledConcentrationControlCoefficientIds(h))
    assert ucc[0][1][0].startswith('ucc(') and len(ucc) == 5
    assert take_list(L.getUnscaledElasticityCoefficientIds(h))[0][1][0] == 'uec(J1, %s)' % r.getFloatingSpeciesIds()[0]


def test_integrator_and_solver_settings():
    c = Capi(CHAIN)
    h = c.h
    names = [take_text(L.getRegisteredIntegratorName(i)) for i in range(L.getNumRegisteredIntegrators())]
    assert names == ['cvode', 'gillespie', 'rk4', 'rk45', 'euler']
    assert take_text(L.getRegisteredIntegratorHint(3)) == 'Internal RK45 ODE solver'
    assert not L.getRegisteredIntegratorName(5) and L.hasError()
    assert [take_text(L.getRegisteredSteadyStateSolverName(i)) for i in range(L.getNumRegisteredSteadyStateSolvers())] == \
        ['nleq1', 'nleq2', 'newton', 'newton_linesearch']
    assert take_text(L.getCurrentIntegratorName(h)) == 'cvode' and L.getNumInstantiatedIntegrators(h) == 1
    keys = take_strings(L.getListOfCurrentIntegratorParameterNames(h))
    assert keys[:2] == ['relative_tolerance', 'absolute_tolerance'] and L.getNumberOfCurrentIntegratorParameters(h) == len(keys)
    assert take_text(L.getCurrentIntegratorNthParameterName(h, 0)) == 'relative_tolerance'
    assert L.getCurrentIntegratorNthParameterType(h, 0) == 8 and L.getCurrentIntegratorParameterType(h, b'stiff') == 2
    assert L.getCurrentIntegratorParameterType(h, b'maximum_num_steps') == 3
    assert L.setCurrentIntegratorParameterInt(h, b'maximum_num_steps', 5000) == 1
    assert L.getCurrentIntegratorParameterInt(h, b'maximum_num_steps') == 5000
    assert take_text(L.getCurrentIntegratorParameterString(h, b'maximum_num_steps')) == '5000'
    assert L.setCurrentIntegratorParameterBoolean(h, b'stiff', 0) == 1 and L.getCurrentIntegratorParameterBoolean(h, b'stiff') == 0
    assert L.setCurrentIntegratorParameterDouble(h, b'relative_tolerance', 1e-9) == 1
    assert L.resetCurrentIntegratorParameters(h) == 1 and L.getCurrentIntegratorParameterDouble(h, b'relative_tolerance') == 1e-6
    assert L.setCurrentIntegratorParameterDouble(h, b'nope', 1.0) == -1 and "no setting 'nope'" in err()
    assert L.getCurrentIntegratorParameterDouble(h, b'nope') == -1.0
    assert L.getCurrentIntegratorParameterType(h, b'variable_step_size') == 2
    assert L.setCurrentIntegratorParameterBoolean(h, b'variable_step_size', 1) == 1
    assert L.getCurrentIntegratorParameterBoolean(h, b'variable_step_size') == 1
    for key in (b'initial_time_step', b'minimum_time_step', b'maximum_time_step'):
        assert L.getCurrentIntegratorParameterType(h, key) == 8
        assert L.setCurrentIntegratorParameterDouble(h, key, 0.02) == 1
        assert L.getCurrentIntegratorParameterDouble(h, key) == 0.02
    assert L.setCurrentIntegratorIndividualTolerance(h, b'k1', 1e-9) == -1 and 'neither a floating species nor has a rate rule' in err()
    assert L.setCurrentIntegrator(h, b'lsoda') == -1 and L.hasError()
    assert L.setCurrentIntegrator(h, b'gillespie') == 1
    assert L.setCurrentIntegratorParameterUInt(h, b'seed', 1234) == 1 and L.getCurrentIntegratorParameterUInt(h, b'seed') == 1234
    s = C.c_int64()
    assert L.getSeed(h, C.byref(s)) and s.value == 1234 and L.setSeed(h, 99) and L.getSeed(h, C.byref(s)) and s.value == 99
    assert L.setCurrentSteadyStateSolver(h, b'newton') == 1 and take_text(L.getCurrentSteadyStateSolverName(h)) == 'newton'
    assert L.getCurrentSteadyStateSolverParameterBoolean(h, b'allow_presimulation') == 1
    assert L.setCurrentSteadyStateSolverParameterDouble(h, b'presimulation_time', 10) == 1
    assert L.getCurrentSteadyStateSolverParameterDouble(h, b'presimulation_time') == 10
    assert 'presimulation_times' in take_strings(L.getListOfCurrentSteadyStateSolverParameterNames(h))
    assert L.resetCurrentSteadyStateSolverParameters(h) == 1 and L.getCurrentSteadyStateSolverParameterDouble(h, b'presimulation_time') == 5
    # configuration XML round trip
    L.setCurrentIntegrator(h, b'cvode')
    L.setCurrentIntegratorParameterDouble(h, b'relative_tolerance', 1e-9)
    x = take_text(L.getConfigurationXML(h))
    assert '<section name="Integration" method="cvode"' in x and 'name="relative_tolerance" value="1e-09"' in x
    c2 = Capi(CHAIN)
    assert L.setConfigurationXML(c2.h, x.encode()), err()
    assert L.getCurrentIntegratorParameterDouble(c2.h, b'relative_tolerance') == 1e-9
    assert take_text(L.getCurrentSteadyStateSolverName(c2.h)) == 'newton'


def test_gillespie_variants():
    sbml = rn.antimony_to_sbml('J1: S1 -> S2; k1*S1; S1 = 100; S2 = 0; k1 = 0.5')
    c, r = pair(sbml)
    h = c.h
    assert L.setSeed(h, 7)
    y1, cols = take_cdata(L.gillespieOnGridEx(h, 0, 5, 11))
    assert cols == ['time', '[S1]', '[S2]'] and y1.shape == (11, 3) and np.all(y1[:, 1] + y1[:, 2] == 100)
    assert take_text(L.getCurrentIntegratorName(h)) == 'gillespie'  # stays current (as upstream)
    r.setIntegrator('gillespie'); r.setIntegratorSetting('variable_step_size', False); r.setSeed(7)
    assert np.array_equal(y1, np.asarray(r.simulate(0, 5, 11)))  # same seed stream as the core
    assert L.reset(h) and L.setSeed(h, 7)
    assert np.array_equal(y1, take_cdata(L.gillespieOnGrid(h))[0])
    L.reset(h)
    yv, _ = take_cdata(L.gillespieEx(h, 0, 5))  # variable steps: a row per firing
    assert yv.shape[0] > 50 and yv[-1, 0] == 5 and np.all(np.diff(yv[:, 0]) >= 0)
    L.reset(h)
    assert take_cdata(L.gillespie(h))[0].shape[0] > 50
    m, _ = take_cdata(L.gillespieMeanOnGridEx(h, 0, 5, 6, 3000))
    np.testing.assert_allclose(m[:, 1], 100 * np.exp(-0.5 * m[:, 0]), rtol=0.02)
    m2, _, sd = take_cdata(L.gillespieMeanSDOnGridEx(h, 0, 5, 6, 3000), weights=True)
    p = np.exp(-0.5 * m2[:, 0])
    np.testing.assert_allclose(sd[1:, 1], np.sqrt(100 * p * (1 - p))[1:], rtol=0.06)
    m0, _ = take_cdata(L.gillespieMeanOnGridEx(h, 0, 5, 6, 0))
    assert np.all(m0 == 0)


def test_reset_family_and_state(tmp_path):
    c, r = pair(CHAIN)
    h = c.h
    c.sim(0, 2, 5)
    c.setValue('k1', 0.9)
    assert L.reset(h)
    assert dval(L.getValue, h, b'k1') == 0.9 and L._getTime(h) == 0 and dval(L.getValue, h, b'S1') == 0
    assert L.resetParameter(h) and dval(L.getValue, h, b'k1') == 0.3
    c.setValue('k1', 0.9); c.setValue('X', 5)
    assert L.resetAll(h) and dval(L.getValue, h, b'k1') == 0.3 and dval(L.getValue, h, b'X') == 2
    c.setValue('X', 5)
    assert L.resetToOrigin(h) and dval(L.getValue, h, b'X') == 2
    c.sim(0, 2, 3)
    c.setValue('k1', 0.7)
    p = str(tmp_path / 'state.txt')
    assert L.saveState(h, p.encode())
    c2 = Capi(CHAIN)
    assert L.loadState(c2.h, p.encode()), err()
    assert L._getTime(c2.h) == 2 and dval(L.getValue, c2.h, b'k1') == 0.7
    np.testing.assert_allclose(c2.sim(2, 4, 5), c.sim(2, 4, 5), rtol=1e-9)
    assert not L.loadState(c2.h, b'/nonexistent/state') and L.hasError()


def test_editing_and_sbml():
    c, r = pair(CHAIN)
    h = c.h
    assert L.addParameter(h, b'kx', 0.5) and dval(L.getValue, h, b'kx') == 0.5
    assert L.addParameterNoRegen(h, b'ky', 0.25)
    assert not L.getValue(h, b'ky', C.byref(C.c_double()))  # not compiled until regenerate()
    assert L.addRateRuleNoRegen(h, b'ky', b'-ky') and L.regenerate(h)
    assert dval(L.getValue, h, b'ky') == 0.25 and L._getNumRateRules(h) == 1
    assert L.addReaction(h, b'J9', cstrs(['S3']), 1, cstrs(['S1']), 1, b'kx*S3')
    assert take_text(L.getKineticLaw(h, b'J9')) == 'kx * S3'
    assert L.setKineticLaw(h, b'J9', b'kx*S3*S3') and 'S3' in take_text(L.getKineticLaw(h, b'J9'))
    assert L.addReactionNoRegen(h, b'J10', cstrs(['S2']), 1, cstrs([]), 0, b'0.1*S2') and L.regenerate(h)
    assert L.getNumberOfReactions(h) == 6
    assert L.removeReaction(h, b'J10') and L.getNumberOfReactions(h) == 5
    comp = c.getCompartmentIds()[0]
    assert L.addSpeciesConcentration(h, b'Z', comp.encode(), 1.5, False, False) and dval(L.getValue, h, b'[Z]') == 1.5
    assert L.addEvent(h, b'E1', True, b'time > 1') and L.addEventAssignment(h, b'E1', b'Z', b'7')
    assert L.validateCurrentSBML(h), err()
    sb = take_text(L.getSBML(h))
    assert 'id="J9"' in sb and 'id="E1"' in sb
    # same edits on the core give the same model
    for args in [('add_parameter', 'kx', 0.5), ('add_parameter', 'ky', 0.25), ('add_rate_rule', 'ky', '-ky'),
                 ('add_reaction', 'J9', ['S3'], ['S1'], 'kx*S3'), ('set_kinetic_law', 'J9', 'kx*S3*S3'),
                 ('add_species_concentration', 'Z', comp, 1.5, False, False, ''), ('add_event', 'E1', True, 'time > 1'),
                 ('add_event_assignment', 'E1', 'Z', '7')]:
        r.load(rn.sbml_edit(args[0], r.getSBML(), *args[1:]))
    L.resetToOrigin(h)
    r.timeCourseSelections = c.timeCourseSelections
    np.testing.assert_allclose(c.sim(0, 3, 7), np.asarray(r.simulate(0, 3, 7)), rtol=1e-9)
    cur = take_text(L.getCurrentSBML(h))
    assert 'id="Z"' in cur
    pp = take_text(L.getParamPromotedSBML(h, take_text(L.getSBML(h)).encode()))
    assert '<sbml' in pp
    assert not L.addParameter(h, b'kx', 1.0) and L.hasError()  # duplicate id
    assert not L.removeReaction(h, b'nope') and L.hasError()


def test_selections_rate_rules_and_plugin_helpers():
    sbml = rn.antimony_to_sbml(T.ANT['raterule'])
    c, r = pair(sbml)
    h = c.h
    assert L._getNumRateRules(h) == 1 and L._getNumIndFloatingSpecies(h) == 2
    np.testing.assert_allclose(c._v('getRateRuleValues'), [1.0])
    n = L._getStateVector(h)
    assert n == len(r.getRatesOfChangeIds()) == 3
    buf = (C.c_double * n)()
    L._getStateVectorRate(h, 0.0, buf)
    np.testing.assert_allclose(list(buf), r.getRatesOfChange(), rtol=1e-12)
    assert c.getRatesOfChangeIds() == ["q'", "S1'", "S2'"]
    np.testing.assert_allclose(c.getRatesOfChange(), r.getRatesOfChange(), rtol=1e-12)  # aligned with the ids


def test_error_paths_and_ownership():
    assert not L.loadSBML(None, b'<sbml/>') and 'Failed to cast to a valid RoadRunner handle' in err()
    h = L.createRRInstance()
    try:
        assert not L.isModelLoaded(h) and not L.clearModel(h)
        assert L.getNumberOfFloatingSpecies(h) == -1 and 'A model needs to be loaded' in err()
        assert not L.simulate(h) and not L.getFloatingSpeciesConcentrations(h)
        assert not L.loadSBML(h, b'<sbml>broken') and err().startswith('RoadRunner exception: ')
        assert not L.loadSBMLFromFile(h, b'/nonexistent/m.xml') and err() == 'The file /nonexistent/m.xml was not found'
        assert L.loadSBML(h, CHAIN.encode()) and L.isModelLoaded(h)
        assert not L.getValue(h, b'nope', C.byref(C.c_double())) and L.hasError()
        assert not L.setValue(h, b'J1', 1.0) and L.hasError()
        assert not L.getValue(h, b'S1', None) and 'NULL output' in err()
        assert not L.getConservationMatrix(h) and L.hasError()  # moiety analysis off
        assert not L.oneStep(h, 0, 1, None)
        assert L.clearModel(h) and not L.isModelLoaded(h)
    finally:
        assert L.freeRRInstance(h)
    assert not L.freeRRInstance(None)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q']))

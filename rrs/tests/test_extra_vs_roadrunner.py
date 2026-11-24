"""crate::rr::extra against libroadrunner's recorded answers (keys 'extra/...' in reference/<model>.npz, made
by make_reference.py --extra): K / extended / reactant / product stoichiometry matrices, stoichiometry ids
and values, setStoichiometry, getIds by SelectionRecord type, frequency response, internalOneStep (rk45),
piecewise triggers, events. SBML conversion, root evaluation and exact frequency
response are also checked independently.

    PYTHONPATH=<built roadrunner_nitro>:. python -m pytest rrs/tests/test_extra_vs_roadrunner.py
"""
import os
from xml.etree import ElementTree as ET
import numpy as np
import pytest
import roadrunner_nitro as rn
from test_core_vs_roadrunner import MODELS, IDS, REF, Rec, Sub, ANT, cmp

ALL, BIG = -1, (787492, 787584, 65536, 786980, 787520)
SMALL = BIG + (-262145, -524289, -262659, 787072, 8192, 16384, 6144, 262184)


class Core:
    def __init__(self, sbml): self.r = rn.RoadRunner(sbml)
    def cma(self, on): self.r.conservedMoietyAnalysis = on
    def __getattr__(self, n): return getattr(self.r, n)
    def sids(self): return self.r.getStoichiometryIds()
    def stoich(self, i): return self.r.getStoichiometry(i)
    def ids(self, t, model): return self.r.getIds(t, not model)
    def npw(self): return self.r.getNumPiecewiseTriggers()
    def events(self): return self.r.getEventIds()
    def sim(self, a, b, n): return np.asarray(self.r.simulate(a, b, n))


def sc_extra(r, rec, name, cma):
    for n in ('getKMatrix', 'getExtendedStoichiometryMatrix', 'getReactantsStoichiometryMatrix', 'getProductsStoichiometryMatrix'):
        result = rec(n, lambda: getattr(r, n)(), mirror=('',))
        if rec.rec.ref is not None and result is not None:
            _, rows, columns = result
            key = rec.prefix + n
            assert rows == rec.rec.ref[key + '/rownames'].tolist()
            assert columns == rec.rec.ref[key + '/colnames'].tolist()
    ids = rec('sids', lambda: list(r.sids()))
    rec('svals', lambda: [r.stoich(i) for i in range(len(ids) if ids is not None else 0)])
    small = len(r.getFloatingSpeciesIds()) <= 20
    for t in SMALL if small else BIG:
        for model in (False, True):
            rec('ids%d:%d' % (model, t), lambda: list(r.ids(t, model)))
    if not cma:
        rec('npw', r.npw)
        rec('events', lambda: list(r.events()))
    if name == 'chain':
        fs = r.getFloatingSpeciesIds()
        for p in ('k0', 'k1'):
            for db, hz in ((False, False), (True, True)):
                rec('fr:%s:%d' % (p, db), lambda: r.getFrequencyResponse(0.01, 3, 7, p, fs[-1], db, hz), 1e-5, 1e-9, scaled=True, mirror=('',))


def sc_setstoich(r, rec, name):
    if not r.sids():
        return
    r.setStoichiometry(0, 2.5)
    rec('val', lambda: r.stoich(0))
    # The reference uses CVODE's default tolerances on the time-varying coefficient.
    rec('sim', lambda: r.sim(0, 2, 5), 2e-5 if name == 'raterule' else 1e-5, 1e-9, scaled=True)


def run_extra(make, sbml, name, rec):
    for cma in (False, True):
        r = make(sbml)
        sr = Sub(rec, 'extra/%d/' % cma)
        if sr('cma', lambda: r.cma(cma) or 0.0, mirror=('',)) is None:
            continue
        sc_extra(r, sr, name, cma)
    if name in ('chain', 'moieties', 'events', 'raterule'):
        sc_setstoich(make(sbml), Sub(rec, 'extra/setstoich/'), name)


@pytest.mark.parametrize('name,sbml', MODELS, ids=IDS)
def test_extra_against_libroadrunner(name, sbml):
    with np.load(os.path.join(REF, name + '.npz')) as z:
        ref = dict(z)
    run_extra(Core, sbml, name, Rec(ref))


@pytest.mark.parametrize('name,sbml', MODELS, ids=IDS)
def test_sbml_level_version_round_trip(name, sbml):
    a = rn.RoadRunner(sbml)
    a.simulate(0, 1, 3)
    a.setValue(a.getGlobalParameterIds()[0], 0.37) if a.getGlobalParameterIds() else None
    for cur in (False, True):
        for lv in ((2, 1), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2)):
            method = a.getCurrentSBML if cur else a.getSBML
            try:
                text = method(*lv)
            except ValueError as e:  # refused: content the target can't express
                assert name != 'chain' and 'converting' in str(e)
                continue
            root = ET.fromstring(text)
            assert (int(root.get('level')), int(root.get('version'))) == lv
            b, c = rn.RoadRunner(text), rn.RoadRunner(method())
            for f in ('getFloatingSpeciesIds', 'getBoundarySpeciesIds', 'getGlobalParameterIds', 'getReactionIds', 'getCompartmentIds'):
                assert getattr(b, f)() == getattr(c, f)(), f
            cmp(c.getGlobalParameterValues(), b.getGlobalParameterValues(), 1e-12, 1e-15, False)
            cmp(c.simulate(0, 2, 5), b.simulate(0, 2, 5), 1e-6, 1e-10, True)
    with pytest.raises(ValueError):
        a.getSBML(4, 1)


def test_internal_one_step_and_more():
    r = rn.RoadRunner(rn.antimony_to_sbml('J0: S1 -> S2; k1*S1; S1 = 10; k1 = 0.3'))
    t = r.internalOneStep(0, 1.0)
    assert 0 < t <= 1 and r.getCurrentTime() == t
    np.testing.assert_allclose(r.getValue('S1'), 10 * np.exp(-0.3 * t), rtol=1e-6)
    next_time = r.internalOneStep(t, 1.0, reset=False)
    assert next_time - t > t
    np.testing.assert_allclose(r.getValue('S1'), 10 * np.exp(-0.3 * next_time), rtol=1e-8)
    first_time = next_time
    r.setIntegrator('rk45')  # libroadrunner: 0.010281535624713385, S1 9.96920291374347
    t = r.internalOneStep(0, 1.0)
    r2 = rn.RoadRunner(rn.antimony_to_sbml('J0: S1 -> S2; k1*S1; S1 = 10; k1 = 0.3'))
    r2.setIntegrator('rk45')
    np.testing.assert_allclose(r2.internalOneStep(0, 1.0), 0.010281535624713385, rtol=1e-12)
    assert 0 < t < 1
    np.testing.assert_allclose(r.getValue('S1'), 10 * np.exp(-0.3 * (t + first_time)), rtol=1e-6)
    np.testing.assert_allclose(r2.getValue('S1'), 9.96920291374347, rtol=1e-12)
    with np.load(os.path.join(REF, 'chain.npz')) as ref:
        supported = ref['extra/0/internal_step_support']
        for n, accepts in zip(('cvode', 'rk45', 'gillespie', 'rk4', 'euler'), supported):
            r2.setIntegrator(n)
            if accepts:
                assert r2.internalOneStep(0, 1.0) > 0
            else:
                with pytest.raises(ValueError, match='no internal'):
                    r2.internalOneStep(0, 1.0)
    # frequency response, exactly: H = k0 X / (iw + k1) for S1 of the chain's first step
    c = rn.RoadRunner(rn.antimony_to_sbml(ANT['chain']))
    f = c.getFrequencyResponse(0.1, 2, 5, 'k0', 'S1', False, False)
    w = np.logspace(-1, 1, 5)
    h = 2 / (1j * w + 0.3)
    np.testing.assert_allclose(f, np.c_[w, abs(h), np.degrees(np.angle(h)) % 360], rtol=1e-12)
    with pytest.raises(ValueError):
        c.getFrequencyResponse(0.1, 1, 2, 'k0', 'J1', False, False)  # libroadrunner: zeros
    e = rn.RoadRunner(rn.antimony_to_sbml(ANT['events']))
    assert e.getEventRoots() == [-1.0] and e.getEventRoots(3.0) == [1.0]
    assert 0 < e.getRandom() < 1



def test_facade_matrices_and_selection_types():
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner, SelectionRecord

    r = RoadRunner(rn.antimony_to_sbml(ANT['chain']))
    with np.load(os.path.join(REF, 'chain.npz')) as ref:
        prefix = 'extra/0/'
        for name in ('getKMatrix', 'getExtendedStoichiometryMatrix',
                     'getReactantsStoichiometryMatrix', 'getProductsStoichiometryMatrix'):
            matrix = getattr(r, name)()
            key = prefix + name
            assert matrix.rownames == ref[key + '/rownames'].tolist()
            assert matrix.colnames == ref[key + '/colnames'].tolist()
            np.testing.assert_allclose(matrix, ref[key])
        extended = r.getExtendedStoichiometryMatrix()
        reactants = r.getReactantsStoichiometryMatrix(boundary=True)
        products = r.getProductsStoichiometryMatrix(boundary=True)
        assert reactants.rownames == products.rownames == extended.rownames
        np.testing.assert_allclose(extended, np.asarray(reactants) + np.asarray(products))
        assert extended['X'][0] == -1
        assert r.getSupportedIdTypes() == ref[prefix + 'supported']
        assert r.model.getSupportedIdTypes() == ref[prefix + 'model_supported']
        assert r.getIds() == ref[prefix + 'ids_all'].tolist()
        t = r.internalOneStep(currentTime=0, stepSize=0.1)
        assert 0 < t <= 0.1 and r.getCurrentTime() == t
        response = r.getFrequencyResponse(0.1, 2, 5, 'k0', 'S1', False, False)
        assert not hasattr(response, 'colnames') and ref[prefix + 'frequency_columns'].size == 0
        assert r.getIds(SelectionRecord.FLOATING_AMOUNT) == ['S1', 'S2', 'S3']
        for key, value in ref[prefix + 'selection_flags']:
            assert getattr(SelectionRecord, key) == int(value), key


def test_facade_stoichiometry_and_settings():
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner

    r = RoadRunner(rn.antimony_to_sbml(ANT['chain']))
    m = r.model
    with np.load(os.path.join(REF, 'chain.npz')) as ref:
        prefix = 'extra/0/'
        ids = m.getStoichiometryIds()
        assert ids == ref[prefix + 'sids'].tolist()
        assert [m.getStoichiometryId(i) for i in range(len(ids))] == ids
        assert [m.getStoichiometryIndex(s) for s in [*ids, 'missing']] == ref[prefix + 'stoich_index'].tolist()
        sp, rx = r.getFloatingSpeciesIds(), r.getReactionIds()
        assert [[m.getStoichiometryIndex(s, x) for x in rx] for s in sp] == ref[prefix + 'stoich_pair_index'].tolist()
        np.testing.assert_allclose([m.getInitStoichiometry(i) for i in range(len(ids))], ref[prefix + 'initial_stoich'])
        for getter, key in ((m.getStoichiometry, 'stoich_pair'), (m.getInitStoichiometry, 'initial_stoich_pair')):
            values = [[getter(i, j) for j in range(len(rx))] for i in range(len(sp))]
            np.testing.assert_allclose(values, ref[prefix + key])
        assert np.isnan(m.getInitStoichiometry(-1))
        assert np.isnan(m.getInitStoichiometry(0, 100))
        original_sbml = r.getSBML()
        result = [m.setStoichiometry(0, 2.5), m.getStoichiometry(0), m.getInitStoichiometry(0)]
        assert r.getSBML() == original_sbml
        assert RoadRunner(r.getCurrentSBML()).model.getStoichiometry(0) == 2.5
        for reset in (r.reset, r.resetAll, r.resetToOrigin):
            reset()
            result.append(m.getStoichiometry(0))
        np.testing.assert_allclose(result, ref[prefix + 'set_stoich_reset'])
        assert m.setStoichiometry(0, 1, -2.5) == 1
        assert m.getStoichiometry(0, 1) == -2.5
        assert m.getInitStoichiometry(0, 1) == -1
        assert m.setStoichiometry(0, 3, 1.0) == 0
        assert r.getSteadyStateThreshold() == ref[prefix + 'threshold']
        r.setSteadyStateThreshold(0.0001)
        assert r.steadyStateThresh == 0.0001
        r.steadyStateThresh = 0.001
        assert r.getSteadyStateThreshold() == 0.001
        r.setSeed(42)
        first = m.getRandom()
        r.setSeed(42)
        assert m.getRandom() == first
        r.simulate(0, 1, 2)
        r.resetSeed()
        assert r.getSeed() != 42
        np.testing.assert_allclose([r.getCurrentTime(), r['S1']], ref[prefix + 'reset_seed_state'])


def test_exact_frequency_response_with_dynamic_stoichiometry():
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner

    sbml = rn.antimony_to_sbml("J0: -> S; p; J1: S -> ; k*S; S = 6; p = 3; k = 0.5; a := p/2")
    sbml = sbml.replace('<speciesReference species="S" stoichiometry="1" constant="true"/>',
                        '<speciesReference id="a" species="S" constant="false"/>', 1)
    # Move the assignment rule from a global parameter onto the named species reference.
    import re
    sbml = re.sub(r'^.*<parameter id="a"[^>]*' + chr(10), '', sbml, flags=re.MULTILINE)
    r = RoadRunner(sbml)
    model = r.model
    assert model.getStoichiometryId(0) == 'a'
    assert model.getStoichiometryIndex('a') == 0
    assert model.getInitStoichiometry(0) == 1.5
    with pytest.raises(ValueError, match='assignment rule'):
        model.setStoichiometry(0, 2)
    r['p'] = 4
    assert model.getStoichiometry(0) == 2
    # dS/dt = p^2/2 - k*S, hence d f/d p = p at the current steady state.
    omega = np.logspace(-1, 1, 5)
    h = 4 / (1j * omega + 0.5)
    result = r.getFrequencyResponse(startFrequency=0.1, numberOfDecades=2, numberOfPoints=5,
                                    parameterName='p', variableName='S', useDB=True, useHz=True)
    np.testing.assert_allclose(result, np.c_[omega / (2 * np.pi), 20 * np.log10(np.abs(h)),
                                             np.degrees(np.angle(h)) % 360], rtol=1e-12)


def test_facade_event_and_piecewise_roots():
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner

    r = RoadRunner(rn.antimony_to_sbml("""
        compartment C = 3
        S in C; S = 1
        k := piecewise(2, time < 1, 3, S > 2, 4)
        J: S -> ; k*S
        E: at time >= 2 && S > 2: S = 0
    """))
    assert r.getNumEvents() == 1
    assert r.getEventIds() == ['E']
    assert r.model.getEventIndex('E') == 0 and r.model.getEventIndex('missing') == -1
    assert r.getNumPiecewiseTriggers() == 3
    np.testing.assert_array_equal(r.getEventRoots(2), [-1])
    gdot = np.zeros(1)
    assert r.model.getEventRoots(2, np.array([9.0]), gdot) is None
    np.testing.assert_array_equal(gdot, [1])
    np.testing.assert_array_equal(r.getPiecewiseTriggerRoots(0, [9]), [1, -1, -1])
    np.testing.assert_array_equal(r.getPiecewiseTriggerRoots(2, [9]), [-1, 1, -1])
    np.testing.assert_array_equal(r.getPiecewiseTriggerRoots(2, [3]), [-1, -1, 1])
    assert r['S'] == 3 and r.getCurrentTime() == 0


@pytest.mark.parametrize('level,version', [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2)])
def test_facade_sbml_conversion(level, version):
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner

    r = RoadRunner(rn.antimony_to_sbml("J: S -> ; k*S; S = 3; k = 0.25"))
    r['S'] = 5
    for export, initial in [(r.getSBML, 3), (r.getCurrentSBML, 5)]:
        text = export(level=level, version=version)
        root = ET.fromstring(text)
        assert (root.get('level'), root.get('version')) == (str(level), str(version))
        assert root.tag.startswith('{http://www.sbml.org/sbml/level')
        converted = RoadRunner(text)
        assert converted['S'] == initial
        assert converted['k'] == 0.25



def test_sbml_conversion_rejects_unsupported_features():
    from tellurium_nitro.roadrunner.extended_roadrunner import RoadRunner

    r = RoadRunner(rn.antimony_to_sbml(ANT['rules']))
    with pytest.raises(ValueError, match='initial assignments'):
        r.getSBML(2, 1)
    assert RoadRunner(r.getCurrentSBML(2, 1))['S2'] == r['S2']

    r = RoadRunner(rn.antimony_to_sbml("J: S -> ; max(1, S); S = 3"))
    with pytest.raises(ValueError, match='Level 3 Version 2'):
        r.getSBML(3, 1)

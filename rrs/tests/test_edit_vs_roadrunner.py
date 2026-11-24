"""Replay SBML edits against the outcomes in reference/edit.npz."""

import glob
import os
import re
import signal

import antimony
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'reference', 'edit.npz')
SNAPSHOT_IDS = ('FloatingSpecies', 'BoundarySpecies', 'GlobalParameter', 'Compartment', 'Reaction', 'Event')
SNAPSHOT_VALUES = ('FloatingSpeciesConcentrations', 'BoundarySpeciesConcentrations',
                   'GlobalParameterValues', 'CompartmentVolumes')
# Keep integration error below the comparison tolerance on both engines.
INTEGRATOR_TOLERANCES = dict(relative_tolerance=1e-12, absolute_tolerance=1e-15)
# The original comparison had no outcome for these libroadrunner segfaults.
REFERENCE_SEGFAULTS = {
    f'{model}/34-add_initial_assignment'
    for model in ('00014-sbml-l2v2.xml', '00022-sbml-l2v1.xml', '00033-sbml-l2v1.xml',
                  '10000-sbml-l2v1.xml', '10001-sbml-l2v1.xml', 'repressilator-sbml-l2v1.xml')
} | {
    f'{model}/35-add_initial_assignment'
    for model in ('00033-sbml-l2v1.xml', '10000-sbml-l2v1.xml', '10001-sbml-l2v1.xml',
                  'repressilator-sbml-l2v1.xml')
}

ANT = ['''function f(x, y)
  x*y
end
model a
compartment C = 2
J0: 2 S1 -> S2; f(k, S1)*exp(-t1*time)/(1 + S1^2)
J1: S2 -> ; k2*S2
k = 0.1; k2 = 0.2; S1 = 10; S2 = 0; t1 = 0.5
E1: at time > 1 && S1 < 5: k = 2
end''', '''model b
compartment V = 1.5
species A in V, B in V, $X in V
A = 5; B = 1; X = 2
J0: X -> A; k0*X
J1: A -> B; kf*A - kr*B
J2: B -> ; kd*B
k0 = 0.3; kf = 0.5; kr = 0.1; kd = 0.2; q = 0
p := 2*A
q' = 0.1*B
E0: at time > 4: kf = 1
end''']


def models():
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, '../../benchmarks/sosbench/*/*.xml'))):
        out.append((os.path.basename(f), open(f).read()))
    for src in ANT:
        assert antimony.loadAntimonyString(src) >= 0, antimony.getLastError()
        out.append(('ant:' + antimony.getMainModuleName(), antimony.getSBMLString(antimony.getMainModuleName())))
    # a local parameter for getParamPromotedSBML / removeReaction cleanup
    s = out[-1][1].replace('<ci> kd </ci>', '<apply><times/><ci> kd </ci><ci> kl </ci></apply>', 1)
    j = s.index('</math>', s.index('<reaction id="J2"')) + len('</math>')
    out.append(('ant:b+local', s[:j] + '<listOfLocalParameters><localParameter id="kl" value="3"/>'
                '<localParameter id="km" value="4"/></listOfLocalParameters>' + s[j:]))
    return out


def cases(sbml, ids):
    fl, _, gp, cp, rx, _ = ids
    level = 3 if 'level="3"' in sbml else 2
    a, b = fl[0], fl[-1]
    p = gp[0] if gp else 'k_none'
    rules = re.findall(r"<(?:assignment|rate)Rule[^>]*variable=\"(\w+)\"", sbml)
    con = ' constant="true"' if level == 3 else ''
    xml = (f'<reaction id="Jx" reversible="false"><listOfReactants><speciesReference species="{a}" stoichiometry="1"{con}/>'
           f'</listOfReactants><listOfProducts><speciesReference species="{b}" stoichiometry="2"{con}/></listOfProducts>'
           f'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><cn> 0.1 </cn><ci> {a} </ci></apply>'
           f'</math></kineticLaw></reaction>')
    ev = [('add_event', ('Ex', True, 'time > 2'))]
    one = lambda op, *args: [(op, args)]
    c = [
        one('add_species_amount', 'Sx', cp[0], 3.0, False, False, ''),
        one('add_species_amount', 'Sx', cp[0], 3.0, True, True, 'mole'),
        one('add_species_concentration', 'Sx', cp[0], 2.0, False, False, ''),
        one('add_species_amount', a, cp[0], 1.0, False, False, ''),
        one('add_species_concentration', 'Sx', 'nocomp', 1.0, False, False, ''),
        one('remove_species', a), one('remove_species', b), one('remove_species', 'nope'),
        one('set_boundary', a, True), one('set_has_only_substance_units', a, True),
        one('set_init_amount', a, 7.5), one('set_init_concentration', b, 0.25), one('set_init_amount', 'nope', 1.0),
        one('set_constant', p, False), one('set_constant', cp[0], False), one('set_constant', 'nope', True),
        one('add_reaction_from_sbml', xml),
        one('add_reaction', 'Jx', ['2' + a], [b], f'0.1*{a}'),
        one('add_reaction', 'Jx', [a], ['3.5' + b], f'0.1*{a}*{fl[len(fl) // 2]}/(1+{b})'),
        one('add_reaction', 'Jx', ['nope'], [b], '1'), one('add_reaction', rx[0], [a], [b], '1'),
        one('remove_reaction', rx[0], False), one('remove_reaction', rx[-1], True), one('remove_reaction', 'nope', False),
        one('set_reversible', rx[0], True), one('set_kinetic_law', rx[0], f'0.5*{a}'), one('set_kinetic_law', 'nope', '1'),
        one('add_parameter', 'px', 2.5), one('add_parameter', a, 1.0),
        one('remove_parameter', p), one('remove_parameter', 'nope'),
        one('add_compartment', 'cx', 2.0), one('add_compartment', cp[0], 2.0), one('remove_compartment', cp[0]),
        one('add_initial_assignment', a, '2*3'), one('add_initial_assignment', p, f'{a}*2'),
        one('add_initial_assignment', 'nope', '1'),
        one('add_assignment_rule', 'px0', '1') , [('add_parameter', ('px', 1.0)), ('add_assignment_rule', ('px', f'time*0.1 + {a}'))],
        one('add_assignment_rule', a, '1'), one('add_rate_rule', p, '0.01'),
        [('set_boundary', (a, True)), ('add_rate_rule', (a, '-0.1'))],
        ev, ev + one('add_event_assignment', 'Ex', p, '1.5'), ev + one('add_event_assignment', 'Ex', a, f'{a}/2'),
        ev + one('add_event_assignment', 'Ex', 'nope', '1'),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('add_delay', 'Ex', '0.5'),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('add_priority', 'Ex', '1'),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('set_persistent', 'Ex', True),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('set_trigger_initial_value', 'Ex', True),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('add_trigger', 'Ex', 'time > 3'),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('remove_event_assignments', 'Ex', p),
        ev + one('remove_event_assignments', 'Ex', 'nope'),
        ev + one('add_event_assignment', 'Ex', p, '1.5') + one('remove_event', 'Ex'), one('remove_event', 'nope'),
        ev + one('add_event_assignment', 'Ex', p, '3') + one('remove_parameter', p),
    ]
    c += [one('remove_rules', v) for v in rules] + [one('remove_rules', 'nope')]
    c += [one('get_kinetic_law', rid) for rid in rx]
    c += [one('get_param_promoted_sbml'), one('validate_current_sbml')]
    return c


def keyed_cases(name, sbml, ids):
    for i, steps in enumerate(cases(sbml, ids)):
        ops = '+'.join(op for op, _ in steps)
        yield f'{name}/{i:02d}-{ops}', steps


def model_ids(model):
    return tuple(list(getattr(model, f'get{kind}Ids')()) for kind in SNAPSHOT_IDS)


def snapshot(r, model):
    ids = model_ids(model)
    init = np.concatenate([getattr(model, f'get{kind}')() for kind in SNAPSHOT_VALUES])
    return dict(ids=ids, init=init, sim=np.asarray(r.simulate(0, 10, 51)))


def replay_cases():
    with np.load(REF, allow_pickle=True) as ref:
        for name, sbml in models():
            for key, steps in keyed_cases(name, sbml, ref[f'{name}/ids']):
                expected = ref[key].item()
                if 'signal' in expected:
                    assert key in REFERENCE_SEGFAULTS and expected['signal'] == signal.SIGSEGV, key
                    continue
                yield pytest.param(sbml, steps, expected, id=key)


def assert_close(actual, expected):
    scale = max(1.0, np.nanmax(np.abs(expected), initial=0))
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-9 * scale)


def pytest_generate_tests(metafunc):
    if metafunc.function is test_edit:
        metafunc.parametrize('sbml, steps, expected', list(replay_cases()))


def test_edit(sbml, steps, expected):
    import roadrunner_nitro as nitro

    op, args = steps[0]
    if op == 'get_kinetic_law':
        rid, = args

        def rates(formula):
            # Keep the original reaction's local parameters and model functions in scope.
            edited = nitro.sbml_edit('set_kinetic_law', sbml, rid, formula)
            return np.asarray(nitro.RoadRunner(edited).simulate(0, 1, 6, [rid]))

        formula = nitro.sbml_edit(op, sbml, rid)
        assert_close(rates(formula), rates(expected['formula']))
        return
    if op == 'validate_current_sbml':
        if expected['valid']:
            nitro.sbml_edit(op, sbml)
        else:
            with pytest.raises(ValueError):
                nitro.sbml_edit(op, sbml)
        return

    def apply():
        edited = sbml
        for op, args in steps:
            edited = nitro.sbml_edit(op, edited, *args)
        r = nitro.RoadRunner(edited)
        for setting, value in INTEGRATOR_TOLERANCES.items():
            r.setIntegratorSetting(setting, value)
        return snapshot(r, r)

    if 'error' in expected:
        with pytest.raises(ValueError):
            apply()
    else:
        actual = apply()
        assert actual['ids'] == expected['ids']
        assert_close(actual['init'], expected['init'])
        assert_close(actual['sim'], expected['sim'])

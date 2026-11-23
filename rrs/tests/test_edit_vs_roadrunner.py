"""Model editing (roadrunner_nitro.sbml_edit, rrs/src/sbml/edit.rs) vs libroadrunner's own editing.

For every model and every edit sequence: apply it with libroadrunner (r.addSpecies... then
r.getSBML()) and with sbml_edit, load both results into a fresh libroadrunner, and compare
ids, initial values and a short simulation (rtol 1e-6). Both failing (edit or load) counts as a match.
Run with the venv python, roadrunner_nitro built from this tree first on PYTHONPATH.
"""
import glob, os, pickle, re, sys
import numpy as np
import antimony, roadrunner
import roadrunner_nitro as nitro

roadrunner.Logger.setLevel(roadrunner.Logger.LOG_FATAL)
HERE = os.path.dirname(os.path.abspath(__file__))

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


def camel(op):
    special = {'add_reaction_from_sbml': 'addReaction', 'add_reaction': 'addReaction',
               'validate_current_sbml': 'validateCurrentSBML'}
    if op in special:
        return special[op]
    a, *b = op.split('_')
    return a + ''.join(x.title() for x in b)


def rr_apply(sbml, steps):
    r = roadrunner.RoadRunner(sbml)
    for op, args in steps:
        if op == 'remove_reaction':
            r.removeReaction(args[0], args[1])
        elif op == 'add_reaction_from_sbml':  # the Python wrapper's addReaction(*args) hides this overload
            roadrunner._roadrunner.RoadRunner_addReaction(r, args[0], True)
        else:
            getattr(r, camel(op))(*args)
    return r.getSBML()  # the edited document (getCurrentSBML bakes in current values, dropping initial assignments)


def ours_apply(sbml, steps):
    for op, args in steps:
        sbml = nitro.sbml_edit(op, sbml, *args)
    return sbml


def snapshot(sbml):
    r = roadrunner.RoadRunner(sbml)
    ids = (r.model.getFloatingSpeciesIds(), r.model.getBoundarySpeciesIds(), r.model.getGlobalParameterIds(),
           r.model.getCompartmentIds(), r.model.getReactionIds(), list(r.model.getEventIds()))
    init = np.concatenate([r.model.getFloatingSpeciesConcentrations(), r.model.getBoundarySpeciesConcentrations(),
                           r.model.getGlobalParameterValues(), r.model.getCompartmentVolumes()])
    sim = np.asarray(r.simulate(0, 10, 51))
    return ids, init, sim


def close(a, b):
    return a.shape == b.shape and np.allclose(a, b, rtol=1e-6, atol=1e-9 * max(1.0, np.nanmax(np.abs(a), initial=0)), equal_nan=True)


def result(fn, sbml, steps):
    """snapshot(fn(sbml, steps)) in a forked child (libroadrunner segfaults on some edits): the snapshot,
    'ERR: ...' if it raised, CRASH if the child died."""
    rd, wr = os.pipe()
    pid = os.fork()
    if pid == 0:
        os.close(rd)
        try:
            out = snapshot(fn(sbml, steps))
        except Exception as e:
            out = 'ERR: ' + str(e).splitlines()[0][:100] if str(e) else 'ERR'
        with os.fdopen(wr, 'wb') as f:
            pickle.dump(out, f)
        os._exit(0)
    os.close(wr)
    with os.fdopen(rd, 'rb') as f:
        data = f.read()
    os.waitpid(pid, 0)
    return pickle.loads(data) if data else CRASH


CRASH = 'CRASH'


def show(x, other):
    """An error, or what differs from `other`: ids not in it, else 'values'."""
    if isinstance(x, str) or isinstance(other, str):
        return x if isinstance(x, str) else 'loaded fine'
    d = [sorted(set(a) - set(b)) for a, b in zip(x[0], other[0])]
    return f'ids {d}, order {[a == b for a, b in zip(x[0], other[0])]}' if x[0] != other[0] else 'same ids, values differ'


def cases(sbml):
    r = roadrunner.RoadRunner(sbml)
    fl, bd, gp = r.model.getFloatingSpeciesIds(), r.model.getBoundarySpeciesIds(), r.model.getGlobalParameterIds()
    cp, rx = r.model.getCompartmentIds(), r.model.getReactionIds()
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
    return c, rx


def main():
    n = bad = both_err = 0
    tally, crashed = {}, []
    for name, sbml in models():
        cs, rx = cases(sbml)
        for steps in cs:
            x, y = result(rr_apply, sbml, steps), result(ours_apply, sbml, steps)
            if x == CRASH:  # libroadrunner died: nothing to compare with
                crashed.append((name, steps[-1][0], 'ours: ' + (y if isinstance(y, str) else 'loads and simulates')))
                continue
            if isinstance(x, str) or isinstance(y, str):
                ok = isinstance(x, str) and isinstance(y, str)
                both_err += ok
            else:
                ok = x[0] == y[0] and close(x[1], y[1]) and close(x[2], y[2])
            n += 1
            key = '+'.join(op for op, _ in steps)
            t = tally.setdefault(key, [0, 0]); t[0] += ok; t[1] += 1
            if not ok:
                bad += 1
                print(f'MISMATCH {name} {[(op, tuple(str(v)[:40] for v in a)) for op, a in steps]}\n  roadrunner: {show(x, y)}\n  ours:       {show(y, x)}')
        # getKineticLaw: our formula, read back by libsbml, prints as libsbml's own
        r = roadrunner.RoadRunner(sbml)
        for rid in rx:
            ours = nitro.sbml_edit('get_kinetic_law', sbml, rid)
            r2 = roadrunner.RoadRunner(sbml)
            r2.setKineticLaw(rid, ours, False)
            ok = r2.getKineticLaw(rid) == r.getKineticLaw(rid)
            n += 1; bad += not ok
            t = tally.setdefault('get_kinetic_law', [0, 0]); t[0] += ok; t[1] += 1
            if not ok:
                print(f'MISMATCH {name} getKineticLaw {rid}: {r.getKineticLaw(rid)!r} vs ours {ours!r}')
        # getParamPromotedSBML
        x, y = result(lambda s, _: roadrunner.RoadRunner.getParamPromotedSBML(s), sbml, None), \
            result(lambda s, _: nitro.sbml_edit('get_param_promoted_sbml', s), sbml, None)
        ok = not isinstance(x, str) and not isinstance(y, str) and x[0] == y[0] and close(x[1], y[1]) and close(x[2], y[2])
        n += 1; bad += not ok
        t = tally.setdefault('get_param_promoted_sbml', [0, 0]); t[0] += ok; t[1] += 1
        if not ok:
            print(f'MISMATCH {name} getParamPromotedSBML\n  {x if isinstance(x, str) else x[0]}\n  {y if isinstance(y, str) else y[0]}')
        # validateCurrentSBML: pass / fail only
        try:
            roadrunner.RoadRunner(sbml).validateCurrentSBML(); vx = True
        except Exception:
            vx = False
        try:
            nitro.sbml_edit('validate_current_sbml', sbml); vy = True
        except ValueError:
            vy = False
        t = tally.setdefault('validate_current_sbml', [0, 0]); t[0] += vx == vy; t[1] += 1
        n += 1; bad += vx != vy
        if vx != vy:
            print(f'MISMATCH {name} validateCurrentSBML: roadrunner {vx}, ours {vy}')
    for k, (g, t) in sorted(tally.items()):
        print(f'{g:3d}/{t:<3d} {k}')
    for c in crashed:
        print('libroadrunner crashed (skipped):', *c)
    print(f'{n - bad}/{n} match ({both_err} of them: both raise)')
    return bad


if __name__ == '__main__':
    print(nitro.__file__)
    sys.exit(1 if main() else 0)

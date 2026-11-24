"""Record the missing libroadrunner edit outcomes in reference/edit.npz.

Run with the Python environment that has libroadrunner installed.
"""
import os
import pickle

import numpy as np
import roadrunner

from test_edit_vs_roadrunner import INTEGRATOR_TOLERANCES, REF, keyed_cases, model_ids, models, snapshot


def camel(op):
    a, *b = op.split('_')
    return a + ''.join(x.title() for x in b)


def record(sbml, steps):
    op, args = steps[0]
    if op == 'get_param_promoted_sbml':
        edited = roadrunner.RoadRunner.getParamPromotedSBML(sbml)
    else:
        r = roadrunner.RoadRunner(sbml)
        if op == 'get_kinetic_law':
            return dict(formula=r.getKineticLaw(*args))
        if op == 'validate_current_sbml':
            r.validateCurrentSBML()
            return dict(valid=True)
        for op, args in steps:
            if op == 'add_reaction_from_sbml':
                # The Python wrapper hides this overload.
                roadrunner._roadrunner.RoadRunner_addReaction(r, args[0], True)
            else:
                getattr(r, camel(op))(*args)
        edited = r.getSBML()
    r = roadrunner.RoadRunner(edited)
    for setting, value in INTEGRATOR_TOLERANCES.items():
        setattr(r.integrator, setting, value)
    return snapshot(r, r.model)


def fork_record(sbml, steps):
    """Run one case in a child because libroadrunner can segfault."""
    rd, wr = os.pipe()
    pid = os.fork()
    if pid == 0:
        os.close(rd)
        try:
            result = record(sbml, steps)
        except Exception as exc:
            result = dict(error=str(exc))
        with os.fdopen(wr, 'wb') as f:
            pickle.dump(result, f)
        os._exit(0)
    os.close(wr)
    with os.fdopen(rd, 'rb') as f:
        data = f.read()
    _, status = os.waitpid(pid, 0)
    if os.WIFSIGNALED(status):
        return dict(signal=os.WTERMSIG(status))
    if status or not data:
        raise RuntimeError(f'Edit recorder exited with status {status}')
    return pickle.loads(data)


def main():
    roadrunner.Logger.setLevel(roadrunner.Logger.LOG_FATAL)
    ref = {}
    if os.path.exists(REF):
        with np.load(REF, allow_pickle=True) as z:
            ref = dict(z)
    for name, sbml in models():
        ids_key = f'{name}/ids'
        changed = ids_key not in ref
        if changed:
            r = roadrunner.RoadRunner(sbml)
            ref[ids_key] = np.array(model_ids(r.model), dtype=object)
        added = 0
        for key, steps in keyed_cases(name, sbml, ref[ids_key]):
            if key not in ref:
                result = fork_record(sbml, steps)
                if steps[0][0] == 'validate_current_sbml' and 'error' in result:
                    result = dict(valid=False)
                if 'signal' in result:
                    print(f'{key}: libroadrunner died with signal {result["signal"]}', flush=True)
                ref[key] = np.array(result, dtype=object)
                added += 1
        if added or changed:
            tmp = REF + '.tmp.npz'
            np.savez_compressed(tmp, **ref)
            os.replace(tmp, REF)
        print(f'{name}: {added} new cases', flush=True)


if __name__ == '__main__':
    main()

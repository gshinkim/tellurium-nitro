"""Record libroadrunner's answers for test_turbo.py into reference/turbo.npz.

    python tellurium_nitro/tests/make_reference.py      (needs libroadrunner and python-libsbml)
"""
import os
import sys

import libsbml
import numpy as np
import roadrunner

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_turbo as T  # noqa: E402

roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)


def load(sbml):
    return roadrunner.RoadRunner(sbml)


def edit(ant, f):
    """SBML of `ant` after f(libsbml model) (features Antimony cannot write)."""
    d = libsbml.readSBMLFromString(T.te.antimonyToSBML(ant))
    f(d.getModel())
    return libsbml.writeSBMLToString(d)


def cf(m, model_level, constant=True):
    p = m.createParameter()
    p.setId('cf'), p.setValue(2.0), p.setConstant(constant)
    (m if model_level else m.getSpecies('S2')).setConversionFactor('cf')


def cf_event(m):
    cf(m, True, False)
    ea = m.getEvent(0).createEventAssignment()
    ea.setVariable('cf'), ea.setMath(libsbml.parseL3Formula('3'))


def vstoich(m, how):
    sr = m.getReaction('J1').getProduct(0)
    sr.setId('n2'), sr.setConstant(False)
    if how == 'event':
        e = m.createEvent()
        e.setUseValuesFromTriggerTime(True)
        t = e.createTrigger()
        t.setMath(libsbml.parseL3Formula('time > 3')), t.setPersistent(True), t.setInitialValue(True)
        ea = e.createEventAssignment()
        ea.setVariable('n2'), ea.setMath(libsbml.parseL3Formula('3'))
    else:
        rule = m.createAssignmentRule() if how == 'rule' else m.createRateRule()
        rule.setVariable('n2'), rule.setMath(libsbml.parseL3Formula('1 + time/5' if how == 'rule' else '0.2'))


def stoich_math():
    d = libsbml.readSBMLFromString(T.te.antimonyToSBML("J1: S1 -> S2; k*S1; k=1; S1=10; S2=0; n = 2"))
    d.setLevelAndVersion(2, 4, False)
    d.getModel().getReaction('J1').getProduct(0).createStoichiometryMath().setMath(libsbml.parseL3Formula('n + time/10'))
    return libsbml.writeSBMLToString(d)


SBML_FEATURES = {
    'conversion factor (species)': edit(T.BASE, lambda m: cf(m, False)),
    'conversion factor (model)': edit(T.BASE, lambda m: cf(m, True)),
    'conversion factor (event-changed)': edit(T.BASE + "; at time>2: k=2", cf_event),
    'stoichiometry: assignment rule': edit(T.BASE, lambda m: vstoich(m, 'rule')),
    'stoichiometry: rate rule': edit(T.BASE, lambda m: vstoich(m, 'rate')),
    'stoichiometry: event': edit(T.BASE, lambda m: vstoich(m, 'event')),
    'stoichiometryMath (L2)': stoich_math(),
}

if __name__ == '__main__':
    out = {}
    for name, (f, _, _) in T.WORKFLOWS.items():
        out.update(T.flatten(name, f(load)))
        print(name, flush=True)
    for name in list(T.FEATURES) + list(SBML_FEATURES):
        sbml = SBML_FEATURES.get(name) or T.feature_sbml(name)
        extra = [e if isinstance(e, tuple) else (e, e) for e in T.FEATURES.get(name, (None, []))[1]]
        out.update(T.flatten('feature/' + name, T.f_feature(load, sbml, extra)))
        if name in SBML_FEATURES:
            out['feature/%s/sbml' % name] = np.array(sbml)
    out.update(T.flatten('events', T.f_events(load)))
    out.update(T.flatten('rate_mca', T.f_rate_mca(load)))
    os.makedirs(os.path.dirname(T.REF), exist_ok=True)
    np.savez_compressed(T.REF, **out)
    print(len(out), 'entries')

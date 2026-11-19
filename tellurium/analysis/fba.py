"""Flux balance analysis of SBML models with the fbc package (versions 1 and 2): a linear program,
solved with SciPy's HiGHS. roadrunner does not do FBA; this is a small standalone solver."""
import os

import libsbml
import numpy as np


def fba(model):
    """Optimal steady-state fluxes of an fbc model (SBML text or file path).

    Maximizes or minimizes the active fbc objective subject to N v = 0 over the non-boundary,
    non-constant species, and to the flux bounds (fbc v2 bound parameters, or fbc v1 fluxBounds;
    without a lower bound an irreversible reaction is >= 0). Parameter and stoichiometry values
    include initial assignments.

    :returns: dict reaction id -> flux, plus the objective id -> objective value; all NaN if the
        program is infeasible or unbounded.
    """
    from scipy.optimize import linprog
    from .rrs_bridge import from_sbml, _drop_empty_math
    text, base = model, None
    if not text.lstrip().startswith('<') and os.path.isfile(text):
        base = os.path.dirname(os.path.abspath(text))
        with open(text, encoding='utf-8') as f:
            text = f.read()
    info = {}
    from_sbml(text, info, base)  # initial values as the simulator computes them
    vals = info['values']
    m = libsbml.readSBMLFromString(text).getModel()
    _drop_empty_math(m)
    fb = m.getPlugin('fbc')
    if fb is None:
        raise ValueError("the model has no fbc package")

    def value(x):
        return float(vals[x]) if x in vals else m.getElementBySId(x).getValue()

    rxns = list(m.getListOfReactions())
    rows = [s.getId() for s in m.getListOfSpecies() if not s.getBoundaryCondition() and not s.getConstant()]
    ri = {s: i for i, s in enumerate(rows)}
    N = np.zeros((len(rows), len(rxns)))
    lb, ub = np.full(len(rxns), -np.inf), np.full(len(rxns), np.inf)
    for j, rx in enumerate(rxns):
        for sg, lst in ((-1, rx.getListOfReactants()), (1, rx.getListOfProducts())):
            for sr in lst:
                if sr.getSpecies() in ri:
                    n = value(sr.getId()) if sr.isSetId() and sr.getId() in vals else (
                        sr.getStoichiometry() if sr.isSetStoichiometry() else 1.0)
                    N[ri[sr.getSpecies()], j] += sg * n
        rp = rx.getPlugin('fbc')
        if rp is not None and rp.isSetLowerFluxBound():
            lb[j] = value(rp.getLowerFluxBound())
        if rp is not None and rp.isSetUpperFluxBound():
            ub[j] = value(rp.getUpperFluxBound())
    col = {rx.getId(): j for j, rx in enumerate(rxns)}
    for b in fb.getListOfFluxBounds():  # fbc version 1
        j, op, v = col[b.getReaction()], b.getOperation(), b.getValue()
        if op in ('lessEqual', 'less'):
            ub[j] = min(ub[j], v)
        elif op in ('greaterEqual', 'greater'):
            lb[j] = max(lb[j], v)
        elif op == 'equal':
            lb[j], ub[j] = max(lb[j], v), min(ub[j], v)
    for j, rx in enumerate(rxns):  # no lower bound given: an irreversible reaction runs forward only
        if lb[j] == -np.inf and not rx.getReversible():
            lb[j] = 0.0
    obj = fb.getActiveObjective() or fb.getObjective(0)
    c = np.zeros(len(rxns))
    for fo in obj.getListOfFluxObjectives():
        c[col[fo.getReaction()]] += fo.getCoefficient()
    sign = -1.0 if obj.getType() == 'maximize' else 1.0
    res = linprog(sign * c, A_eq=N if len(rows) else None, b_eq=np.zeros(len(rows)) if len(rows) else None,
                  bounds=list(zip(lb, ub)), method='highs')
    x = res.x if res.status == 0 else np.full(len(rxns), np.nan)  # infeasible/unbounded: no optimum
    out = {rx.getId(): float(v) for rx, v in zip(rxns, x)}
    out[obj.getId()] = float(c @ x)
    return out

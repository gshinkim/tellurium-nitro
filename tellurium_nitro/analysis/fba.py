"""Flux balance analysis of SBML models with the fbc package (versions 1 and 2): a linear program,
solved with SciPy's HiGHS. roadrunner does not do FBA; this is a small standalone solver."""
import os

import xml.etree.ElementTree as ET

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
    from .rrs_bridge import from_sbml
    text, base = model, None
    if not text.lstrip().startswith('<') and os.path.isfile(text):
        base = os.path.dirname(os.path.abspath(text))
        with open(text, encoding='utf-8') as f:
            text = f.read()
    info = {}
    from_sbml(text, info, base)  # initial values as the simulator computes them
    vals = info['values']
    root = ET.fromstring(text)
    m = _child(root, 'model')
    fb = _find(m, 'listOfObjectives')
    if fb is None:
        raise ValueError("the model has no fbc package")
    elems = {e.get('id'): e for e in m.iter() if e.get('id')}

    def value(x):
        return float(vals[x]) if x in vals else float(elems[x].get('value'))

    rxns = _find(m, 'listOfReactions', 'reaction')
    ids = [rx.get('id') for rx in rxns]
    rows = [s.get('id') for s in _find(m, 'listOfSpecies', 'species')
            if s.get('boundaryCondition') != 'true' and s.get('constant') != 'true']
    ri = {s: i for i, s in enumerate(rows)}
    N = np.zeros((len(rows), len(rxns)))
    lb, ub = np.full(len(rxns), -np.inf), np.full(len(rxns), np.inf)
    for j, rx in enumerate(rxns):
        for sg, lst in ((-1, 'listOfReactants'), (1, 'listOfProducts')):
            for sr in _find(rx, lst, 'speciesReference'):
                if sr.get('species') in ri:
                    n = value(sr.get('id')) if sr.get('id') in vals else float(sr.get('stoichiometry', 1.0))
                    N[ri[sr.get('species')], j] += sg * n
        if _attr(rx, 'lowerFluxBound') is not None:
            lb[j] = value(_attr(rx, 'lowerFluxBound'))
        if _attr(rx, 'upperFluxBound') is not None:
            ub[j] = value(_attr(rx, 'upperFluxBound'))
    col = {r: j for j, r in enumerate(ids)}
    for b in _find(m, 'listOfFluxBounds', 'fluxBound'):  # fbc version 1
        j, op, v = col[_attr(b, 'reaction')], _attr(b, 'operation'), float(_attr(b, 'value'))
        if op in ('lessEqual', 'less'):
            ub[j] = min(ub[j], v)
        elif op in ('greaterEqual', 'greater'):
            lb[j] = max(lb[j], v)
        elif op == 'equal':
            lb[j], ub[j] = max(lb[j], v), min(ub[j], v)
    for j, rx in enumerate(rxns):  # no lower bound given: an irreversible reaction runs forward only
        if lb[j] == -np.inf and rx.get('reversible', 'true') == 'false':
            lb[j] = 0.0
    objs = [o for o in fb if _local(o.tag) == 'objective']
    act = _attr(fb, 'activeObjective')
    obj = next((o for o in objs if _attr(o, 'id') == act), objs[0])
    c = np.zeros(len(rxns))
    for fo in _find(obj, 'listOfFluxObjectives', 'fluxObjective'):
        c[col[_attr(fo, 'reaction')]] += float(_attr(fo, 'coefficient'))
    sign = -1.0 if _attr(obj, 'type') == 'maximize' else 1.0
    res = linprog(sign * c, A_eq=N if len(rows) else None, b_eq=np.zeros(len(rows)) if len(rows) else None,
                  bounds=list(zip(lb, ub)), method='highs')
    x = res.x if res.status == 0 else np.full(len(rxns), np.nan)  # infeasible/unbounded: no optimum
    out = {r: float(v) for r, v in zip(ids, x)}
    out[_attr(obj, 'id')] = float(c @ x)
    return out


def _local(tag):
    return tag.rsplit('}', 1)[-1]


def _attr(e, name):  # attribute by local name (fbc attributes carry the fbc namespace)
    return next((v for k, v in e.attrib.items() if _local(k) == name), None)


def _child(e, name):
    return next((c for c in e if _local(c.tag) == name), None)


def _find(e, lst, item=None):
    """Child list ``lst`` of e (None if missing); with ``item``, its ``item`` children ([] if missing)."""
    l = _child(e, lst)
    if item is None:
        return l
    return [] if l is None else [c for c in l if _local(c.tag) == item]

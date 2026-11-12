"""Build a ``roadrunner_rs.Model`` (Rust steady-state/MCA/simulation engine) from a roadrunner model.

    m = from_roadrunner(r)
    res = m.mca_batch(['k1'], np.array([[0.5], [1.0]]))
"""
import re

import libsbml
import numpy as np
import roadrunner
import roadrunner_rs


def _names(ast, out):
    if ast.isName():
        out.add(ast.getName())
    for i in range(ast.getNumChildren()):
        _names(ast.getChild(i), out)
    return out


_OUT = libsbml.L3ParserSettings()
_OUT.setParseUnits(False)  # no '1 dimensionless' in formulas


def _l3(ast):
    return libsbml.formulaToL3StringWithSettings(ast, _OUT)


def _name(n):
    ast = libsbml.ASTNode(libsbml.AST_NAME)  # plain name (parseL3Formula('time') gives a csymbol)
    ast.setName(n)
    return ast


def from_roadrunner(r, info=None):
    """Compile roadrunner model ``r`` (not modified) into a ``roadrunner_rs.Model``.

    Covers what roadrunner loads (it refuses delay(), algebraic rules and fast reactions), except
    distrib random draws and the rare V := f(species in V) or rateOf of a piecewise assignment rule
    (NotImplementedError): rate and
    assignment rules on any symbol, events, function definitions, avogadro/rateOf/time csymbols,
    reaction ids in math, reactions without kinetic law (rate 0), conversion factors, variable
    stoichiometry and variable compartments.

    State is concentrations. Assignment-rule species/compartments/parameters are outputs (in
    output_ids, not state_ids). A species whose d[S]/dt is not stoich * rates / constant volume
    (variable compartment or stoichiometry, non-constant conversion factor) gets a zero stoichiometry
    row and d[S]/dt = cf * sum_j n_j v_j / V - [S] V'/V as a rate rule; amounts are kept when an
    event resizes its compartment, as in roadrunner.

    ``info`` (dict), if given, receives: comp_of (species -> compartment), changing (ids set by rules or
    events), amount_params (non-floating
    species whose parameter value is an amount, not a concentration), general (species whose d[S]/dt
    is a rate rule in Rust although they have reactions: variable compartment/stoichiometry/cf).
    """
    sbml = r.getCurrentSBML()
    rr = r  # read values from the loaded model: a fresh RoadRunner(sbml) is a full re-load (~10 s at 1000 reactions)
    doc = libsbml.readSBMLFromString(sbml)
    m = doc.getModel()
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    odoc = libsbml.readSBMLFromString(r.getSBML())  # getCurrentSBML() drops stoichiometryMath

    rules, rate_rules = {}, {}
    for rule in m.getListOfRules():
        (rate_rules if rule.isRate() else rules)[rule.getVariable()] = rule.getMath()
    targets = {ea.getVariable() for ev in m.getListOfEvents() for ea in ev.getListOfEventAssignments()}
    changing = set(rules) | set(rate_rules) | targets

    if m.getNumReactions():
        Nm = rr.getFullStoichiometryMatrix()
        if sorted(Nm.colnames) != sorted(x.getId() for x in m.getListOfReactions()):
            # roadrunner's structure is stale after some model edits: re-load from the current SBML
            rr = roadrunner.roadrunner.RoadRunner(sbml)  # plain libroadrunner, not tellurium's facade
            Nm = rr.getFullStoichiometryMatrix()
        rows, reactions = list(Nm.rownames), list(Nm.colnames)
    else:  # roadrunner has no stoichiometry matrix then
        rows, reactions = list(rr.model.getFloatingSpeciesIds()), []
        Nm = np.zeros((len(rows), 0))
    conc = dict(zip(rr.model.getFloatingSpeciesIds(), rr.model.getFloatingSpeciesConcentrations()))
    sp = {s.getId(): s for s in m.getListOfSpecies()}
    species = [s for s in rows if s not in rules]  # assignment-rule species: outputs only
    floating = set(species)
    comp_of = {s: sp[s].getCompartment() for s in sp}
    hosu = {s for s in sp if sp[s].getHasOnlySubstanceUnits()}
    vol = {s: rr.getValue(comp_of[s]) for s in sp}

    params = {}
    for x in list(m.getListOfParameters()) + list(m.getListOfCompartments()):
        if x.getId() not in rules:
            params[x.getId()] = rr.getValue(x.getId())
    for s in sp:  # boundary / non-floating species: parameters in their SBML units
        if s not in rules and s not in floating:
            params[s] = rr.getValue('[%s]' % s) * (vol[s] if s in hosu else 1.0)

    # variable stoichiometry: (species, reaction) -> [(sign, stoichiometry AST)] over all references
    def refs_of(model):
        return [(sg, rxn.getId(), sr) for rxn in model.getListOfReactions()
                for sg, lst in ((-1, rxn.getListOfReactants()), (1, rxn.getListOfProducts())) for sr in lst]
    refs = refs_of(m)
    smath = {(j, sr.getSpecies(), sg, i): sr.getStoichiometryMath().getMath()
             for i, (sg, j, sr) in enumerate(refs_of(odoc.getModel()))
             if sr.isSetStoichiometryMath() and sr.getStoichiometryMath().isSetMath()}
    smath = {i: smath.get((j, sr.getSpecies(), sg, i)) for i, (sg, j, sr) in enumerate(refs)}
    for _, _, sr in refs:
        if sr.isSetId() and sr.getId() not in rules:
            params[sr.getId()] = rr.getValue(sr.getId()) if sr.getId() in changing else sr.getStoichiometry()
    var_ref = {(sr.getSpecies(), j) for i, (_, j, sr) in enumerate(refs)
               if (sr.isSetId() and sr.getId() in changing) or smath[i] is not None}
    vstoich = {}
    for i, (sg, j, sr) in enumerate(refs):
        if (sr.getSpecies(), j) in var_ref:
            if smath[i] is not None:
                a = smath[i]
            else:
                a = _name(sr.getId()) if sr.isSetId() else libsbml.parseL3Formula(repr(sr.getStoichiometry()))
            vstoich.setdefault((sr.getSpecies(), j), []).append((sg, a))

    def cf_of(s):
        return (sp[s].getConversionFactor() if sp[s].isSetConversionFactor() else m.getConversionFactor()) or None

    cvar = {c.getId() for c in m.getListOfCompartments() if c.getId() in changing}
    # boundary species in a variable compartment keep their amount (parameter); [S] = amount / V
    bamt = {s for s in sp if s not in floating and s not in rules and s not in rate_rules and s not in hosu
            and comp_of[s] in cvar}
    for s in bamt:
        params[s] *= vol[s]
    general = {s for s in species if comp_of[s] in cvar or any(k[0] == s for k in vstoich)
               or cf_of(s) in changing}

    local, kl_math = {}, {}
    for rxn in m.getListOfReactions():
        rid, kl = rxn.getId(), rxn.getKineticLaw()
        kl_math[rid] = kl.getMath() if kl is not None and kl.getMath() is not None else None
        local[rid] = []
        if kl is not None:
            for p in list(kl.getListOfParameters()) + list(kl.getListOfLocalParameters()):
                local[rid].append((p.getId(), '%s__%s' % (rid, p.getId())))
                params[local[rid][-1][1]] = p.getValue()

    def kinetic(rid):
        if kl_math[rid] is None:  # no kinetic law: rate 0, as roadrunner
            return libsbml.parseL3Formula('0')
        ast = kl_math[rid].deepCopy()
        for old, new in local[rid]:
            ast.renameSIdRefs(old, new)
        return expand(ast)

    def expand(ast):
        """Copy with function calls inlined, csymbols resolved and reaction ids -> kinetic laws."""
        t = ast.getType()
        if t == libsbml.AST_NAME_AVOGADRO:
            return libsbml.parseL3Formula(repr(ast.getReal()))
        if t == libsbml.AST_NAME_TIME:
            return _name('time')
        if t == libsbml.AST_FUNCTION_RATE_OF:
            return _name('__rateof__' + ast.getChild(0).getName())
        if libsbml.AST_DISTRIB_FUNCTION_NORMAL <= t <= libsbml.AST_DISTRIB_FUNCTION_RAYLEIGH:
            # ponytail: roadrunner draws a fresh random number per evaluation; add an RNG op if needed
            raise NotImplementedError("distrib random function %s()" % ast.getName())
        if t == libsbml.AST_NAME and ast.getName() in kl_math:
            return kinetic(ast.getName())
        out = ast.deepCopy()
        for i in range(ast.getNumChildren()):
            out.replaceChild(i, expand(ast.getChild(i)), True)
        if t == libsbml.AST_FUNCTION and ast.getName() in fdefs:
            lam = fdefs[ast.getName()]
            n = lam.getNumChildren() - 1
            body = lam.getChild(n).deepCopy()
            for i in range(n):  # via temporary names: an argument may mention another parameter's name
                body.renameSIdRefs(lam.getChild(i).getName(), '__arg%d__' % i)
            for i in range(n):
                body.replaceArgument('__arg%d__' % i, out.getChild(i))
            return expand(body)
        return out

    def inline(ast):
        """Substitute assignment rules and amounts ([S] * V) of hasOnlySubstanceUnits species."""
        for _ in range(1000):
            todo = _names(ast, set()) & (set(rules) | (hosu & floating) | bamt)
            if not todo:
                for v in (hosu & floating) | bamt:
                    ast.renameSIdRefs('__conc__' + v, v)
                return ast
            for v in todo:
                new = expand(rules[v]) if v in rules else libsbml.parseL3Formula(
                    '(__conc__%s %s %s)' % (v, '/' if v in bamt else '*', comp_of[v]))
                if ast.isName() and ast.getName() == v:  # replaceIDWithFunction skips the root
                    ast = new
                else:
                    ast.replaceIDWithFunction(v, new)
        raise ValueError("cyclic assignment rules")

    known = floating | set(params) | {'time'}

    def formula(math, what):
        ast = inline(expand(math))
        unknown = {n for n in _names(ast, set()) if not n.startswith('__rateof__')} - known
        if unknown:
            raise ValueError("%s references %s" % (what, sorted(unknown)))
        return re.sub(r'\b__rateof__(\w+)', lambda g: '(%s)' % rate_of(g.group(1)), _l3(ast))

    memo, busy = {}, set()

    def lazy(key, f):
        if key not in memo:
            if key in busy:
                raise NotImplementedError("implicit rate of change through %s" % (key,))  # e.g. V := f([S in V])
            busy.add(key)
            memo[key] = f()
            busy.discard(key)
        return memo[key]

    def rate(j):
        return lazy(('v', j), lambda: formula(kinetic(reactions[j]), "kinetic law of %s" % reactions[j]))

    def rate_of(x):
        """d x/dt in x's SBML units (rateOf csymbol)."""
        def f():
            if x in rate_rules:
                return formula(rate_rules[x], "rate rule for %s" % x)
            if x in floating:
                if x not in hosu:
                    return dconc(x)
                V = formula(_name(comp_of[x]), x)
                return '(%s) * %s + %s * (%s)' % (dconc(x), V, x, rate_of(comp_of[x]))
            if x in rules:  # total time derivative of the rule expression
                e = inline(expand(rules[x]))
                ds = [(e.derivative(v), v) for v in _names(e, set())]
                if any(d is None for d, _ in ds):
                    raise NotImplementedError("rateOf(%s): libsbml cannot differentiate its assignment rule" % x)
                terms = [_l3(d) + ('' if v == 'time' else ' * __rateof__' + v) for d, v in ds]
                return formula(libsbml.parseL3Formula(' + '.join(terms) or '0'), "rateOf(%s)" % x)
            return '0'
        return lazy(('d', x), f)

    def dconc(s):
        """d[S]/dt (concentration / time) of floating species s."""
        if s in rate_rules:
            return conc_of(s, formula(rate_rules[s], "rate rule for %s" % s))
        i = rows.index(s)
        terms = []
        for j, rid in enumerate(reactions):
            if (s, rid) in vstoich:
                n = ' + '.join('%d * (%s)' % (sg, formula(a, "stoichiometry of %s in %s" % (s, rid)))
                               for sg, a in vstoich[(s, rid)])
            elif Nm[i, j] != 0:
                n = repr(float(Nm[i, j]))
            else:
                continue
            terms.append('(%s) * (%s)' % (n, rate(j)))
        flux = ' + '.join(terms) or '0'
        if cf_of(s):
            flux = '(%s) * (%s)' % (formula(_name(cf_of(s)), "conversion factor of %s" % s), flux)
        V = formula(_name(comp_of[s]), "compartment of %s" % s)
        dV = rate_of(comp_of[s])
        return '(%s) / (%s)' % (flux, V) + ('' if dV == '0' else ' - %s * (%s) / (%s)' % (s, dV, V))

    def conc_of(var, f):  # value for species var in its SBML units -> state units (concentration)
        if var in hosu and (var in floating or var in rules):
            return '(%s) / (%s)' % (f, formula(_name(comp_of[var]), "compartment of %s" % var))
        return f

    rates = [rate(j) for j in range(len(reactions))]
    stoich = np.array(Nm, dtype=float)
    for i, s in enumerate(rows):
        if cf_of(s) and s not in general:  # constant conversion factor: scale the stoichiometry row
            stoich[i] *= params[cf_of(s)]
    stoich = stoich[[rows.index(s) for s in species]]
    rate_list = []
    for v in rate_rules:
        if v not in floating:
            rate_list.append((v, formula(rate_rules[v], "rate rule for %s" % v)))
    for i, s in enumerate(species):
        if s in rate_rules or s in general:
            stoich[i] = 0.0
            rate_list.append((s, dconc(s)))
    out_rules = [(v, conc_of(v, formula(f, "assignment rule for %s" % v))) for v, f in rules.items()]

    events = []
    for ev in m.getListOfEvents():
        trig = ev.getTrigger()
        what = "event %s" % (ev.getId() or len(events))
        asg = []
        for ea in ev.getListOfEventAssignments():
            v = ea.getVariable()
            f = formula(ea.getMath(), what)
            asg.append((v, '(%s) * (%s)' % (f, formula(_name(comp_of[v]), what)) if v in bamt else conc_of(v, f)))
        # a compartment resized directly or through its assignment rule keeps its species' amounts
        # ponytail: with useValuesFromTriggerTime the amount is also taken at trigger time
        new = dict(asg)
        pat = re.compile(r'\b(%s)\b' % '|'.join(map(re.escape, new)))
        for c in cvar:
            old = formula(_name(c), what)
            cnew = pat.sub(lambda g: '(%s)' % new[g.group(1)], old) if new else old
            if cnew != old:  # an assigned concentration is converted with the old volume, as roadrunner
                for s in species:
                    if comp_of[s] == c:
                        new[s] = '(%s) * (%s) / (%s)' % (new.get(s, s), old, cnew)
        asg = list(new.items())
        events.append((formula(trig.getMath(), what), asg,
                       formula(ev.getDelay().getMath(), what) if ev.isSetDelay() else None,
                       formula(ev.getPriority().getMath(), what) if ev.isSetPriority() else None,
                       bool(trig.getInitialValue()), bool(trig.getPersistent()),
                       bool(ev.getUseValuesFromTriggerTime())))

    if info is not None:
        info.update(comp_of=comp_of, general=general, changing=changing,
                    amount_params={s for s in sp if s in params and (s in hosu or s in bamt)})
    pn = list(params)
    return roadrunner_rs.Model(species, [conc[s] for s in species], [vol[s] for s in species], reactions,
                               rates, np.ascontiguousarray(stoich), pn,
                               [float(params[p]) for p in pn], rate_rules=rate_list,
                               assignment_rules=out_rules, events=events)

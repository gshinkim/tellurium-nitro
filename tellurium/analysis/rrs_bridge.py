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
    odoc = libsbml.readSBMLFromString(r.getSBML())  # getCurrentSBML() drops stoichiometryMath
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

    def val(x):  # value in SBML math units (species: concentration, amount if hasOnlySubstanceUnits)
        if x in sp:
            v = rr.getValue('[%s]' % x)
            return v * rr.getValue(sp[x].getCompartment()) if sp[x].getHasOnlySubstanceUnits() else v
        return rr.getValue(x)

    def sref(sr, changing):
        return rr.getValue(sr.getId()) if sr.getId() in changing else sr.getStoichiometry()
    return _compile(m, odoc.getModel(), rows, reactions, Nm, conc, val, sref, info)


def _compile(m, om, rows, reactions, Nm, conc, val, sref, info):
    """Rust model from SBML model ``m`` (``om``: same model with stoichiometryMath), stoichiometry ``Nm``
    (rows: floating species incl. rule ones, columns: reactions), floating species concentrations
    ``conc`` and values ``val(id)`` / ``sref(speciesReference, changing)`` in SBML units."""
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    rules, rate_rules = {}, {}
    for rule in m.getListOfRules():
        (rate_rules if rule.isRate() else rules)[rule.getVariable()] = rule.getMath()
    targets = {ea.getVariable() for ev in m.getListOfEvents() for ea in ev.getListOfEventAssignments()}
    changing = set(rules) | set(rate_rules) | targets

    sp = {s.getId(): s for s in m.getListOfSpecies()}
    species = [s for s in rows if s not in rules]  # assignment-rule species: outputs only
    floating = set(species)
    comp_of = {s: sp[s].getCompartment() for s in sp}
    hosu = {s for s in sp if sp[s].getHasOnlySubstanceUnits()}
    vol = {s: val(comp_of[s]) for s in sp}

    params = {}
    for x in list(m.getListOfParameters()) + list(m.getListOfCompartments()):
        if x.getId() not in rules:
            params[x.getId()] = val(x.getId())
    for s in sp:  # boundary / non-floating species: parameters in their SBML units
        if s not in rules and s not in floating:
            params[s] = val(s)

    # variable stoichiometry: (species, reaction) -> [(sign, stoichiometry AST)] over all references
    def refs_of(model):
        return [(sg, rxn.getId(), sr) for rxn in model.getListOfReactions()
                for sg, lst in ((-1, rxn.getListOfReactants()), (1, rxn.getListOfProducts())) for sr in lst]
    refs = refs_of(m)
    smath = {(j, sr.getSpecies(), sg, i): sr.getStoichiometryMath().getMath()
             for i, (sg, j, sr) in enumerate(refs_of(om))
             if sr.isSetStoichiometryMath() and sr.getStoichiometryMath().isSetMath()}
    smath = {i: smath.get((j, sr.getSpecies(), sg, i)) for i, (sg, j, sr) in enumerate(refs)}
    for _, _, sr in refs:
        if sr.isSetId() and sr.getId() not in rules:
            params[sr.getId()] = sref(sr, changing)
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


# ---------------------------------------------------------------------------------------------
# Without roadrunner: initial values and stoichiometry from the SBML alone
# ---------------------------------------------------------------------------------------------
_A = libsbml
_MATH = {_A.AST_FUNCTION_ABS: np.abs, _A.AST_FUNCTION_EXP: np.exp, _A.AST_FUNCTION_LN: np.log,
         _A.AST_FUNCTION_FLOOR: np.floor, _A.AST_FUNCTION_CEILING: np.ceil, _A.AST_FUNCTION_SIN: np.sin,
         _A.AST_FUNCTION_COS: np.cos, _A.AST_FUNCTION_TAN: np.tan, _A.AST_FUNCTION_ARCSIN: np.arcsin,
         _A.AST_FUNCTION_ARCCOS: np.arccos, _A.AST_FUNCTION_ARCTAN: np.arctan, _A.AST_FUNCTION_SINH: np.sinh,
         _A.AST_FUNCTION_COSH: np.cosh, _A.AST_FUNCTION_TANH: np.tanh}
_REL = {_A.AST_RELATIONAL_EQ: np.equal, _A.AST_RELATIONAL_NEQ: np.not_equal, _A.AST_RELATIONAL_GT: np.greater,
        _A.AST_RELATIONAL_LT: np.less, _A.AST_RELATIONAL_GEQ: np.greater_equal, _A.AST_RELATIONAL_LEQ: np.less_equal}


def evaluate(ast, get, fdefs):
    """Value of a libsbml AST at time 0; names through get(name). NotImplementedError if unsupported."""
    t, n = ast.getType(), ast.getNumChildren()

    def e(i):
        return evaluate(ast.getChild(i), get, fdefs)
    with np.errstate(all='ignore'):
        if ast.isNumber() or t in (_A.AST_CONSTANT_E, _A.AST_CONSTANT_PI):
            return np.float64(ast.getValue())
        if t in (_A.AST_CONSTANT_TRUE, _A.AST_CONSTANT_FALSE):
            return np.float64(t == _A.AST_CONSTANT_TRUE)
        if t == _A.AST_NAME:
            return np.float64(get(ast.getName()))
        if t == _A.AST_NAME_TIME:
            return np.float64(0.0)
        if t == _A.AST_NAME_AVOGADRO:
            return np.float64(ast.getReal())
        if t == _A.AST_PLUS:
            return sum((e(i) for i in range(n)), np.float64(0.0))
        if t == _A.AST_MINUS:
            return -e(0) if n == 1 else e(0) - e(1)
        if t == _A.AST_TIMES:
            return np.prod([e(i) for i in range(n)], dtype=float)
        if t == _A.AST_DIVIDE:
            return e(0) / e(1)
        if t in (_A.AST_POWER, _A.AST_FUNCTION_POWER):
            return np.power(e(0), e(1))
        if t == _A.AST_FUNCTION_ROOT:
            return np.sqrt(e(0)) if n == 1 else np.power(e(1), 1.0 / e(0))
        if t == _A.AST_FUNCTION_LOG:
            return np.log10(e(0)) if n == 1 else np.log(e(1)) / np.log(e(0))
        if t in _MATH and n == 1:
            return _MATH[t](e(0))
        if t in _REL:
            return np.float64(all(_REL[t](e(i), e(i + 1)) for i in range(n - 1)))
        if t == _A.AST_LOGICAL_AND:
            return np.float64(all(e(i) for i in range(n)))
        if t == _A.AST_LOGICAL_OR:
            return np.float64(any(e(i) for i in range(n)))
        if t == _A.AST_LOGICAL_XOR:
            return np.float64(sum(bool(e(i)) for i in range(n)) % 2)
        if t == _A.AST_LOGICAL_NOT:
            return np.float64(not e(0))
        if t == _A.AST_FUNCTION_PIECEWISE:
            for i in range(0, n - 1, 2):
                if e(i + 1):
                    return e(i)
            return e(n - 1) if n % 2 else np.float64(np.nan)
        if t == _A.AST_FUNCTION and ast.getName() in fdefs:
            lam = fdefs[ast.getName()]
            args = {lam.getChild(i).getName(): e(i) for i in range(lam.getNumChildren() - 1)}
            return evaluate(lam.getChild(lam.getNumChildren() - 1), lambda x: args[x] if x in args else get(x), fdefs)
    raise NotImplementedError("initial value of %s" % _l3(ast))


def _rr_order(xs, ruled):  # roadrunner lists rule-determined symbols after the others
    return [x for x in xs if x not in ruled] + [x for x in xs if x in ruled]


def from_sbml(sbml, info=None):
    """``from_roadrunner`` without roadrunner: compile SBML text into a ``roadrunner_rs.Model`` at the
    SBML initial values (initial assignments, assignment rules, amounts / concentrations as roadrunner).

    ``info`` also receives the model in roadrunner's terms: ids (floating, boundary, global parameter,
    compartment, reaction, rate-rule ids in roadrunner's order), rows/N (full stoichiometry), ia
    (initial assignments), fdefs, values (initial values, SBML units), doc (libsbml document).
    NotImplementedError / ValueError: not supported (load roadrunner instead)."""
    doc = libsbml.readSBMLFromString(sbml)
    m = doc.getModel()
    if m is None or doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) or doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL):
        raise ValueError("SBML errors")
    ns = doc.getNamespaces()
    if any(p in ns.getURI(i) for i in range(ns.getNumNamespaces()) for p in ('/comp/', '/distrib/', '/arrays/')):
        raise NotImplementedError("SBML package")  # roadrunner flattens comp models
    if any(r.isAlgebraic() for r in m.getListOfRules()) or any(x.getFast() for x in m.getListOfReactions()) \
            or 'symbols/delay' in sbml:
        raise NotImplementedError("algebraic rule, fast reaction or delay")
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    ia = {x.getSymbol(): x.getMath() for x in m.getListOfInitialAssignments()}
    ar = {x.getVariable(): x.getMath() for x in m.getListOfRules() if x.isAssignment()}
    rrules = [x.getVariable() for x in m.getListOfRules() if x.isRate()]
    sp = {s.getId(): s for s in m.getListOfSpecies()}
    comps = {c.getId(): c for c in m.getListOfCompartments()}
    pars = {p.getId(): p for p in m.getListOfParameters()}
    refs = {sr.getId(): sr for rx in m.getListOfReactions()
            for lst in (rx.getListOfReactants(), rx.getListOfProducts()) for sr in lst if sr.isSetId()}
    vals, busy = {}, set()

    def get(x):  # initial value in SBML units, as roadrunner computes it
        if x in vals:
            return vals[x]
        if x in busy:
            raise ValueError("cyclic initial values through %s" % x)
        busy.add(x)
        if x in ar or x in ia:
            v = evaluate(ar[x] if x in ar else ia[x], get, fdefs)
        elif x in comps and comps[x].isSetSize():
            v = comps[x].getSize()
        elif x in pars and pars[x].isSetValue():
            v = pars[x].getValue()
        elif x in sp and sp[x].isSetInitialConcentration():
            v = sp[x].getInitialConcentration()
            if sp[x].getHasOnlySubstanceUnits():
                v *= get(sp[x].getCompartment())
        elif x in sp and sp[x].isSetInitialAmount():
            v = sp[x].getInitialAmount()
            if not sp[x].getHasOnlySubstanceUnits():
                v /= get(sp[x].getCompartment())
        elif x in refs and not refs[x].isSetStoichiometryMath() and refs[x].isSetStoichiometry():
            v = refs[x].getStoichiometry()
        elif x in comps or x in sp:  # roadrunner's defaults
            v = 1.0 if x in comps else 0.0
        else:  # unset parameter (roadrunner refuses it), reaction ids in initial assignments
            raise NotImplementedError("initial value of %s" % x)
        busy.discard(x)
        vals[x] = float(v)
        return vals[x]

    def stoich(sr):
        if sr.isSetId():
            return get(sr.getId())
        if sr.isSetStoichiometryMath():
            return float(evaluate(sr.getStoichiometryMath().getMath(), get, fdefs))
        if sr.isSetStoichiometry() or m.getLevel() < 3:
            return sr.getStoichiometry()
        raise NotImplementedError("unset stoichiometry")

    rows = [s for s in sp if not sp[s].getBoundaryCondition()]
    reactions = [x.getId() for x in m.getListOfReactions()]
    N = np.zeros((len(rows), len(reactions)))
    ri = {s: i for i, s in enumerate(rows)}
    for j, rx in enumerate(m.getListOfReactions()):
        for sg, lst in ((-1, rx.getListOfReactants()), (1, rx.getListOfProducts())):
            for sr in lst:
                if sr.getSpecies() in ri:
                    N[ri[sr.getSpecies()], j] += sg * stoich(sr)
    for x in list(comps) + list(pars) + list(sp):
        get(x)
    conc = {s: vals[s] / vals[sp[s].getCompartment()] if sp[s].getHasOnlySubstanceUnits() else vals[s] for s in rows}
    out = _compile(m, m, rows, reactions, N, conc, get, lambda sr, changing: get(sr.getId()), info)
    if info is not None:
        ruled = set(ar) | set(rrules)
        info.update(doc=doc, ia=ia, fdefs=fdefs, values=vals, rows=rows, N=N, rate_rule_ids=rrules,
                    ids=(_rr_order(rows, ruled), _rr_order([s for s in sp if s not in ri], ruled),
                         _rr_order(list(pars), ruled), _rr_order(list(comps), ruled), reactions))
    return out

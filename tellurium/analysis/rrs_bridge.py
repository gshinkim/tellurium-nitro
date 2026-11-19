"""Build a ``roadrunner_rs.Model`` (Rust steady-state/MCA/simulation engine) from a roadrunner model.

    m = from_roadrunner(r)
    res = m.mca_batch(['k1'], np.array([[0.5], [1.0]]))
"""
import math
import os
import re

import libsbml
import numpy as np
import roadrunner
import roadrunner_rs


def _type(ast):
    """ast.getType() as an int. Importing libsedml after libsbml makes libsbml's SWIG binding return
    the enum as a pointer object; read the value through it then."""
    t = ast.getType()
    if isinstance(t, int):
        return t
    import ctypes
    return ctypes.c_int.from_address(int(t)).value


def _names(ast, out):
    if ast.isName():
        out.add(ast.getName())
    for i in range(ast.getNumChildren()):
        _names(ast.getChild(i), out)
    return out


_OUT = libsbml.L3ParserSettings()
_OUT.setParseUnits(False)  # no '1 dimensionless' in formulas


def _l3(ast):
    """Infix text; constants spelled so no model id can shadow them (a parameter may be called INF or pi)."""
    def plain(a):
        if a.isReal() and not np.isfinite(a.getReal()) or a.isConstant() and a.getName() in ('pi', 'exponentiale', 'true', 'false'):
            return True
        return any(plain(a.getChild(i)) for i in range(a.getNumChildren()))
    if not plain(ast):
        return libsbml.formulaToL3StringWithSettings(ast, _OUT)
    ast = ast.deepCopy()

    def fix(a):
        for i in range(a.getNumChildren()):
            c = a.getChild(i)
            v = None
            if c.isReal() and not np.isfinite(c.getReal()):
                x = c.getReal()
                v = '(0/0)' if np.isnan(x) else '(1/0)' if x > 0 else '(-1/0)'
            elif c.isConstant() and c.getName() in ('pi', 'exponentiale', 'true', 'false'):
                v = repr({'pi': math.pi, 'exponentiale': math.e, 'true': 1.0, 'false': 0.0}[c.getName()])
            if v is None:
                fix(c)
            else:
                a.replaceChild(i, libsbml.parseL3Formula(v), True)
    root = libsbml.ASTNode(libsbml.AST_PLUS)  # a parent, so the root itself can be replaced
    root.addChild(ast)
    fix(root)
    return libsbml.formulaToL3StringWithSettings(root.getChild(0), _OUT)


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
    m = _drop_empty_math(doc.getModel())
    odoc = libsbml.readSBMLFromString(r.getSBML())  # getCurrentSBML() drops stoichiometryMath
    _drop_empty_math(odoc.getModel())
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


def _compile(m, om, rows, reactions, Nm, conc, val, sref, info, dae=(), pre=None):
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

    local, kl_math, delays = {}, {}, {}
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
        t = _type(ast)
        if t == libsbml.AST_NAME_AVOGADRO:
            return libsbml.parseL3Formula(repr(ast.getReal()))
        if t == libsbml.AST_NAME_TIME:
            return _name('time')
        if t == libsbml.AST_FUNCTION_RATE_OF:
            return _name('__rateof__' + ast.getChild(0).getName())
        if t == libsbml.AST_FUNCTION_DELAY:  # a parameter the engine refills with the delayed value
            key = (formula(ast.getChild(0), "delay()"), formula(ast.getChild(1), "delay()"),
                   formula(before_start(ast.getChild(0)), "delay()"))
            if key not in delays:
                delays[key] = '__delay_%d__' % len(delays)
                params[delays[key]] = 0.0
                known.add(delays[key])
            return _name(delays[key])
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

    ias = {a.getSymbol(): a.getMath() for a in m.getListOfInitialAssignments()}

    def before_start(ast):
        """delay() before the run: everything holds its initial value (rateOf 0), and a symbol with an
        initial assignment is that formula at the earlier time."""
        out = ast.deepCopy()

        def zero_rates(a):
            for i in range(a.getNumChildren()):
                c = a.getChild(i)
                if _type(c) == libsbml.AST_FUNCTION_RATE_OF:
                    a.replaceChild(i, libsbml.parseL3Formula('0'), True)
                else:
                    zero_rates(c)
        if _type(out) == libsbml.AST_FUNCTION_RATE_OF:
            return libsbml.parseL3Formula('0')
        zero_rates(out)
        pre_ = pre or {}  # values before the run that differ from t0 (fast species before equilibrating)
        for x, f in ias.items():
            if x in _names(out, set()) and x not in pre_:
                out = _sub(out, x, _l3(f))
        for x, f in pre_.items():
            if x in _names(out, set()):
                out = _sub(out, x, f)
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
        asg, div = [], {}
        uvftt = ev.getUseValuesFromTriggerTime()
        for ea in ev.getListOfEventAssignments():
            v = ea.getVariable()
            f = formula(ea.getMath(), what)
            if uvftt and ev.isSetDelay() and v in hosu and v in floating:
                asg.append((v, f))  # an amount fixed at trigger time: / volume when it executes
                div[v] = formula(_name(comp_of[v]), what)
                continue
            asg.append((v, '(%s) * (%s)' % (f, formula(_name(comp_of[v]), what)) if v in bamt else conc_of(v, f)))
        # a compartment resized directly or through its assignment rule keeps its species' amounts: the
        # engine scales their concentrations by V_before / V_after when the event executes (an assigned
        # concentration is so converted with the old volume, as roadrunner)
        new = dict(asg)
        pat = re.compile(r'\b(%s)\b' % '|'.join(map(re.escape, new)))
        rescale = []
        for c in cvar:
            old = formula(_name(c), what)
            if new and pat.sub(lambda g: '(%s)' % new[g.group(1)], old) != old:
                held = [x for x in rate_rules if x in sp and x not in floating and x not in hosu]  # conc. states
                rescale += [(s, old) for s in species + held if comp_of[s] == c]
        events.append((formula(trig.getMath(), what), asg,
                       formula(ev.getDelay().getMath(), what) if ev.isSetDelay() else None,
                       formula(ev.getPriority().getMath(), what) if ev.isSetPriority() else None,
                       bool(trig.getInitialValue()), bool(trig.getPersistent()),
                       bool(uvftt), [div.get(v) for v, _ in asg], rescale))

    if info is not None:
        info.update(comp_of=comp_of, general=general, changing=changing,
                    amount_params={s for s in sp if s in params and (s in hosu or s in bamt)})
    pn = list(params)
    return roadrunner_rs.Model(species, [conc[s] for s in species], [vol[s] for s in species], reactions,
                               rates, np.ascontiguousarray(stoich), pn,
                               [float(params[p]) for p in pn], rate_rules=rate_list,
                               assignment_rules=out_rules, events=events,
                               delays=[(n, e, tau, pre) for (e, tau, pre), n in delays.items()],
                               algebraic=sorted(dae))


# ---------------------------------------------------------------------------------------------
# Without roadrunner: initial values and stoichiometry from the SBML alone
# ---------------------------------------------------------------------------------------------
_A = libsbml
_MATH = {_A.AST_FUNCTION_ABS: np.abs, _A.AST_FUNCTION_EXP: np.exp, _A.AST_FUNCTION_LN: np.log,
         _A.AST_FUNCTION_FLOOR: np.floor, _A.AST_FUNCTION_CEILING: np.ceil, _A.AST_FUNCTION_SIN: np.sin,
         _A.AST_FUNCTION_COS: np.cos, _A.AST_FUNCTION_TAN: np.tan, _A.AST_FUNCTION_ARCSIN: np.arcsin,
         _A.AST_FUNCTION_ARCCOS: np.arccos, _A.AST_FUNCTION_ARCTAN: np.arctan, _A.AST_FUNCTION_SINH: np.sinh,
         _A.AST_FUNCTION_COSH: np.cosh, _A.AST_FUNCTION_TANH: np.tanh,
         _A.AST_FUNCTION_SEC: lambda x: 1 / np.cos(x), _A.AST_FUNCTION_CSC: lambda x: 1 / np.sin(x),
         _A.AST_FUNCTION_COT: lambda x: 1 / np.tan(x), _A.AST_FUNCTION_SECH: lambda x: 1 / np.cosh(x),
         _A.AST_FUNCTION_CSCH: lambda x: 1 / np.sinh(x), _A.AST_FUNCTION_COTH: lambda x: 1 / np.tanh(x),
         _A.AST_FUNCTION_ARCSEC: lambda x: np.arccos(1 / x), _A.AST_FUNCTION_ARCCSC: lambda x: np.arcsin(1 / x),
         _A.AST_FUNCTION_ARCCOT: lambda x: np.arctan(1 / x), _A.AST_FUNCTION_ARCSINH: np.arcsinh,
         _A.AST_FUNCTION_ARCCOSH: np.arccosh, _A.AST_FUNCTION_ARCTANH: np.arctanh,
         _A.AST_FUNCTION_ARCSECH: lambda x: np.arccosh(1 / x), _A.AST_FUNCTION_ARCCSCH: lambda x: np.arcsinh(1 / x),
         _A.AST_FUNCTION_ARCCOTH: lambda x: np.arctanh(1 / x),
         _A.AST_FUNCTION_FACTORIAL: lambda x: np.float64(math.factorial(int(round(x)))) if x >= 0 else np.nan}
_REL = {_A.AST_RELATIONAL_EQ: np.equal, _A.AST_RELATIONAL_NEQ: np.not_equal, _A.AST_RELATIONAL_GT: np.greater,
        _A.AST_RELATIONAL_LT: np.less, _A.AST_RELATIONAL_GEQ: np.greater_equal, _A.AST_RELATIONAL_LEQ: np.less_equal}


def evaluate(ast, get, fdefs):
    """Value of a libsbml AST at time 0; names through get(name). NotImplementedError if unsupported."""
    if ast.isConstant():  # by name: libsbml's SWIG binding can return an unusable type for pi
        c = {'pi': np.pi, 'exponentiale': np.e, 'true': 1.0, 'false': 0.0}.get(ast.getName())
        if c is not None:
            return np.float64(c)
    t, n = _type(ast), ast.getNumChildren()

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
        if t in (_A.AST_FUNCTION_QUOTIENT, _A.AST_FUNCTION_REM):  # C truncation, as libroadrunner
            a, b = e(0), e(1)
            q = np.trunc(a / b)
            return q if t == _A.AST_FUNCTION_QUOTIENT else a - b * q
        if t == _A.AST_LOGICAL_IMPLIES:
            return np.float64(not e(0) or bool(e(1)))
        if t in (_A.AST_FUNCTION_MIN, _A.AST_FUNCTION_MAX) and n:
            return (min if t == _A.AST_FUNCTION_MIN else max)(e(i) for i in range(n))
        if t == _A.AST_FUNCTION_RATE_OF:
            return np.float64(get('__rateof__' + ast.getChild(0).getName()))
        if t == _A.AST_FUNCTION_DELAY:  # before the start, every value is its initial one
            return e(0)
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


def _drop_empty_math(m):
    """SBML L3V2 lets math be missing; such an element has no effect: drop rules, initial and event
    assignments without math, events whose trigger has none, and delays/priorities/stoichiometryMath without it."""
    for lst in (m.getListOfRules(), m.getListOfInitialAssignments(), m.getListOfConstraints()):
        for i in reversed(range(lst.size())):
            if lst.get(i).getMath() is None:
                lst.remove(i)
    evs = m.getListOfEvents()
    for i in reversed(range(evs.size())):
        ev = evs.get(i)
        if ev.getTrigger() is None or ev.getTrigger().getMath() is None:
            evs.remove(i)
            continue
        if ev.isSetDelay() and ev.getDelay().getMath() is None:
            ev.unsetDelay()
        if ev.isSetPriority() and ev.getPriority().getMath() is None:
            ev.unsetPriority()
        eas = ev.getListOfEventAssignments()
        for j in reversed(range(eas.size())):
            if eas.get(j).getMath() is None:
                eas.remove(j)
    for rx in m.getListOfReactions():
        for lst in (rx.getListOfReactants(), rx.getListOfProducts()):
            for sr in lst:
                if sr.isSetStoichiometryMath() and sr.getStoichiometryMath().getMath() is None:
                    sr.unsetStoichiometryMath()
    return m


def _rr_order(xs, ruled):  # roadrunner lists rule-determined symbols after the others
    return [x for x in xs if x not in ruled] + [x for x in xs if x in ruled]


def _flatten_comp(doc, base=None):
    """Flatten a hierarchical (comp) model in place, as roadrunner does; ``base``: directory that
    external model definitions are relative to."""
    if base:
        doc.setLocationURI('file://' + os.path.join(os.path.abspath(base), ''))
    props = libsbml.ConversionProperties()
    props.addOption('flatten comp', True)
    if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS:
        errs = [doc.getError(i).getMessage().strip() for i in range(doc.getNumErrors())]
        raise NotImplementedError("comp model could not be flattened: %s" % '; '.join(errs[:3]))


def _sub(ast, x, f):
    """Copy of ast with symbol x replaced by the formula f (an L3 string)."""
    out = ast.deepCopy()
    new = libsbml.parseL3Formula('(%s)' % f)
    if out.isName() and out.getName() == x:
        return new
    out.replaceIDWithFunction(x, new)
    return out


def _free_symbols(m):
    """Symbols an algebraic rule may determine: not constant, not set by another rule, not a floating
    species that reactions change."""
    taken = {r.getVariable() for r in m.getListOfRules() if not r.isAlgebraic()} | \
        {a.getVariable() for ev in m.getListOfEvents() for a in ev.getListOfEventAssignments()}
    moved = {sr.getSpecies() for rx in m.getListOfReactions()
             for lst in (rx.getListOfReactants(), rx.getListOfProducts()) for sr in lst}
    free = set()
    for x in m.getListOfParameters():
        if not x.getConstant():
            free.add(x.getId())
    for x in m.getListOfCompartments():
        if not x.getConstant():
            free.add(x.getId())
    for x in m.getListOfSpecies():
        if not x.getConstant() and (x.getBoundaryCondition() or x.getId() not in moved):
            free.add(x.getId())
    for rx in m.getListOfReactions():
        for lst in (rx.getListOfReactants(), rx.getListOfProducts()):
            for sr in lst:
                if sr.isSetId() and not sr.getConstant():
                    free.add(sr.getId())
    return free - taken


def _left_null(rows, cols, a):
    """Exact basis (Fractions) of {y : y a = 0} for the rows x cols matrix a (list of lists)."""
    from fractions import Fraction
    t = [[Fraction(a[i][j]) for i in range(rows)] for j in range(cols)]  # a^T: solve a^T y = 0
    piv, r = [], 0
    for c in range(rows):
        p = next((i for i in range(r, cols) if t[i][c] != 0), None)
        if p is None:
            continue
        t[r], t[p] = t[p], t[r]
        t[r] = [x / t[r][c] for x in t[r]]
        for i in range(cols):
            if i != r and t[i][c] != 0:
                t[i] = [x - t[i][c] * y for x, y in zip(t[i], t[r])]
        piv.append(c)
        r += 1
    basis = []
    for f in (c for c in range(rows) if c not in piv):
        y = [Fraction(0)] * rows
        y[f] = Fraction(1)
        for i, c in enumerate(piv):
            y[c] = -t[i][f]
        basis.append(y)
    return basis, len(piv)


def _fast_to_algebraic(m, pre=None):  # returns the fast species, which the new algebraic rules determine
    """SBML fast reactions are at equilibrium at every instant. Rewrite them as: a rate-rule parameter
    per conserved total of the fast species' amounts (driven by the slow reactions), algebraic rules
    v_fast = 0 and sum(l_i n_i) = total, and the fast species as boundary species (rule-set)."""
    fast = [rx for rx in m.getListOfReactions() if rx.getFast()]
    if not fast:
        return set()
    # a fast reaction whose rate is 0 for any values (e.g. k * A * ln(1)) never moves anything: drop it
    rng = np.random.default_rng(2)
    fd = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    rnd = {}
    dead = set()
    for rx in fast:
        kl = rx.getKineticLaw()
        if kl is None or kl.getMath() is None:
            continue
        try:
            vals = [evaluate(kl.getMath(), lambda x: rnd.setdefault(x, rng.uniform(0.5, 2.0)), fd) for _ in range(2)]
        except NotImplementedError:
            continue
        if all(v == 0 for v in vals):
            dead.add(rx.getId())
    if dead:
        lst = m.getListOfReactions()
        for i in reversed([i for i in range(lst.size()) if lst.get(i).getId() in dead]):
            lst.remove(i)
        fast = [rx for rx in m.getListOfReactions() if rx.getFast()]
        if not fast:
            return set()
    sp = {s.getId(): s for s in m.getListOfSpecies()}
    fs = []
    for rx in fast:
        for lst in (rx.getListOfReactants(), rx.getListOfProducts()):
            for sr in lst:
                x = sr.getSpecies()
                if not sp[x].getBoundaryCondition() and not sp[x].getConstant() and x not in fs:
                    fs.append(x)
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    ias = {a.getSymbol(): a for a in m.getListOfInitialAssignments()}
    moving = {r.getVariable() for r in m.getListOfRules()} | \
        {a.getVariable() for ev in m.getListOfEvents() for a in ev.getListOfEventAssignments()}

    def const(x):  # value of a symbol fixed for the whole run (stoichiometries, conversion factors)
        if x in moving:
            raise NotImplementedError("fast reaction with stoichiometry or conversion factor %s that changes" % x)
        if x in ias:
            return float(evaluate(ias[x].getMath(), const, fdefs))
        e = m.getElementBySId(x)
        if isinstance(e, libsbml.SpeciesReference):
            return e.getStoichiometry() if e.isSetStoichiometry() else 1.0
        if isinstance(e, libsbml.Parameter) and e.isSetValue():
            return e.getValue()
        if isinstance(e, libsbml.Compartment) and e.isSetSize():
            return e.getSize()
        raise NotImplementedError("fast reaction: no constant value for %s" % x)
    if any(sr.isSetStoichiometryMath() for rx in m.getListOfReactions()
           for lst in (rx.getListOfReactants(), rx.getListOfProducts()) for sr in lst if sr.getSpecies() in fs):
        raise NotImplementedError("fast reaction with stoichiometryMath")

    def n(rx, x):  # net stoichiometry of species x in reaction rx
        return sum(sg * (const(sr.getId()) if sr.isSetId() else sr.getStoichiometry())
                   for sg, lst in ((-1, rx.getListOfReactants()), (1, rx.getListOfProducts()))
                   for sr in lst if sr.getSpecies() == x)

    def cf(x):  # conversion factor id of species x, or None
        return sp[x].getConversionFactor() if sp[x].isSetConversionFactor() else (m.getConversionFactor() or None)
    for x in fs:
        if cf(x):
            const(cf(x))
    basis, rank = _left_null(len(fs), len(fast), [[n(rx, x) for rx in fast] for x in fs])
    if rank < len(fast):
        raise NotImplementedError("dependent fast reactions")
    # amounts over the conversion factor: reactions change these by n * rate, so totals are conserved
    amt = {x: (x if sp[x].getHasOnlySubstanceUnits() else '%s * %s' % (x, sp[x].getCompartment()))
           + (' / %s' % cf(x) if cf(x) else '') for x in fs}

    def amt0(x):  # x's initial amount (over its conversion factor), from its own initial value
        s, c = sp[x], sp[x].getCompartment()
        if x in ias:
            v = '(%s)' % _l3(ias[x].getMath())
            v = v if s.getHasOnlySubstanceUnits() else '%s * %s' % (v, c)
        elif s.isSetInitialAmount():
            v = repr(s.getInitialAmount())
        else:
            v = '%s * %s' % (repr(s.getInitialConcentration() if s.isSetInitialConcentration() else 0.0), c)
        return '(%s) / %s' % (v, cf(x)) if cf(x) else v
    slow = [rx for rx in m.getListOfReactions() if not rx.getFast()]
    for k, y in enumerate(basis):
        tid = '__fast_total_%d' % k
        while m.getElementBySId(tid) is not None:
            tid += '_'
        p = m.createParameter()
        p.setId(tid)
        p.setConstant(False)
        tot = ' + '.join('(%s) * (%s)' % (float(c), amt[x]) for c, x in zip(y, fs) if c)
        ia = m.createInitialAssignment()
        ia.setSymbol(tid)
        ia.setMath(libsbml.parseL3Formula(' + '.join('(%s) * (%s)' % (float(c), amt0(x)) for c, x in zip(y, fs) if c)))
        d = ' + '.join('(%s) * %s' % (float(c * n(rx, x)), rx.getId()) for rx in slow for c, x in zip(y, fs) if c and n(rx, x))
        rr = m.createRateRule()
        rr.setVariable(tid)
        rr.setMath(libsbml.parseL3Formula(d or '0'))
        ar = m.createAlgebraicRule()
        ar.setMath(libsbml.parseL3Formula('%s - %s' % (tot, tid)))
    for rx in fast:
        ar = m.createAlgebraicRule()
        kl = rx.getKineticLaw()
        if kl is None or kl.getMath() is None:
            raise NotImplementedError("fast reaction %s without kinetic law" % rx.getId())
        math = kl.getMath().deepCopy()
        for q in list(kl.getListOfParameters()) + list(kl.getListOfLocalParameters()):
            math = _sub(math, q.getId(), repr(q.getValue()))
        ar.setMath(math)
    # starting guess where an equilibrium stays a DAE constraint (e.g. k*A*B = 0 has two roots): every
    # fast reaction run to completion in the direction of its initial rate
    def conc_from(x, u):  # species value (SBML units) from an amount over the conversion factor
        a = '(%s) * %s' % (u, cf(x)) if cf(x) else u
        return a if sp[x].getHasOnlySubstanceUnits() else '(%s) / %s' % (a, sp[x].getCompartment())
    conc0 = {x: conc_from(x, amt0(x)) for x in fs}
    if pre is not None:  # delay() before the run sees the values before they equilibrate
        pre.update(conc0)
    ext = {}
    for rx in fast:
        v0 = rx.getKineticLaw().getMath().deepCopy()
        for q in list(rx.getKineticLaw().getListOfParameters()) + list(rx.getKineticLaw().getListOfLocalParameters()):
            v0 = _sub(v0, q.getId(), repr(q.getValue()))
        for x in fs:
            v0 = _sub(v0, x, conc0[x])
        v0 = _l3(v0)
        lim = lambda sg: ['(%s) / %s' % (amt0(x), abs(n(rx, x))) for x in fs if sg * n(rx, x) < 0]
        fwd, bwd = lim(1), lim(-1)
        ext[rx.getId()] = 'piecewise(%s, (%s) > 0, %s, (%s) < 0, 0)' % (
            'min(%s)' % ', '.join(fwd) if fwd else '0', v0, '-min(%s)' % ', '.join(bwd) if bwd else '0', v0)
    for x in fs:
        sp[x].setBoundaryCondition(True)
        if x in ias:  # now set by a rule
            m.getListOfInitialAssignments().remove(x)
        guess = '%s + %s' % (amt0(x), ' + '.join('(%s) * (%s)' % (n(rx, x), ext[rx.getId()]) for rx in fast if n(rx, x)) or '0')
        ia = m.createInitialAssignment()
        ia.setSymbol(x)
        ia.setMath(libsbml.parseL3Formula(conc_from(x, guess)))
    ids = {rx.getId() for rx in fast}
    named = [sr for rx in fast for l in (rx.getListOfReactants(), rx.getListOfProducts()) for sr in l if sr.isSetId()]
    lst = m.getListOfReactions()
    keep = [(sr.getId(), sr.getStoichiometry() if sr.isSetStoichiometry() else 1.0) for sr in named]
    for i in reversed([i for i in range(lst.size()) if lst.get(i).getId() in ids]):
        lst.remove(i)
    for sid, v in keep:  # a named stoichiometry other math may use: kept as a constant parameter
        p = m.createParameter()
        p.setId(sid)
        p.setValue(v)
        p.setConstant(sid not in moving)
    return set(fs)


def _algebraic_to_assignment(m, prefer=()):
    """Replace algebraic rules by assignment rules: match each rule to a symbol it determines, then
    eliminate in turn. A rule affine in its symbol (checked numerically at random values) is solved
    for it; any other stays a constraint: its symbol becomes an algebraic state (a rate rule holding
    the residual, flagged in the returned set) and the model a DAE (diffsol)."""
    alg = [r for r in m.getListOfRules() if r.isAlgebraic()]
    if not alg:
        return set()
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    free = _free_symbols(m)
    cand = [sorted(_names(r.getMath(), set()) & free, key=lambda x: (x not in prefer, x)) for r in alg]
    match = {}  # symbol -> rule index (maximum bipartite matching, augmenting paths)

    def assign(i, seen):
        for x in cand[i]:
            if x not in seen:
                seen.add(x)
                if x not in match or assign(match[x], seen):
                    match[x] = i
                    return True
        return False
    for i in range(len(alg)):
        if not assign(i, set()):
            raise NotImplementedError("algebraic rule %s determines no free symbol" % _l3(alg[i].getMath()))
    var = {i: x for x, i in match.items()}
    rng = np.random.default_rng(1)
    vals = {}

    def val(n):  # random positive stand-ins for every other symbol
        if n not in vals:
            vals[n] = rng.uniform(0.5, 2.0)
        return vals[n]
    g = [r.getMath().deepCopy() for r in alg]
    sol, dae = [], []
    for i in range(len(alg)):
        z = var[i]
        g0, g1 = _l3(_sub(g[i], z, '0')), _l3(_sub(g[i], z, '1'))
        f = '-(%s) / ((%s) - (%s))' % (g0, g1, g0)
        linear = True
        for zt in (2.0, -1.3, 3.7):  # affine in z?
            ex = evaluate(g[i], lambda n: zt if n == z else val(n), fdefs)
            lin = evaluate(libsbml.parseL3Formula('(%s) + (%s) * ((%s) - (%s))' % (g0, zt, g1, g0)), val, fdefs)
            linear = linear and bool(np.isclose(ex, lin, rtol=1e-9, atol=1e-12))
        if not linear:  # a true constraint: z stays a state of the DAE
            dae.append((z, _l3(g[i])))
            continue
        sol.append((z, f))
        for k in range(i + 1, len(alg)):  # eliminate z from the rules still to solve
            g[k] = _sub(g[k], z, f)
    lst = m.getListOfRules()
    for i in reversed([i for i in range(lst.size()) if lst.get(i).isAlgebraic()]):
        lst.remove(i)
    for z, f in sol:
        r = m.createAssignmentRule()
        r.setVariable(z)
        r.setMath(libsbml.parseL3Formula(f))
    sp = {x.getId(): x for x in m.getListOfSpecies()}
    for z, f in dae:
        r = m.createRateRule()  # its "rate" is the residual; the engine holds it at 0
        r.setVariable(z)
        r.setMath(libsbml.parseL3Formula(f))
        if z in sp:  # a rule-set species may not take part in reactions as a non-boundary one
            sp[z].setBoundaryCondition(True)
    return {z for z, _ in dae}


def _delayed_params_to_states(m):
    """delay() reads the run's state history: a parameter (or compartment, boundary species) that only
    events change, used inside delay(), becomes a state with rate 0 so its history is kept."""
    found = set()

    def walk(a, inside):
        if a is None:
            return
        if inside and a.isName():
            found.add(a.getName())
        d = _type(a) == libsbml.AST_FUNCTION_DELAY
        for i in range(a.getNumChildren()):
            walk(a.getChild(i), inside or (d and i == 0))
    for x in list(m.getListOfRules()) + list(m.getListOfInitialAssignments()) + list(m.getListOfConstraints()):
        walk(x.getMath(), False)
    for rx in m.getListOfReactions():
        if rx.getKineticLaw() is not None:
            walk(rx.getKineticLaw().getMath(), False)
    for ev in m.getListOfEvents():
        for x in [ev.getTrigger(), ev.getDelay(), ev.getPriority()] + list(ev.getListOfEventAssignments()):
            if x is not None:
                walk(x.getMath(), False)
    ruled = {r.getVariable() for r in m.getListOfRules()}
    targets = {a.getVariable() for ev in m.getListOfEvents() for a in ev.getListOfEventAssignments()}
    for x in found & targets - ruled:
        e = m.getElementBySId(x)
        if isinstance(e, (libsbml.Parameter, libsbml.Compartment)) or isinstance(e, libsbml.Species) and e.getBoundaryCondition():
            r = m.createRateRule()
            r.setVariable(x)
            r.setMath(libsbml.parseL3Formula('0'))


def from_sbml(sbml, info=None, base=None):
    """``from_roadrunner`` without roadrunner: compile SBML text into a ``roadrunner_rs.Model`` at the
    SBML initial values (initial assignments, assignment rules, amounts / concentrations as roadrunner).

    ``info`` also receives the model in roadrunner's terms: ids (floating, boundary, global parameter,
    compartment, reaction, rate-rule ids in roadrunner's order), rows/N (full stoichiometry), ia
    (initial assignments), fdefs, values (initial values, SBML units), doc (libsbml document).
    NotImplementedError / ValueError: not supported (load roadrunner instead). ``base``: the file's
    directory, for comp external model definitions."""
    doc = libsbml.readSBMLFromString(sbml)
    if doc.getModel() is not None and doc.getPlugin('comp') is not None and doc.isPackageEnabled('comp'):
        _flatten_comp(doc, base)
        sbml = libsbml.writeSBMLToString(doc)
    m = doc.getModel()
    if m is not None:
        _drop_empty_math(m)
    if m is None or doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) or doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL):
        raise ValueError("SBML errors")
    ns = doc.getNamespaces()
    if any(p in ns.getURI(i) for i in range(ns.getNumNamespaces()) for p in ('/distrib/', '/arrays/')):
        raise NotImplementedError("SBML package")
    pre = {}
    dae = _algebraic_to_assignment(m, prefer=_fast_to_algebraic(m, pre))
    _delayed_params_to_states(m)
    fdefs = {f.getId(): f.getMath() for f in m.getListOfFunctionDefinitions() if f.getMath() is not None}
    ia = {x.getSymbol(): x.getMath() for x in m.getListOfInitialAssignments()}
    ar = {x.getVariable(): x.getMath() for x in m.getListOfRules() if x.isAssignment()}
    rrules = [x.getVariable() for x in m.getListOfRules() if x.isRate()]
    rr_math = {x.getVariable(): x.getMath() for x in m.getListOfRules() if x.isRate()}
    rxns = {x.getId(): x for x in m.getListOfReactions()}
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
        elif x in refs and refs[x].isSetStoichiometryMath():  # SBML L2
            v = evaluate(refs[x].getStoichiometryMath().getMath(), get, fdefs)
        elif x in refs:
            v = refs[x].getStoichiometry() if refs[x].isSetStoichiometry() else 1.0  # unset: 1, as roadrunner
        elif x in rxns:  # a reaction id: its rate at t0
            kl = rxns[x].getKineticLaw()
            loc = {p.getId(): p.getValue() for p in (list(kl.getListOfParameters()) + list(kl.getListOfLocalParameters()))} \
                if kl is not None else {}
            v = evaluate(kl.getMath(), lambda y: loc[y] if y in loc else get(y), fdefs) \
                if kl is not None and kl.getMath() is not None else 0.0
        elif x.startswith('__rateof__'):
            v = rate0(x[10:])
        elif x in comps or x in sp:  # roadrunner's defaults
            v = 1.0 if x in comps else 0.0
        elif x in dae:  # an algebraic state without a value: a starting guess; the solver makes it consistent
            v = 0.0
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
        return sr.getStoichiometry() if sr.isSetStoichiometry() else 1.0  # unset: 1, as roadrunner

    def rate0(x):  # rateOf(x) at t0, in x's SBML units
        if x in rr_math:
            return evaluate(rr_math[x], get, fdefs)
        if x in sp and not sp[x].getBoundaryCondition() and x not in ar:
            flux = 0.0
            for rx in m.getListOfReactions():
                for sg, lst in ((-1, rx.getListOfReactants()), (1, rx.getListOfProducts())):
                    for sr in lst:
                        if sr.getSpecies() == x:
                            flux += sg * stoich(sr) * get(rx.getId())
            cf = sp[x].getConversionFactor() if sp[x].isSetConversionFactor() else m.getConversionFactor()
            flux *= get(cf) if cf else 1.0
            return flux if sp[x].getHasOnlySubstanceUnits() else flux / get(sp[x].getCompartment())
        if x in ar:
            raise NotImplementedError("initial rateOf(%s) of an assignment rule" % x)
        return 0.0

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
    out = _compile(m, m, rows, reactions, N, conc, get, lambda sr, changing: get(sr.getId()), info, dae, pre)
    if info is not None:
        ruled = set(ar) | set(rrules)
        info.update(doc=doc, ia=ia, fdefs=fdefs, values=vals, rows=rows, N=N, rate_rule_ids=rrules,
                    ids=(_rr_order(rows, ruled), _rr_order([s for s in sp if s not in ri], ruled),
                         _rr_order(list(pars), ruled), _rr_order(list(comps), ruled), reactions))
    return out

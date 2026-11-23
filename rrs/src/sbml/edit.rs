//! Model editing as SBML text -> SBML text: the "Model editing methods" of libroadrunner's C API
//! (rrc_api.h), with RoadRunner::addSpecies & co's semantics (rrRoadRunner.cpp). Formulas are
//! libsbml L3 infix. No regeneration here: the caller recompiles the returned SBML.
use crate::expr::{lex, Tok};
use sbml_io::*;

type R<T = String> = Result<T, String>;

/// Read, edit the model (given level, version), write.
fn edit(sbml: &str, f: impl FnOnce(&mut Model, u32, u32) -> R<()>) -> R {
    let mut doc = read_doc(sbml)?;
    let (l, v) = (doc.level, doc.version);
    f(doc.model.as_mut().ok_or("the document has no model")?, l, v)?;
    Ok(write(&doc))
}

fn read_doc(sbml: &str) -> R<Document> {
    read(sbml).map(|r| r.document).map_err(|d| d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n"))
}

// ---------------------------------------------------------------- ids and math helpers

/// libsbml Model::getElementBySId (what checkID rejects). ponytail: unit definitions and local
/// parameters left out, as libsbml's lookup does.
fn has_sid(m: &Model, id: &str) -> bool {
    m.function_definitions.iter().any(|x| x.id == id)
        || m.compartments.iter().any(|x| x.id == id)
        || m.species.iter().any(|x| x.id == id)
        || m.parameters.iter().any(|x| x.id == id)
        || m.events.iter().any(|x| x.id.as_deref() == Some(id))
        || m.reactions.iter().any(|r| r.id == id || sref(r, id))
}

fn sref(r: &Reaction, id: &str) -> bool {
    r.reactants.iter().chain(&r.products).any(|s| s.id.as_deref() == Some(id))
        || r.modifiers.iter().any(|s| s.id.as_deref() == Some(id))
}

fn check_id(m: &Model, f: &str, id: &str) -> R<()> {
    if id.is_empty() || !id.chars().next().unwrap().is_alphabetic() && !id.starts_with('_')
        || !id.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
        return Err(format!("Roadrunner::{f} failed: invalid id '{id}'."));
    }
    if has_sid(m, id) { return Err(format!("Roadrunner::{f} failed, identifier {id} already existed in the model")); }
    Ok(())
}

/// A compartment, species, parameter or species reference (addInitialAssignment, addEventAssignment).
fn is_symbol(m: &Model, id: &str) -> bool {
    m.compartments.iter().any(|x| x.id == id) || m.species.iter().any(|x| x.id == id)
        || m.parameters.iter().any(|x| x.id == id) || m.reactions.iter().any(|r| sref(r, id))
}

/// libroadrunner's hasVariable: any ci / function name equal to `id`.
fn uses(m: &Option<Math>, id: &str) -> bool {
    let mut hit = false;
    if let Some(m) = m { m.visit_idents(&mut |x| hit |= x == id); }
    hit
}

fn vars(m: &Option<Math>, out: &mut Vec<String>) {
    if let Some(m) = m { m.visit_idents(&mut |x| out.push(x.to_string())); }
}

fn species<'a>(m: &'a mut Model, f: &str, id: &str) -> R<&'a mut Species> {
    m.species.iter_mut().find(|s| s.id == id).ok_or(format!("Roadrunner::{f} failed, no species with ID {id} existed in the model"))
}
fn reaction<'a>(m: &'a mut Model, f: &str, id: &str) -> R<&'a mut Reaction> {
    m.reactions.iter_mut().find(|s| s.id == id).ok_or(format!("Roadrunner::{f} failed, no reaction with ID {id} existed in the model"))
}
fn event<'a>(m: &'a mut Model, f: &str, id: &str) -> R<&'a mut Event> {
    m.events.iter_mut().find(|s| s.id.as_deref() == Some(id)).ok_or(format!("Roadrunner::{f} failed, no event {id} existed in the model"))
}
/// Remove the first element matching `p`; whether there was one.
fn take<T>(v: &mut Vec<T>, p: impl Fn(&T) -> bool) -> bool {
    match v.iter().position(p) { Some(i) => { v.remove(i); true } None => false }
}

// ---------------------------------------------------------------- infix <-> Math

/// libsbml L3 infix (SBML_parseL3Formula's syntax) -> MathML tree. Lexer shared with expr.rs;
/// expr::parse itself folds into an evaluation DAG (no user function calls), so not reusable here.
pub fn parse_formula(src: &str) -> R<Math> {
    let mut p = P { t: lex(src)?, i: 0 };
    let r = p.or()?;
    if p.i != p.t.len() { return Err(format!("trailing input in '{src}'")); }
    Ok(r)
}

struct P { t: Vec<Tok>, i: usize }

/// a op b, flattened into an n-ary apply when a already is one (libsbml does the same).
fn nary(o: Op, a: Math, b: Math) -> Math {
    match a {
        Math::Apply(x, mut v) if x == o => { v.push(b); Math::Apply(o, v) }
        a => Math::Apply(o, vec![a, b]),
    }
}

impl P {
    fn eat(&mut self, s: &str) -> bool {
        let m = matches!(self.t.get(self.i), Some(Tok::Sym(x)) if *x == s);
        if m { self.i += 1; }
        m
    }
    fn or(&mut self) -> R<Math> {
        let mut a = self.and()?;
        while self.eat("||") { a = nary(Op::Or, a, self.and()?); }
        Ok(a)
    }
    fn and(&mut self) -> R<Math> {
        let mut a = self.cmp()?;
        while self.eat("&&") { a = nary(Op::And, a, self.cmp()?); }
        Ok(a)
    }
    fn cmp(&mut self) -> R<Math> {
        let a = self.sum()?;
        let ops = [("<=", Op::Leq), (">=", Op::Geq), ("==", Op::Eq), ("!=", Op::Neq), ("<", Op::Lt), (">", Op::Gt)];
        let Some(&(_, o)) = ops.iter().find(|(s, _)| self.eat(s)) else { return Ok(a) };
        let mut v = vec![a, self.sum()?];
        while self.eat(o_sym(o)) { v.push(self.sum()?); } // a < b < c: one n-ary relation
        Ok(Math::Apply(o, v))
    }
    fn sum(&mut self) -> R<Math> {
        let mut a = self.prod()?;
        loop {
            if self.eat("+") { a = nary(Op::Plus, a, self.prod()?); }
            else if self.eat("-") { a = Math::Apply(Op::Minus, vec![a, self.prod()?]); }
            else { return Ok(a); }
        }
    }
    fn prod(&mut self) -> R<Math> {
        let mut a = self.unary()?;
        loop {
            if self.eat("*") { a = nary(Op::Times, a, self.unary()?); }
            else if self.eat("/") { a = Math::Apply(Op::Divide, vec![a, self.unary()?]); }
            else { return Ok(a); }
        }
    }
    fn unary(&mut self) -> R<Math> {
        if self.eat("-") { return Ok(Math::Apply(Op::Minus, vec![self.unary()?])); }
        if self.eat("+") { return self.unary(); }
        if self.eat("!") { return Ok(Math::Apply(Op::Not, vec![self.unary()?])); }
        let a = self.atom()?;
        if self.eat("^") { return Ok(Math::Apply(Op::Power, vec![a, self.unary()?])); }
        Ok(a)
    }
    fn atom(&mut self) -> R<Math> {
        let tok = self.t.get(self.i).cloned().ok_or("unexpected end of formula")?;
        self.i += 1;
        match tok {
            Tok::Num(x) => Ok(if x.fract() == 0.0 && x.abs() < 9e15 { Math::int(x as i64) } else { Math::real(x) }),
            Tok::Sym("(") => { let a = self.or()?; if self.eat(")") { Ok(a) } else { Err("expected ')'".into()) } }
            Tok::Id(f) if self.eat("(") => {
                let mut a = vec![];
                if !self.eat(")") {
                    loop { a.push(self.or()?); if self.eat(")") { break; } if !self.eat(",") { return Err("expected ','".into()); } }
                }
                call(&f, a)
            }
            Tok::Id(n) => Ok(match n.as_str() {
                "time" => Math::Symbol(Symbol::new(Symbol::TIME, "time")),
                "avogadro" => Math::Symbol(Symbol::new(Symbol::AVOGADRO, "avogadro")),
                "true" | "false" | "pi" | "exponentiale" => Math::Constant(Constant::from_tag(&n).unwrap()),
                "infinity" | "INF" | "inf" | "Infinity" => Math::Constant(Constant::Infinity),
                "notanumber" | "NaN" | "nan" => Math::Constant(Constant::NotANumber),
                _ => Math::Ident(n),
            }),
            Tok::Sym(s) => Err(format!("unexpected '{s}'")),
        }
    }
}

fn o_sym(o: Op) -> &'static str {
    match o { Op::Leq => "<=", Op::Geq => ">=", Op::Eq => "==", Op::Neq => "!=", Op::Lt => "<", _ => ">" }
}

fn call(f: &str, mut a: Vec<Math>) -> R<Math> {
    let n = a.len();
    let bx = |m: Math| Box::new(m);
    Ok(match (f, n) {
        ("log10", 1) | ("log", 1) => Math::Log { base: Some(bx(Math::int(10))), arg: bx(a.remove(0)) },
        ("log", 2) => { let b = a.remove(0); Math::Log { base: Some(bx(b)), arg: bx(a.remove(0)) } }
        ("sqrt", 1) | ("root", 1) => Math::Root { degree: None, arg: bx(a.remove(0)) },
        ("root", 2) => { let d = a.remove(0); Math::Root { degree: Some(bx(d)), arg: bx(a.remove(0)) } }
        ("piecewise", _) => {
            let otherwise = (n % 2 == 1).then(|| bx(a.pop().unwrap()));
            let mut it = a.into_iter();
            let mut pieces = vec![];
            while let (Some(v), Some(c)) = (it.next(), it.next()) { pieces.push((v, c)); }
            Math::Piecewise { pieces, otherwise }
        }
        ("delay", _) => Math::CallSymbol(Symbol::new(Symbol::DELAY, "delay"), a),
        ("rateOf", _) => Math::CallSymbol(Symbol::new(Symbol::RATE_OF, "rateOf"), a),
        _ => {
            let tag = match f {
                "pow" => "power", "ceil" => "ceiling", "asin" => "arcsin", "acos" => "arccos", "atan" => "arctan",
                "asinh" => "arcsinh", "acosh" => "arccosh", "atanh" => "arctanh", f => f,
            };
            match Op::from_tag(tag) { Some(o) => Math::Apply(o, a), None => Math::Call(f.to_string(), a) }
        }
    })
}

/// MathML tree -> L3 infix with minimal parentheses (SBML_formulaToL3String's syntax; libsbml's
/// exact spacing and grouping are not copied, parse_formula and libsbml read both back the same).
pub fn formula(m: &Math) -> String { fmt(m, 0) }

fn fmt(m: &Math, min: u8) -> String {
    let all = |a: &[Math], p: u8, s: &str| a.iter().map(|x| fmt(x, p)).collect::<Vec<_>>().join(s);
    let call = |f: &str, a: &[Math]| format!("{f}({})", all(a, 0, ", "));
    let (s, p) = match m {
        Math::Number { value, .. } => (match value {
            Number::Integer(i) => i.to_string(),
            Number::ENotation(x, e) => format!("{x}e{e}"),
            v => format!("{:?}", v.value()),
        }, 8),
        Math::Ident(x) => (x.clone(), 8),
        Math::Symbol(s) => (if s.url == Symbol::TIME { "time" } else { "avogadro" }.into(), 8),
        Math::Constant(c) => (match c { Constant::Infinity => "INF", Constant::NotANumber => "NaN", c => c.tag() }.into(), 8),
        Math::Apply(o, a) => match (o, a.len()) {
            // later operands of the same precedence in parentheses, so the tree reads back unchanged
            (Op::Plus, 2..) => (format!("{} + {}", fmt(&a[0], 4), all(&a[1..], 5, " + ")), 4),
            (Op::Times, 2..) => (format!("{} * {}", fmt(&a[0], 5), all(&a[1..], 6, " * ")), 5),
            (Op::Plus | Op::Times, 1) => return fmt(&a[0], min),
            (Op::Minus, 2) => (format!("{} - {}", fmt(&a[0], 4), fmt(&a[1], 5)), 4),
            (Op::Divide, 2) => (format!("{} / {}", fmt(&a[0], 5), fmt(&a[1], 6)), 5),
            (Op::Minus, 1) => (format!("-{}", fmt(&a[0], 6)), 6),
            (Op::Not, 1) => (format!("!{}", fmt(&a[0], 6)), 6),
            (Op::Power, 2) => (format!("{}^{}", fmt(&a[0], 8), fmt(&a[1], 6)), 7),
            (Op::Lt | Op::Gt | Op::Leq | Op::Geq | Op::Eq | Op::Neq, 2..) => (all(a, 4, &format!(" {} ", o_sym(*o))), 3),
            (Op::And, 2..) => (all(a, 3, " && "), 2),
            (Op::Or, 2..) => (all(a, 2, " || "), 1),
            (Op::Ceiling, _) => (call("ceil", a), 8),
            _ => (call(o.tag(), a), 8),
        },
        Math::Call(f, a) => (call(f, a), 8),
        Math::CallSymbol(s, a) => (call(if s.url == Symbol::DELAY { "delay" } else if s.url == Symbol::RATE_OF { "rateOf" } else { &s.name }, a), 8),
        Math::Root { degree: None, arg } => (format!("sqrt({})", fmt(arg, 0)), 8),
        Math::Root { degree: Some(d), arg } => (format!("root({}, {})", fmt(d, 0), fmt(arg, 0)), 8),
        Math::Log { base, arg } => match base.as_deref() {
            None | Some(Math::Number { value: Number::Integer(10), .. }) => (format!("log10({})", fmt(arg, 0)), 8),
            Some(b) => (format!("log({}, {})", fmt(b, 0), fmt(arg, 0)), 8),
        },
        Math::Piecewise { pieces, otherwise } => {
            let mut v: Vec<String> = pieces.iter().flat_map(|(x, c)| [fmt(x, 0), fmt(c, 0)]).collect();
            v.extend(otherwise.iter().map(|o| fmt(o, 0)));
            (format!("piecewise({})", v.join(", ")), 8)
        }
        Math::Lambda { args, body } => (format!("lambda({}, {})", args.join(", "), fmt(body, 0)), 8),
        Math::Semantics { math, .. } => return fmt(math, min),
        Math::Raw(x) => (x.clone(), 8),
    };
    if p < min { format!("({s})") } else { s }
}

fn rename(m: &mut Math, from: &str, to: &str) {
    match m {
        Math::Ident(x) if x == from => *x = to.to_string(),
        Math::Apply(_, a) | Math::Call(_, a) | Math::CallSymbol(_, a) => a.iter_mut().for_each(|x| rename(x, from, to)),
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            if let Some(q) = q { rename(q, from, to); }
            rename(arg, from, to);
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, c) in pieces { rename(v, from, to); rename(c, from, to); }
            if let Some(o) = otherwise { rename(o, from, to); }
        }
        Math::Semantics { math, .. } => rename(math, from, to),
        _ => {}
    }
}

// ---------------------------------------------------------------- removal cascade

/// RoadRunner::removeVariable: drop or trim everything whose math refers to `sid`.
fn remove_variable(m: &mut Model, level: u32, sid: &str) {
    let smath = |s: &SpeciesReference| s.stoichiometry_math.as_ref().is_some_and(|x| uses(&Some(x.clone()), sid));
    m.reactions.retain(|r| level != 2 || !r.reactants.iter().chain(&r.products).any(smath));
    for r in &mut m.reactions {
        if r.kinetic_law.as_ref().is_some_and(|k| uses(&k.math, sid)) { r.kinetic_law = None; }
    }
    m.rules.retain(|r| r.variable() != Some(sid) && !uses(&r.math().cloned(), sid));
    m.function_definitions.retain(|f| !uses(&f.math, sid));
    m.constraints.retain(|c| !uses(&c.math, sid));
    m.initial_assignments.retain(|a| a.symbol != sid && !uses(&a.math, sid));
    m.events.retain(|e| !e.trigger.as_ref().is_some_and(|t| uses(&t.math, sid)));
    for e in &mut m.events {
        if e.priority.as_ref().is_some_and(|p| uses(&p.math, sid)) { e.priority = None; }
        if e.delay.as_ref().is_some_and(|p| uses(&p.math, sid)) { e.delay = None; }
        e.assignments.retain(|a| a.variable != sid && !uses(&a.math, sid));
    }
}

fn remove_species_in(m: &mut Model, level: u32, sid: &str) -> R<()> {
    if !take(&mut m.species, |s| s.id == sid) {
        return Err(format!("Roadrunner::removeSpecies failed, no species with ID {sid} existed in the model"));
    }
    m.reactions.retain(|r| !r.reactants.iter().chain(&r.products).any(|s| s.species == sid) && !r.modifiers.iter().any(|s| s.species == sid));
    remove_variable(m, level, sid);
    Ok(())
}

fn remove_parameter_in(m: &mut Model, level: u32, pid: &str) -> R<()> {
    if !take(&mut m.parameters, |p| p.id == pid) {
        return Err(format!("Roadrunner::removeParameter failed, no parameter with ID {pid} existed in the model"));
    }
    remove_variable(m, level, pid);
    Ok(())
}

/// RoadRunner::isParameterUsed (reactions, rules, initial assignments; not events, as upstream).
fn parameter_used(m: &Model, level: u32, sid: &str) -> bool {
    m.reactions.iter().any(|r| {
        (level == 2 && r.reactants.iter().chain(&r.products).any(|s| uses(&s.stoichiometry_math, sid)))
            || r.kinetic_law.as_ref().is_some_and(|k| uses(&k.math, sid))
    }) || m.rules.iter().any(|r| r.variable() == Some(sid) || uses(&r.math().cloned(), sid))
        || m.initial_assignments.iter().any(|a| a.symbol == sid || uses(&a.math, sid))
}

// ---------------------------------------------------------------- species

const UNIT_KINDS: &[&str] = &["ampere", "avogadro", "becquerel", "candela", "celsius", "coulomb", "dimensionless",
    "farad", "gram", "gray", "henry", "hertz", "item", "joule", "katal", "kelvin", "kilogram", "liter", "litre", "lumen",
    "lux", "meter", "metre", "mole", "newton", "ohm", "pascal", "radian", "second", "siemens", "sievert", "steradian",
    "tesla", "volt", "watt", "weber"];

fn add_species(sbml: &str, f: &str, id: &str, comp: &str, amount: Option<f64>, conc: Option<f64>, hosu: bool, boundary: bool, units: &str) -> R {
    edit(sbml, |m, _, _| {
        check_id(m, f, id)?;
        // ponytail: upstream skips this check in the NoRegen variants (compartment added later); we always check
        if !m.compartments.iter().any(|c| c.id == comp) {
            return Err(format!("Roadrunner::{f} failed, no compartment {comp} existed in the model"));
        }
        let ok_unit = !units.is_empty() && (m.unit_definitions.iter().any(|u| u.id == units) || UNIT_KINDS.contains(&units));
        m.species.push(Species {
            id: id.into(), compartment: comp.into(), initial_amount: amount, initial_concentration: conc,
            has_only_substance_units: Some(hosu), boundary_condition: Some(boundary), constant: Some(false),
            substance_units: ok_unit.then(|| units.into()), ..Default::default()
        });
        Ok(())
    })
}

/// addSpeciesAmount
pub fn add_species_amount(sbml: &str, id: &str, comp: &str, amount: f64, hosu: bool, boundary: bool, units: &str) -> R {
    add_species(sbml, "addSpeciesAmount", id, comp, Some(amount), None, hosu, boundary, units)
}

/// addSpeciesConcentration
pub fn add_species_concentration(sbml: &str, id: &str, comp: &str, conc: f64, hosu: bool, boundary: bool, units: &str) -> R {
    add_species(sbml, "addSpeciesConcentration", id, comp, None, Some(conc), hosu, boundary, units)
}

/// removeSpecies: also every reaction it takes part in, and all math using it.
pub fn remove_species(sbml: &str, id: &str) -> R { edit(sbml, |m, l, _| remove_species_in(m, l, id)) }

/// setBoundary
pub fn set_boundary(sbml: &str, id: &str, b: bool) -> R {
    edit(sbml, |m, _, _| { species(m, "setBoundary", id)?.boundary_condition = Some(b); Ok(()) })
}

/// setHasOnlySubstanceUnits
pub fn set_has_only_substance_units(sbml: &str, id: &str, b: bool) -> R {
    edit(sbml, |m, _, _| { species(m, "setHasOnlySubstanceUnits", id)?.has_only_substance_units = Some(b); Ok(()) })
}

/// setInitAmount / setInitConcentration: the SBML value, minus any initial assignment.
fn set_init(sbml: &str, f: &str, id: &str, v: f64, conc: bool) -> R {
    edit(sbml, |m, _, _| {
        // upstream sets the SBML value, then fails setting init() on the regenerated model
        if m.rules.iter().any(|r| matches!(r, Rule::Assignment { variable, .. } if variable == id)) {
            return Err(format!("Could not set value for {id}, as it is defined by an assignment rule, and can not be set independently."));
        }
        let s = species(m, f, id)?;
        (s.initial_amount, s.initial_concentration) = if conc { (None, Some(v)) } else { (Some(v), None) };
        take(&mut m.initial_assignments, |a| a.symbol == id);
        Ok(())
    })
}
/// setInitAmount
pub fn set_init_amount(sbml: &str, id: &str, v: f64) -> R { set_init(sbml, "setInitAmount", id, v, false) }
/// setInitConcentration
pub fn set_init_concentration(sbml: &str, id: &str, v: f64) -> R { set_init(sbml, "setInitConcentration", id, v, true) }

/// setConstant: species, else parameter, else compartment.
pub fn set_constant(sbml: &str, id: &str, b: bool) -> R {
    edit(sbml, |m, _, _| {
        if let Some(s) = m.species.iter_mut().find(|x| x.id == id) { s.constant = Some(b) }
        else if let Some(p) = m.parameters.iter_mut().find(|x| x.id == id) { p.constant = Some(b) }
        else if let Some(c) = m.compartments.iter_mut().find(|x| x.id == id) { c.constant = Some(b) }
        else { return Err(format!("Roadrunner::setConstant failed, no species/ parameter/ compartment with ID {id} existed in the model")) }
        Ok(())
    })
}

// ---------------------------------------------------------------- reactions

/// addReactionFromSBML: a `<reaction>` element, read at the document's level and version.
pub fn add_reaction_from_sbml(sbml: &str, xml: &str) -> R {
    edit(sbml, |m, l, v| {
        let x = xml.trim_start();
        let x = if x.starts_with("<?xml") { &x[x.find("?>").map_or(0, |i| i + 2)..] } else { x };
        let wrap = format!(r#"<sbml xmlns="{}" level="{l}" version="{v}"><model><listOfReactions>{x}</listOfReactions></model></sbml>"#,
            ns::core(l, v).unwrap_or(""));
        let r = read_doc(&wrap)?.model.and_then(|mut m| m.reactions.pop()).ok_or("no <reaction> in the given SBML")?;
        check_id(m, "addReaction", &r.id)?; // ponytail: upstream only fails later, at regeneration
        m.reactions.push(r);
        Ok(())
    })
}

/// strtod-style "2.5S1" -> (2.5, "S1"); stoichiometry 1 without a number.
fn parse_species(s: &str) -> (f64, &str) {
    let s = s.trim();
    let b = s.as_bytes();
    let mut i = 0;
    let digits = |i: &mut usize| { let st = *i; while *i < b.len() && b[*i].is_ascii_digit() { *i += 1; } *i > st };
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') { i += 1; }
    let mut num = digits(&mut i);
    if i < b.len() && b[i] == b'.' { i += 1; num |= digits(&mut i); }
    if !num { return (1.0, s); }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') { j += 1; }
        if digits(&mut j) { i = j; }
    }
    (s[..i].parse().unwrap_or(1.0), s[i..].trim_start())
}

/// addReaction: reactants/products as "S1" or "2 S1" / "2S1" (stoichiometry first). Irreversible,
/// as upstream. Species in the rate law that are not in the reactant/product strings become modifiers
/// (string comparison and repeats as upstream).
pub fn add_reaction(sbml: &str, id: &str, reactants: &[&str], products: &[&str], law: &str) -> R {
    edit(sbml, |m, l, v| {
        check_id(m, "addReaction", id)?;
        // fast is gone in L3V2
        let mut r = Reaction { id: id.into(), reversible: Some(false), fast: ((l, v) != (3, 2)).then_some(false), ..Default::default() };
        for (list, out) in [(reactants, &mut r.reactants), (products, &mut r.products)] {
            for s in list {
                let (st, sid) = parse_species(s);
                if sid.is_empty() { return Err(format!("Roadrunner::addReaction failed, invalid stoichiomety value and species ID {s}")); }
                if !m.species.iter().any(|x| x.id == sid) {
                    return Err(format!("Roadrunner::addReaction failed, no species with ID {s} existed in the model"));
                }
                out.push(SpeciesReference { species: sid.into(), stoichiometry: Some(st), constant: (l == 3).then_some(true), ..Default::default() });
            }
        }
        let math = parse_formula(law)?;
        let mut ids = vec![];
        math.visit_idents(&mut |x| if m.species.iter().any(|s| s.id == x) { ids.push(x.to_string()) });
        for s in ids {
            if !products.contains(&s.as_str()) && !reactants.contains(&s.as_str()) {
                r.modifiers.push(ModifierSpeciesReference { species: s, ..Default::default() });
            }
        }
        r.kinetic_law = Some(KineticLaw { math: Some(math), ..Default::default() });
        m.reactions.push(r);
        Ok(())
    })
}

/// removeReaction; with `delete_unused`, also the global parameters its rate law / stoichiometries
/// used that nothing else (reactions, rules, initial assignments) uses.
pub fn remove_reaction(sbml: &str, id: &str, delete_unused: bool) -> R {
    edit(sbml, |m, l, _| {
        let i = m.reactions.iter().position(|r| r.id == id)
            .ok_or(format!("Roadrunner::removeReaction failed, no reaction with ID {id} existed in the model"))?;
        let r = m.reactions.remove(i);
        if delete_unused {
            let mut check = vec![];
            vars(&r.kinetic_law.and_then(|k| k.math), &mut check);
            for s in r.reactants.iter().chain(&r.products) {
                check.extend(s.id.clone());
                vars(&s.stoichiometry_math, &mut check);
            }
            check.sort();
            check.dedup();
            for p in check {
                if m.parameters.iter().any(|x| x.id == p) && !parameter_used(m, l, &p) { remove_parameter_in(m, l, &p)?; }
            }
        }
        Ok(())
    })
}

/// setReversible
pub fn set_reversible(sbml: &str, id: &str, b: bool) -> R {
    edit(sbml, |m, _, _| { reaction(m, "setReversible", id)?.reversible = Some(b); Ok(()) })
}

/// setKineticLaw: new math (local parameters kept); species it uses that are not reactants or
/// products are added as modifiers (repeats as upstream).
pub fn set_kinetic_law(sbml: &str, id: &str, law: &str) -> R {
    edit(sbml, |m, _, _| {
        let math = parse_formula(law)?;
        let mut ids = vec![];
        math.visit_idents(&mut |x| if m.species.iter().any(|s| s.id == x) { ids.push(x.to_string()) });
        let r = reaction(m, "setKineticLaw", id)?;
        for s in ids {
            if !r.reactants.iter().chain(&r.products).any(|x| x.species == s) {
                r.modifiers.push(ModifierSpeciesReference { species: s, ..Default::default() });
            }
        }
        r.kinetic_law.get_or_insert_with(Default::default).math = Some(math);
        Ok(())
    })
}

/// getKineticLaw: the rate law as infix ("" if none).
pub fn get_kinetic_law(sbml: &str, id: &str) -> R {
    let doc = read_doc(sbml)?;
    let r = doc.model.as_ref().and_then(|m| m.reactions.iter().find(|r| r.id == id))
        .ok_or(format!("Roadrunner::getKineticLaw failed, no reaction with ID {id} exists in the model"))?;
    Ok(r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()).map(formula).unwrap_or_default())
}

// ---------------------------------------------------------------- parameters, compartments

/// addParameter (constant="false", as upstream)
pub fn add_parameter(sbml: &str, id: &str, value: f64) -> R {
    edit(sbml, |m, _, _| {
        check_id(m, "addParameter", id)?;
        m.parameters.push(Parameter { id: id.into(), value: Some(value), constant: Some(false), ..Default::default() });
        Ok(())
    })
}

/// removeParameter: also all math using it.
pub fn remove_parameter(sbml: &str, id: &str) -> R { edit(sbml, |m, l, _| remove_parameter_in(m, l, id)) }

/// addCompartment (constant="false", as upstream)
pub fn add_compartment(sbml: &str, id: &str, size: f64) -> R {
    edit(sbml, |m, _, _| {
        check_id(m, "addCompartment", id)?;
        m.compartments.push(Compartment { id: id.into(), size: Some(size), constant: Some(false), ..Default::default() });
        Ok(())
    })
}

/// removeCompartment: also its species (with their cascade) and all math using it.
pub fn remove_compartment(sbml: &str, id: &str) -> R {
    edit(sbml, |m, l, _| {
        if !take(&mut m.compartments, |c| c.id == id) {
            return Err(format!("Roadrunner::removeCompartment failed, no compartment with ID {id} existed in the model"));
        }
        let gone: Vec<String> = m.species.iter().filter(|s| s.compartment == id).map(|s| s.id.clone()).collect();
        for s in gone { remove_species_in(m, l, &s)?; }
        remove_variable(m, l, id);
        Ok(())
    })
}

// ---------------------------------------------------------------- rules, initial assignments

/// checkAddRule: the target must be a species (boundary, or boundary unset), parameter,
/// compartment or species reference; it is made non-constant; no other rule may set it.
fn check_add_rule(m: &mut Model, l: u32, vid: &str) -> R<()> {
    if let Some(s) = m.species.iter_mut().find(|x| x.id == vid) {
        s.constant = Some(false);
        if s.boundary_condition == Some(false) || s.boundary_condition.is_none() && l == 2 { // L2: unset is false
            return Err(format!("Unable to add rule because the species with ID {vid} is not a boundary species, and therefore may only change due to reactions.  To change this, the boundary species flag can be set with the function 'setBoundary(id, value)' or when added with 'addSpecies'."));
        }
        s.boundary_condition = Some(true);
    } else if let Some(p) = m.parameters.iter_mut().find(|x| x.id == vid) {
        p.constant = Some(false);
    } else if let Some(c) = m.compartments.iter_mut().find(|x| x.id == vid) {
        c.constant = Some(false);
    } else if let Some(s) = m.reactions.iter_mut().flat_map(|r| r.reactants.iter_mut().chain(r.products.iter_mut())).find(|s| s.id.as_deref() == Some(vid)) {
        if l == 3 { s.constant = Some(false); }
    } else if m.reactions.iter().any(|r| r.id == vid) {
        return Err(format!("Unable to add rule because the variable with ID {vid} is a Reaction, which does not have mathematical meaning."));
    } else {
        return Err(format!("Unable to add rule because no variable with ID {vid} was found in the model."));
    }
    if m.rules.iter().any(|r| r.variable() == Some(vid)) {
        return Err(format!("Unable to add rule because the variable {vid} already has a rule existing in the model."));
    }
    Ok(())
}

/// addAssignmentRule
pub fn add_assignment_rule(sbml: &str, vid: &str, f: &str) -> R {
    edit(sbml, |m, l, _| {
        check_add_rule(m, l, vid)?;
        m.rules.push(Rule::Assignment { variable: vid.into(), math: Some(parse_formula(f)?), meta: Default::default() });
        Ok(())
    })
}

/// addRateRule
pub fn add_rate_rule(sbml: &str, vid: &str, f: &str) -> R {
    edit(sbml, |m, l, _| {
        check_add_rule(m, l, vid)?;
        m.rules.push(Rule::Rate { variable: vid.into(), math: Some(parse_formula(f)?), meta: Default::default() });
        Ok(())
    })
}

/// removeRules: the rule setting `vid`. ponytail: upstream's useInitialValueAsCurrent only touches
/// the live model's state, nothing in the SBML, so it is the core's business.
pub fn remove_rules(sbml: &str, vid: &str) -> R {
    edit(sbml, |m, _, _| {
        if take(&mut m.rules, |r| r.variable() == Some(vid)) { Ok(()) }
        else { Err(format!("Roadrunner::removeRules failed, no rules for variable {vid} existed in the model")) }
    })
}

/// addInitialAssignment
pub fn add_initial_assignment(sbml: &str, vid: &str, f: &str) -> R {
    edit(sbml, |m, l, v| {
        if (l, v) == (2, 1) { // upstream crashes here (libsbml has no initial assignments in L2V1)
            return Err("Roadrunner::addInitialAssignment failed, SBML Level 2 Version 1 has no initial assignments".into());
        }
        if !is_symbol(m, vid) {
            return Err(format!("Roadrunner::addInitialAssignment failed, no symbol with ID {vid} existed in the model"));
        }
        if m.rules.iter().any(|r| matches!(r, Rule::Assignment { variable, .. } if variable == vid)) {
            return Err(format!("Roadrunner::addInitialAssignment failed, symbol {vid} already has an assignment rule existing in the model"));
        }
        if m.initial_assignments.iter().any(|a| a.symbol == vid) {
            return Err(format!("Roadrunner::addInitialAssignment failed, symbol {vid} already has an initial assignment existing in the model"));
        }
        m.initial_assignments.push(InitialAssignment { symbol: vid.into(), math: Some(parse_formula(f)?), meta: Default::default() });
        Ok(())
    })
}

// ---------------------------------------------------------------- events

/// addEvent: trigger persistent="false" initialValue="false" in Level 3, as upstream.
pub fn add_event(sbml: &str, id: &str, use_values_from_trigger_time: bool, trigger: &str) -> R {
    edit(sbml, |m, l, v| {
        check_id(m, "addEvent", id)?;
        let math = parse_formula(trigger).map_err(|e| format!("Roadrunner::addEvent failed, an error occurred in parsing the trigger formula: {e}"))?;
        m.events.push(Event {
            id: Some(id.into()),
            use_values_from_trigger_time: (l == 3 || v >= 4).then_some(use_values_from_trigger_time),
            trigger: Some(Trigger { math: Some(math), persistent: (l == 3).then_some(false), initial_value: (l == 3).then_some(false), ..Default::default() }),
            ..Default::default()
        });
        Ok(())
    })
}

/// addTrigger: replaces the event's trigger. Upstream sets persistent / initialValue only for
/// Level > 3, i.e. never; kept as is.
pub fn add_trigger(sbml: &str, eid: &str, trigger: &str) -> R {
    edit(sbml, |m, _, _| {
        let math = parse_formula(trigger).map_err(|e| format!("Roadrunner::addTrigger failed, an error occurred in parsing the trigger formula: {e}"))?;
        event(m, "addTrigger", eid)?.trigger = Some(Trigger { math: Some(math), ..Default::default() });
        Ok(())
    })
}

fn trigger_flag(sbml: &str, f: &str, attr: &str, eid: &str, set: impl FnOnce(&mut Trigger)) -> R {
    edit(sbml, |m, l, _| {
        if l < 3 {
            return Err(format!("Roadrunner::{f} failed, current SBML level and version does not have {attr} attribute in the trigger"));
        }
        let t = event(m, f, eid)?.trigger.as_mut().ok_or(format!("Roadrunner::{f} failed, given event {eid} does not have a trigger"))?;
        set(t);
        Ok(())
    })
}

/// setPersistent
pub fn set_persistent(sbml: &str, eid: &str, b: bool) -> R {
    trigger_flag(sbml, "setPersistent", "peresistent", eid, |t| t.persistent = Some(b))
}

/// setTriggerInitialValue
pub fn set_trigger_initial_value(sbml: &str, eid: &str, b: bool) -> R {
    trigger_flag(sbml, "setTriggerInitialValue", "initialValue", eid, |t| t.initial_value = Some(b))
}

/// addPriority (Level 3 only)
pub fn add_priority(sbml: &str, eid: &str, f: &str) -> R {
    edit(sbml, |m, l, _| {
        let e = event(m, "addPriority", eid)?;
        if l < 3 { return Err("Roadrunner::addPriority failed, current SBML level and version does not support Priority in event".into()); }
        e.priority = Some(Priority { math: Some(parse_formula(f)?), ..Default::default() });
        Ok(())
    })
}

/// addDelay
pub fn add_delay(sbml: &str, eid: &str, f: &str) -> R {
    edit(sbml, |m, _, _| {
        let math = parse_formula(f)?;
        event(m, "addDelay", eid)?.delay = Some(Delay { math: Some(math), ..Default::default() });
        Ok(())
    })
}

/// addEventAssignment
pub fn add_event_assignment(sbml: &str, eid: &str, vid: &str, f: &str) -> R {
    edit(sbml, |m, _, _| {
        event(m, "addEventAssignment", eid)?;
        if !is_symbol(m, vid) {
            return Err(format!("Roadrunner::addEventAssignment failed, no variable with ID {vid} existed in the model"));
        }
        if m.rules.iter().any(|r| matches!(r, Rule::Assignment { variable, .. } if variable == vid)) {
            return Err(format!("Roadrunner::addEventAssignment failed, variable {vid} already has a rule existing in the model"));
        }
        let math = parse_formula(f)?;
        let e = event(m, "addEventAssignment", eid)?;
        if e.assignments.iter().any(|a| a.variable == vid) {
            return Err(format!("Roadrunner::addEventAssignment failed, variable {vid} already has an assignment in the event {eid}"));
        }
        e.assignments.push(EventAssignment { variable: vid.into(), math: Some(math), meta: Default::default() });
        Ok(())
    })
}

/// removeEventAssignments: the (first) assignment to `vid` in event `eid`.
pub fn remove_event_assignments(sbml: &str, eid: &str, vid: &str) -> R {
    edit(sbml, |m, _, _| {
        let e = m.events.iter_mut().find(|e| e.id.as_deref() == Some(eid))
            .ok_or(format!("Roadrunner::removeEventAssignment failed, no event with ID {eid} existed in the model"))?;
        if take(&mut e.assignments, |a| a.variable == vid) { Ok(()) }
        else { Err(format!("Roadrunner::removeEventAssignment failed, no event assignment for variable {vid} existed in the event {eid}")) }
    })
}

/// removeEvent
pub fn remove_event(sbml: &str, eid: &str) -> R {
    edit(sbml, |m, _, _| {
        if take(&mut m.events, |e| e.id.as_deref() == Some(eid)) { Ok(()) }
        else { Err(format!("Roadrunner::removeEvent failed, no event with ID {eid} existed in the model")) }
    })
}

// ---------------------------------------------------------------- whole document

/// validateCurrentSBML: Err with the problems found. ponytail: sbml-io's structural checks (ids,
/// references), not libsbml's full general / identifier / MathML / overdetermined validation.
pub fn validate_current_sbml(sbml: &str) -> R<()> {
    let issues: Vec<String> = validate(&read_doc(sbml)?).iter().map(|i| i.to_string()).collect();
    if issues.is_empty() { Ok(()) } else { Err(issues.join("\n")) }
}

/// getParamPromotedSBML (libsbml's SBMLLocalParameterConverter): every local parameter becomes a
/// constant global `<reaction id>_<parameter id>`. Takes SBML text or a file path, as upstream.
/// ponytail: no renaming if `<rid>_<pid>` is already taken.
pub fn get_param_promoted_sbml(sbml: &str) -> R {
    let text = if sbml.trim_start().starts_with('<') { sbml.to_string() } else { std::fs::read_to_string(sbml).map_err(|e| e.to_string())? };
    edit(&text, |m, _, _| {
        for r in &mut m.reactions {
            let Some(k) = &mut r.kinetic_law else { continue };
            for p in std::mem::take(&mut k.parameters).into_iter().rev() { // libsbml promotes the last first
                let id = format!("{}_{}", r.id, p.id);
                if let Some(x) = &mut k.math { rename(x, &p.id, &id); }
                m.parameters.push(Parameter { id, name: p.name, value: p.value, units: p.units, constant: Some(true), meta: p.meta });
            }
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const M: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  <model id="m">
    <listOfCompartments><compartment id="C" size="1" constant="true"/></listOfCompartments>
    <listOfSpecies>
      <species id="S1" compartment="C" initialConcentration="10" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="S2" compartment="C" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="0.1" constant="true"/>
      <parameter id="k2" value="0.2" constant="true"/>
    </listOfParameters>
    <listOfReactions>
      <reaction id="J0" reversible="false">
        <listOfReactants><speciesReference species="S1" stoichiometry="1" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="S2" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k</ci><ci>S1</ci><ci>kl</ci></apply></math>
          <listOfLocalParameters><localParameter id="kl" value="3"/></listOfLocalParameters>
        </kineticLaw>
      </reaction>
      <reaction id="J1" reversible="false">
        <listOfReactants><speciesReference species="S2" stoichiometry="1" constant="true"/></listOfReactants>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k2</ci><ci>S2</ci></apply></math></kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"#;

    fn model(s: &str) -> Model { read(s).unwrap().document.model.unwrap() }

    #[test]
    fn formulas_round_trip() {
        for f in ["k * S1 / (1 + S1^2) - 1e-05 * S2", "a - (b - c)", "a / (b * c)", "-a^2 + (-a)^2", "(a^b)^c", "a^b^c",
                  "f(a, b) * exp(-t * time)", "piecewise(1, a < b, 0)", "a < b && (c || !d)", "log(2, x) + log10(x) + sqrt(x) + root(3, x)",
                  "delay(x, 1) + ceil(x) + pi + INF"] {
            let m = parse_formula(f).unwrap();
            assert_eq!(parse_formula(&formula(&m)).unwrap(), m, "{f} -> {}", formula(&m));
        }
        assert_eq!(formula(&parse_formula("a+b+c").unwrap()), "a + b + c");
        assert_eq!(formula(&parse_formula("(a+b)*c").unwrap()), "(a + b) * c");
        assert!(parse_formula("a +").is_err());
    }

    #[test]
    fn species_edits() {
        let s = add_species_amount(M, "S3", "C", 5.0, false, true, "mole").unwrap();
        let sp = model(&s).species.pop().unwrap();
        assert_eq!((sp.initial_amount, sp.boundary_condition, sp.substance_units.as_deref()), (Some(5.0), Some(true), Some("mole")));
        assert!(add_species_amount(M, "S1", "C", 1.0, false, false, "").unwrap_err().contains("already existed"));
        assert!(add_species_concentration(M, "S3", "X", 1.0, false, false, "").unwrap_err().contains("no compartment"));
        let m = model(&remove_species(M, "S1").unwrap());
        assert_eq!(m.reactions.iter().map(|r| r.id.as_str()).collect::<Vec<_>>(), ["J1"]);
        let m = model(&set_init_amount(M, "S2", 4.0).unwrap());
        assert_eq!((m.species[1].initial_amount, m.species[1].initial_concentration), (Some(4.0), None));
        assert!(set_constant(M, "nope", true).is_err());
        assert!(add_assignment_rule(M, "S1", "k").unwrap_err().contains("not a boundary species"));
    }

    #[test]
    fn reaction_edits() {
        let s = add_reaction(M, "J2", &["2S1", "S2"], &["3 S2"], "k*S1*S2").unwrap();
        let r = model(&s).reactions.pop().unwrap();
        assert_eq!((r.reactants[0].species.as_str(), r.reactants[0].stoichiometry, r.products[0].stoichiometry), ("S1", Some(2.0), Some(3.0)));
        assert_eq!(r.modifiers.len(), 1); // "S1" is not literally in the reactant strings: upstream adds it
        assert_eq!(get_kinetic_law(&s, "J2").unwrap(), "k * S1 * S2");
        assert!(add_reaction(M, "J3", &["X"], &[], "1").is_err());
        let m = model(&remove_reaction(M, "J1", true).unwrap());
        assert!(!m.parameters.iter().any(|p| p.id == "k2") && m.parameters.iter().any(|p| p.id == "k"));
        let m = model(&remove_parameter(M, "k2").unwrap());
        assert!(m.reactions[1].kinetic_law.is_none());
        let m = model(&get_param_promoted_sbml(M).unwrap());
        assert_eq!(m.parameters.last().unwrap().id, "J0_kl");
        assert_eq!(get_kinetic_law(&write(&Document { model: Some(m), ..Document::new(3, 1) }), "J0").unwrap(), "k * S1 * J0_kl");
        let s = add_reaction_from_sbml(M, r#"<reaction id="J9" reversible="false"><listOfReactants><speciesReference species="S2" stoichiometry="1" constant="true"/></listOfReactants></reaction>"#).unwrap();
        assert_eq!(model(&s).reactions.len(), 3);
    }

    #[test]
    fn event_edits() {
        let s = add_event(M, "E1", true, "time > 1").unwrap();
        let s = add_event_assignment(&s, "E1", "k", "2*k").unwrap();
        let s = add_delay(&s, "E1", "0.5").unwrap();
        let s = set_persistent(&s, "E1", true).unwrap();
        assert!(add_event_assignment(&s, "E1", "k", "1").unwrap_err().contains("already has an assignment"));
        let e = &model(&s).events[0];
        assert_eq!((e.trigger.as_ref().unwrap().persistent, e.assignments.len(), e.delay.is_some()), (Some(true), 1, true));
        // removing k drops the assignment that sets it
        assert!(model(&remove_parameter(&s, "k").unwrap()).events[0].assignments.is_empty());
        assert!(model(&remove_event(&s, "E1").unwrap()).events.is_empty());
        assert!(remove_event(M, "E1").is_err());
        assert!(validate_current_sbml(M).is_ok());
    }
}

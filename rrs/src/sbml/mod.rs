//! SBML text -> engine `Model`, without roadrunner or libsbml (port of rrs_bridge.from_sbml).
//!
//! Values are the SBML initial values as roadrunner computes them (initial assignments, assignment
//! rules, amounts / concentrations). Fast reactions and algebraic rules are rewritten first
//! (`transform`); the compiler then writes every rate, rule and event as an infix formula over
//! species concentrations and parameters, as `Model::new` takes them.
mod math;
mod transform;
pub mod edit;

use crate::model::{EventDef, Model, Rules};
use math::*;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Valid SBML the engine does not cover (Python: NotImplementedError).
    NotImplemented(String),
    /// Invalid or inconsistent SBML (Python: ValueError).
    Invalid(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { Error::NotImplemented(s) | Error::Invalid(s) => f.write_str(s) }
    }
}

pub(crate) fn nimp(s: String) -> Error { Error::NotImplemented(s) }
pub(crate) fn inval(s: String) -> Error { Error::Invalid(s) }

/// What the Python side needs besides the model (rrs_bridge's `info`).
#[derive(Debug, Default)]
pub struct Info {
    /// (floating, boundary, global parameter, compartment, reaction) ids in roadrunner's order.
    pub ids: [Vec<String>; 5],
    pub comp_of: HashMap<String, String>,
    /// Species whose d[S]/dt is a rate rule although they have reactions.
    pub general: HashSet<String>,
    /// Ids set by rules or events.
    pub changing: HashSet<String>,
    /// Non-floating species whose parameter value is an amount.
    pub amount_params: HashSet<String>,
    pub rate_rule_ids: Vec<String>,
    /// Full stoichiometry matrix rows (non-boundary species) x reactions, row-major.
    pub rows: Vec<String>,
    pub n: Vec<f64>,
    /// Initial assignments (symbol, infix) and function definitions (id, "lambda(args, body)").
    pub ia: Vec<(String, String)>,
    pub fdefs: Vec<(String, String)>,
    /// Initial values in SBML units.
    pub values: HashMap<String, f64>,
    // the model after the rewrites, as a libsbml document would show it
    pub n_events: usize,
    pub rule_vars: Vec<String>,
    pub hosu: HashSet<String>,
    /// Conversion factor ids: the model's, then the species' (those set).
    pub conversion_factors: Vec<String>,
    /// Species references: (species, reaction, -1 reactant / 1 product, id, has stoichiometryMath).
    pub srefs: Vec<(String, String, i32, Option<String>, bool)>,
}

// ---------------------------------------------------------------- the model, after reading

pub(crate) struct Comp { id: String, size: Option<f64>, constant: bool }
pub(crate) struct Sp {
    id: String, comp: String, amt: Option<f64>, conc: Option<f64>, hosu: bool, boundary: bool,
    constant: bool, cf: Option<String>,
}
pub(crate) struct Par { id: String, value: Option<f64>, constant: bool }
pub(crate) struct SRef { id: Option<String>, species: String, stoich: Option<f64>, smath: Option<E>, constant: bool }
pub(crate) struct Rx {
    id: String, fast: bool, reactants: Vec<SRef>, products: Vec<SRef>,
    /// (math, local parameters)
    kl: Option<(Option<E>, Vec<(String, f64)>)>,
}
#[derive(PartialEq, Clone, Copy)]
pub(crate) enum RK { Assign, Rate, Alg }
pub(crate) struct Rule { kind: RK, var: String, math: E }
pub(crate) struct Ev {
    id: Option<String>, uvftt: bool, trigger: E, init: bool, persistent: bool,
    delay: Option<E>, priority: Option<E>, asg: Vec<(String, E)>,
}
#[derive(Default)]
pub(crate) struct M {
    cf: Option<String>,
    fdefs: Vec<(String, Vec<String>, E)>,
    comps: Vec<Comp>,
    species: Vec<Sp>,
    params: Vec<Par>,
    ias: Vec<(String, E)>,
    rules: Vec<Rule>,
    constraints: Vec<E>,
    reactions: Vec<Rx>,
    events: Vec<Ev>,
}

impl Rx {
    fn refs(&self) -> impl Iterator<Item = (i32, &SRef)> {
        self.reactants.iter().map(|r| (-1, r)).chain(self.products.iter().map(|r| (1, r)))
    }
}

impl M {
    fn fdefs(&self) -> Fdefs { self.fdefs.iter().map(|(f, a, b)| (f.clone(), (a.clone(), b.clone()))).collect() }
    fn sids(&self) -> HashSet<String> {
        let mut s: HashSet<String> = self.comps.iter().map(|c| c.id.clone()).chain(self.species.iter().map(|x| x.id.clone()))
            .chain(self.params.iter().map(|p| p.id.clone())).chain(self.fdefs.iter().map(|f| f.0.clone()))
            .chain(self.events.iter().filter_map(|e| e.id.clone())).collect();
        for r in &self.reactions {
            s.insert(r.id.clone());
            s.extend(r.refs().filter_map(|(_, sr)| sr.id.clone()));
            if let Some(k) = &r.kl { s.extend(k.1.iter().map(|p| p.0.clone())); }
        }
        s
    }

    /// From sbml-io's model, with libsbml's defaults; elements whose math is missing are dropped
    /// (SBML L3V2 lets math be absent; such an element has no effect).
    fn read(m: &sbml_io::Model) -> Result<M, Error> {
        let c = |x: &Option<sbml_io::Math>| x.as_ref().map(conv).transpose();
        let mut fdefs = vec![];
        for f in &m.function_definitions {
            let mut x = f.math.as_ref();
            while let Some(sbml_io::Math::Semantics { math, .. }) = x { x = Some(math); }
            match x {
                Some(sbml_io::Math::Lambda { args, body }) => fdefs.push((f.id.clone(), args.clone(), conv(body)?)),
                Some(b) => fdefs.push((f.id.clone(), vec![], conv(b)?)),
                None => {}
            }
        }
        let sref = |r: &sbml_io::SpeciesReference| -> Result<SRef, Error> {
            Ok(SRef { id: r.id.clone(), species: r.species.clone(), stoich: r.stoichiometry, smath: c(&r.stoichiometry_math)?,
                      constant: r.constant.unwrap_or(false) })
        };
        let mut reactions = vec![];
        for r in &m.reactions {
            reactions.push(Rx {
                id: r.id.clone(), fast: r.fast.unwrap_or(false),
                reactants: r.reactants.iter().map(sref).collect::<Result<_, _>>()?,
                products: r.products.iter().map(sref).collect::<Result<_, _>>()?,
                kl: match &r.kinetic_law {
                    Some(k) => Some((c(&k.math)?, k.parameters.iter().map(|p| (p.id.clone(), p.value.unwrap_or(f64::NAN))).collect())),
                    None => None,
                },
            });
        }
        let mut rules = vec![];
        for r in &m.rules {
            let Some(x) = c(&r.math().cloned())? else { continue };
            let (kind, var) = match r {
                sbml_io::Rule::Assignment { variable, .. } => (RK::Assign, variable.clone()),
                sbml_io::Rule::Rate { variable, .. } => (RK::Rate, variable.clone()),
                sbml_io::Rule::Algebraic { .. } => (RK::Alg, String::new()),
            };
            rules.push(Rule { kind, var, math: x });
        }
        let mut events = vec![];
        for e in &m.events {
            let Some(t) = e.trigger.as_ref() else { continue };
            let Some(trigger) = c(&t.math)? else { continue };
            let mut asg = vec![];
            for a in &e.assignments {
                if let Some(x) = c(&a.math)? { asg.push((a.variable.clone(), x)); }
            }
            events.push(Ev {
                id: e.id.clone(), uvftt: e.use_values_from_trigger_time.unwrap_or(true), trigger,
                init: t.initial_value.unwrap_or(true), persistent: t.persistent.unwrap_or(true),
                delay: match &e.delay { Some(d) => c(&d.math)?, None => None },
                priority: match &e.priority { Some(p) => c(&p.math)?, None => None }, asg,
            });
        }
        let mut ias = vec![];
        for a in &m.initial_assignments {
            if let Some(x) = c(&a.math)? { ias.push((a.symbol.clone(), x)); }
        }
        let mut constraints = vec![];
        for k in &m.constraints {
            if let Some(x) = c(&k.math)? { constraints.push(x); }
        }
        Ok(M {
            cf: m.conversion_factor.clone().filter(|s| !s.is_empty()),
            fdefs,
            comps: m.compartments.iter().map(|x| Comp { id: x.id.clone(), size: x.size, constant: x.constant.unwrap_or(true) }).collect(),
            species: m.species.iter().map(|x| Sp {
                id: x.id.clone(), comp: x.compartment.clone(), amt: x.initial_amount, conc: x.initial_concentration,
                hosu: x.has_only_substance_units.unwrap_or(false), // a constant species never changes: roadrunner lists it with the boundary species
                boundary: x.boundary_condition.unwrap_or(false) || x.constant.unwrap_or(false),
                constant: x.constant.unwrap_or(false), cf: x.conversion_factor.clone().filter(|s| !s.is_empty()),
            }).collect(),
            params: m.parameters.iter().map(|x| Par { id: x.id.clone(), value: x.value, constant: x.constant.unwrap_or(true) }).collect(),
            ias, rules, constraints, reactions, events,
        })
    }
}

/// roadrunner lists rule-determined symbols after the others.
fn rr_order(xs: &[String], ruled: &HashSet<String>) -> Vec<String> {
    xs.iter().filter(|x| !ruled.contains(*x)).chain(xs.iter().filter(|x| ruled.contains(*x))).cloned().collect()
}

// ---------------------------------------------------------------- initial values

/// Initial values in SBML units, as roadrunner computes them (memoized).
struct Init<'a> {
    m: &'a M,
    fd: &'a Fdefs,
    dae: &'a HashSet<String>,
    vals: HashMap<String, f64>,
    busy: HashSet<String>,
}

impl<'a> Init<'a> {
    fn get(&mut self, x: &str) -> Result<f64, Error> {
        if let Some(&v) = self.vals.get(x) { return Ok(v); }
        if !self.busy.insert(x.to_string()) { return Err(inval(format!("cyclic initial values through {x}"))); }
        let (m, fd) = (self.m, self.fd);
        let ar = m.rules.iter().find(|r| r.kind == RK::Assign && r.var == x).map(|r| &r.math);
        let ia = m.ias.iter().find(|a| a.0 == x).map(|a| &a.1);
        let sp = m.species.iter().find(|s| s.id == x);
        let comp = m.comps.iter().find(|c| c.id == x);
        let par = m.params.iter().find(|p| p.id == x);
        let sr = m.reactions.iter().flat_map(|r| r.refs()).map(|p| p.1).find(|r| r.id.as_deref() == Some(x));
        let rx = m.reactions.iter().find(|r| r.id == x);
        let v = if let Some(e) = ar.or(ia) {
            eval(e, &mut |y| self.get(y), fd)?
        } else if let Some(s) = comp.and_then(|c| c.size) {
            s
        } else if let Some(v) = par.and_then(|p| p.value) {
            v
        } else if let Some(s) = sp.filter(|s| s.conc.is_some()) {
            s.conc.unwrap() * if s.hosu { self.get(&s.comp)? } else { 1.0 }
        } else if let Some(s) = sp.filter(|s| s.amt.is_some()) {
            s.amt.unwrap() / if s.hosu { 1.0 } else { self.get(&s.comp)? }
        } else if let Some(e) = sr.and_then(|r| r.smath.as_ref()) { // SBML L2
            eval(e, &mut |y| self.get(y), fd)?
        } else if let Some(r) = sr {
            r.stoich.unwrap_or(1.0) // unset: 1, as roadrunner
        } else if let Some(r) = rx { // a reaction id: its rate at t0
            match &r.kl {
                Some((Some(k), loc)) => eval(k, &mut |y| match loc.iter().find(|p| p.0 == y) {
                    Some(p) => Ok(p.1), None => self.get(y) }, fd)?,
                _ => 0.0,
            }
        } else if let Some(y) = x.strip_prefix("__rateof__") {
            self.rate0(y)?
        } else if comp.is_some() || sp.is_some() { // roadrunner's defaults
            if comp.is_some() { 1.0 } else { 0.0 }
        } else if self.dae.contains(x) { // an algebraic state without a value: a starting guess
            0.0
        } else { // unset parameter (roadrunner refuses it), reaction ids in initial assignments
            return Err(nimp(format!("initial value of {x}")));
        };
        self.busy.remove(x);
        self.vals.insert(x.to_string(), v);
        Ok(v)
    }

    fn stoich(&mut self, sr: &SRef) -> Result<f64, Error> {
        if let Some(i) = &sr.id { return self.get(i); }
        if let Some(e) = &sr.smath { let fd = self.fd; return eval(e, &mut |y| self.get(y), fd); }
        Ok(sr.stoich.unwrap_or(1.0)) // unset: 1, as roadrunner
    }

    /// rateOf(x) at t0, in x's SBML units.
    fn rate0(&mut self, x: &str) -> Result<f64, Error> {
        let m = self.m;
        if let Some(r) = m.rules.iter().find(|r| r.kind == RK::Rate && r.var == x) { let fd = self.fd; return eval(&r.math, &mut |y| self.get(y), fd); }
        let ar = m.rules.iter().any(|r| r.kind == RK::Assign && r.var == x);
        if let Some(s) = m.species.iter().find(|s| s.id == x && !s.boundary && !ar) {
            let mut flux = 0.0;
            for rx in &m.reactions {
                for (sg, sr) in rx.refs() {
                    if sr.species == x { flux += sg as f64 * self.stoich(sr)? * self.get(&rx.id)?; }
                }
            }
            if let Some(cf) = s.cf.as_ref().or(m.cf.as_ref()) { flux *= self.get(cf)?; }
            return Ok(if s.hosu { flux } else { flux / self.get(&s.comp)? });
        }
        if ar { return Err(nimp(format!("initial rateOf({x}) of an assignment rule"))); }
        Ok(0.0)
    }
}

// ---------------------------------------------------------------- entry point

/// Compile SBML text into a `Model` at the SBML initial values; `Info` describes it in roadrunner's terms.
pub fn from_sbml(sbml: &str) -> Result<(Model, Info), Error> {
    let doc = sbml_io::read(sbml).map_err(|e| inval(format!("SBML errors: {}", e.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("; "))))?.document;
    if doc.namespaces.iter().any(|(_, u)| u.contains("/comp/")) { return Err(nimp("comp: flatten first".into())); }
    if doc.namespaces.iter().any(|(_, u)| u.contains("/arrays/")) { return Err(nimp("SBML package".into())); }
    let Some(sm) = doc.model.as_ref() else { return Err(inval("SBML errors: no model".into())) };
    let mut m = M::read(sm)?;
    let mut pre = vec![];
    let prefer = transform::fast_to_algebraic(&mut m, &mut pre)?;
    let dae = transform::algebraic_to_assignment(&mut m, &prefer)?;
    transform::delayed_params_to_states(&mut m);

    let fd = m.fdefs();
    let mut init = Init { m: &m, fd: &fd, dae: &dae, vals: HashMap::new(), busy: HashSet::new() };
    let rows: Vec<String> = m.species.iter().filter(|s| !s.boundary).map(|s| s.id.clone()).collect();
    let reactions: Vec<String> = m.reactions.iter().map(|r| r.id.clone()).collect();
    let ri: HashMap<&str, usize> = rows.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let nr = reactions.len();
    let mut nm = vec![0.0; rows.len() * nr];
    for (j, rx) in m.reactions.iter().enumerate() {
        for (sg, sr) in rx.refs() {
            if let Some(&i) = ri.get(sr.species.as_str()) { nm[i * nr + j] += sg as f64 * init.stoich(sr)?; }
        }
    }
    for x in m.comps.iter().map(|c| &c.id).chain(m.params.iter().map(|p| &p.id)).chain(m.species.iter().map(|s| &s.id)) {
        init.get(x)?;
    }
    let mut conc = HashMap::new();
    for s in m.species.iter().filter(|s| !s.boundary) {
        let v = init.vals[&s.id];
        conc.insert(s.id.clone(), if s.hosu { v / init.vals[&s.comp] } else { v });
    }
    let mut info = Info::default();
    let model = compile(&m, &rows, &reactions, &nm, &conc, &mut init, &mut info, &dae, &pre)?;

    let ruled: HashSet<String> = m.rules.iter().filter(|r| r.kind != RK::Alg).map(|r| r.var.clone()).collect();
    let ids = |v: Vec<String>| rr_order(&v, &ruled);
    info.ids = [ids(rows.clone()), ids(m.species.iter().filter(|s| !ri.contains_key(s.id.as_str())).map(|s| s.id.clone()).collect()),
                ids(m.params.iter().map(|p| p.id.clone()).collect()), ids(m.comps.iter().map(|c| c.id.clone()).collect()), reactions];
    info.rate_rule_ids = m.rules.iter().filter(|r| r.kind == RK::Rate).map(|r| r.var.clone()).collect();
    info.ia = m.ias.iter().map(|(x, e)| (x.clone(), s(e))).collect();
    info.fdefs = m.fdefs.iter().map(|(f, a, b)| {
        let mut parts = a.clone();
        parts.push(s(b));
        (f.clone(), format!("lambda({})", parts.join(", ")))
    }).collect();
    info.values = init.vals;
    info.rows = rows;
    info.n = nm;
    info.n_events = m.events.len();
    info.rule_vars = m.rules.iter().filter(|r| r.kind != RK::Alg).map(|r| r.var.clone()).collect();
    info.hosu = m.species.iter().filter(|s| s.hosu).map(|s| s.id.clone()).collect();
    info.conversion_factors = m.cf.iter().chain(m.species.iter().filter_map(|s| s.cf.as_ref())).cloned().collect();
    info.srefs = m.reactions.iter().flat_map(|rx| rx.refs().map(move |(sg, sr)|
        (sr.species.clone(), rx.id.clone(), sg, sr.id.clone(), sr.smath.is_some()))).collect();
    Ok((model, info))
}

// ---------------------------------------------------------------- compiler

/// Insertion-ordered parameter values.
#[derive(Default)]
struct Params { names: Vec<String>, vals: HashMap<String, f64> }
impl Params {
    fn set(&mut self, k: &str, v: f64) {
        if self.vals.insert(k.to_string(), v).is_none() { self.names.push(k.to_string()); }
    }
}

struct C<'a> {
    fd: Fdefs,
    rules: HashMap<String, E>,
    rule_order: Vec<String>,
    rate_rules: HashMap<String, E>,
    rate_order: Vec<String>,
    comp_of: HashMap<String, String>,
    hosu: HashSet<String>,
    floating: HashSet<String>,
    bamt: HashSet<String>,
    kl_math: HashMap<String, Option<E>>,
    local: HashMap<String, Vec<(String, String)>>,
    ias: Vec<(String, E)>,
    pre: &'a [(String, E)],
    params: Params,
    known: HashSet<String>,
    delays: Vec<((String, String, String), String)>,
    memo: HashMap<(char, String), String>,
    busy: HashSet<(char, String)>,
    rows: &'a [String],
    reactions: &'a [String],
    nm: &'a [f64],
    vstoich: HashMap<(String, String), Vec<(i32, E)>>,
    cf: HashMap<String, String>,
}

impl<'a> C<'a> {
    fn kinetic(&mut self, rid: &str) -> Result<E, Error> {
        let Some(k) = self.kl_math[rid].clone() else { return Ok(num(0.0)) }; // no kinetic law: rate 0, as roadrunner
        let ren: HashMap<String, E> = self.local[rid].iter().map(|(o, n)| (o.clone(), name(n))).collect();
        self.expand(&k.subst(&ren))
    }

    /// Copy with function calls inlined, csymbols resolved and reaction ids -> kinetic laws.
    fn expand(&mut self, e: &E) -> Result<E, Error> {
        match e {
            E::Avogadro => return Ok(num(AVOGADRO)),
            E::Time => return Ok(name("__time__")), // the clock, even when a model id is named 'time'
            E::Op(f, a) if f == "rateOf" => return match a.first() {
                Some(E::Name(x)) => Ok(name(&format!("__rateof__{x}"))),
                _ => Err(nimp("rateOf of an expression".into())),
            },
            E::Op(f, a) if f == "delay" => { // a parameter the engine refills with the delayed value
                if a.len() != 2 { return Err(inval("delay() takes 2 arguments".into())); }
                let bs = self.before_start(&a[0]);
                let key = (self.formula(&a[0], "delay()")?, self.formula(&a[1], "delay()")?, self.formula(&bs, "delay()")?);
                if let Some((_, n)) = self.delays.iter().find(|d| d.0 == key) { return Ok(name(n)); }
                let n = format!("__delay_{}__", self.delays.len());
                self.params.set(&n, 0.0);
                self.known.insert(n.clone());
                self.delays.push((key, n.clone()));
                return Ok(name(&n));
            }
            E::Name(n) if self.kl_math.contains_key(n) => return self.kinetic(&n.clone()),
            _ => {}
        }
        let kids = e.kids().iter().map(|k| self.expand(k)).collect::<Result<Vec<_>, _>>()?;
        Ok(match e {
            E::Call(f, _) => {
                let Some((args, body)) = self.fd.get(f).cloned() else { return Err(inval(format!("undefined function {f}"))) };
                if args.len() != kids.len() { return Err(inval(format!("{f}() takes {} arguments", args.len()))); }
                let b = body.subst(&args.into_iter().zip(kids).collect());
                self.expand(&b)?
            }
            E::Op(f, _) => E::Op(f.clone(), kids),
            _ => e.clone(),
        })
    }

    /// delay() before the run: everything holds its initial value (rateOf 0), and a symbol with an
    /// initial assignment is that formula at the earlier time.
    fn before_start(&self, e: &E) -> E {
        fn zero_rates(e: &E) -> E {
            match e {
                E::Op(f, _) if f == "rateOf" => num(0.0),
                E::Op(f, a) => E::Op(f.clone(), a.iter().map(zero_rates).collect()),
                E::Call(f, a) => E::Call(f.clone(), a.iter().map(zero_rates).collect()),
                _ => e.clone(),
            }
        }
        let mut out = zero_rates(e);
        for (x, f) in &self.ias {
            if out.names().contains(x) && !self.pre.iter().any(|p| &p.0 == x) { out = out.sub1(x, f); }
        }
        for (x, f) in self.pre {
            if out.names().contains(x) { out = out.sub1(x, f); }
        }
        out
    }

    /// Substitute assignment rules and amounts ([S] * V) of hasOnlySubstanceUnits species.
    fn inline(&mut self, e: E) -> Result<E, Error> {
        let mut e = e;
        for _ in 0..1000 {
            let todo: Vec<String> = e.names().into_iter().filter(|v| self.rules.contains_key(v)
                || (self.hosu.contains(v) && self.floating.contains(v)) || self.bamt.contains(v)).collect();
            if todo.is_empty() {
                let back: HashMap<String, E> = self.hosu.iter().filter(|v| self.floating.contains(*v)).chain(&self.bamt)
                    .map(|v| (format!("__conc__{v}"), name(v))).collect();
                return Ok(e.subst(&back));
            }
            let mut rep = HashMap::new();
            for v in todo {
                let new = if let Some(r) = self.rules.get(&v).cloned() { self.expand(&r)? } else {
                    let c = name(&format!("__conc__{v}"));
                    let comp = name(&self.comp_of[&v]);
                    if self.bamt.contains(&v) { div(c, comp) } else { mul(c, comp) }
                };
                rep.insert(v, new);
            }
            e = e.subst(&rep);
        }
        Err(inval("cyclic assignment rules".into()))
    }

    fn formula(&mut self, e: &E, what: &str) -> Result<String, Error> {
        let x = self.expand(e)?;
        let ast = self.inline(x)?;
        let names = ast.names();
        let mut unknown: Vec<&String> = names.iter().filter(|n| !n.starts_with("__rateof__") && !self.known.contains(*n)).collect();
        if !unknown.is_empty() {
            unknown.sort();
            return Err(inval(format!("{what} references {unknown:?}")));
        }
        let mut rep = HashMap::new();
        for n in &names {
            if let Some(y) = n.strip_prefix("__rateof__") { rep.insert(n.clone(), self.rate_of(y)?); }
        }
        Ok(l3(&ast, &rep))
    }

    fn lazy(&mut self, key: (char, String), f: impl FnOnce(&mut Self) -> Result<String, Error>) -> Result<String, Error> {
        if let Some(v) = self.memo.get(&key) { return Ok(v.clone()); }
        if !self.busy.insert(key.clone()) {
            return Err(nimp(format!("implicit rate of change through {key:?}"))); // e.g. V := f([S in V])
        }
        let v = f(self)?;
        self.busy.remove(&key);
        self.memo.insert(key, v.clone());
        Ok(v)
    }

    fn rate(&mut self, j: usize) -> Result<String, Error> {
        let rid = self.reactions[j].clone();
        self.lazy(('v', j.to_string()), |c| { let k = c.kinetic(&rid)?; c.formula(&k, &format!("kinetic law of {rid}")) })
    }

    fn comp_formula(&mut self, s: &str, what: &str) -> Result<String, Error> {
        let c = name(&self.comp_of[s]);
        self.formula(&c, what)
    }

    /// d x/dt in x's SBML units (rateOf csymbol).
    fn rate_of(&mut self, x: &str) -> Result<String, Error> {
        let x = x.to_string();
        self.lazy(('d', x.clone()), move |c| {
            if let Some(r) = c.rate_rules.get(&x).cloned() { return c.formula(&r, &format!("rate rule for {x}")); }
            if c.floating.contains(&x) {
                if !c.hosu.contains(&x) { return c.dconc(&x); }
                let v = c.comp_formula(&x, &x)?;
                let comp = c.comp_of[&x].clone();
                return Ok(format!("({}) * {} + {} * ({})", c.dconc(&x)?, v, x, c.rate_of(&comp)?));
            }
            if let Some(r) = c.rules.get(&x).cloned() { // total time derivative of the rule expression
                let ex = c.expand(&r)?;
                let e = c.inline(ex)?;
                let terms: Vec<E> = derivatives(&e)?.into_iter().filter(|(_, d)| *d != num(0.0))
                    .map(|(v, d)| if v == "__time__" { d } else { mul(d, name(&format!("__rateof__{v}"))) }).collect();
                let t = if terms.is_empty() { num(0.0) } else { sum(terms) };
                return c.formula(&t, &format!("rateOf({x})"));
            }
            Ok("0".into())
        })
    }

    /// d[S]/dt (concentration / time) of floating species s.
    fn dconc(&mut self, s: &str) -> Result<String, Error> {
        if let Some(r) = self.rate_rules.get(s).cloned() {
            let f = self.formula(&r, &format!("rate rule for {s}"))?;
            return self.conc_of(s, f);
        }
        let i = self.rows.iter().position(|r| r == s).unwrap();
        let nr = self.reactions.len();
        let mut terms = vec![];
        for j in 0..nr {
            let rid = self.reactions[j].clone();
            let n = if let Some(v) = self.vstoich.get(&(s.to_string(), rid.clone())).cloned() {
                let mut p = vec![];
                for (sg, a) in v { p.push(format!("{sg} * ({})", self.formula(&a, &format!("stoichiometry of {s} in {rid}"))?)); }
                p.join(" + ")
            } else if self.nm[i * nr + j] != 0.0 {
                fnum(self.nm[i * nr + j])
            } else {
                continue;
            };
            terms.push(format!("({n}) * ({})", self.rate(j)?));
        }
        let mut flux = if terms.is_empty() { "0".to_string() } else { terms.join(" + ") };
        if let Some(cf) = self.cf.get(s).cloned() {
            flux = format!("({}) * ({flux})", self.formula(&name(&cf), &format!("conversion factor of {s}"))?);
        }
        let v = self.comp_formula(s, &format!("compartment of {s}"))?;
        let comp = self.comp_of[s].clone();
        let dv = self.rate_of(&comp)?;
        Ok(format!("({flux}) / ({v})") + &if dv == "0" { String::new() } else { format!(" - {s} * ({dv}) / ({v})") })
    }

    /// value for species var in its SBML units -> state units (concentration)
    fn conc_of(&mut self, var: &str, f: String) -> Result<String, Error> {
        if self.hosu.contains(var) && (self.floating.contains(var) || self.rules.contains_key(var)) {
            return Ok(format!("({f}) / ({})", self.comp_formula(var, &format!("compartment of {var}"))?));
        }
        Ok(f)
    }
}

/// Words (ids) of a formula.
fn words(s: &str) -> impl Iterator<Item = &str> {
    s.split(|c: char| !(c.is_alphanumeric() || c == '_')).filter(|w| !w.is_empty())
}

#[allow(clippy::too_many_arguments)]
fn compile(m: &M, rows: &[String], reactions: &[String], nm: &[f64], conc: &HashMap<String, f64>, init: &mut Init,
           info: &mut Info, dae: &HashSet<String>, pre: &[(String, E)]) -> Result<Model, Error> {
    let mut rules = HashMap::new();
    let mut rule_order = vec![];
    let mut rate_rules = HashMap::new();
    let mut rate_order = vec![];
    for r in &m.rules {
        match r.kind {
            RK::Assign => { if rules.insert(r.var.clone(), r.math.clone()).is_none() { rule_order.push(r.var.clone()); } }
            RK::Rate => { if rate_rules.insert(r.var.clone(), r.math.clone()).is_none() { rate_order.push(r.var.clone()); } }
            RK::Alg => {}
        }
    }
    let targets: HashSet<String> = m.events.iter().flat_map(|e| e.asg.iter().map(|a| a.0.clone())).collect();
    let changing: HashSet<String> = rules.keys().chain(rate_rules.keys()).cloned().chain(targets).collect();
    let sp: HashMap<&str, &Sp> = m.species.iter().map(|s| (s.id.as_str(), s)).collect();
    let species: Vec<String> = rows.iter().filter(|s| !rules.contains_key(*s)).cloned().collect(); // rule species: outputs only
    let floating: HashSet<String> = species.iter().cloned().collect();
    let comp_of: HashMap<String, String> = m.species.iter().map(|s| (s.id.clone(), s.comp.clone())).collect();
    let hosu: HashSet<String> = m.species.iter().filter(|s| s.hosu).map(|s| s.id.clone()).collect();
    let mut vol = HashMap::new();
    for s in &m.species { vol.insert(s.id.clone(), init.get(&s.comp)?); }

    let mut params = Params::default();
    for x in m.params.iter().map(|p| &p.id).chain(m.comps.iter().map(|c| &c.id)) {
        if !rules.contains_key(x) { params.set(x, init.get(x)?); }
    }
    for s in &m.species { // boundary / non-floating species: parameters in their SBML units
        if !rules.contains_key(&s.id) && !floating.contains(&s.id) { params.set(&s.id, init.get(&s.id)?); }
    }
    // variable stoichiometry: (species, reaction) -> [(sign, stoichiometry)] over all references
    let refs: Vec<(i32, &str, &SRef)> = m.reactions.iter().flat_map(|r| r.refs().map(move |(sg, sr)| (sg, r.id.as_str(), sr))).collect();
    for (_, _, sr) in &refs {
        if let Some(i) = sr.id.as_ref().filter(|i| !rules.contains_key(*i)) { params.set(i, init.get(i)?); }
    }
    let var_ref: HashSet<(String, String)> = refs.iter().filter(|(_, _, sr)| sr.id.as_ref().is_some_and(|i| changing.contains(i)) || sr.smath.is_some())
        .map(|(_, j, sr)| (sr.species.clone(), j.to_string())).collect();
    let mut vstoich: HashMap<(String, String), Vec<(i32, E)>> = HashMap::new();
    for (sg, j, sr) in &refs {
        let k = (sr.species.clone(), j.to_string());
        if var_ref.contains(&k) {
            let a = if let Some(e) = &sr.smath { e.clone() } else if let Some(i) = &sr.id { name(i) } else { num(sr.stoich.unwrap_or(1.0)) };
            vstoich.entry(k).or_default().push((*sg, a));
        }
    }
    let cf: HashMap<String, String> = m.species.iter().filter_map(|s| s.cf.clone().or_else(|| m.cf.clone()).map(|c| (s.id.clone(), c))).collect();
    let cvar: HashSet<&str> = m.comps.iter().map(|c| c.id.as_str()).filter(|c| changing.contains(*c)).collect();
    // boundary species in a variable compartment keep their amount (parameter); [S] = amount / V
    let bamt: HashSet<String> = m.species.iter().filter(|s| !floating.contains(&s.id) && !rules.contains_key(&s.id)
        && !rate_rules.contains_key(&s.id) && !s.hosu && cvar.contains(s.comp.as_str())).map(|s| s.id.clone()).collect();
    for s in &bamt { let v = params.vals[s] * vol[s]; params.set(s, v); }
    let general: HashSet<String> = species.iter().filter(|s| cvar.contains(comp_of[*s].as_str())
        || vstoich.keys().any(|k| &k.0 == *s) || cf.get(*s).is_some_and(|c| changing.contains(c))).cloned().collect();

    let mut kl_math = HashMap::new();
    let mut local = HashMap::new();
    for rx in &m.reactions {
        kl_math.insert(rx.id.clone(), rx.kl.as_ref().and_then(|k| k.0.clone()));
        let mut l = vec![];
        if let Some((_, ps)) = &rx.kl {
            for (p, v) in ps {
                let n = format!("{}__{p}", rx.id);
                params.set(&n, *v);
                l.push((p.clone(), n));
            }
        }
        local.insert(rx.id.clone(), l);
    }
    let known: HashSet<String> = floating.iter().cloned().chain(params.names.iter().cloned())
        .chain(["time".to_string(), "__time__".to_string()]).collect();
    let mut c = C {
        fd: m.fdefs(), rules, rule_order, rate_rules, rate_order, comp_of, hosu, floating, bamt, kl_math, local,
        ias: m.ias.clone(), pre, params, known, delays: vec![], memo: HashMap::new(), busy: HashSet::new(),
        rows, reactions, nm, vstoich, cf,
    };

    let rates = (0..reactions.len()).map(|j| c.rate(j)).collect::<Result<Vec<_>, _>>()?;
    let nr = reactions.len();
    let mut stoich: Vec<f64> = vec![];
    for s in &species {
        let i = rows.iter().position(|r| r == s).unwrap();
        let mut row = nm[i * nr..(i + 1) * nr].to_vec();
        if let Some(f) = c.cf.get(s).filter(|_| !general.contains(s)) { // constant conversion factor: scale the row
            let k = *c.params.vals.get(f).ok_or_else(|| inval(format!("conversion factor {f} has no value")))?;
            row.iter_mut().for_each(|x| *x *= k);
        }
        if c.rate_rules.contains_key(s) || general.contains(s) { row.iter_mut().for_each(|x| *x = 0.0); }
        stoich.extend(row);
    }
    let mut rate_list = vec![];
    for v in c.rate_order.clone() {
        if !c.floating.contains(&v) {
            let r = c.rate_rules[&v].clone();
            rate_list.push((v.clone(), c.formula(&r, &format!("rate rule for {v}"))?));
        }
    }
    for s in &species {
        if c.rate_rules.contains_key(s) || general.contains(s) { rate_list.push((s.clone(), c.dconc(s)?)); }
    }
    let mut out_rules = vec![];
    for v in c.rule_order.clone() {
        let r = c.rules[&v].clone();
        let f = c.formula(&r, &format!("assignment rule for {v}"))?;
        out_rules.push((v.clone(), c.conc_of(&v, f)?));
    }

    let mut events = vec![];
    for ev in &m.events {
        let what = format!("event {}", ev.id.clone().unwrap_or_else(|| events.len().to_string()));
        let (mut asg, mut div) = (vec![], vec![]);
        for (v, e) in &ev.asg {
            let f = c.formula(e, &what)?;
            if ev.uvftt && ev.delay.is_some() && c.hosu.contains(v) && c.floating.contains(v) {
                asg.push((v.clone(), f)); // an amount fixed at trigger time: / volume when it executes
                div.push(Some(c.comp_formula(v, &what)?));
                continue;
            }
            let f = if c.bamt.contains(v) { format!("({f}) * ({})", c.comp_formula(v, &what)?) } else { c.conc_of(v, f)? };
            asg.push((v.clone(), f));
            div.push(None);
        }
        // a compartment resized directly or through its assignment rule keeps its species' amounts: the
        // engine scales their concentrations by V_before / V_after when the event executes
        let new: HashSet<&str> = asg.iter().map(|a| a.0.as_str()).collect();
        let mut rescale = vec![];
        for comp in m.comps.iter().map(|x| &x.id).filter(|x| cvar.contains(x.as_str())) {
            let old = c.formula(&name(comp), &what)?;
            if words(&old).any(|w| new.contains(w)) {
                let held = c.rate_order.iter().filter(|x| sp.contains_key(x.as_str()) && !c.floating.contains(*x) && !c.hosu.contains(*x));
                rescale.extend(species.iter().chain(held).filter(|s| &c.comp_of[*s] == comp).map(|s| (s.clone(), old.clone())));
            }
        }
        let f = |c: &mut C, e: &Option<E>| e.as_ref().map(|e| c.formula(e, &what)).transpose();
        events.push(EventDef {
            trigger: c.formula(&ev.trigger, &what)?, assign: asg, delay: f(&mut c, &ev.delay)?, priority: f(&mut c, &ev.priority)?,
            initial_value: ev.init, persistent: ev.persistent, use_trigger_values: ev.uvftt, div, rescale,
        });
    }

    info.amount_params = m.species.iter().filter(|s| c.params.vals.contains_key(&s.id) && (c.hosu.contains(&s.id) || c.bamt.contains(&s.id)))
        .map(|s| s.id.clone()).collect();
    info.comp_of = c.comp_of.clone();
    info.general = general;
    info.changing = changing;
    let mut alg: Vec<String> = dae.iter().cloned().collect();
    alg.sort();
    let pvals = c.params.names.iter().map(|p| c.params.vals[p]).collect();
    let x = Rules {
        rate: rate_list, assign: out_rules, events,
        delays: c.delays.iter().map(|((e, tau, pre), n)| (n.clone(), e.clone(), tau.clone(), pre.clone())).collect(),
        algebraic: alg,
    };
    let init_c = species.iter().map(|s| conc[s]).collect();
    let vols = species.iter().map(|s| vol[s]).collect();
    Model::new(species, init_c, vols, reactions.to_vec(), rates, stoich, c.params.names.clone(), pvals, x).map_err(inval)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc(body: &str) -> String {
        format!(r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2"><model id="m">{body}</model></sbml>"#)
    }
    const MATH: &str = r#"xmlns="http://www.w3.org/1998/Math/MathML""#;

    fn at(m: &Model, t: f64, ids: &[&str]) -> Vec<f64> {
        let outs = m.output_regs(&ids.iter().map(|s| s.to_string()).collect::<Vec<_>>()).unwrap();
        let times = [0.0, t];
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let mut out = vec![0.0; 2 * outs.len()];
        assert!(m.simulate(&mut w, &times, 1e-10, 1e-12, &outs, &mut out));
        out[outs.len()..].to_vec()
    }

    #[test]
    fn hosu_species_amount_rate() {
        // S has only substance units in V = 2: dS/dt = -k S (amount), so [S] = S / V
        let s = doc(&format!(r#"<listOfCompartments><compartment id="V" size="2" constant="true"/></listOfCompartments>
<listOfSpecies><species id="S" compartment="V" initialAmount="10" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters><parameter id="k" value="0.5" constant="true"/></listOfParameters>
<listOfReactions><reaction id="J" reversible="false"><listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
<kineticLaw><math {MATH}><apply><times/><ci>k</ci><ci>S</ci></apply></math></kineticLaw></reaction></listOfReactions>"#));
        let (m, info) = from_sbml(&s).unwrap();
        assert!(info.hosu.contains("S"));
        assert_eq!(info.values["S"], 10.0);
        let v = at(&m, 1.0, &["S", "J"]);
        assert!((v[0] - 5.0 * (-0.5f64).exp()).abs() < 1e-7, "{v:?}"); // concentration
        assert!((v[1] - 0.5 * 10.0 * (-0.5f64).exp()).abs() < 1e-6, "{v:?}"); // k * amount
    }

    #[test]
    fn constant_species_is_boundary() {
        // constant="true" without boundaryCondition (BioModels 33): roadrunner lists it as a boundary species
        let s = doc(&format!(r#"<listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
<species id="E" compartment="c" initialConcentration="2" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="true"/></listOfSpecies>
<listOfReactions><reaction id="J" reversible="false"><listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
<listOfModifiers><modifierSpeciesReference species="E"/></listOfModifiers>
<kineticLaw><math {MATH}><apply><times/><ci>E</ci><ci>S</ci></apply></math></kineticLaw></reaction></listOfReactions>"#));
        let (_, info) = from_sbml(&s).unwrap();
        assert_eq!(info.ids[0], vec!["S".to_string()]);
        assert_eq!(info.ids[1], vec!["E".to_string()]);
    }

    #[test]
    fn fast_reaction_equilibrium() {
        // A <-> B fast (kf A - kr B), A + B = 3: equilibrium B = 2 A -> A = 1, B = 2
        let s = doc(&format!(r#"<listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" initialConcentration="3" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
<species id="B" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfReactions><reaction id="F" reversible="true" fast="true"><listOfReactants><speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>
<listOfProducts><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>
<kineticLaw><math {MATH}><apply><minus/><apply><times/><cn>2</cn><ci>A</ci></apply><ci>B</ci></apply></math></kineticLaw></reaction></listOfReactions>"#))
            .replace("level3/version2/core\" level=\"3\" version=\"2\"", "level3/version1/core\" level=\"3\" version=\"1\"");
        let (m, info) = from_sbml(&s).unwrap();
        assert!(info.rule_vars.contains(&"A".to_string()) && info.ids[4].is_empty());
        let v = at(&m, 1.0, &["A", "B"]);
        assert!((v[0] - 1.0).abs() < 1e-8 && (v[1] - 2.0).abs() < 1e-8, "{v:?}");
    }

    #[test]
    fn algebraic_rule_solved() {
        // x' = 1, x + y - 5 = 0 -> y = 5 - x
        let s = doc(&format!(r#"<listOfParameters><parameter id="x" value="0" constant="false"/><parameter id="y" value="5" constant="false"/></listOfParameters>
<listOfRules><rateRule variable="x"><math {MATH}><cn>1</cn></math></rateRule>
<algebraicRule><math {MATH}><apply><minus/><apply><plus/><ci>x</ci><ci>y</ci></apply><cn>5</cn></apply></math></algebraicRule></listOfRules>"#));
        let (m, info) = from_sbml(&s).unwrap();
        assert!(info.rule_vars.contains(&"y".to_string()));
        let v = at(&m, 2.0, &["x", "y"]);
        assert!((v[0] - 2.0).abs() < 1e-8 && (v[1] - 3.0).abs() < 1e-8, "{v:?}");
    }

    #[test]
    fn delay_and_rateof() {
        // p' = 1 from 0; q = delay(p, 1) (0 before t = 1); r' = rateOf(z), z = 3 p (an assignment rule)
        let s = doc(&format!(r#"<listOfParameters><parameter id="p" value="0" constant="false"/><parameter id="q" constant="false"/>
<parameter id="r" value="0" constant="false"/><parameter id="z" constant="false"/></listOfParameters>
<listOfRules><rateRule variable="p"><math {MATH}><cn>1</cn></math></rateRule>
<assignmentRule variable="q"><math {MATH}><apply><csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/delay">delay</csymbol><ci>p</ci><cn>1</cn></apply></math></assignmentRule>
<assignmentRule variable="z"><math {MATH}><apply><times/><cn>3</cn><ci>p</ci></apply></math></assignmentRule>
<rateRule variable="r"><math {MATH}><apply><csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/rateOf">rateOf</csymbol><ci>z</ci></apply></math></rateRule></listOfRules>"#));
        let (m, _) = from_sbml(&s).unwrap();
        let v = at(&m, 3.0, &["p", "q", "r"]);
        assert!((v[0] - 3.0).abs() < 1e-8 && (v[1] - 2.0).abs() < 1e-6 && (v[2] - 9.0).abs() < 1e-6, "{v:?}");
    }

    #[test]
    fn comp_needs_flattening() {
        let s = doc("").replace("<sbml ", r#"<sbml xmlns:comp="http://www.sbml.org/sbml/level3/version1/comp/version1" comp:required="true" "#);
        assert_eq!(from_sbml(&s).err(), Some(Error::NotImplemented("comp: flatten first".into())));
    }
}

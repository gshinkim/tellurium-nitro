//! SBML math as a small tree: conversion from sbml-io MathML, infix printing (the syntax
//! `expr::parse` and libsbml's L3 parser read), evaluation at t0, substitution and derivatives.
use super::{nimp, Error};
use crate::expr::{parse, Dag, Op as D};
use sbml_io::{Constant, Math, Symbol};
use std::collections::{HashMap, HashSet};

/// libsbml's value of the avogadro csymbol.
pub const AVOGADRO: f64 = 6.02214179e23;

#[derive(Clone, Debug, PartialEq)]
pub enum E {
    Num(f64),
    Name(String),
    Time,
    Avogadro,
    /// Built-in: MathML operator names ("plus", "lt", "piecewise", "root" [degree, x], "log" [base, x],
    /// "rateOf", "delay", distrib names).
    Op(String, Vec<E>),
    /// Call of a function definition.
    Call(String, Vec<E>),
}

pub type Fdefs = HashMap<String, (Vec<String>, E)>;

pub fn num(x: f64) -> E { E::Num(x) }
pub fn name(s: &str) -> E { E::Name(s.to_string()) }
pub fn op(f: &str, a: Vec<E>) -> E { E::Op(f.to_string(), a) }
pub fn add(a: E, b: E) -> E { op("plus", vec![a, b]) }
pub fn sub(a: E, b: E) -> E { op("minus", vec![a, b]) }
pub fn mul(a: E, b: E) -> E { op("times", vec![a, b]) }
pub fn div(a: E, b: E) -> E { op("divide", vec![a, b]) }
pub fn neg(a: E) -> E { op("minus", vec![a]) }
/// a1 + a2 + ... (0 if empty)
pub fn sum(a: Vec<E>) -> E { if a.len() == 1 { a.into_iter().next().unwrap() } else { op("plus", a) } }

impl E {
    pub fn kids(&self) -> &[E] {
        match self { E::Op(_, a) | E::Call(_, a) => a, _ => &[] }
    }

    /// Ids referenced (ci names; not function names or csymbols).
    pub fn names(&self) -> HashSet<String> {
        let mut out = HashSet::new();
        self.names_into(&mut out);
        out
    }
    pub fn names_into(&self, out: &mut HashSet<String>) {
        if let E::Name(n) = self { out.insert(n.clone()); }
        for k in self.kids() { k.names_into(out); }
    }

    /// Copy with every name in `m` replaced at once by its expression.
    pub fn subst(&self, m: &HashMap<String, E>) -> E {
        match self {
            E::Name(n) => m.get(n).cloned().unwrap_or_else(|| self.clone()),
            E::Op(f, a) => E::Op(f.clone(), a.iter().map(|x| x.subst(m)).collect()),
            E::Call(f, a) => E::Call(f.clone(), a.iter().map(|x| x.subst(m)).collect()),
            _ => self.clone(),
        }
    }
    pub fn sub1(&self, x: &str, f: &E) -> E { self.subst(&HashMap::from([(x.to_string(), f.clone())])) }
}

/// sbml-io MathML -> E (constants become numbers, semantics unwrapped).
pub fn conv(m: &Math) -> Result<E, Error> {
    let all = |a: &[Math]| a.iter().map(conv).collect::<Result<Vec<_>, _>>();
    Ok(match m {
        Math::Number { value, .. } => E::Num(value.value()),
        Math::Ident(s) => E::Name(s.clone()),
        Math::Symbol(s) if s.url == Symbol::TIME => E::Time,
        Math::Symbol(s) if s.url == Symbol::AVOGADRO => E::Avogadro,
        Math::Symbol(s) => return Err(nimp(format!("csymbol {}", s.url))),
        Math::Constant(c) => E::Num(match c {
            Constant::True => 1.0, Constant::False => 0.0, Constant::Pi => std::f64::consts::PI,
            Constant::ExponentialE => std::f64::consts::E, Constant::NotANumber => f64::NAN,
            Constant::Infinity => f64::INFINITY,
        }),
        Math::Apply(o, a) => E::Op(o.tag().to_string(), all(a)?),
        Math::Call(f, a) => E::Call(f.clone(), all(a)?),
        Math::CallSymbol(s, a) => {
            let f = if s.url == Symbol::DELAY { "delay" } else if s.url == Symbol::RATE_OF { "rateOf" } else {
                let last = s.url.rsplit('/').next().unwrap_or("");
                match crate::rand::NAMES.iter().find(|(n, _)| *n == last) {
                    Some((n, _)) if s.url.contains("distrib") => n,
                    _ => return Err(nimp(format!("csymbol {}", s.url))),
                }
            };
            E::Op(f.to_string(), all(a)?)
        }
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            let mut a = vec![];
            if let Some(q) = q { a.push(conv(q)?); }
            a.push(conv(arg)?);
            E::Op(if matches!(m, Math::Root { .. }) { "root" } else { "log" }.into(), a)
        }
        Math::Piecewise { pieces, otherwise } => {
            let mut a = vec![];
            for (v, c) in pieces { a.push(conv(v)?); a.push(conv(c)?); }
            if let Some(o) = otherwise { a.push(conv(o)?); }
            E::Op("piecewise".into(), a)
        }
        Math::Semantics { math, .. } => conv(math)?,
        Math::Lambda { .. } => return Err(super::inval("lambda outside a function definition".into())),
        Math::Raw(x) => return Err(super::inval(format!("unsupported MathML {}", &x[..x.len().min(80)]))),
    })
}

/// Number text that round-trips; non-finite as (1/0) etc. so no model id can shadow it.
pub fn fnum(x: f64) -> String {
    if x.is_nan() { "(0/0)".into() } else if x.is_infinite() { if x > 0.0 { "(1/0)".into() } else { "(-1/0)".into() } }
    else if x.is_sign_negative() { format!("(-{:?})", -x) } else { format!("{x:?}") }
}

/// Infix text; names in `rep` are printed as "(text)".
pub fn l3(e: &E, rep: &HashMap<String, String>) -> String {
    let p = |x: &E| l3(x, rep);
    let join = |a: &[E], s: &str| a.iter().map(|x| l3(x, rep)).collect::<Vec<_>>().join(s);
    match e {
        E::Num(x) => fnum(*x),
        E::Name(n) => rep.get(n).map(|r| format!("({r})")).unwrap_or_else(|| n.clone()),
        E::Time => "time".into(),
        E::Avogadro => "avogadro".into(),
        E::Call(f, a) => format!("{f}({})", join(a, ", ")),
        E::Op(f, a) => match (f.as_str(), a.len()) {
            ("plus", 0) => "0".into(),
            ("times", 0) | ("and", 0) => "1".into(),
            ("or", 0) => "0".into(),
            ("plus", _) => format!("({})", join(a, " + ")),
            ("times", _) => format!("({})", join(a, " * ")),
            ("and", _) => format!("({})", join(a, " && ")),
            ("or", _) => format!("({})", join(a, " || ")),
            ("minus", 1) => format!("(-{})", p(&a[0])),
            ("minus", 2) => format!("({} - {})", p(&a[0]), p(&a[1])),
            ("divide", 2) => format!("({} / {})", p(&a[0]), p(&a[1])),
            ("power", 2) => format!("({} ^ {})", p(&a[0]), p(&a[1])),
            ("not", 1) => format!("(!{})", p(&a[0])),
            ("lt" | "gt" | "leq" | "geq" | "eq" | "neq", n) => {
                if n < 2 { return "1".into(); }
                let s = match f.as_str() { "lt" => "<", "gt" => ">", "leq" => "<=", "geq" => ">=", "eq" => "==", _ => "!=" };
                let w: Vec<String> = a.windows(2).map(|w| format!("{} {s} {}", p(&w[0]), p(&w[1]))).collect();
                format!("({})", w.join(" && "))
            }
            _ => format!("{f}({})", join(a, ", ")),
        },
    }
}

pub fn s(e: &E) -> String { l3(e, &HashMap::new()) }

fn draw(f: &str, a: &[f64]) -> f64 {
    let k = crate::rand::NAMES.iter().position(|(n, _)| *n == f).unwrap();
    let one = matches!(f, "cauchy" | "laplace") && (a.len() == 1 || a.len() == 3);
    let np = if one { 1 } else { crate::rand::NAMES[k].1 };
    let (lo, hi) = if a.len() == np + 2 { (a[np], a[np + 1]) } else { (f64::NEG_INFINITY, f64::INFINITY) };
    let at = |i: usize| a.get(i).copied().unwrap_or(f64::NAN);
    let (x, y) = if one { (0.0, at(0)) } else if np == 1 { (at(0), f64::NAN) } else { (at(0), at(1)) };
    crate::rand::draw(k as u8, x, y, lo, hi)
}

/// Value at time 0; names through `get`. NotImplemented if unsupported.
pub fn eval(e: &E, get: &mut dyn FnMut(&str) -> Result<f64, Error>, fd: &Fdefs) -> Result<f64, Error> {
    let b = |x: bool| if x { 1.0 } else { 0.0 };
    let (f, a) = match e {
        E::Num(x) => return Ok(*x),
        E::Name(n) => return get(n),
        E::Time => return Ok(0.0),
        E::Avogadro => return Ok(AVOGADRO),
        E::Call(f, a) => {
            let Some((args, body)) = fd.get(f) else { return Err(nimp(format!("initial value of {}", s(e)))) };
            let v = a.iter().map(|x| eval(x, get, fd)).collect::<Result<Vec<_>, _>>()?;
            return eval(body, &mut |x| match args.iter().position(|y| y == x) { Some(i) => Ok(v[i]), None => get(x) }, fd);
        }
        E::Op(f, a) => (f.as_str(), a),
    };
    if f == "rateOf" {
        return match a.first() { Some(E::Name(x)) => get(&format!("__rateof__{x}")), _ => Err(nimp(format!("initial value of {}", s(e)))) };
    }
    if f == "delay" { return eval(&a[0], get, fd); } // before the start, every value is its initial one
    if f == "piecewise" {
        let n = a.len();
        for i in (0..n.saturating_sub(1)).step_by(2) {
            if eval(&a[i + 1], get, fd)? != 0.0 { return eval(&a[i], get, fd); }
        }
        return if n % 2 == 1 { eval(&a[n - 1], get, fd) } else { Ok(f64::NAN) };
    }
    let v = a.iter().map(|x| eval(x, get, fd)).collect::<Result<Vec<_>, _>>()?;
    let n = v.len();
    let one = |g: fn(f64) -> f64| if n == 1 { Ok(g(v[0])) } else { Err(nimp(format!("initial value of {}", s(e)))) };
    let rel = |g: fn(f64, f64) -> bool| Ok(b(v.windows(2).all(|w| g(w[0], w[1]))));
    match f {
        "plus" => Ok(v.iter().sum()),
        "minus" if n == 1 => Ok(-v[0]),
        "minus" if n == 2 => Ok(v[0] - v[1]),
        "times" => Ok(v.iter().product()),
        "divide" if n == 2 => Ok(v[0] / v[1]),
        "power" if n == 2 => Ok(v[0].powf(v[1])),
        "root" if n == 1 => Ok(v[0].sqrt()),
        "root" if n == 2 => Ok(v[1].powf(1.0 / v[0])),
        "log" if n == 1 => Ok(v[0].log10()),
        "log" if n == 2 => Ok(v[1].ln() / v[0].ln()),
        "abs" => one(f64::abs), "exp" => one(f64::exp), "ln" => one(f64::ln), "floor" => one(f64::floor),
        "ceiling" => one(f64::ceil), "sin" => one(f64::sin), "cos" => one(f64::cos), "tan" => one(f64::tan),
        "arcsin" => one(f64::asin), "arccos" => one(f64::acos), "arctan" => one(f64::atan),
        "sinh" => one(f64::sinh), "cosh" => one(f64::cosh), "tanh" => one(f64::tanh),
        "sec" => one(|x| 1.0 / x.cos()), "csc" => one(|x| 1.0 / x.sin()), "cot" => one(|x| 1.0 / x.tan()),
        "sech" => one(|x| 1.0 / x.cosh()), "csch" => one(|x| 1.0 / x.sinh()), "coth" => one(|x| 1.0 / x.tanh()),
        "arcsec" => one(|x| (1.0 / x).acos()), "arccsc" => one(|x| (1.0 / x).asin()), "arccot" => one(|x| (1.0 / x).atan()),
        "arcsinh" => one(f64::asinh), "arccosh" => one(f64::acosh), "arctanh" => one(f64::atanh),
        "arcsech" => one(|x| (1.0 / x).acosh()), "arccsch" => one(|x| (1.0 / x).asinh()), "arccoth" => one(|x| (1.0 / x).atanh()),
        "factorial" => one(crate::expr::factorial),
        "eq" => rel(|x, y| x == y), "neq" => rel(|x, y| x != y), "gt" => rel(|x, y| x > y),
        "lt" => rel(|x, y| x < y), "geq" => rel(|x, y| x >= y), "leq" => rel(|x, y| x <= y),
        "and" => Ok(b(v.iter().all(|&x| x != 0.0))),
        "or" => Ok(b(v.iter().any(|&x| x != 0.0))),
        "xor" => Ok((v.iter().filter(|&&x| x != 0.0).count() % 2) as f64),
        "not" if n == 1 => Ok(b(v[0] == 0.0)),
        "quotient" | "rem" if n == 2 => { // C truncation, as libroadrunner
            let q = (v[0] / v[1]).trunc();
            Ok(if f == "quotient" { q } else { v[0] - v[1] * q })
        }
        "implies" if n == 2 => Ok(b(v[0] == 0.0 || v[1] != 0.0)),
        "min" if n > 0 => Ok(v.iter().copied().fold(f64::INFINITY, |x, y| if y < x || y.is_nan() { y } else { x })),
        "max" if n > 0 => Ok(v.iter().copied().fold(f64::NEG_INFINITY, |x, y| if y > x || y.is_nan() { y } else { x })),
        _ if crate::rand::NAMES.iter().any(|(nm, _)| *nm == f) => Ok(draw(f, &v)),
        _ => Err(nimp(format!("initial value of {}", s(e)))),
    }
}

/// d e / d x for every name x in e (expr.rs's symbolic derivative), as (x, derivative).
pub fn derivatives(e: &E) -> Result<Vec<(String, E)>, Error> {
    let mut xs: Vec<String> = e.names().into_iter().collect();
    xs.sort();
    let idx: HashMap<&str, u32> = xs.iter().enumerate().map(|(i, x)| (x.as_str(), i as u32)).collect();
    let mut dag = Dag::default();
    let root = parse(&s(e), &mut dag, |n| idx.get(n).map(|&k| D::Var(k)))
        .map_err(|m| nimp(format!("cannot differentiate: {m}")))?;
    let mut out = vec![];
    for (k, x) in xs.iter().enumerate() {
        let d = dag.diff(root, D::Var(k as u32));
        out.push((x.clone(), from_dag(&dag, d, &xs)?));
    }
    Ok(out)
}

fn from_dag(dag: &Dag, n: u32, xs: &[String]) -> Result<E, Error> {
    let r = |m: u32| from_dag(dag, m, xs);
    let u = |f: &str, a: u32| -> Result<E, Error> { Ok(op(f, vec![r(a)?])) };
    let bi = |f: &str, a: u32, b: u32| -> Result<E, Error> { Ok(op(f, vec![r(a)?, r(b)?])) };
    use D::*;
    Ok(match dag.nodes[n as usize] {
        Const(x) => num(f64::from_bits(x)),
        Var(k) => name(&xs[k as usize]),
        Par(_) | Time => return Err(nimp("derivative".into())),
        Neg(a) => u("minus", a)?, Add(a, b) => bi("plus", a, b)?, Sub(a, b) => bi("minus", a, b)?,
        Mul(a, b) => bi("times", a, b)?, Div(a, b) => bi("divide", a, b)?, Pow(a, b) => bi("power", a, b)?,
        PowI(a, k) => op("power", vec![r(a)?, num(k as f64)]),
        Exp(a) => u("exp", a)?, Ln(a) => u("ln", a)?, Log10(a) => u("log", a)?, Sqrt(a) => u("root", a)?,
        Abs(a) => u("abs", a)?, Floor(a) => u("floor", a)?, Ceil(a) => u("ceiling", a)?,
        Sign(a) => { let x = r(a)?; op("piecewise", vec![num(1.0), op("gt", vec![x.clone(), num(0.0)]), num(-1.0), op("lt", vec![x, num(0.0)]), num(0.0)]) }
        Sin(a) => u("sin", a)?, Cos(a) => u("cos", a)?, Tan(a) => u("tan", a)?,
        Asin(a) => u("arcsin", a)?, Acos(a) => u("arccos", a)?, Atan(a) => u("arctan", a)?, Fact(a) => u("factorial", a)?,
        Lt(a, b) => bi("lt", a, b)?, Le(a, b) => bi("leq", a, b)?, Eq(a, b) => bi("eq", a, b)?, Ne(a, b) => bi("neq", a, b)?,
        And(a, b) => bi("and", a, b)?, Or(a, b) => bi("or", a, b)?, Xor(a, b) => bi("xor", a, b)?, Not(a) => u("not", a)?,
        Sel(c, x, y) => op("piecewise", vec![r(x)?, r(c)?, r(y)?]),
        Min(a, b) => bi("min", a, b)?, Max(a, b) => bi("max", a, b)?,
        Rand(..) => return Err(nimp("derivative of a random draw".into())),
    })
}

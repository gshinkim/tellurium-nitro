//! Dimensional checks for conversions through SBML versions requiring consistent units.
use super::convert::{for_each_math, MathContext};
use sbml_io::{Math, Model, Op, Rule, Symbol};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug, Default, PartialEq)]
struct Dimension(BTreeMap<String, f64>);

impl Dimension {
    fn base(kind: &str) -> Self {
        let (kind, exponent) = match kind {
            "dimensionless" | "radian" | "steradian" => return Self::default(),
            "meter" | "length" => ("metre", 1.0),
            "area" => ("metre", 2.0),
            "liter" | "litre" | "volume" => ("metre", 3.0),
            "substance" => ("mole", 1.0),
            "time" => ("second", 1.0),
            "gram" => ("kilogram", 1.0),
            "hertz" | "becquerel" => ("second", -1.0),
            _ => (kind, 1.0),
        };
        Self(BTreeMap::from([(kind.to_owned(), exponent)]))
    }

    fn combine(mut self, other: &Self, power: f64) -> Self {
        for (kind, exponent) in &other.0 {
            *self.0.entry(kind.clone()).or_default() += power * exponent;
        }
        self.0.retain(|_, exponent| *exponent != 0.0);
        self
    }
}

struct Scope {
    ids: HashMap<String, Option<Dimension>>,
    units: HashMap<String, Dimension>,
    time: Option<Dimension>,
    implicit_time: Option<Dimension>,
}

impl Scope {
    fn unit(&self, name: &str) -> Dimension {
        self.units.get(name).cloned().unwrap_or_else(|| Dimension::base(name))
    }

    fn expression(&self, math: &Math, locals: &HashMap<&str, Option<Dimension>>) -> Result<Option<Dimension>, ()> {
        let infer = |math| self.expression(math, locals);
        let scalar = || Some(Dimension::default());
        let compatible = |args: &[Option<Dimension>]| {
            let known = args.iter().flatten().next().cloned();
            if args.iter().all(Option::is_some) && args.iter().any(|u| *u != known) {
                return Err(());
            }
            Ok(known)
        };
        Ok(match math {
            Math::Number { units, .. } => units.as_ref().map(|u| self.unit(u)),
            Math::Ident(id) => locals.get(id.as_str()).or_else(|| self.ids.get(id)).cloned().flatten(),
            Math::Symbol(s) if s.url == Symbol::TIME => self.time.clone(),
            Math::Apply(op, args) => {
                let args_units: Vec<_> = args.iter().map(infer).collect::<Result<_, _>>()?;
                match op {
                    Op::Plus | Op::Minus | Op::Min | Op::Max => compatible(&args_units)?,
                    Op::Eq | Op::Neq | Op::Lt | Op::Gt | Op::Leq | Op::Geq => {
                        let comparison_units: Vec<_> = args.iter().zip(args_units).map(|(arg, units)| {
                            if matches!(arg, Math::Symbol(s) if s.url == Symbol::TIME) {
                                units.or_else(|| self.implicit_time.clone())
                            } else { units }
                        }).collect();
                        compatible(&comparison_units)?;
                        scalar()
                    }
                    Op::Times | Op::Divide if args_units.iter().all(Option::is_some) => {
                        Some(args_units.iter().flatten().enumerate().fold(Dimension::default(), |u, (i, x)| {
                            u.combine(x, if *op == Op::Divide && i > 0 { -1.0 } else { 1.0 })
                        }))
                    }
                    Op::Power => match (args_units.first().cloned().flatten(), args.get(1)) {
                        (Some(base), Some(Math::Number { value, .. })) => Some(Dimension::default().combine(&base, value.value())),
                        _ => None,
                    },
                    Op::Exp | Op::Ln | Op::Sin | Op::Cos | Op::Tan | Op::Sec | Op::Csc | Op::Cot |
                    Op::Sinh | Op::Cosh | Op::Tanh | Op::Sech | Op::Csch | Op::Coth |
                    Op::ArcSin | Op::ArcCos | Op::ArcTan | Op::ArcSec | Op::ArcCsc | Op::ArcCot |
                    Op::ArcSinh | Op::ArcCosh | Op::ArcTanh | Op::ArcSech | Op::ArcCsch | Op::ArcCoth => {
                        if args_units.iter().flatten().any(|u| !u.0.is_empty()) { return Err(()); }
                        scalar()
                    }
                    Op::Abs | Op::Floor | Op::Ceiling => args_units.into_iter().next().flatten(),
                    _ => None,
                }
            }
            Math::Log { base, arg } => {
                let base = base.as_ref().map(|b| infer(b)).transpose()?.flatten();
                if [infer(arg)?, base].into_iter().flatten().any(|u| !u.0.is_empty()) { return Err(()); }
                scalar()
            }
            Math::Root { degree, arg } => {
                let degree = match degree.as_deref() {
                    None => 2.0,
                    Some(Math::Number { value, .. }) => value.value(),
                    Some(expr) => { infer(expr)?; return Ok(None); }
                };
                infer(arg)?.map(|u| Dimension::default().combine(&u, degree.recip()))
            }
            Math::CallSymbol(s, args) if s.url == Symbol::DELAY => {
                for a in args.iter().skip(1) { infer(a)?; }
                args.first().map(infer).transpose()?.flatten()
            }
            Math::Call(_, args) | Math::CallSymbol(_, args) => {
                for a in args { infer(a)?; }
                None
            }
            Math::Piecewise { pieces, otherwise } => {
                let mut values = vec![];
                for (v, condition) in pieces {
                    infer(condition)?;
                    values.push(infer(v)?);
                }
                if let Some(v) = otherwise { values.push(infer(v)?); }
                compatible(&values)?
            }
            Math::Semantics { math, .. } => infer(math)?,
            _ => None,
        })
    }
}

pub(super) fn is_consistent(m: &mut Model, from_level: u32, to_level: u32) -> bool {
    let units = m.unit_definitions.iter().map(|d| {
        let unit = d.units.iter().fold(Dimension::default(), |acc, u| {
            acc.combine(&Dimension::base(&u.kind), u.exponent.unwrap_or(1.0))
        });
        (d.id.clone(), unit)
    }).collect();
    let mut scope = Scope { ids: HashMap::new(), units, time: None, implicit_time: None };
    let declared = |name: Option<&str>, default: &str, scope: &Scope| {
        name.or((from_level == 2).then_some(default)).map(|u| scope.unit(u))
    };
    scope.time = declared(m.time_units.as_deref(), "time", &scope);
    let default_substance = m.substance_units.as_deref().or((from_level == 2).then_some("substance"));
    for p in &m.parameters {
        scope.ids.insert(p.id.clone(), p.units.as_ref().map(|u| scope.unit(u)));
    }
    for c in &m.compartments {
        let units = match c.spatial_dimensions.unwrap_or(3.0) {
            3.0 => declared(m.volume_units.as_deref(), "volume", &scope),
            2.0 => declared(m.area_units.as_deref(), "area", &scope),
            1.0 => declared(m.length_units.as_deref(), "length", &scope),
            _ => Some(Dimension::default()),
        };
        scope.ids.insert(c.id.clone(), c.units.as_ref().map(|u| scope.unit(u)).or(units));
    }
    for s in &m.species {
        let substance = s.substance_units.as_deref().or(default_substance).map(|u| scope.unit(u));
        let units = if s.has_only_substance_units == Some(true) {
            substance
        } else {
            let volume = s.spatial_size_units.as_ref().map(|u| scope.unit(u))
                .or_else(|| scope.ids.get(&s.compartment).cloned().flatten());
            substance.zip(volume).map(|(s, v)| s.combine(&v, -1.0))
        };
        scope.ids.insert(s.id.clone(), units);
    }
    let flux = m.extent_units.as_deref().or(default_substance).map(|u| scope.unit(u))
        .zip(scope.time.clone()).map(|(s, t)| s.combine(&t, -1.0));
    let default_time = scope.unit("time");
    for r in &m.reactions {
        let rate = flux.clone().or_else(|| {
            (to_level == 2).then(|| scope.unit("substance").combine(&default_time, -1.0))
        });
        scope.ids.insert(r.id.clone(), rate);
    }

    let compare = |given: Option<Dimension>, expected: Option<Dimension>| {
        match (given, expected) {
            (Some(g), Some(e)) => g == e,
            _ => true,
        }
    };
    let globals = HashMap::new();
    for r in &m.rules {
        if let (Some(id), Some(math)) = (r.variable(), r.math()) {
            let mut expected = scope.ids.get(id).cloned().flatten();
            if matches!(r, Rule::Rate { .. }) {
                expected = expected.zip(scope.time.clone()).map(|(d, t)| d.combine(&t, -1.0));
            }
            if !scope.expression(math, &globals).is_ok_and(|u| compare(u, expected)) { return false; }
        }
    }
    for r in &m.reactions {
        if let Some(k) = &r.kinetic_law {
            let locals = k.parameters.iter().map(|p| (p.id.as_str(), p.units.as_ref().map(|u| scope.unit(u)))).collect();
            if let Some(math) = &k.math {
                if !scope.expression(math, &locals).is_ok_and(|u| compare(u, flux.clone())) { return false; }
            }
        }
    }
    for_each_math(m, |expr, locals, context| {
        if context == MathContext::Function || !locals.is_empty() { return Ok(()); }
        // A direct trigger comparison uses the default Level 2 clock units.
        scope.implicit_time = (context == MathContext::Trigger && to_level == 2).then(|| default_time.clone());
        scope.expression(expr, &globals).map(|_| ())
    }).is_ok()
}

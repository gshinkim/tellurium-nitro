//! Convert the SBML core between Level 2 and Level 3.

use super::{Init, M};
use sbml_io::{Constant, Document, Math, Model, Op, Rule, Symbol};
use std::collections::{HashMap, HashSet};

/// Convert core constructs without changing the model's equations. Reject constructs
/// for which this conversion has no equivalent in the target.
pub(crate) fn convert(doc: &mut Document, level: u32, version: u32) -> Result<(), String> {
    let bad = |what: &str| Err(format!("Error converting SBML to level {level} version {version}: {what}"));
    if !matches!((level, version), (2, 1..=5) | (3, 1..=2)) {
        return bad("unsupported level/version");
    }
    let (from_level, from_version) = (doc.level, doc.version);
    let l3v2 = (level, version) == (3, 2);
    let Some(m) = doc.model.as_mut() else {
        (doc.level, doc.version) = (level, version);
        return Ok(());
    };
    let typed = !m.compartment_types.is_empty() || !m.species_types.is_empty();
    if typed && (level == 3 || matches!(version, 1 | 5)) {
        return bad("compartment or species types");
    }
    let l2v1 = (level, version) == (2, 1);
    if !l3v2 && (from_level, from_version) == (3, 2) {
        remove_empty_math(m);
    }
    if l2v1 && !m.constraints.is_empty() {
        return bad("constraints require Level 2 Version 2 or later");
    }
    if !l3v2 && m.reactions.iter().any(|r| r.reactants.is_empty() && r.products.is_empty()) {
        return bad("reactions without reactants and products require Level 3 Version 2");
    }
    if l3v2 && m.reactions.iter().any(|r| r.fast == Some(true)) {
        return bad("fast reactions require an older SBML level or version");
    }
    if (level, version) != (2, 1) &&
        m.unit_definitions.iter().flat_map(|u| &u.units).any(|u| u.offset.is_some_and(|x| x != 0.0)) {
        return bad("nonzero unit offsets require Level 2 Version 1");
    }

    let refs: HashSet<_> = if level == 2 && from_level == 3 {
        species_references(m).filter_map(|sr| sr.id.clone()).collect()
    } else { HashSet::new() };
    if level == 2 {
        if m.conversion_factor.is_some() || m.species.iter().any(|s| s.conversion_factor.is_some()) {
            return bad("conversion factors");
        }
        if m.compartments.iter().any(|c| c.spatial_dimensions.is_some_and(|d| d.fract() != 0.0 || !(0.0..=3.0).contains(&d))) ||
            m.unit_definitions.iter().flat_map(|u| &u.units).any(|u| u.exponent.is_some_and(|e| e.fract() != 0.0)) {
            return bad("dimensions or unit exponents outside the Level 2 range");
        }
        if from_level == 3 && m.extent_units.is_some() && m.extent_units != m.substance_units {
            return bad("separate extent and substance units");
        }
        if m.compartments.iter().any(|c| c.spatial_dimensions == Some(0.0) &&
            (c.size.is_some() || m.species.iter().any(|s| s.compartment == c.id && s.initial_concentration.is_some()))) {
            return bad("size or initial concentration in a zero-dimensional compartment");
        }
        if m.events.iter().any(|e| e.assignments.is_empty() ||
            e.assignments.iter().any(|a| refs.contains(&a.variable))) {
            return bad("events without assignments or with species reference assignments");
        }
    }
    if l3v2 && m.reactions.iter().any(|r| r.kinetic_law.as_ref().is_some_and(|k| {
        k.parameters.iter().any(|p| r.reactants.iter().chain(&r.products).any(|sr| sr.species == p.id)
            || r.modifiers.iter().any(|sr| sr.species == p.id))
    })) {
        return bad("local parameter shadows a reaction species");
    }
    if l2v1 {
        let mut assigned: HashSet<_> = m.rules.iter().filter_map(|r|
            if let Rule::Assignment { variable, .. } = r {
                (!refs.contains(variable)).then_some(variable.as_str())
            } else { None }).collect();
        for r in &m.rules {
            if let Rule::Assignment { variable, math: Some(expr), .. } = r {
                if refs.contains(variable) { continue; }
                let mut invalid = false;
                expr.visit_idents(&mut |id| invalid |= assigned.contains(id));
                if invalid { return bad("assignment rule dependencies must precede their use"); }
                assigned.remove(variable.as_str());
            }
        }
    }
    let reactions: HashSet<_> = m.reactions.iter().map(|r| r.id.clone()).collect();

    if level == 2 && version <= 3 || from_level == 2 && from_version <= 3 {
        if !super::units::is_consistent(m, from_level, level) {
            return bad("inconsistent units");
        }
    }
    let functions: HashMap<_, _> = m.function_definitions.iter().filter_map(|f| {
        f.math.as_ref().map(|expr| (f.id.clone(), expr.clone()))
    }).collect();
    for_each_math(m, |expr, locals, context| {
        if !l3v2 && context != MathContext::Function {
            let expected = match context {
                MathContext::Trigger | MathContext::Constraint => Some(true),
                MathContext::Initial if l2v1 || uses_new_math(expr) => None,
                _ => Some(false),
            };
            let kind = math_type(expr, &functions, &HashMap::new());
            if kind.is_err() || expected.is_some_and(|x| kind == Ok(Some(!x))) {
                return bad("MathML uses boolean and numeric expressions interchangeably");
            }
        }
        visit_math(expr, &mut |node| {
            match node {
                Math::Apply(Op::Min | Op::Max | Op::Rem | Op::Quotient | Op::Implies, _)
                    if !l3v2 && context != MathContext::Initial => {
                    return bad("MathML operator requires Level 3 Version 2");
                }
                Math::CallSymbol(s, _) if s.url == Symbol::RATE_OF && !l3v2 => {
                    return bad("rateOf requires Level 3 Version 2");
                }
                Math::Ident(id) if !locals.contains(id) => {
                    if refs.contains(id) {
                        return bad("species reference used as a MathML symbol");
                    }
                    if l2v1 && reactions.contains(id) && context != MathContext::Initial {
                        return bad("reaction rates used as MathML symbols require Level 2 Version 2 or later");
                    }
                }
                Math::Symbol(s) if level == 2 && s.url == Symbol::AVOGADRO => return bad("Avogadro csymbol"),
                Math::Number { units, .. } if level == 2 => { *units = None; }
                _ => {}
            }
            Ok(())
        })
    })?;
    if l2v1 || !l3v2 && (from_level, from_version) == (3, 2) {
        fold_initial_assignments(m, l2v1).map_err(|e| format!(
            "Error converting SBML to level {level} version {version}: {e}"))?;
    }
    if level == 2 && from_level == 3 {
        lower_stoichiometry(m);
    } else if level == 3 && from_level == 2 {
        raise_stoichiometry(m)?;
    }

    for r in &mut m.reactions {
        if level == 3 {
            r.reversible.get_or_insert(true);
            r.fast = if l3v2 { None } else { Some(r.fast.unwrap_or(false)) };
        } else {
            r.compartment = None;
        }
        if let Some(k) = &mut r.kinetic_law {
            if (level, version) != (2, 1) { (k.time_units, k.substance_units) = (None, None); }
        }
        for sr in r.reactants.iter_mut().chain(&mut r.products) {
            if level == 3 && from_level == 2 {
                if sr.constant.unwrap_or(true) { sr.stoichiometry.get_or_insert(1.0); }
                sr.constant.get_or_insert(true);
            }
            if level == 2 {
                sr.constant = None;
                if version == 1 { (sr.id, sr.name) = (None, None); }
            }
        }
        if (level, version) == (2, 1) {
            for sr in &mut r.modifiers { (sr.id, sr.name) = (None, None); }
        }
    }

    for c in &mut m.compartments {
        if level == 3 && from_level == 2 {
            c.constant.get_or_insert(true);
            c.spatial_dimensions.get_or_insert(3.0);
            c.outside = None;
        }
        if level == 2 && from_level == 3 && c.units.is_none() {
            c.units = match c.spatial_dimensions.unwrap_or(3.0) as i32 {
                1 => m.length_units.clone(), 2 => m.area_units.clone(), 3 => m.volume_units.clone(), _ => None,
            };
        }
    }
    for s in &mut m.species {
        if level == 3 && from_level == 2 {
            s.has_only_substance_units.get_or_insert(false);
            s.boundary_condition.get_or_insert(false);
            s.constant.get_or_insert(false);
        }
        if level == 3 || version >= 3 { (s.charge, s.spatial_size_units) = (None, None); }
        if level == 2 && from_level == 3 && s.substance_units.is_none() {
            s.substance_units.clone_from(&m.substance_units);
        }
    }
    if level == 3 && from_level == 2 {
        for p in &mut m.parameters { p.constant.get_or_insert(true); }
    }

    for e in &mut m.events {
        if level == 3 {
            e.use_values_from_trigger_time.get_or_insert(true);
            e.time_units = None;
            if let Some(t) = &mut e.trigger {
                t.initial_value.get_or_insert(true);
                t.persistent.get_or_insert(true);
            }
        } else {
            if e.priority.is_some() { return bad("event priorities"); }
            if let Some(t) = &mut e.trigger {
                if t.initial_value == Some(false) || t.persistent == Some(false) {
                    return bad("non-persistent triggers or an initial trigger value of false");
                }
                (t.initial_value, t.persistent) = (None, None);
            }
            if version < 4 {
                if e.use_values_from_trigger_time == Some(false) { return bad("useValuesFromTriggerTime=false"); }
                e.use_values_from_trigger_time = None;
            } else {
                e.use_values_from_trigger_time.get_or_insert(true);
            }
            if version >= 3 { e.time_units = None; }
        }
    }
    if level == 2 {
        if let Some(units) = m.time_units.take() {
            if units != "time" {
                let mut time = if let Some(u) = m.unit_definitions.iter().find(|u| u.id == units) { u.clone() }
                    else { sbml_io::UnitDefinition { units: vec![sbml_io::Unit { kind: units, ..Default::default() }], ..Default::default() } };
                time.id = "time".into();
                m.unit_definitions.retain(|u| u.id != "time");
                m.unit_definitions.push(time);
            }
        }
        (m.substance_units, m.volume_units, m.area_units, m.length_units, m.extent_units) = (None, None, None, None, None);
    }
    if (level, version) != (2, 1) {
        for u in m.unit_definitions.iter_mut().flat_map(|d| &mut d.units) { u.offset = None; }
    }
    (doc.level, doc.version) = (level, version);
    Ok(())
}

pub(crate) fn has_level_three_packages(doc: &Document) -> bool {
    doc.namespaces.iter().any(|(prefix, url)| prefix.is_some() &&
        url.starts_with("http://www.sbml.org/sbml/level3/") && !url.ends_with("/core"))
}

fn remove_empty_math(m: &mut Model) {
    m.function_definitions.retain(|x| x.math.is_some());
    m.initial_assignments.retain(|x| x.math.is_some());
    m.constraints.retain(|x| x.math.is_some());
    m.rules.retain(|x| x.math().is_some());
    for r in &mut m.reactions {
        if r.kinetic_law.as_ref().is_some_and(|x| x.math.is_none()) { r.kinetic_law = None; }
    }
    m.events.retain(|x| x.trigger.as_ref().is_some_and(|t| t.math.is_some()));
    for e in &mut m.events {
        if e.priority.as_ref().is_some_and(|x| x.math.is_none()) { e.priority = None; }
        if e.delay.as_ref().is_some_and(|x| x.math.is_none()) { e.delay = None; }
        e.assignments.retain(|x| x.math.is_some());
    }
}

fn uses_new_math(expr: &mut Math) -> bool {
    let mut found = false;
    visit_math(expr, &mut |node| {
        found |= matches!(node, Math::Apply(Op::Min | Op::Max | Op::Rem | Op::Quotient | Op::Implies, _));
        Ok::<(), std::convert::Infallible>(())
    }).unwrap();
    found
}

fn fold_initial_assignments(m: &mut Model, all: bool) -> Result<(), super::Error> {
    let symbols: HashSet<_> = m.initial_assignments.iter_mut()
        .filter_map(|a| a.math.as_mut().is_some_and(|expr| all || uses_new_math(expr))
            .then(|| a.symbol.clone())).collect();
    if symbols.is_empty() { return Ok(()); }
    let source = M::read(m)?;
    let fd = source.fdefs();
    let mut init = Init { m: &source, fd: &fd, dae: &HashSet::new(), vals: HashMap::new(), busy: HashSet::new() };
    for id in &symbols {
        if init.get(id)?.is_nan() {
            return Err(super::inval(format!("initial assignment {id} has no numeric value")));
        }
    }
    if all { preserve_delay_history(m)?; }
    let get = |id: &str| symbols.contains(id).then(|| init.vals[id]);
    for p in &mut m.parameters {
        if let Some(value) = get(&p.id) { p.value = Some(value); }
    }
    for c in &mut m.compartments {
        if let Some(value) = get(&c.id) { c.size = Some(value); }
    }
    for s in &mut m.species {
        if let Some(value) = get(&s.id) {
            if s.has_only_substance_units.unwrap_or(false) {
                (s.initial_amount, s.initial_concentration) = (Some(value), None);
            } else {
                (s.initial_amount, s.initial_concentration) = (None, Some(value));
            }
        }
    }
    for sr in species_references(m) {
        if let Some(value) = sr.id.as_deref().and_then(get) { sr.stoichiometry = Some(value); }
    }
    m.initial_assignments.retain(|a| !symbols.contains(&a.symbol));
    Ok(())
}

// Level 2 Version 1 stores the evaluated initial attributes. The original
// initial expressions still provide the history before time zero in delay().
fn preserve_delay_history(m: &mut Model) -> Result<(), super::Error> {
    let initial: HashMap<_, _> = m.initial_assignments.iter()
        .filter_map(|a| a.math.as_ref().map(|e| (a.symbol.clone(), e.clone()))).collect();
    for_each_math(m, |expr, locals, _| history_math(expr, &initial, locals))
}

fn history_math(expr: &mut Math, initial: &HashMap<String, Math>, locals: &HashSet<String>) -> Result<(), super::Error> {
    visit_children(expr, &mut |child| history_math(child, initial, locals))?;
    let Math::CallSymbol(symbol, args) = expr else { return Ok(()) };
    if symbol.url != Symbol::DELAY { return Ok(()); }
    let [formula, lag] = args.as_slice() else { return Ok(()) };
    let mut history = formula.clone();
    visit_math(&mut history, &mut |node| {
        if let Math::Ident(id) = node {
            if let Some(value) = initial.get(id).filter(|_| !locals.contains(id)) {
                *node = value.clone();
            }
        }
        Ok(())
    })?;
    if &history != formula {
        let time = Math::Symbol(Symbol { url: Symbol::TIME.into(), name: "time".into() });
        let before = Math::Apply(Op::Lt, vec![time, lag.clone()]);
        let past = Math::CallSymbol(symbol.clone(), vec![history, lag.clone()]);
        *expr = Math::Piecewise { pieces: vec![(past, before)], otherwise: Some(Box::new(expr.clone())) };
    }
    Ok(())
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum MathContext { Function, Initial, Constraint, Trigger, Numeric }

pub(super) fn for_each_math<E>(m: &mut Model, mut visit: impl FnMut(&mut Math, &HashSet<String>, MathContext) -> Result<(), E>) -> Result<(), E> {
    let empty = HashSet::new();
    for f in &mut m.function_definitions {
        if let Some(expr) = &mut f.math {
            let args = lambda(expr).map_or_else(HashSet::new, |(ids, _)| ids.iter().cloned().collect());
            visit(expr, &args, MathContext::Function)?;
        }
    }
    for r in &mut m.rules {
        let expr = match r {
            Rule::Assignment { math, .. } | Rule::Rate { math, .. } | Rule::Algebraic { math, .. } => math,
        };
        if let Some(expr) = expr { visit(expr, &empty, MathContext::Numeric)?; }
    }
    for a in &mut m.initial_assignments {
        if let Some(expr) = &mut a.math { visit(expr, &empty, MathContext::Initial)?; }
    }
    for c in &mut m.constraints {
        if let Some(expr) = &mut c.math { visit(expr, &empty, MathContext::Constraint)?; }
    }
    for r in &mut m.reactions {
        for sr in r.reactants.iter_mut().chain(&mut r.products) {
            if let Some(expr) = &mut sr.stoichiometry_math { visit(expr, &empty, MathContext::Numeric)?; }
        }
        if let Some(k) = &mut r.kinetic_law {
            if let Some(expr) = &mut k.math {
                visit(expr, &k.parameters.iter().map(|p| p.id.clone()).collect(), MathContext::Numeric)?;
            }
        }
    }
    for e in &mut m.events {
        if let Some(expr) = e.trigger.as_mut().and_then(|t| t.math.as_mut()) {
            visit(expr, &empty, MathContext::Trigger)?;
        }
        let expressions = e.delay.iter_mut().filter_map(|d| d.math.as_mut())
            .chain(e.priority.iter_mut().filter_map(|p| p.math.as_mut()))
            .chain(e.assignments.iter_mut().filter_map(|a| a.math.as_mut()));
        for expr in expressions { visit(expr, &empty, MathContext::Numeric)?; }
    }
    Ok(())
}

fn lambda(expr: &Math) -> Option<(&[String], &Math)> {
    match expr {
        Math::Lambda { args, body } => Some((args, body)),
        Math::Semantics { math, .. } => lambda(math),
        _ => None,
    }
}

// An unknown lambda argument or empty branch imposes no constraint.
fn math_type(expr: &Math, functions: &HashMap<String, Math>, bindings: &HashMap<String, Option<bool>>) -> Result<Option<bool>, ()> {
    let infer = |x| math_type(x, functions, bindings);
    let require = |value: Option<bool>, expected| if value.is_some_and(|x| x != expected) { Err(()) } else { Ok(()) };
    Ok(match expr {
        Math::Ident(id) => bindings.get(id).copied().unwrap_or(Some(false)),
        Math::Constant(Constant::True | Constant::False) => Some(true),
        Math::Apply(op, args) => {
            let logical = matches!(op, Op::And | Op::Or | Op::Xor | Op::Not | Op::Implies);
            let relation = matches!(op, Op::Eq | Op::Neq | Op::Gt | Op::Lt | Op::Geq | Op::Leq);
            let equality = matches!(op, Op::Eq | Op::Neq);
            let mut first = None;
            for x in args {
                let kind = infer(x)?;
                if !equality { require(kind, logical)?; }
                if equality && first.is_some() { require(kind, first.unwrap())?; }
                first = first.or(kind);
            }
            Some(logical || relation)
        }
        Math::Call(id, args) => {
            let Some((params, body)) = functions.get(id).and_then(lambda) else { return Ok(None) };
            let arguments = params.iter().zip(args).map(|(p, x)| Ok((p.clone(), infer(x)?)))
                .collect::<Result<HashMap<_, _>, ()>>()?;
            math_type(body, functions, &arguments)?
        }
        Math::Piecewise { pieces, otherwise } => {
            let mut first = None;
            for (v, c) in pieces {
                require(infer(c)?, true)?;
                let kind = infer(v)?;
                if let Some(expected) = first { require(kind, expected)?; }
                first = first.or(kind);
            }
            if let Some(other) = otherwise {
                let kind = infer(other)?;
                if let Some(expected) = first { require(kind, expected)?; }
                first = first.or(kind);
            }
            first
        }
        Math::CallSymbol(_, args) => {
            for x in args { require(infer(x)?, false)?; }
            Some(false)
        }
        Math::Root { degree, arg } | Math::Log { base: degree, arg } => {
            require(infer(arg)?, false)?;
            if let Some(x) = degree { require(infer(x)?, false)?; }
            Some(false)
        }
        Math::Semantics { math, .. } => infer(math)?,
        _ => Some(false),
    })
}

fn species_references(m: &mut Model) -> impl Iterator<Item = &mut sbml_io::SpeciesReference> {
    m.reactions.iter_mut().flat_map(|r| r.reactants.iter_mut().chain(&mut r.products))
}

fn raise_stoichiometry(m: &mut Model) -> Result<(), String> {
    let mut ids = M::read(m).map_err(|e| e.to_string())?.sids();
    let mut rules = vec![];
    let mut next = 0;
    for sr in species_references(m) {
        if let Some(expr) = sr.stoichiometry_math.take() {
            let id = sr.id.get_or_insert_with(|| loop {
                let id = format!("stoichiometry_{next}");
                next += 1;
                if ids.insert(id.clone()) { break id; }
            });
            rules.push(Rule::Assignment { variable: id.clone(), math: Some(expr), meta: Default::default() });
            sr.constant = Some(false);
            sr.stoichiometry = None;
        }
    }
    m.rules.extend(rules);
    Ok(())
}

const STOICHIOMETRY_PARAMETER: &str = "https://tellurium-nitro.org/ns/sbml/stoichiometry";

fn lower_stoichiometry(m: &mut Model) {
    let ids: HashSet<_> = species_references(m).filter_map(|sr| sr.id.clone()).collect();
    let dynamic: HashSet<_> = m.rules.iter().filter_map(|r| match r {
        Rule::Rate { variable, .. } => Some(variable.clone()), _ => None,
    }).collect();
    let mut assignments = HashMap::new();
    m.rules.retain_mut(|r| {
        if let Rule::Assignment { variable, math, .. } = r {
            if ids.contains(variable) {
                if let Some(expr) = math.take() { assignments.insert(variable.clone(), expr); }
                return false;
            }
        }
        true
    });
    let initial: HashSet<_> = m.initial_assignments.iter().map(|a| a.symbol.clone()).collect();
    let mut params = vec![];
    for sr in species_references(m) {
        let Some(id) = sr.id.clone() else { continue };
        let value = sr.stoichiometry.unwrap_or(1.0);
        if dynamic.contains(&id) || initial.contains(&id) {
            // Level 2 represents a changing coefficient through a parameter.
            // Annotate the parameter so Nitro keeps it in the stoichiometry ID space.
            let mut p = sbml_io::Parameter {
                id: id.clone(), value: Some(value), constant: sr.constant, ..Default::default()
            };
            p.meta.annotation = Some(format!(r#"<stoichiometryParameter xmlns="{STOICHIOMETRY_PARAMETER}"/>"#));
            params.push(p);
            sr.stoichiometry_math = Some(Math::Ident(id));
            sr.id = None;
            sr.stoichiometry = None;
        } else if let Some(expr) = assignments.remove(&id) {
            sr.stoichiometry_math = Some(expr);
            sr.stoichiometry = None;
        }
    }
    m.parameters.extend(params);
}

pub(super) fn is_stoichiometry_parameter(p: &sbml_io::Parameter) -> bool {
    p.meta.annotation.as_ref().is_some_and(|a| a.contains(STOICHIOMETRY_PARAMETER))
}

fn visit_math<E>(e: &mut Math, visit: &mut impl FnMut(&mut Math) -> Result<(), E>) -> Result<(), E> {
    visit(e)?;
    visit_children(e, &mut |child| visit_math(child, visit))
}

fn visit_children<E>(e: &mut Math, visit: &mut impl FnMut(&mut Math) -> Result<(), E>) -> Result<(), E> {
    match e {
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            for arg in args { visit(arg)?; }
        }
        Math::Root { degree: base, arg } | Math::Log { base, arg } => {
            if let Some(base) = base { visit(base)?; }
            visit(arg)?;
        }
        Math::Piecewise { pieces, otherwise } => {
            for (value, condition) in pieces {
                visit(value)?;
                visit(condition)?;
            }
            if let Some(value) = otherwise { visit(value)?; }
        }
        Math::Semantics { math: body, .. } | Math::Lambda { body, .. } => visit(body)?,
        _ => {}
    }
    Ok(())
}

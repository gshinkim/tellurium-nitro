//! Structural checks: ids are unique, and references between elements
//! resolve. Not libSBML's full set of consistency rules.

use std::collections::{HashMap, HashSet};

use crate::layout::{AdditionalGlyph, Layout};
use crate::math::Math;
use crate::model::*;

/// A problem `validate` found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Issue {
    pub message: String,
}

impl std::fmt::Display for Issue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

/// Checks that ids are unique and that references resolve.
pub fn validate(doc: &Document) -> Vec<Issue> {
    let mut v = V::default();
    if let Some(m) = &doc.model {
        v.model(m, doc.level);
    }
    v.issues
}

/// Base units, and in Level 2 the built-in units.
const BASE_UNITS: &[&str] = &[
    "ampere",
    "avogadro",
    "becquerel",
    "candela",
    "celsius",
    "coulomb",
    "dimensionless",
    "farad",
    "gram",
    "gray",
    "henry",
    "hertz",
    "item",
    "joule",
    "katal",
    "kelvin",
    "kilogram",
    "liter",
    "litre",
    "lumen",
    "lux",
    "meter",
    "metre",
    "mole",
    "newton",
    "ohm",
    "pascal",
    "radian",
    "second",
    "siemens",
    "sievert",
    "steradian",
    "tesla",
    "volt",
    "watt",
    "weber",
];
const L2_BUILT_IN_UNITS: &[&str] = &["substance", "volume", "area", "length", "time"];

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Function,
    CompartmentType,
    SpeciesType,
    Compartment,
    Species,
    Parameter,
    Reaction,
    SpeciesReference,
    Event,
}

#[derive(Default)]
struct V {
    issues: Vec<Issue>,
}

impl V {
    fn issue(&mut self, message: String) {
        self.issues.push(Issue { message });
    }

    fn model(&mut self, m: &Model, level: u32) {
        // The model-wide id namespace.
        fn add<'m>(v: &mut V, id: &'m str, kind: Kind, ids: &mut HashMap<&'m str, Kind>) {
            if !id.is_empty() && ids.insert(id, kind).is_some() {
                v.issue(format!("the id \"{id}\" is used more than once"));
            }
        }
        let mut ids: HashMap<&str, Kind> = HashMap::new();
        for f in &m.function_definitions {
            add(self, &f.id, Kind::Function, &mut ids);
        }
        for t in &m.compartment_types {
            add(self, &t.id, Kind::CompartmentType, &mut ids);
        }
        for t in &m.species_types {
            add(self, &t.id, Kind::SpeciesType, &mut ids);
        }
        for c in &m.compartments {
            add(self, &c.id, Kind::Compartment, &mut ids);
        }
        for s in &m.species {
            add(self, &s.id, Kind::Species, &mut ids);
        }
        for p in &m.parameters {
            add(self, &p.id, Kind::Parameter, &mut ids);
        }
        for r in &m.reactions {
            add(self, &r.id, Kind::Reaction, &mut ids);
            for sr in r.reactants.iter().chain(&r.products) {
                if let Some(id) = &sr.id {
                    add(self, id, Kind::SpeciesReference, &mut ids);
                }
            }
            for mr in &r.modifiers {
                if let Some(id) = &mr.id {
                    add(self, id, Kind::SpeciesReference, &mut ids);
                }
            }
        }
        for e in &m.events {
            if let Some(id) = &e.id {
                add(self, id, Kind::Event, &mut ids);
            }
        }

        let mut unit_ids = HashSet::new();
        for u in &m.unit_definitions {
            if !unit_ids.insert(u.id.as_str()) {
                self.issue(format!(
                    "the unit definition \"{}\" is defined more than once",
                    u.id
                ));
            }
            for unit in &u.units {
                if !BASE_UNITS.contains(&unit.kind.as_str()) {
                    self.issue(format!(
                        "unit definition \"{}\": \"{}\" isn't a base unit",
                        u.id, unit.kind
                    ));
                }
            }
        }
        let units_ok = |u: &str| {
            unit_ids.contains(u)
                || BASE_UNITS.contains(&u)
                || (level == 2 && L2_BUILT_IN_UNITS.contains(&u))
        };
        let check_units = |v: &mut V, what: &str, u: &Option<String>| {
            if let Some(u) = u
                && !units_ok(u)
            {
                v.issue(format!("{what}: the units \"{u}\" aren't defined"));
            }
        };
        for (what, u) in [
            ("model substanceUnits", &m.substance_units),
            ("model timeUnits", &m.time_units),
            ("model volumeUnits", &m.volume_units),
            ("model areaUnits", &m.area_units),
            ("model lengthUnits", &m.length_units),
            ("model extentUnits", &m.extent_units),
        ] {
            check_units(self, what, u);
        }

        let refs = |v: &mut V, what: &str, id: &str, kinds: &[Kind]| {
            if !ids.get(id).is_some_and(|k| kinds.contains(k)) {
                v.issue(format!("{what} refers to \"{id}\", which isn't defined"));
            }
        };
        let opt_refs = |v: &mut V, what: &str, id: &Option<String>, kinds: &[Kind]| {
            if let Some(id) = id {
                refs(v, what, id, kinds);
            }
        };
        let variables = [
            Kind::Compartment,
            Kind::Species,
            Kind::Parameter,
            Kind::SpeciesReference,
        ];

        if let Some(cf) = &m.conversion_factor {
            refs(self, "the model's conversionFactor", cf, &[Kind::Parameter]);
        }
        for c in &m.compartments {
            let what = format!("compartment \"{}\"", c.id);
            opt_refs(self, &what, &c.outside, &[Kind::Compartment]);
            opt_refs(self, &what, &c.compartment_type, &[Kind::CompartmentType]);
            check_units(self, &what, &c.units);
        }
        for s in &m.species {
            let what = format!("species \"{}\"", s.id);
            refs(self, &what, &s.compartment, &[Kind::Compartment]);
            opt_refs(self, &what, &s.species_type, &[Kind::SpeciesType]);
            opt_refs(self, &what, &s.conversion_factor, &[Kind::Parameter]);
            check_units(self, &what, &s.substance_units);
            check_units(self, &what, &s.spatial_size_units);
        }
        for p in &m.parameters {
            check_units(self, &format!("parameter \"{}\"", p.id), &p.units);
        }

        // Ids that may appear in math, besides function names.
        let math_ok = |id: &str| ids.get(id).is_some_and(|k| *k != Kind::Function);
        let check_math = |v: &mut V, what: &str, math: &Option<Math>, local: &[&str]| {
            if let Some(math) = math {
                v.math(
                    what,
                    math,
                    &math_ok,
                    &|f| ids.get(f) == Some(&Kind::Function),
                    local,
                );
            }
        };

        for f in &m.function_definitions {
            let what = format!("function definition \"{}\"", f.id);
            match &f.math {
                Some(Math::Lambda { args, body }) => {
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    // A function body sees only its arguments (and other functions).
                    self.math(
                        &what,
                        body,
                        &|_| false,
                        &|f| ids.get(f) == Some(&Kind::Function),
                        &args,
                    );
                }
                Some(Math::Semantics { .. }) | None => {}
                Some(_) => self.issue(format!("{what}: its math isn't a lambda")),
            }
        }
        for ia in &m.initial_assignments {
            let what = format!("initial assignment to \"{}\"", ia.symbol);
            refs(self, &what, &ia.symbol, &variables);
            check_math(self, &what, &ia.math, &[]);
        }
        for r in &m.rules {
            let what = match r.variable() {
                Some(var) => {
                    let what = format!("rule for \"{var}\"");
                    refs(self, &what, var, &variables);
                    what
                }
                None => "algebraic rule".to_string(),
            };
            check_math(self, &what, &r.math().cloned(), &[]);
        }
        for (i, c) in m.constraints.iter().enumerate() {
            check_math(self, &format!("constraint {}", i + 1), &c.math, &[]);
        }
        for r in &m.reactions {
            let what = format!("reaction \"{}\"", r.id);
            opt_refs(self, &what, &r.compartment, &[Kind::Compartment]);
            for sr in r.reactants.iter().chain(&r.products) {
                refs(self, &what, &sr.species, &[Kind::Species]);
                check_math(self, &what, &sr.stoichiometry_math, &[]);
            }
            for mr in &r.modifiers {
                refs(self, &what, &mr.species, &[Kind::Species]);
            }
            if let Some(kl) = &r.kinetic_law {
                let mut local = HashSet::new();
                for p in &kl.parameters {
                    if !local.insert(p.id.as_str()) {
                        self.issue(format!(
                            "{what}: the local parameter \"{}\" is defined more than once",
                            p.id
                        ));
                    }
                    check_units(
                        self,
                        &format!("{what}, local parameter \"{}\"", p.id),
                        &p.units,
                    );
                }
                let local: Vec<&str> = local.into_iter().collect();
                check_math(self, &what, &kl.math, &local);
            }
        }
        for (i, e) in m.events.iter().enumerate() {
            let what = match &e.id {
                Some(id) => format!("event \"{id}\""),
                None => format!("event {}", i + 1),
            };
            if let Some(t) = &e.trigger {
                check_math(self, &what, &t.math, &[]);
            }
            if let Some(d) = &e.delay {
                check_math(self, &what, &d.math, &[]);
            }
            if let Some(p) = &e.priority {
                check_math(self, &what, &p.math, &[]);
            }
            for ea in &e.assignments {
                refs(self, &what, &ea.variable, &variables);
                check_math(self, &what, &ea.math, &[]);
            }
        }

        let mut layout_ids = HashSet::new();
        for l in &m.layouts {
            if !layout_ids.insert(l.id.as_str()) {
                self.issue(format!("the layout \"{}\" is defined more than once", l.id));
            }
            self.layout(l, &ids);
        }
    }

    fn math(
        &mut self,
        what: &str,
        math: &Math,
        id_ok: &dyn Fn(&str) -> bool,
        function_ok: &dyn Fn(&str) -> bool,
        local: &[&str],
    ) {
        let mut bad = Vec::new();
        let walk = |m: &Math, bad: &mut Vec<String>| {
            collect(m, &mut |name, is_call| {
                let ok = if is_call {
                    function_ok(name)
                } else {
                    local.contains(&name) || id_ok(name)
                };
                if !ok && !bad.iter().any(|b| b == name) {
                    bad.push(name.to_string());
                }
            })
        };
        walk(math, &mut bad);
        for name in bad {
            self.issue(format!(
                "{what}: the math refers to \"{name}\", which isn't defined"
            ));
        }
    }

    fn layout(&mut self, l: &Layout, ids: &HashMap<&str, Kind>) {
        let what = format!("layout \"{}\"", l.id);
        // Every glyph's id, with what kind of glyph it is.
        fn additional<'l>(g: &'l AdditionalGlyph, out: &mut Vec<(&'l str, &'static str)>) {
            match g {
                AdditionalGlyph::General(g) => {
                    out.push((&g.id, "general"));
                    out.extend(
                        g.reference_glyphs
                            .iter()
                            .map(|r| (r.id.as_str(), "reference")),
                    );
                    for s in &g.sub_glyphs {
                        additional(s, out);
                    }
                }
                AdditionalGlyph::Species(g) => out.push((&g.id, "species")),
                AdditionalGlyph::Reaction(g) => out.push((&g.id, "reaction")),
                g => out.push((g.id(), "other")),
            }
        }
        let mut all: Vec<(&str, &'static str)> = Vec::new();
        all.extend(
            l.compartment_glyphs
                .iter()
                .map(|g| (g.id.as_str(), "compartment")),
        );
        all.extend(l.species_glyphs.iter().map(|g| (g.id.as_str(), "species")));
        for g in &l.reaction_glyphs {
            all.push((&g.id, "reaction"));
            all.extend(
                g.species_reference_glyphs
                    .iter()
                    .map(|s| (s.id.as_str(), "species reference")),
            );
        }
        all.extend(l.text_glyphs.iter().map(|g| (g.id.as_str(), "text")));
        for g in &l.additional {
            additional(g, &mut all);
        }
        let mut glyphs: HashMap<&str, &str> = HashMap::new();
        for (id, kind) in all {
            if glyphs.insert(id, kind).is_some() {
                self.issue(format!(
                    "{what}: the glyph id \"{id}\" is used more than once"
                ));
            }
        }

        let model_ref = |v: &mut V, id: &Option<String>, kinds: &[Kind]| {
            if let Some(id) = id
                && !ids.get(id.as_str()).is_some_and(|k| kinds.contains(k))
            {
                v.issue(format!(
                    "{what}: a glyph refers to \"{id}\", which isn't in the model"
                ));
            }
        };
        for g in &l.compartment_glyphs {
            model_ref(self, &g.compartment, &[Kind::Compartment]);
        }
        for g in &l.species_glyphs {
            model_ref(self, &g.species, &[Kind::Species]);
        }
        for g in &l.reaction_glyphs {
            model_ref(self, &g.reaction, &[Kind::Reaction]);
            for s in &g.species_reference_glyphs {
                model_ref(self, &s.species_reference, &[Kind::SpeciesReference]);
                if glyphs.get(s.species_glyph.as_str()) != Some(&"species") {
                    self.issue(format!(
                        "{what}: species reference glyph \"{}\" refers to \"{}\", which isn't a species glyph",
                        s.id, s.species_glyph
                    ));
                }
            }
        }
        for g in &l.text_glyphs {
            if let Some(go) = &g.graphical_object
                && !glyphs.contains_key(go.as_str())
            {
                self.issue(format!(
                    "{what}: text glyph \"{}\" refers to \"{go}\", which isn't a glyph",
                    g.id
                ));
            }
            if let Some(o) = &g.origin_of_text
                && !ids.contains_key(o.as_str())
            {
                self.issue(format!(
                    "{what}: text glyph \"{}\" refers to \"{o}\", which isn't in the model",
                    g.id
                ));
            }
        }
    }
}

/// Calls `f(name, is_call)` for every identifier in `m`, leaving out lambda
/// arguments inside their lambda.
fn collect(m: &Math, f: &mut dyn FnMut(&str, bool)) {
    match m {
        Math::Ident(id) => f(id, false),
        Math::Call(name, args) => {
            f(name, true);
            args.iter().for_each(|a| collect(a, f));
        }
        Math::Lambda { args, body } => collect(body, &mut |name, call| {
            if call || !args.iter().any(|a| a == name) {
                f(name, call)
            }
        }),
        Math::Number { .. } | Math::Symbol(_) | Math::Constant(_) | Math::Raw(_) => {}
        Math::Apply(_, args) | Math::CallSymbol(_, args) => args.iter().for_each(|a| collect(a, f)),
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            if let Some(q) = q {
                collect(q, f);
            }
            collect(arg, f);
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, c) in pieces {
                collect(v, f);
                collect(c, f);
            }
            if let Some(o) = otherwise {
                collect(o, f);
            }
        }
        Math::Semantics { math, .. } => collect(math, f),
    }
}

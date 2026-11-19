//! Reading SBML: XML is parsed with `roxmltree`, then walked into the data
//! model. Anything not understood is kept as raw XML, cut from the source.

mod layout;
mod math;
mod render;

use roxmltree::{Node, ParsingOptions};

use crate::model::*;
use crate::ns;

/// A problem found while reading, with its place in the source (1-based).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub line: u32,
    pub column: u32,
    pub message: String,
}

impl std::fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}: {}", self.line, self.column, self.message)
    }
}

/// A document read, with anything worth knowing that didn't stop the read.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub document: Document,
    pub warnings: Vec<Diagnostic>,
}

/// Reads an SBML document. Fails on XML that isn't well formed, on levels
/// and versions other than Level 2 Versions 1–5 and Level 3 Versions 1–2,
/// and on missing or malformed attributes `sbml-io` needs.
pub fn read(xml: &str) -> Result<Read, Vec<Diagnostic>> {
    let opts = ParsingOptions {
        allow_dtd: true,
        ..ParsingOptions::default()
    };
    let doc = roxmltree::Document::parse_with_options(xml, opts).map_err(|e| {
        let p = e.pos();
        vec![Diagnostic {
            line: p.row,
            column: p.col,
            message: e.to_string(),
        }]
    })?;
    let root = doc.root_element();
    let mut r = Reader {
        src: xml,
        doc: &doc,
        level: 0,
        version: 0,
        core: "",
        root_ns: root
            .namespaces()
            .filter(|n| n.name() != Some("xml"))
            .map(|n| (n.name().map(str::to_string), n.uri().to_string()))
            .collect(),
        errors: Vec::new(),
        warnings: Vec::new(),
    };
    if root.tag_name().name() != "sbml" {
        r.error(root, "the root element isn't <sbml>");
        return Err(r.errors);
    }
    let mut e = E::new(root);
    let level = r.uint(&mut e, "level");
    let version = r.uint(&mut e, "version");
    let (Some(level), Some(version)) = (level, version) else {
        if r.errors.is_empty() {
            r.error(root, "<sbml> needs both level and version");
        }
        return Err(r.errors);
    };
    let Some(core) = ns::core(level, version) else {
        r.error(
            root,
            &format!("SBML Level {level} Version {version} isn't supported"),
        );
        return Err(r.errors);
    };
    (r.level, r.version, r.core) = (level, version, core);
    if root.tag_name().namespace() != Some(core) {
        r.error(
            root,
            &format!(
                "<sbml> is in namespace {:?}, but Level {level} Version {version} uses {core}",
                root.tag_name().namespace().unwrap_or("")
            ),
        );
        return Err(r.errors);
    }

    let mut model = None;
    let meta = r.meta(e, None, |r, c| match r.core_name(c) {
        Some("model") => {
            if model.is_some() {
                r.warn(c, "a second <model> is kept as raw XML");
                return false;
            }
            model = Some(r.model(c));
            true
        }
        _ => false,
    });
    let document = Document {
        level,
        version,
        namespaces: r.root_ns.clone(),
        model,
        meta,
    };
    if r.errors.is_empty() {
        Ok(Read {
            document,
            warnings: r.warnings,
        })
    } else {
        Err(r.errors)
    }
}

/// An element whose attributes are taken one by one, so that what's left can
/// be kept as `Extra`.
struct E<'a, 'i> {
    n: Node<'a, 'i>,
    used: Vec<bool>,
}

impl<'a, 'i> E<'a, 'i> {
    fn new(n: Node<'a, 'i>) -> Self {
        E {
            n,
            used: vec![false; n.attributes().len()],
        }
    }

    /// Takes an attribute in the given namespace (`None`: no prefix).
    fn take_in(&mut self, namespace: Option<&str>, name: &str) -> Option<&'a str> {
        let i = self
            .n
            .attributes()
            .position(|a| a.name() == name && a.namespace() == namespace)?;
        self.used[i] = true;
        Some(self.n.attributes().nth(i).unwrap().value())
    }

    fn take(&mut self, name: &str) -> Option<&'a str> {
        self.take_in(None, name)
    }

    /// Takes an attribute of a package element, which may be written with the
    /// package's prefix or without one.
    fn take_pkg(&mut self, pkg: Option<&str>, name: &str) -> Option<&'a str> {
        pkg.and_then(|p| self.take_in(Some(p), name))
            .or_else(|| self.take(name))
    }

    fn take_string(&mut self, name: &str) -> Option<String> {
        self.take(name).map(str::to_string)
    }
}

struct Reader<'a, 'i> {
    src: &'i str,
    doc: &'a roxmltree::Document<'i>,
    level: u32,
    version: u32,
    core: &'static str,
    root_ns: Vec<(Option<String>, String)>,
    errors: Vec<Diagnostic>,
    warnings: Vec<Diagnostic>,
}

impl<'a, 'i: 'a> Reader<'a, 'i> {
    fn diag(&self, n: Node, message: &str) -> Diagnostic {
        let p = self.doc.text_pos_at(n.range().start);
        Diagnostic {
            line: p.row,
            column: p.col,
            message: message.to_string(),
        }
    }

    fn error(&mut self, n: Node, message: &str) {
        let d = self.diag(n, message);
        self.errors.push(d);
    }

    /// `base`, or `base_2`, `base_3`, …: the first that isn't an attribute
    /// value anywhere in the document, so it clashes with no id.
    fn unused_id(&self, base: &str) -> String {
        let used = |id: &str| {
            self.doc
                .descendants()
                .any(|n| n.attributes().any(|a| a.value() == id))
        };
        std::iter::once(base.to_owned())
            .chain((2..).map(|n| format!("{base}_{n}")))
            .find(|id| !used(id))
            .expect("some suffix is free")
    }

    fn warn(&mut self, n: Node, message: &str) {
        let d = self.diag(n, message);
        self.warnings.push(d);
    }

    /// The local name of a core element.
    fn core_name<'n>(&self, n: Node<'n, 'i>) -> Option<&'n str> {
        (n.is_element() && n.tag_name().namespace() == Some(self.core)).then(|| n.tag_name().name())
    }

    /// The local name of an element in namespace `ns`.
    fn name_in<'n>(&self, n: Node<'n, 'i>, ns: &str) -> Option<&'n str> {
        (n.is_element() && n.tag_name().namespace() == Some(ns)).then(|| n.tag_name().name())
    }

    /// An element as raw XML. Namespaces it inherits from elements other than
    /// `<sbml>` are declared on it, so that it stands on its own.
    fn raw(&self, n: Node) -> String {
        let text = &self.src[n.range()];
        if !n.is_element() {
            return text.to_string();
        }
        let parent_ns: Vec<_> = n
            .parent_element()
            .map(|p| p.namespaces().map(|x| (x.name(), x.uri())).collect())
            .unwrap_or_default();
        let mut add = String::new();
        for x in n.namespaces() {
            let (prefix, uri) = (x.name(), x.uri());
            if prefix == Some("xml") {
                continue;
            }
            let inherited =
                parent_ns.contains(&(prefix, uri)) && !declares(text, prefix) && uses(text, prefix);
            let on_root = self
                .root_ns
                .iter()
                .any(|(p, u)| p.as_deref() == prefix && u == uri);
            if inherited && !on_root {
                match prefix {
                    Some(p) => add.push_str(&format!(" xmlns:{p}=\"{}\"", escape_attr(uri))),
                    None => add.push_str(&format!(" xmlns=\"{}\"", escape_attr(uri))),
                }
            }
        }
        if add.is_empty() {
            return text.to_string();
        }
        // After `<name`.
        let name_end = text[1..]
            .find(|c: char| c.is_whitespace() || c == '>' || c == '/')
            .map_or(text.len(), |i| i + 1);
        format!("{}{}{}", &text[..name_end], add, &text[name_end..])
    }

    /// The content of an element as raw XML.
    fn inner(&self, n: Node) -> String {
        n.children().map(|c| self.raw(c)).collect()
    }

    /// Reads `metaid`, `sboTerm`, `<notes>` and `<annotation>`, gives every
    /// other child element to `child` (which returns whether it took it), and
    /// keeps whatever is left as `Extra`. `pkg` is the namespace of a package
    /// element, whose attributes may carry the package's prefix.
    fn meta(
        &mut self,
        mut e: E<'a, 'i>,
        pkg: Option<&str>,
        mut child: impl FnMut(&mut Self, Node<'a, 'i>) -> bool,
    ) -> Meta {
        let mut meta = Meta {
            metaid: e.take_pkg(pkg, "metaid").map(str::to_string),
            sbo_term: e.take_pkg(pkg, "sboTerm").map(str::to_string),
            ..Meta::default()
        };
        for c in e.n.children().filter(|c| c.is_element()) {
            let own = self.core_name(c).is_some()
                || pkg.is_some_and(|p| c.tag_name().namespace() == Some(p));
            let name = c.tag_name().name();
            if own && name == "notes" && meta.notes.is_none() {
                meta.notes = Some(self.inner(c));
            } else if child(self, c) {
                // Taken.
            } else if own && name == "annotation" && meta.annotation.is_none() {
                meta.annotation = Some(self.inner(c));
            } else {
                meta.extra.elements.push(self.raw(c));
            }
        }
        meta.extra.attributes =
            e.n.attributes()
                .zip(&e.used)
                .filter(|(_, used)| !**used)
                .map(|(a, _)| Attribute {
                    name: self.src[a.range_qname()].to_string(),
                    namespace: a.namespace().map(str::to_string),
                    value: a.value().to_string(),
                })
                .collect();
        meta
    }

    /// The items of a `listOf…` element, read by `item`. Anything else in the
    /// list is reported and dropped.
    fn list<T>(
        &mut self,
        n: Node<'a, 'i>,
        ns: &str,
        names: &[&str],
        mut item: impl FnMut(&mut Self, Node<'a, 'i>) -> T,
    ) -> Vec<T> {
        let mut out = Vec::new();
        let list_name = n.tag_name().name();
        for c in n.children().filter(|c| c.is_element()) {
            if self
                .name_in(c, ns)
                .is_some_and(|name| names.contains(&name))
            {
                out.push(item(self, c));
            } else {
                let msg = format!("<{}> in <{list_name}> is dropped", c.tag_name().name());
                self.warn(c, &msg);
            }
        }
        if n.attributes().len() > 0 {
            let msg = format!("the attributes of <{list_name}> are dropped");
            self.warn(n, &msg);
        }
        out
    }

    /// A required attribute; reports it if it's missing.
    fn req(&mut self, e: &mut E<'a, 'i>, name: &str) -> String {
        self.req_pkg(e, None, name)
    }

    fn req_pkg(&mut self, e: &mut E<'a, 'i>, pkg: Option<&str>, name: &str) -> String {
        match e.take_pkg(pkg, name) {
            Some(v) => v.to_string(),
            None => {
                let msg = format!("<{}> needs a {name}", e.n.tag_name().name());
                self.error(e.n, &msg);
                String::new()
            }
        }
    }

    fn parsed<T>(
        &mut self,
        e: &mut E<'a, 'i>,
        pkg: Option<&str>,
        name: &str,
        what: &str,
        parse: impl Fn(&str) -> Option<T>,
    ) -> Option<T> {
        let v = e.take_pkg(pkg, name)?;
        let parsed = parse(v.trim());
        if parsed.is_none() {
            let msg = format!("{name}=\"{v}\" on <{}> isn't {what}", e.n.tag_name().name());
            self.error(e.n, &msg);
        }
        parsed
    }

    fn double(&mut self, e: &mut E<'a, 'i>, name: &str) -> Option<f64> {
        self.parsed(e, None, name, "a number", parse_double)
    }

    fn double_pkg(&mut self, e: &mut E<'a, 'i>, pkg: Option<&str>, name: &str) -> Option<f64> {
        self.parsed(e, pkg, name, "a number", parse_double)
    }

    fn boolean(&mut self, e: &mut E<'a, 'i>, name: &str) -> Option<bool> {
        self.parsed(e, None, name, "true or false", |v| match v {
            "true" | "1" => Some(true),
            "false" | "0" => Some(false),
            _ => None,
        })
    }

    fn int(&mut self, e: &mut E<'a, 'i>, name: &str) -> Option<i32> {
        self.parsed(e, None, name, "an integer", |v| v.parse().ok())
    }

    fn uint(&mut self, e: &mut E<'a, 'i>, name: &str) -> Option<u32> {
        self.parsed(e, None, name, "a whole number", |v| v.parse().ok())
    }

    /// The `<math>` child of an element, if any.
    fn math_child(&mut self, c: Node<'a, 'i>, slot: &mut Option<crate::Math>) -> bool {
        if self.name_in(c, ns::MATHML) == Some("math") {
            *slot = self.math(c);
            true
        } else {
            false
        }
    }

    fn model(&mut self, n: Node<'a, 'i>) -> Model {
        let mut e = E::new(n);
        let mut m = Model {
            id: e.take_string("id"),
            name: e.take_string("name"),
            substance_units: e.take_string("substanceUnits"),
            time_units: e.take_string("timeUnits"),
            volume_units: e.take_string("volumeUnits"),
            area_units: e.take_string("areaUnits"),
            length_units: e.take_string("lengthUnits"),
            extent_units: e.take_string("extentUnits"),
            conversion_factor: e.take_string("conversionFactor"),
            ..Model::default()
        };
        let level = self.level;
        let core = self.core;
        let mut l2_annotation = None;
        let mut layouts = Vec::new();
        let mut layouts_meta = Meta::default();
        let mut global_render = None;
        let meta = self.meta(e, None, |r, c| {
            if level == 2 && r.core_name(c) == Some("annotation") {
                l2_annotation = Some(r.l2_model_annotation(
                    c,
                    &mut layouts,
                    &mut layouts_meta,
                    &mut global_render,
                ));
                return true;
            }
            if level == 3 && r.name_in(c, ns::LAYOUT) == Some("listOfLayouts") {
                let (ls, meta, global) = r.layout_list(c, ns::LAYOUT);
                layouts.extend(ls);
                (layouts_meta, global_render) = (meta, global);
                return true;
            }
            let Some(name) = r.core_name(c) else {
                return false;
            };
            match name {
                "listOfFunctionDefinitions" => {
                    m.function_definitions =
                        r.list(c, core, &["functionDefinition"], Self::function_definition)
                }
                "listOfUnitDefinitions" => {
                    m.unit_definitions = r.list(c, core, &["unitDefinition"], Self::unit_definition)
                }
                "listOfCompartmentTypes" => {
                    m.compartment_types = r.list(c, core, &["compartmentType"], |r, c| {
                        let mut e = E::new(c);
                        CompartmentType {
                            id: r.req(&mut e, "id"),
                            name: e.take_string("name"),
                            meta: r.meta(e, None, |_, _| false),
                        }
                    })
                }
                "listOfSpeciesTypes" => {
                    m.species_types = r.list(c, core, &["speciesType"], |r, c| {
                        let mut e = E::new(c);
                        SpeciesType {
                            id: r.req(&mut e, "id"),
                            name: e.take_string("name"),
                            meta: r.meta(e, None, |_, _| false),
                        }
                    })
                }
                "listOfCompartments" => {
                    m.compartments = r.list(c, core, &["compartment"], Self::compartment)
                }
                "listOfSpecies" => m.species = r.list(c, core, &["species"], Self::species),
                "listOfParameters" => {
                    m.parameters = r.list(c, core, &["parameter"], Self::parameter)
                }
                "listOfInitialAssignments" => {
                    m.initial_assignments = r.list(c, core, &["initialAssignment"], |r, c| {
                        let mut e = E::new(c);
                        let symbol = r.req(&mut e, "symbol");
                        let mut math = None;
                        let meta = r.meta(e, None, |r, c| r.math_child(c, &mut math));
                        InitialAssignment { symbol, math, meta }
                    })
                }
                "listOfRules" => {
                    m.rules = r.list(
                        c,
                        core,
                        &["algebraicRule", "assignmentRule", "rateRule"],
                        Self::rule,
                    )
                }
                "listOfConstraints" => {
                    m.constraints = r.list(c, core, &["constraint"], Self::constraint)
                }
                "listOfReactions" => m.reactions = r.list(c, core, &["reaction"], Self::reaction),
                "listOfEvents" => m.events = r.list(c, core, &["event"], Self::event),
                _ => return false,
            }
            true
        });
        m.meta = meta;
        if l2_annotation.is_some() {
            m.meta.annotation = l2_annotation.flatten();
        }
        m.layouts = layouts;
        m.layouts_meta = layouts_meta;
        m.global_render = global_render;
        m
    }

    /// A Level 2 model annotation: layouts are taken out of it, the rest is
    /// kept. Returns `None` if nothing is left.
    fn l2_model_annotation(
        &mut self,
        n: Node<'a, 'i>,
        layouts: &mut Vec<crate::Layout>,
        layouts_meta: &mut Meta,
        global_render: &mut Option<crate::render::RenderList>,
    ) -> Option<String> {
        let mut rest = String::new();
        let mut took = false;
        let is_layouts = |c: Option<Node>| {
            c.is_some_and(|c| {
                c.is_element()
                    && c.tag_name().namespace() == Some(ns::LAYOUT_L2)
                    && c.tag_name().name() == "listOfLayouts"
            })
        };
        for c in n.children() {
            // The whitespace before the layouts goes with them.
            if c.is_text()
                && c.text().is_some_and(|t| t.trim().is_empty())
                && is_layouts(c.next_sibling())
            {
                continue;
            }
            if is_layouts(Some(c)) {
                let (ls, meta, global) = self.layout_list(c, ns::LAYOUT_L2);
                layouts.extend(ls);
                (*layouts_meta, *global_render) = (meta, global);
                took = true;
            } else {
                rest.push_str(&self.raw(c));
            }
        }
        if took && rest.trim().is_empty() {
            None
        } else {
            Some(rest)
        }
    }

    fn function_definition(&mut self, n: Node<'a, 'i>) -> FunctionDefinition {
        let mut e = E::new(n);
        let id = self.req(&mut e, "id");
        let name = e.take_string("name");
        let mut math = None;
        let meta = self.meta(e, None, |r, c| r.math_child(c, &mut math));
        FunctionDefinition {
            id,
            name,
            math,
            meta,
        }
    }

    fn unit_definition(&mut self, n: Node<'a, 'i>) -> UnitDefinition {
        let mut e = E::new(n);
        let id = self.req(&mut e, "id");
        let name = e.take_string("name");
        let core = self.core;
        let mut units = Vec::new();
        let meta = self.meta(e, None, |r, c| {
            if r.core_name(c) == Some("listOfUnits") {
                units = r.list(c, core, &["unit"], |r, c| {
                    let mut e = E::new(c);
                    Unit {
                        kind: r.req(&mut e, "kind"),
                        exponent: r.double(&mut e, "exponent"),
                        scale: r.int(&mut e, "scale"),
                        multiplier: r.double(&mut e, "multiplier"),
                        offset: r.double(&mut e, "offset"),
                        meta: r.meta(e, None, |_, _| false),
                    }
                });
                true
            } else {
                false
            }
        });
        UnitDefinition {
            id,
            name,
            units,
            meta,
        }
    }

    fn compartment(&mut self, n: Node<'a, 'i>) -> Compartment {
        let mut e = E::new(n);
        Compartment {
            id: self.req(&mut e, "id"),
            name: e.take_string("name"),
            spatial_dimensions: self.double(&mut e, "spatialDimensions"),
            size: self.double(&mut e, "size"),
            units: e.take_string("units"),
            outside: e.take_string("outside"),
            constant: self.boolean(&mut e, "constant"),
            compartment_type: e.take_string("compartmentType"),
            meta: self.meta(e, None, |_, _| false),
        }
    }

    fn species(&mut self, n: Node<'a, 'i>) -> Species {
        let mut e = E::new(n);
        Species {
            id: self.req(&mut e, "id"),
            name: e.take_string("name"),
            compartment: self.req(&mut e, "compartment"),
            initial_amount: self.double(&mut e, "initialAmount"),
            initial_concentration: self.double(&mut e, "initialConcentration"),
            substance_units: e.take_string("substanceUnits"),
            spatial_size_units: e.take_string("spatialSizeUnits"),
            has_only_substance_units: self.boolean(&mut e, "hasOnlySubstanceUnits"),
            boundary_condition: self.boolean(&mut e, "boundaryCondition"),
            charge: self.int(&mut e, "charge"),
            constant: self.boolean(&mut e, "constant"),
            conversion_factor: e.take_string("conversionFactor"),
            species_type: e.take_string("speciesType"),
            meta: self.meta(e, None, |_, _| false),
        }
    }

    fn parameter(&mut self, n: Node<'a, 'i>) -> Parameter {
        let mut e = E::new(n);
        Parameter {
            id: self.req(&mut e, "id"),
            name: e.take_string("name"),
            value: self.double(&mut e, "value"),
            units: e.take_string("units"),
            constant: self.boolean(&mut e, "constant"),
            meta: self.meta(e, None, |_, _| false),
        }
    }

    fn rule(&mut self, n: Node<'a, 'i>) -> Rule {
        let mut e = E::new(n);
        let kind = n.tag_name().name();
        let variable = (kind != "algebraicRule").then(|| self.req(&mut e, "variable"));
        let mut math = None;
        let meta = self.meta(e, None, |r, c| r.math_child(c, &mut math));
        match (kind, variable) {
            ("assignmentRule", Some(variable)) => Rule::Assignment {
                variable,
                math,
                meta,
            },
            ("rateRule", Some(variable)) => Rule::Rate {
                variable,
                math,
                meta,
            },
            _ => Rule::Algebraic { math, meta },
        }
    }

    fn constraint(&mut self, n: Node<'a, 'i>) -> Constraint {
        let e = E::new(n);
        let mut math = None;
        let mut message = None;
        let meta = self.meta(e, None, |r, c| {
            if r.core_name(c) == Some("message") {
                message = Some(r.inner(c));
                true
            } else {
                r.math_child(c, &mut math)
            }
        });
        Constraint {
            math,
            message,
            meta,
        }
    }

    fn reaction(&mut self, n: Node<'a, 'i>) -> Reaction {
        let mut e = E::new(n);
        let mut rx = Reaction {
            id: self.req(&mut e, "id"),
            name: e.take_string("name"),
            reversible: self.boolean(&mut e, "reversible"),
            fast: self.boolean(&mut e, "fast"),
            compartment: e.take_string("compartment"),
            ..Reaction::default()
        };
        let core = self.core;
        rx.meta = self.meta(e, None, |r, c| {
            match r.core_name(c) {
                Some("listOfReactants") => {
                    rx.reactants = r.list(c, core, &["speciesReference"], Self::species_reference)
                }
                Some("listOfProducts") => {
                    rx.products = r.list(c, core, &["speciesReference"], Self::species_reference)
                }
                Some("listOfModifiers") => {
                    rx.modifiers = r.list(c, core, &["modifierSpeciesReference"], |r, c| {
                        let mut e = E::new(c);
                        ModifierSpeciesReference {
                            id: e.take_string("id"),
                            name: e.take_string("name"),
                            species: r.req(&mut e, "species"),
                            meta: r.meta(e, None, |_, _| false),
                        }
                    })
                }
                Some("kineticLaw") => rx.kinetic_law = Some(r.kinetic_law(c)),
                _ => return false,
            }
            true
        });
        rx
    }

    fn species_reference(&mut self, n: Node<'a, 'i>) -> SpeciesReference {
        let mut e = E::new(n);
        let mut sr = SpeciesReference {
            id: e.take_string("id"),
            name: e.take_string("name"),
            species: self.req(&mut e, "species"),
            stoichiometry: self.double(&mut e, "stoichiometry"),
            constant: self.boolean(&mut e, "constant"),
            ..SpeciesReference::default()
        };
        sr.meta = self.meta(e, None, |r, c| {
            if r.core_name(c) != Some("stoichiometryMath") {
                return false;
            }
            let math = c
                .children()
                .find(|m| r.name_in(*m, ns::MATHML) == Some("math"));
            if let Some(m) = math {
                sr.stoichiometry_math = r.math(m);
            }
            true
        });
        sr
    }

    fn kinetic_law(&mut self, n: Node<'a, 'i>) -> KineticLaw {
        let mut e = E::new(n);
        let mut kl = KineticLaw {
            time_units: e.take_string("timeUnits"),
            substance_units: e.take_string("substanceUnits"),
            ..KineticLaw::default()
        };
        let core = self.core;
        kl.meta = self.meta(e, None, |r, c| {
            let item = match r.core_name(c) {
                Some("listOfLocalParameters") => "localParameter",
                Some("listOfParameters") => "parameter",
                _ => return r.math_child(c, &mut kl.math),
            };
            kl.parameters = r.list(c, core, &[item], |r, c| {
                let mut e = E::new(c);
                LocalParameter {
                    id: r.req(&mut e, "id"),
                    name: e.take_string("name"),
                    value: r.double(&mut e, "value"),
                    units: e.take_string("units"),
                    meta: r.meta(e, None, |_, _| false),
                }
            });
            true
        });
        kl
    }

    fn event(&mut self, n: Node<'a, 'i>) -> Event {
        let mut e = E::new(n);
        let mut ev = Event {
            id: e.take_string("id"),
            name: e.take_string("name"),
            use_values_from_trigger_time: self.boolean(&mut e, "useValuesFromTriggerTime"),
            time_units: e.take_string("timeUnits"),
            ..Event::default()
        };
        let core = self.core;
        ev.meta = self.meta(e, None, |r, c| {
            match r.core_name(c) {
                Some("trigger") => {
                    let mut e = E::new(c);
                    let mut t = Trigger {
                        initial_value: r.boolean(&mut e, "initialValue"),
                        persistent: r.boolean(&mut e, "persistent"),
                        ..Trigger::default()
                    };
                    t.meta = r.meta(e, None, |r, c| r.math_child(c, &mut t.math));
                    ev.trigger = Some(t);
                }
                Some("delay") => {
                    let mut math = None;
                    let meta = r.meta(E::new(c), None, |r, c| r.math_child(c, &mut math));
                    ev.delay = Some(Delay { math, meta });
                }
                Some("priority") => {
                    let mut math = None;
                    let meta = r.meta(E::new(c), None, |r, c| r.math_child(c, &mut math));
                    ev.priority = Some(Priority { math, meta });
                }
                Some("listOfEventAssignments") => {
                    ev.assignments = r.list(c, core, &["eventAssignment"], |r, c| {
                        let mut e = E::new(c);
                        let variable = r.req(&mut e, "variable");
                        let mut math = None;
                        let meta = r.meta(e, None, |r, c| r.math_child(c, &mut math));
                        EventAssignment {
                            variable,
                            math,
                            meta,
                        }
                    })
                }
                _ => return false,
            }
            true
        });
        ev
    }
}

/// An SBML double: a number, or `INF`, `-INF` or `NaN`.
pub(crate) fn parse_double(v: &str) -> Option<f64> {
    match v {
        "INF" | "+INF" | "inf" | "Infinity" => Some(f64::INFINITY),
        "-INF" | "-inf" | "-Infinity" => Some(f64::NEG_INFINITY),
        "NaN" | "nan" => Some(f64::NAN),
        _ => v.parse().ok().filter(|x: &f64| x.is_finite()),
    }
}

pub(crate) fn escape_attr(s: &str) -> String {
    crate::write::escape(s, true)
}

/// Whether an element's start tag declares a namespace prefix (`None`: the
/// default namespace).
/// True if raw XML may use a namespace prefix: the default namespace always,
/// a named one if `prefix:` occurs in it.
fn uses(xml: &str, prefix: Option<&str>) -> bool {
    match prefix {
        None => true,
        Some(p) => xml.contains(&format!("{p}:")),
    }
}

fn declares(element: &str, prefix: Option<&str>) -> bool {
    let start_tag = element.split('>').next().unwrap_or("");
    let decl = match prefix {
        Some(p) => format!("xmlns:{p}"),
        None => "xmlns".to_string(),
    };
    start_tag.match_indices(&decl).any(|(i, _)| {
        start_tag[..i].ends_with(char::is_whitespace)
            && start_tag[i + decl.len()..].trim_start().starts_with('=')
    })
}

//! Writing SBML. A document is written at its own level and version; what is
//! `None` is left out, and kept XML is written back after the children
//! `sbml-io` knows.

mod layout;
mod math;
mod render;

use crate::math::Math;
use crate::model::*;
use crate::ns;

/// Writes a document as SBML. The level and version must be ones `read`
/// accepts; others are written with the Level 3 Version 2 namespace.
pub fn write(doc: &Document) -> String {
    let core = ns::core(doc.level, doc.version).unwrap_or(ns::core(3, 2).unwrap());
    let has_layouts = doc.model.as_ref().is_some_and(has_layout_list);
    let has_render = doc
        .model
        .as_ref()
        .is_some_and(|m| m.global_render.is_some() || m.layouts.iter().any(|l| l.render.is_some()));

    // The root's namespaces: as read, with the core namespace for this level
    // and version, and the Layout namespace if it's needed.
    let mut namespaces = doc.namespaces.clone();
    match namespaces.iter_mut().find(|(p, _)| p.is_none()) {
        Some(default) => default.1 = core.to_string(),
        None => namespaces.insert(0, (None, core.to_string())),
    }
    let mut root_attrs = doc.meta.extra.attributes.clone();
    let mut layout_prefix = String::from("layout");
    if doc.level == 3 && has_layouts {
        match namespaces.iter().find(|(_, u)| u == ns::LAYOUT) {
            Some((p, _)) => layout_prefix = p.clone().unwrap(),
            None => {
                while namespaces
                    .iter()
                    .any(|(p, _)| p.as_deref() == Some(&layout_prefix))
                {
                    layout_prefix.push('_');
                }
                namespaces.push((Some(layout_prefix.clone()), ns::LAYOUT.to_string()));
            }
        }
        if !root_attrs
            .iter()
            .any(|a| a.namespace.as_deref() == Some(ns::LAYOUT) && a.name.ends_with(":required"))
        {
            root_attrs.push(Attribute {
                name: format!("{layout_prefix}:required"),
                namespace: Some(ns::LAYOUT.to_string()),
                value: "false".to_string(),
            });
        }
    }

    // The same for Render, in Level 3. In Level 2 its namespace is declared
    // on each list, in the annotation.
    let mut render_prefix = String::from("render");
    if doc.level == 3 && has_layouts && has_render {
        match namespaces.iter().find(|(_, u)| u == ns::RENDER) {
            Some((Some(p), _)) => render_prefix = p.clone(),
            _ => {
                while namespaces
                    .iter()
                    .any(|(p, _)| p.as_deref() == Some(&render_prefix))
                {
                    render_prefix.push('_');
                }
                namespaces.push((Some(render_prefix.clone()), ns::RENDER.to_string()));
            }
        }
        if !root_attrs
            .iter()
            .any(|a| a.namespace.as_deref() == Some(ns::RENDER) && a.name.ends_with(":required"))
        {
            root_attrs.push(Attribute {
                name: format!("{render_prefix}:required"),
                namespace: Some(ns::RENDER.to_string()),
                value: "false".to_string(),
            });
        }
    }

    let mut w = W {
        out: String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"),
        depth: 0,
        level: doc.level,
        core,
        namespaces,
        layout_prefix,
        render_prefix,
    };
    let mut attrs: Vec<(String, String)> = w
        .namespaces
        .iter()
        .map(|(p, u)| match p {
            Some(p) => (format!("xmlns:{p}"), u.clone()),
            None => ("xmlns".to_string(), u.clone()),
        })
        .collect();
    attrs.push(("level".into(), doc.level.to_string()));
    attrs.push(("version".into(), doc.version.to_string()));
    let meta = Meta {
        extra: Extra {
            attributes: root_attrs,
            elements: doc.meta.extra.elements.clone(),
        },
        ..doc.meta.clone()
    };
    w.element("sbml", attrs, &meta, |w| {
        if let Some(m) = &doc.model {
            w.model(m);
        }
    });
    w.out
}

/// Whether a model has a `listOfLayouts` to write: layouts, or the list's
/// own content.
pub(crate) fn has_layout_list(m: &Model) -> bool {
    !m.layouts.is_empty() || m.global_render.is_some() || m.layouts_meta != Meta::default()
}

/// Attributes being gathered for an element: (qualified name, value).
pub(crate) type Attrs = Vec<(String, String)>;

pub(crate) struct W {
    pub(crate) out: String,
    depth: usize,
    pub(crate) level: u32,
    pub(crate) core: &'static str,
    namespaces: Vec<(Option<String>, String)>,
    pub(crate) layout_prefix: String,
    pub(crate) render_prefix: String,
}

/// Adds an attribute if it has a value.
pub(crate) fn attr<T: ToAttr + ?Sized>(attrs: &mut Attrs, name: &str, value: &T) {
    if let Some(v) = value.to_attr() {
        attrs.push((name.to_string(), v));
    }
}

pub(crate) trait ToAttr {
    fn to_attr(&self) -> Option<String>;
}

impl ToAttr for str {
    fn to_attr(&self) -> Option<String> {
        Some(self.to_string())
    }
}

impl ToAttr for String {
    fn to_attr(&self) -> Option<String> {
        Some(self.clone())
    }
}

impl ToAttr for f64 {
    fn to_attr(&self) -> Option<String> {
        Some(fmt_f64(*self))
    }
}

impl ToAttr for bool {
    fn to_attr(&self) -> Option<String> {
        Some(self.to_string())
    }
}

impl ToAttr for i32 {
    fn to_attr(&self) -> Option<String> {
        Some(self.to_string())
    }
}

impl<T: ToAttr> ToAttr for Option<T> {
    fn to_attr(&self) -> Option<String> {
        self.as_ref().and_then(ToAttr::to_attr)
    }
}

impl W {
    fn indent(&mut self) {
        for _ in 0..self.depth {
            self.out.push_str("  ");
        }
    }

    /// Writes an element: its attributes, then `metaid` and `sboTerm` and the
    /// kept attributes from `meta`; then notes and annotation, the children
    /// `body` writes, and the kept elements. An element with no content is
    /// written as `<tag/>`.
    pub(crate) fn element(
        &mut self,
        tag: &str,
        mut attrs: Attrs,
        meta: &Meta,
        body: impl FnOnce(&mut Self),
    ) {
        self.annotated(tag, &mut attrs, meta, None, body);
    }

    /// `element`, with an extra piece of raw XML put first in the annotation
    /// (Level 2 layouts).
    pub(crate) fn annotated(
        &mut self,
        tag: &str,
        attrs: &mut Attrs,
        meta: &Meta,
        annotation_head: Option<String>,
        body: impl FnOnce(&mut Self),
    ) {
        let mut all: Attrs = Vec::new();
        // In a Level 3 package element, metaid and sboTerm are core attributes.
        attr(&mut all, "metaid", &meta.metaid);
        attr(&mut all, "sboTerm", &meta.sbo_term);
        all.append(attrs);
        let mut declared: Vec<String> = Vec::new();
        for a in &meta.extra.attributes {
            // A prefix declared on an element `sbml-io` doesn't keep is
            // declared again here.
            if let (Some(uri), Some((prefix, _))) = (&a.namespace, a.name.split_once(':'))
                && prefix != "xml"
                && !self
                    .namespaces
                    .iter()
                    .any(|(p, u)| p.as_deref() == Some(prefix) && u == uri)
                && !declared.iter().any(|p| p == prefix)
            {
                all.push((format!("xmlns:{prefix}"), uri.clone()));
                declared.push(prefix.to_string());
            }
            all.push((a.name.clone(), a.value.clone()));
        }

        self.indent();
        self.out.push('<');
        self.out.push_str(tag);
        for (k, v) in &all {
            self.out.push(' ');
            self.out.push_str(k);
            self.out.push_str("=\"");
            self.out.push_str(&escape(v, true));
            self.out.push('"');
        }
        let open_end = self.out.len();
        self.out.push_str(">\n");
        self.depth += 1;
        if let Some(notes) = &meta.notes {
            self.raw_element("notes", notes);
        }
        match (&meta.annotation, annotation_head) {
            (None, None) => {}
            (a, head) => {
                let inner = format!("{}{}", head.unwrap_or_default(), a.as_deref().unwrap_or(""));
                self.raw_element("annotation", &inner);
            }
        }
        body(self);
        for e in &meta.extra.elements {
            self.indent();
            self.out.push_str(e);
            self.out.push('\n');
        }
        self.depth -= 1;
        if self.out.len() == open_end + 2 {
            self.out.truncate(open_end);
            self.out.push_str("/>\n");
        } else {
            self.indent();
            self.out.push_str("</");
            self.out.push_str(tag);
            self.out.push_str(">\n");
        }
    }

    /// `<tag>raw</tag>`, the raw content written as it was read.
    fn raw_element(&mut self, tag: &str, raw: &str) {
        self.indent();
        self.out.push_str(&format!("<{tag}>{raw}</{tag}>\n"));
    }

    /// A `listOf…` element, left out if there's nothing in it.
    pub(crate) fn list<T>(&mut self, tag: &str, items: &[T], mut f: impl FnMut(&mut Self, &T)) {
        if items.is_empty() {
            return;
        }
        self.indent();
        self.out.push_str(&format!("<{tag}>\n"));
        self.depth += 1;
        for item in items {
            f(self, item);
        }
        self.depth -= 1;
        self.indent();
        self.out.push_str(&format!("</{tag}>\n"));
    }

    /// A raw line at the current depth.
    pub(crate) fn line(&mut self, s: &str) {
        self.indent();
        self.out.push_str(s);
        self.out.push('\n');
    }

    pub(crate) fn push_depth(&mut self) {
        self.depth += 1;
    }

    pub(crate) fn pop_depth(&mut self) {
        self.depth -= 1;
    }

    fn math_child(&mut self, math: &Option<Math>) {
        if let Some(m) = math {
            self.math(m);
        }
    }

    fn model(&mut self, m: &Model) {
        let mut a = Vec::new();
        attr(&mut a, "id", &m.id);
        attr(&mut a, "name", &m.name);
        attr(&mut a, "substanceUnits", &m.substance_units);
        attr(&mut a, "timeUnits", &m.time_units);
        attr(&mut a, "volumeUnits", &m.volume_units);
        attr(&mut a, "areaUnits", &m.area_units);
        attr(&mut a, "lengthUnits", &m.length_units);
        attr(&mut a, "extentUnits", &m.extent_units);
        attr(&mut a, "conversionFactor", &m.conversion_factor);
        let l2_layouts = (self.level == 2 && has_layout_list(m)).then(|| {
            let mut xml = self.l2_layouts(m);
            if m.meta.annotation.is_none() {
                xml.push('\n');
                xml.push_str(&"  ".repeat(self.depth + 1));
            }
            xml
        });
        self.annotated("model", &mut a, &m.meta, l2_layouts, |w| {
            w.list(
                "listOfFunctionDefinitions",
                &m.function_definitions,
                |w, f| {
                    let mut a = Vec::new();
                    attr(&mut a, "id", &f.id);
                    attr(&mut a, "name", &f.name);
                    w.element("functionDefinition", a, &f.meta, |w| w.math_child(&f.math));
                },
            );
            w.list("listOfUnitDefinitions", &m.unit_definitions, |w, u| {
                let mut a = Vec::new();
                attr(&mut a, "id", &u.id);
                attr(&mut a, "name", &u.name);
                w.element("unitDefinition", a, &u.meta, |w| {
                    w.list("listOfUnits", &u.units, |w, u| {
                        let mut a = Vec::new();
                        attr(&mut a, "kind", &u.kind);
                        attr(&mut a, "exponent", &u.exponent);
                        attr(&mut a, "scale", &u.scale);
                        attr(&mut a, "multiplier", &u.multiplier);
                        attr(&mut a, "offset", &u.offset);
                        w.element("unit", a, &u.meta, |_| {});
                    });
                });
            });
            w.list("listOfCompartmentTypes", &m.compartment_types, |w, t| {
                let mut a = Vec::new();
                attr(&mut a, "id", &t.id);
                attr(&mut a, "name", &t.name);
                w.element("compartmentType", a, &t.meta, |_| {});
            });
            w.list("listOfSpeciesTypes", &m.species_types, |w, t| {
                let mut a = Vec::new();
                attr(&mut a, "id", &t.id);
                attr(&mut a, "name", &t.name);
                w.element("speciesType", a, &t.meta, |_| {});
            });
            w.list("listOfCompartments", &m.compartments, |w, c| {
                let mut a = Vec::new();
                attr(&mut a, "id", &c.id);
                attr(&mut a, "name", &c.name);
                attr(&mut a, "compartmentType", &c.compartment_type);
                attr(&mut a, "spatialDimensions", &c.spatial_dimensions);
                attr(&mut a, "size", &c.size);
                attr(&mut a, "units", &c.units);
                attr(&mut a, "outside", &c.outside);
                attr(&mut a, "constant", &c.constant);
                w.element("compartment", a, &c.meta, |_| {});
            });
            w.list("listOfSpecies", &m.species, |w, s| {
                let mut a = Vec::new();
                attr(&mut a, "id", &s.id);
                attr(&mut a, "name", &s.name);
                attr(&mut a, "speciesType", &s.species_type);
                attr(&mut a, "compartment", &s.compartment);
                attr(&mut a, "initialAmount", &s.initial_amount);
                attr(&mut a, "initialConcentration", &s.initial_concentration);
                attr(&mut a, "substanceUnits", &s.substance_units);
                attr(&mut a, "spatialSizeUnits", &s.spatial_size_units);
                attr(&mut a, "hasOnlySubstanceUnits", &s.has_only_substance_units);
                attr(&mut a, "boundaryCondition", &s.boundary_condition);
                attr(&mut a, "charge", &s.charge);
                attr(&mut a, "constant", &s.constant);
                attr(&mut a, "conversionFactor", &s.conversion_factor);
                w.element("species", a, &s.meta, |_| {});
            });
            w.list("listOfParameters", &m.parameters, |w, p| {
                let mut a = Vec::new();
                attr(&mut a, "id", &p.id);
                attr(&mut a, "name", &p.name);
                attr(&mut a, "value", &p.value);
                attr(&mut a, "units", &p.units);
                attr(&mut a, "constant", &p.constant);
                w.element("parameter", a, &p.meta, |_| {});
            });
            w.list(
                "listOfInitialAssignments",
                &m.initial_assignments,
                |w, ia| {
                    let mut a = Vec::new();
                    attr(&mut a, "symbol", &ia.symbol);
                    w.element("initialAssignment", a, &ia.meta, |w| w.math_child(&ia.math));
                },
            );
            w.list("listOfRules", &m.rules, |w, r| {
                let (tag, variable, math, meta) = match r {
                    Rule::Assignment {
                        variable,
                        math,
                        meta,
                    } => ("assignmentRule", Some(variable), math, meta),
                    Rule::Rate {
                        variable,
                        math,
                        meta,
                    } => ("rateRule", Some(variable), math, meta),
                    Rule::Algebraic { math, meta } => ("algebraicRule", None, math, meta),
                };
                let mut a = Vec::new();
                attr(&mut a, "variable", &variable.cloned());
                w.element(tag, a, meta, |w| w.math_child(math));
            });
            w.list("listOfConstraints", &m.constraints, |w, c| {
                w.element("constraint", Vec::new(), &c.meta, |w| {
                    w.math_child(&c.math);
                    if let Some(msg) = &c.message {
                        w.raw_element("message", msg);
                    }
                });
            });
            w.list("listOfReactions", &m.reactions, |w, r| w.reaction(r));
            w.list("listOfEvents", &m.events, |w, e| w.event(e));
            if w.level == 3 {
                w.l3_layouts(m);
            }
        });
    }

    fn reaction(&mut self, r: &Reaction) {
        let mut a = Vec::new();
        attr(&mut a, "id", &r.id);
        attr(&mut a, "name", &r.name);
        attr(&mut a, "reversible", &r.reversible);
        attr(&mut a, "fast", &r.fast);
        attr(&mut a, "compartment", &r.compartment);
        self.element("reaction", a, &r.meta, |w| {
            let sr = |w: &mut Self, s: &SpeciesReference| {
                let mut a = Vec::new();
                attr(&mut a, "id", &s.id);
                attr(&mut a, "name", &s.name);
                attr(&mut a, "species", &s.species);
                attr(&mut a, "stoichiometry", &s.stoichiometry);
                attr(&mut a, "constant", &s.constant);
                w.element("speciesReference", a, &s.meta, |w| {
                    if let Some(m) = &s.stoichiometry_math {
                        w.line("<stoichiometryMath>");
                        w.push_depth();
                        w.math(m);
                        w.pop_depth();
                        w.line("</stoichiometryMath>");
                    }
                });
            };
            w.list("listOfReactants", &r.reactants, sr);
            w.list("listOfProducts", &r.products, sr);
            w.list("listOfModifiers", &r.modifiers, |w, m| {
                let mut a = Vec::new();
                attr(&mut a, "id", &m.id);
                attr(&mut a, "name", &m.name);
                attr(&mut a, "species", &m.species);
                w.element("modifierSpeciesReference", a, &m.meta, |_| {});
            });
            if let Some(kl) = &r.kinetic_law {
                let mut a = Vec::new();
                attr(&mut a, "timeUnits", &kl.time_units);
                attr(&mut a, "substanceUnits", &kl.substance_units);
                w.element("kineticLaw", a, &kl.meta, |w| {
                    w.math_child(&kl.math);
                    let (list, item) = if w.level >= 3 {
                        ("listOfLocalParameters", "localParameter")
                    } else {
                        ("listOfParameters", "parameter")
                    };
                    w.list(list, &kl.parameters, |w, p| {
                        let mut a = Vec::new();
                        attr(&mut a, "id", &p.id);
                        attr(&mut a, "name", &p.name);
                        attr(&mut a, "value", &p.value);
                        attr(&mut a, "units", &p.units);
                        w.element(item, a, &p.meta, |_| {});
                    });
                });
            }
        });
    }

    fn event(&mut self, e: &Event) {
        let mut a = Vec::new();
        attr(&mut a, "id", &e.id);
        attr(&mut a, "name", &e.name);
        attr(
            &mut a,
            "useValuesFromTriggerTime",
            &e.use_values_from_trigger_time,
        );
        attr(&mut a, "timeUnits", &e.time_units);
        self.element("event", a, &e.meta, |w| {
            if let Some(t) = &e.trigger {
                let mut a = Vec::new();
                attr(&mut a, "initialValue", &t.initial_value);
                attr(&mut a, "persistent", &t.persistent);
                w.element("trigger", a, &t.meta, |w| w.math_child(&t.math));
            }
            if let Some(d) = &e.delay {
                w.element("delay", Vec::new(), &d.meta, |w| w.math_child(&d.math));
            }
            if let Some(p) = &e.priority {
                w.element("priority", Vec::new(), &p.meta, |w| w.math_child(&p.math));
            }
            w.list("listOfEventAssignments", &e.assignments, |w, ea| {
                let mut a = Vec::new();
                attr(&mut a, "variable", &ea.variable);
                w.element("eventAssignment", a, &ea.meta, |w| w.math_child(&ea.math));
            });
        });
    }

    /// Whether a namespace is declared on `<sbml>`.
    pub(crate) fn root_declares(&self, uri: &str) -> bool {
        self.namespaces.iter().any(|(_, u)| u == uri)
    }
}

/// A number as SBML writes it: integers without a decimal point, very large
/// and very small numbers in exponent form, and `INF`, `-INF` and `NaN`.
pub fn fmt_f64(x: f64) -> String {
    if x.is_nan() {
        "NaN".into()
    } else if x.is_infinite() {
        if x > 0.0 { "INF".into() } else { "-INF".into() }
    } else if x == x.trunc() && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else if x.abs() >= 1e15 || x.abs() < 1e-5 {
        format!("{x:e}")
    } else {
        format!("{x}")
    }
}

/// Escapes text for XML; `attr` also escapes double quotes.
pub(crate) fn escape(s: &str, attr: bool) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' if attr => out.push_str("&quot;"),
            '\n' if attr => out.push_str("&#10;"),
            '\t' if attr => out.push_str("&#9;"),
            '\r' => out.push_str("&#13;"),
            c => out.push(c),
        }
    }
    out
}

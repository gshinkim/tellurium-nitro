//! SBML Level 3 Version 2 output (Delphi `uAntSBMLWriter`, Delphi §8), built
//! on `sbml-io`. Only with the `sbml` feature.
//!
//! This follows the Delphi writer's mapping element by element, but builds an
//! [`sbml_io::Document`] rather than writing XML text: `sbml-io` writes it,
//! and checks it with [`sbml_io::validate`]. Like the Delphi writer it reads
//! the result of flattening, never the parse tree (Delphi §2.2.2 seam 2).
//!
//! Not written: the `fbc` package (objectives and flux bounds are modelled,
//! not emitted; Delphi §9.1.5). Annotations are written on the model and on
//! each element a symbol becomes; [`information_loss`] names any that are
//! not, so none is dropped silently. The mapping is described in `SPEC.md`
//! §5.7.
//!
//! ```
//! let loaded = antimony::parse("J1: S1 -> S2; k1*S1\nS1 = 10; S2 = 0; k1 = 0.1");
//! let xml = antimony::sbml::write_sbml(&loaded.document);
//! assert!(xml.contains(r#"<reaction id="J1" reversible="true">"#));
//!
//! // Or as data, to change or check before writing.
//! let flat = antimony::flatten::flatten(&loaded.document);
//! let document = antimony::sbml::to_sbml(&flat);
//! assert!(sbml_io::validate(&document).is_empty());
//! assert_eq!(document.model.unwrap().species.len(), 2);
//! ```

use std::collections::BTreeSet;

use sbml_io::math::{Constant as MathConstant, Math, Number, Op, Symbol as CSymbol};

use crate::diagnostics::{Code, Diagnostics};
use crate::expr::{BinaryOp, Constant, Expr, ExprKind, UnaryOp};
use crate::flatten::flatten;
use crate::model::{
    Constness, ConstraintKind, CvTerm, Document, Metadata, Module, Participant, Qualifier, Symbol,
    SymbolKind, UnitDef,
};
use crate::{messages, number};

const XHTML: &str = "http://www.w3.org/1999/xhtml";
const RDF: &str = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";
const DCTERMS: &str = "http://purl.org/dc/terms/";
const VCARD: &str = "http://www.w3.org/2001/vcard-rdf/3.0#";
const BQBIOL: &str = "http://biomodels.net/biology-qualifiers/";
const BQMODEL: &str = "http://biomodels.net/model-qualifiers/";

/// The compartment invented for species that name none, which SBML requires
/// and Antimony does not. libantimony uses the same id.
pub const DEFAULT_COMPARTMENT: &str = "default_compartment";

/// The SBML base unit kinds a `unit` definition may be built from.
const BASE_UNITS: [&str; 35] = [
    "ampere",
    "avogadro",
    "becquerel",
    "candela",
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
    "litre",
    "lumen",
    "lux",
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
    "liter",
    "meter",
];

/// The main module of a document, flattened, as SBML text.
pub fn write_sbml(document: &Document) -> String {
    sbml_io::write(&to_sbml(&flatten(document)))
}

/// A module, which should be the result of [`flatten`] (or
/// [`flatten_module`](crate::flatten::flatten_module) for one that is not
/// the main one), as SBML text.
pub fn write_sbml_module(module: &Module) -> String {
    sbml_io::write(&to_sbml(module))
}

/// A module, which should be the result of [`flatten`], as an SBML Level 3
/// Version 2 document.
pub fn to_sbml(module: &Module) -> sbml_io::Document {
    let mut document = sbml_io::Document::new(3, 2);
    document.model = Some(Converter { module }.model());
    document
}

/// The annotations of `module` (as given to [`to_sbml`]) that its SBML does
/// not carry, one [`Code::InformationLoss`] warning each: an annotated
/// symbol whose element is not written or has no annotation, and an
/// annotation in a vocabulary this library does not know, which is kept for
/// Antimony only.
///
/// ```
/// let loaded = antimony::parse("S1 = 1\nS1 frobnicates \"http://example.org/x\"");
/// let flat = antimony::flatten::flatten(&loaded.document);
/// let losses = antimony::sbml::information_loss(&flat);
/// assert_eq!(losses.len(), 1);
/// ```
pub fn information_loss(module: &Module) -> Diagnostics {
    let mut losses = Diagnostics::new();
    let mut lose = |what: &str, span| {
        losses.add(
            Code::InformationLoss,
            span,
            messages::annotation_dropped(what),
        );
    };
    let written = to_sbml(module).model.unwrap_or_default();
    let carried = |id: &str| {
        let m = &written;
        m.compartments
            .iter()
            .map(|c| (c.id.as_str(), &c.meta))
            .chain(m.species.iter().map(|s| (s.id.as_str(), &s.meta)))
            .chain(m.parameters.iter().map(|p| (p.id.as_str(), &p.meta)))
            .chain(m.reactions.iter().map(|r| (r.id.as_str(), &r.meta)))
            .chain(
                m.events
                    .iter()
                    .filter_map(|e| Some((e.id.as_deref()?, &e.meta))),
            )
            .chain(
                m.function_definitions
                    .iter()
                    .map(|f| (f.id.as_str(), &f.meta)),
            )
            .any(|(i, meta)| i == id && *meta != sbml_io::Meta::default())
    };
    for symbol in &module.symbols {
        let annotated = has_annotation(&symbol.cv_terms, &symbol.metadata, &symbol.notes)
            || symbol.sbo_term.is_some();
        if annotated && !carried(&symbol.name) {
            lose(&symbol.name, symbol.span);
        }
    }
    for section in &module.unparsed_sections {
        let subject = section
            .text
            .split_whitespace()
            .next()
            .unwrap_or("the model");
        lose(subject, section.span);
    }
    losses
}

// --- text ------------------------------------------------------------------

/// Escapes text for an XML attribute value.
pub fn xml_escape(text: &str) -> String {
    xml_escape_text(text).replace('"', "&quot;")
}

/// Escapes element content, where a quote is an ordinary character.
pub fn xml_escape_text(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// `Some` of a string that is not empty: an empty display name or note is
/// written as nothing, as in Delphi.
fn present(value: &Option<String>) -> Option<String> {
    value.clone().filter(|v| !v.is_empty())
}

/// An RDF description is anchored by `rdf:about="#metaid"`, and Antimony has
/// no syntax for a metaid, so one is made. Prefixed, because an XML ID must
/// be an NCName and an SBML id may begin with a digit.
fn meta_id(id: &str) -> String {
    format!("meta_{id}")
}

fn has_annotation(terms: &[CvTerm], metadata: &Metadata, notes: &Option<String>) -> bool {
    present(notes).is_some() || !terms.is_empty() || !metadata.is_empty()
}

/// `<notes>`' content. SBML notes are XHTML. A note that already is XHTML,
/// as almost every real one is (it came from an SBML file), goes out as
/// markup in a `<body>`; prose is escaped in a `<p>`. The test is the first
/// non-space character.
fn notes_xml(notes: &str) -> String {
    if notes.trim_start().starts_with('<') {
        let mut out = format!("\n<body xmlns=\"{XHTML}\">\n");
        for line in notes.split('\n') {
            out.push_str(line.trim_end());
            out.push('\n');
        }
        out.push_str("</body>\n");
        out
    } else {
        format!("<p xmlns=\"{XHTML}\">{}</p>", xml_escape_text(notes))
    }
}

/// The RDF `<annotation>` content: creators, dates, then CV terms, or `None`
/// when there is nothing to say.
fn annotation_xml(
    metaid: &str,
    terms: &[CvTerm],
    metadata: &Metadata,
    for_model: bool,
) -> Option<String> {
    if terms.is_empty() && metadata.is_empty() {
        return None;
    }
    let mut out = String::new();
    let mut line = |depth: usize, text: &str| {
        out.push_str(&"  ".repeat(depth));
        out.push_str(text);
        out.push('\n');
    };
    line(
        0,
        &format!(
            "<rdf:RDF xmlns:rdf=\"{RDF}\" xmlns:dcterms=\"{DCTERMS}\" xmlns:vCard=\"{VCARD}\" xmlns:bqbiol=\"{BQBIOL}\" xmlns:bqmodel=\"{BQMODEL}\">"
        ),
    );
    line(
        1,
        &format!("<rdf:Description rdf:about=\"#{}\">", xml_escape(metaid)),
    );

    let creators: Vec<_> = metadata.creators.iter().filter(|c| !c.is_empty()).collect();
    if !metadata.creators.is_empty() {
        line(2, "<dcterms:creator>");
        line(3, "<rdf:Bag>");
        for c in creators {
            line(4, "<rdf:li rdf:parseType=\"Resource\">");
            // vCard:N carries the split name. An unsplit `creator.name` is not
            // guessed at: it goes in Family whole.
            let family = present(&c.family_name);
            let given = present(&c.given_name);
            let name = present(&c.name);
            if family.is_some() || given.is_some() || name.is_some() {
                line(5, "<vCard:N rdf:parseType=\"Resource\">");
                if family.is_some() || given.is_some() {
                    if let Some(f) = &family {
                        line(
                            6,
                            &format!("<vCard:Family>{}</vCard:Family>", xml_escape_text(f)),
                        );
                    }
                    if let Some(g) = &given {
                        line(
                            6,
                            &format!("<vCard:Given>{}</vCard:Given>", xml_escape_text(g)),
                        );
                    }
                } else if let Some(n) = &name {
                    line(
                        6,
                        &format!("<vCard:Family>{}</vCard:Family>", xml_escape_text(n)),
                    );
                }
                line(5, "</vCard:N>");
            }
            if let Some(email) = present(&c.email) {
                line(
                    5,
                    &format!("<vCard:EMAIL>{}</vCard:EMAIL>", xml_escape_text(&email)),
                );
            }
            if let Some(org) = present(&c.organization) {
                line(5, "<vCard:ORG rdf:parseType=\"Resource\">");
                line(
                    6,
                    &format!("<vCard:Orgname>{}</vCard:Orgname>", xml_escape_text(&org)),
                );
                line(5, "</vCard:ORG>");
            }
            line(4, "</rdf:li>");
        }
        line(3, "</rdf:Bag>");
        line(2, "</dcterms:creator>");
    }

    let dates = metadata
        .created
        .iter()
        .map(|d| ("created", d))
        .chain(metadata.modified.iter().map(|d| ("modified", d)));
    for (element, date) in dates {
        if date.is_empty() && element == "created" {
            continue;
        }
        line(
            2,
            &format!("<dcterms:{element} rdf:parseType=\"Resource\">"),
        );
        line(
            3,
            &format!("<dcterms:W3CDTF>{}</dcterms:W3CDTF>", xml_escape_text(date)),
        );
        line(2, &format!("</dcterms:{element}>"));
    }

    // Terms sharing a qualifier go in one rdf:Bag: that is what the element
    // means. A model qualifier on an element was warned about when read and
    // has nowhere to go in SBML.
    let mut done: Vec<Qualifier> = Vec::new();
    for term in terms {
        let q = term.qualifier;
        if done.contains(&q) || (q.is_model_only() && !for_model) {
            continue;
        }
        done.push(q);
        // The bqbiol/bqmodel split is by subject: `description` on a model is
        // bqmodel:isDescribedBy, on a species bqbiol:isDescribedBy.
        let prefix = if for_model
            && matches!(
                q,
                Qualifier::ModelIs | Qualifier::IsDerivedFrom | Qualifier::IsDescribedBy
            ) {
            "bqmodel"
        } else {
            "bqbiol"
        };
        let element = format!("{prefix}:{}", q.element());
        line(2, &format!("<{element}>"));
        line(3, "<rdf:Bag>");
        for same in terms.iter().filter(|t| t.qualifier == q) {
            line(
                4,
                &format!("<rdf:li rdf:resource=\"{}\"/>", xml_escape(&same.resource)),
            );
        }
        line(3, "</rdf:Bag>");
        line(2, &format!("</{element}>"));
    }

    line(1, "</rdf:Description>");
    line(0, "</rdf:RDF>");
    Some(format!("\n{out}"))
}

/// An SBO term number as SBML writes it: 247 is `SBO:0000247`.
fn sbo_text(sbo: Option<u32>) -> Option<String> {
    sbo.map(|n| format!("SBO:{n:07}"))
}

/// The `meta` of an element with this id and annotation: a metaid only
/// where there is RDF or a note to anchor. The SBO term is an attribute and
/// needs no metaid.
fn meta(
    id: &str,
    terms: &[CvTerm],
    metadata: &Metadata,
    notes: &Option<String>,
    sbo: Option<u32>,
    for_model: bool,
) -> sbml_io::Meta {
    let sbo_term = sbo_text(sbo);
    if !has_annotation(terms, metadata, notes) {
        return sbml_io::Meta {
            sbo_term,
            ..sbml_io::Meta::default()
        };
    }
    let metaid = meta_id(id);
    sbml_io::Meta {
        notes: present(notes).map(|n| notes_xml(&n)),
        annotation: annotation_xml(&metaid, terms, metadata, for_model),
        metaid: Some(metaid),
        sbo_term,
        ..sbml_io::Meta::default()
    }
}

fn symbol_meta(symbol: &Symbol) -> sbml_io::Meta {
    meta(
        &symbol.name,
        &symbol.cv_terms,
        &symbol.metadata,
        &symbol.notes,
        symbol.sbo_term,
        false,
    )
}

// --- MathML ------------------------------------------------------------------

fn csymbol(url: &str, name: &str) -> CSymbol {
    CSymbol::new(url, name)
}

/// An expression as Content MathML.
pub fn to_math(expr: &Expr) -> Math {
    match &expr.kind {
        ExprKind::Number { value, units } => {
            if value.is_nan() {
                Math::Constant(MathConstant::NotANumber)
            } else if value.is_infinite() {
                let infinity = Math::Constant(MathConstant::Infinity);
                if *value > 0.0 {
                    infinity
                } else {
                    Math::Apply(Op::Minus, vec![infinity])
                }
            } else {
                Math::Number {
                    value: Number::Real(*value),
                    units: units.clone(),
                }
            }
        }
        ExprKind::Constant(c) => Math::Constant(match c {
            Constant::Pi => MathConstant::Pi,
            Constant::ExponentialE => MathConstant::ExponentialE,
            Constant::True => MathConstant::True,
            Constant::False => MathConstant::False,
        }),
        // `time` and `avogadro` are csymbols in SBML, not identifiers.
        ExprKind::Symbol(name) => match name.as_str() {
            "time" => Math::Symbol(csymbol(CSymbol::TIME, "time")),
            "avogadro" => Math::Symbol(csymbol(CSymbol::AVOGADRO, "avogadro")),
            _ => Math::Ident(name.clone()),
        },
        ExprKind::Unary { op, operand } => match op {
            UnaryOp::Minus => Math::Apply(Op::Minus, vec![to_math(operand)]),
            UnaryOp::Not => Math::Apply(Op::Not, vec![to_math(operand)]),
            UnaryOp::Plus => to_math(operand),
        },
        ExprKind::Binary { op, left, right } => {
            let op = match op {
                BinaryOp::Add => Op::Plus,
                BinaryOp::Subtract => Op::Minus,
                BinaryOp::Multiply => Op::Times,
                BinaryOp::Divide => Op::Divide,
                BinaryOp::Modulo => Op::Rem,
                BinaryOp::Power => Op::Power,
                BinaryOp::Equal => Op::Eq,
                BinaryOp::NotEqual => Op::Neq,
                BinaryOp::Less => Op::Lt,
                BinaryOp::LessEqual => Op::Leq,
                BinaryOp::Greater => Op::Gt,
                BinaryOp::GreaterEqual => Op::Geq,
                BinaryOp::And => Op::And,
                BinaryOp::Or => Op::Or,
            };
            Math::Apply(op, vec![to_math(left), to_math(right)])
        }
        // MathML has no ternary operator; a piecewise says the same.
        ExprKind::Conditional {
            condition,
            then,
            otherwise,
        } => Math::Piecewise {
            pieces: vec![(to_math(then), to_math(condition))],
            otherwise: Some(Box::new(to_math(otherwise))),
        },
        ExprKind::Call { name, args } => call(name, args),
    }
}

fn call(name: &str, args: &[Expr]) -> Math {
    let maths = || args.iter().map(to_math).collect::<Vec<_>>();
    match (name, args) {
        // piecewise(v1, c1, v2, c2, …, otherwise).
        ("piecewise", _) => {
            let pieces = args
                .chunks_exact(2)
                .map(|pair| (to_math(&pair[0]), to_math(&pair[1])))
                .collect();
            let otherwise = (args.len() % 2 == 1).then(|| Box::new(to_math(&args[args.len() - 1])));
            Math::Piecewise { pieces, otherwise }
        }
        ("delay", _) => Math::CallSymbol(csymbol(CSymbol::DELAY, "delay"), maths()),
        ("rateOf", _) => Math::CallSymbol(csymbol(CSymbol::RATE_OF, "rateOf"), maths()),
        // MathML has no sqrt or log10: a root of degree 2, a log of base 10.
        ("sqrt", [x]) => Math::Root {
            degree: Some(Box::new(Math::int(2))),
            arg: Box::new(to_math(x)),
        },
        ("log10", [x]) => Math::Log {
            base: Some(Box::new(Math::int(10))),
            arg: Box::new(to_math(x)),
        },
        // root and log take their degree and base as qualifiers, first.
        ("root", [x]) => Math::Root {
            degree: None,
            arg: Box::new(to_math(x)),
        },
        ("root", [n, x]) => Math::Root {
            degree: Some(Box::new(to_math(n))),
            arg: Box::new(to_math(x)),
        },
        ("log", [x]) => Math::Log {
            base: None,
            arg: Box::new(to_math(x)),
        },
        ("log", [b, x]) => Math::Log {
            base: Some(Box::new(to_math(b))),
            arg: Box::new(to_math(x)),
        },
        _ => {
            // The names the two vocabularies spell differently.
            let mathml = match name {
                "ceil" => "ceiling",
                "pow" => "power",
                "asin" => "arcsin",
                "acos" => "arccos",
                "atan" => "arctan",
                other => other,
            };
            match Op::from_tag(mathml) {
                Some(op) => Math::Apply(op, maths()),
                // A call to a function definition.
                None => Math::Call(name.to_string(), maths()),
            }
        }
    }
}

// --- the element mapping -------------------------------------------------------

/// A literal number SBML can carry as an attribute: a number, or a signed
/// one (`-0.1` is a minus over 0.1). A constant expression such as `3/2` is
/// not folded: as an initial assignment it keeps what the author wrote.
fn literal(expr: Option<&Expr>) -> Option<f64> {
    let expr = expr?;
    match &expr.kind {
        ExprKind::Number { value, .. } => Some(*value),
        ExprKind::Unary {
            op: UnaryOp::Minus,
            operand,
        } => literal(Some(operand)).map(|v| -v),
        ExprKind::Unary {
            op: UnaryOp::Plus,
            operand,
        } => literal(Some(operand)),
        _ => None,
    }
}

/// One side of a unit definition as `<unit>`s, or `None` when it is not
/// something SBML can say: juxtaposed base units, `^` exponents and a
/// leading number as multiplier. `metre^2` is one unit with exponent 2, not
/// a kind called "metre^2" (SBML rule 20410).
fn unit_factors(text: &str, sign: f64, into: &mut Vec<sbml_io::Unit>) -> Option<()> {
    let factors: Vec<&str> = text
        .split(' ')
        .map(str::trim)
        .filter(|f| !f.is_empty())
        .collect();
    let mut multiplier = 1.0;
    for (i, factor) in factors.iter().enumerate() {
        if let Some(value) = number::parse(factor) {
            // A trailing number multiplies nothing.
            if i == factors.len() - 1 {
                return None;
            }
            multiplier = value;
            continue;
        }
        let (kind, exponent) = match factor.split_once('^') {
            Some((kind, exponent)) => (kind, number::parse(exponent)?),
            None => (*factor, 1.0),
        };
        let kind = kind.trim();
        if !BASE_UNITS.contains(&kind) {
            return None;
        }
        into.push(sbml_io::Unit {
            kind: kind.to_string(),
            exponent: Some(sign * exponent),
            scale: Some(0),
            multiplier: Some(multiplier),
            ..sbml_io::Unit::default()
        });
        multiplier = 1.0;
    }
    Some(())
}

/// A unit definition. One it cannot represent (a bracketed expression, a
/// reference to another unit, `unit time_unit = time_unit`) is written
/// empty: incomplete rather than wrong.
fn unit_definition(unit: &UnitDef) -> sbml_io::UnitDefinition {
    let mut units = Vec::new();
    if !unit.is_self_referential() {
        let sides: Vec<&str> = unit.definition.trim().split('/').collect();
        let ok = (1..=2).contains(&sides.len())
            && unit_factors(sides[0], 1.0, &mut units).is_some()
            && (sides.len() == 1 || unit_factors(sides[1], -1.0, &mut units).is_some());
        if !ok {
            units.clear();
        }
    }
    sbml_io::UnitDefinition {
        id: unit.name.clone(),
        units,
        meta: meta(
            &unit.name,
            &unit.cv_terms,
            &unit.metadata,
            &None,
            None,
            false,
        ),
        ..sbml_io::UnitDefinition::default()
    }
}

/// An event modifier read as a boolean, or `default` when not given or not
/// computable. The three flags SBML needs from these are required in Level 3,
/// with no default, so they are always written; all three default to true.
fn event_flag(expr: Option<&Expr>, default: bool) -> bool {
    expr.and_then(Expr::evaluate)
        .map_or(default, |value| value != 0.0)
}

struct Converter<'m> {
    module: &'m Module,
}

impl Converter<'_> {
    fn symbols(&self, kind: SymbolKind) -> impl Iterator<Item = &Symbol> {
        self.module.symbols_of(kind)
    }

    fn model(&self) -> sbml_io::Model {
        let m = self.module;
        let id = m.name.clone().unwrap_or_default();
        sbml_io::Model {
            id: m.name.clone(),
            name: present(&m.display_name),
            conversion_factor: m.conversion_factor.clone(),
            function_definitions: self.function_definitions(),
            unit_definitions: m.units.iter().map(unit_definition).collect(),
            compartments: self.compartments(),
            species: self.species(),
            parameters: self.parameters(),
            initial_assignments: self.initial_assignments(),
            rules: self.rules(),
            constraints: self.constraints(),
            reactions: self.reactions(),
            events: self.events(),
            meta: meta(&id, &m.cv_terms, &m.metadata, &m.notes, m.sbo_term, true),
            ..sbml_io::Model::default()
        }
    }

    fn function_definitions(&self) -> Vec<sbml_io::FunctionDefinition> {
        self.module
            .functions
            .iter()
            .map(|f| sbml_io::FunctionDefinition {
                id: f.name.clone(),
                math: f.body.as_ref().map(|body| Math::Lambda {
                    args: f.args.clone(),
                    body: Box::new(to_math(body)),
                }),
                // `mm.sboTerm = 28`, `mm identity "…"` and `mm notes` land
                // on the symbol named after the function.
                meta: self.annotation_meta(&f.name),
                ..sbml_io::FunctionDefinition::default()
            })
            .collect()
    }

    /// Whether some species names no compartment.
    fn needs_default_compartment(&self) -> bool {
        self.symbols(SymbolKind::Species)
            .any(|s| s.compartment.is_none())
    }

    fn compartments(&self) -> Vec<sbml_io::Compartment> {
        let mut out = Vec::new();
        // Unless the model already has one by that name, which a model that
        // libantimony imported from SBML can: two with one id is invalid.
        if self.needs_default_compartment() && self.module.symbol(DEFAULT_COMPARTMENT).is_none() {
            out.push(sbml_io::Compartment {
                id: DEFAULT_COMPARTMENT.into(),
                spatial_dimensions: Some(3.0),
                size: Some(1.0),
                constant: Some(true),
                meta: sbml_io::Meta {
                    // SBO:0000410, implicit compartment: the conversion made it.
                    sbo_term: Some("SBO:0000410".into()),
                    ..sbml_io::Meta::default()
                },
                ..sbml_io::Compartment::default()
            });
        }
        for s in self.symbols(SymbolKind::Compartment) {
            out.push(sbml_io::Compartment {
                id: s.name.clone(),
                name: present(&s.display_name),
                spatial_dimensions: Some(3.0),
                size: literal(s.initial_value.as_ref()),
                units: s.units.clone(),
                constant: Some(self.is_constant(s)),
                meta: symbol_meta(s),
                ..sbml_io::Compartment::default()
            });
        }
        out
    }

    fn species(&self) -> Vec<sbml_io::Species> {
        self.symbols(SymbolKind::Species)
            .map(|s| {
                // `substanceOnly` makes the value an amount.
                let value = literal(s.initial_value.as_ref());
                let (amount, concentration) = if s.is_substance_only {
                    (value, None)
                } else {
                    (None, value)
                };
                sbml_io::Species {
                    id: s.name.clone(),
                    name: present(&s.display_name),
                    compartment: s
                        .compartment
                        .clone()
                        .unwrap_or_else(|| DEFAULT_COMPARTMENT.into()),
                    initial_amount: amount,
                    initial_concentration: concentration,
                    substance_units: s.units.clone(),
                    has_only_substance_units: Some(s.is_substance_only),
                    boundary_condition: Some(s.is_boundary),
                    constant: Some(s.constness == Constness::Const),
                    conversion_factor: s.conversion_factor.clone(),
                    meta: symbol_meta(s),
                    ..sbml_io::Species::default()
                }
            })
            .collect()
    }

    fn parameters(&self) -> Vec<sbml_io::Parameter> {
        self.symbols(SymbolKind::Parameter)
            .map(|s| sbml_io::Parameter {
                id: s.name.clone(),
                name: present(&s.display_name),
                value: literal(s.initial_value.as_ref()),
                units: s.units.clone(),
                constant: Some(self.is_constant(s)),
                meta: symbol_meta(s),
            })
            .collect()
    }

    /// `=` with anything but a literal is an initial assignment; a literal is
    /// the element's own attribute. The same test both ways, so every value
    /// lands in exactly one place.
    fn initial_assignments(&self) -> Vec<sbml_io::InitialAssignment> {
        self.module
            .symbols
            .iter()
            .filter_map(|s| {
                let value = s.initial_value.as_ref()?;
                literal(Some(value))
                    .is_none()
                    .then(|| sbml_io::InitialAssignment {
                        symbol: s.name.clone(),
                        math: Some(to_math(value)),
                        ..sbml_io::InitialAssignment::default()
                    })
            })
            .collect()
    }

    /// Assignment rules, then rate rules, then algebraic rules. A rule with no
    /// expression (`p' = ;`) is not SBML and is left out.
    fn rules(&self) -> Vec<sbml_io::Rule> {
        let mut out = Vec::new();
        for kind in [
            crate::model::RuleKind::Assignment,
            crate::model::RuleKind::Rate,
        ] {
            for s in &self.module.symbols {
                let Some(rule) = s.rule.as_ref().filter(|r| r.kind == kind) else {
                    continue;
                };
                let Some(expr) = &rule.expr else {
                    continue;
                };
                let (variable, math, meta) =
                    (s.name.clone(), Some(to_math(expr)), Default::default());
                out.push(match kind {
                    crate::model::RuleKind::Assignment => sbml_io::Rule::Assignment {
                        variable,
                        math,
                        meta,
                    },
                    crate::model::RuleKind::Rate => sbml_io::Rule::Rate {
                        variable,
                        math,
                        meta,
                    },
                });
            }
        }
        for rule in &self.module.algebraic_rules {
            if let Some(expr) = &rule.expr {
                out.push(sbml_io::Rule::Algebraic {
                    math: Some(to_math(expr)),
                    meta: Default::default(),
                });
            }
        }
        out
    }

    /// A flux bound names a reaction, which core SBML does not allow in a
    /// constraint: it belongs to fbc, which is not written.
    fn is_flux_bound(&self, expr: &Expr) -> bool {
        expr.symbols().iter().any(|name| {
            self.module
                .symbol(name)
                .is_some_and(|s| s.kind == SymbolKind::Reaction)
        })
    }

    fn constraints(&self) -> Vec<sbml_io::Constraint> {
        self.module
            .constraints
            .iter()
            .filter_map(|c| match &c.kind {
                ConstraintKind::Expr(e) if !self.is_flux_bound(e) => Some(sbml_io::Constraint {
                    math: Some(to_math(e)),
                    ..sbml_io::Constraint::default()
                }),
                _ => None,
            })
            .collect()
    }

    /// A compartment's or parameter's `constant`. Constant unless declared
    /// `var` or something changes it: a rule or an event assignment makes it
    /// variable whether or not it was declared so (SBML rules 20903 and
    /// 21213; the Delphi writer read the declaration alone, R30).
    fn is_constant(&self, s: &Symbol) -> bool {
        s.constness != Constness::Var && !self.is_time_varying(&s.name)
    }

    /// Whether something changes `name` after the start: a rule, or an event
    /// assigning to it. That is what `constant="false"` says.
    fn is_time_varying(&self, name: &str) -> bool {
        self.module.symbol(name).is_some_and(|s| s.rule.is_some())
            || self
                .module
                .events
                .iter()
                .any(|e| e.assignments.iter().any(|a| a.target == name))
    }

    fn species_reference(&self, p: &Participant) -> sbml_io::SpeciesReference {
        let mut reference = sbml_io::SpeciesReference {
            species: p.species.clone(),
            ..sbml_io::SpeciesReference::default()
        };
        match &p.stoichiometry {
            None => {
                reference.stoichiometry = Some(1.0);
                reference.constant = Some(true);
            }
            Some(s) => match literal(Some(s)) {
                Some(value) => {
                    reference.stoichiometry = Some(value);
                    reference.constant = Some(true);
                }
                None => {
                    // A named stoichiometry: L3V2's replacement for
                    // StoichiometryMath is an id on the reference. The id is
                    // the quantity, so its value goes here too (case 01434).
                    let name = s.to_string();
                    if let ExprKind::Symbol(symbol) = &s.kind {
                        reference.stoichiometry = self
                            .module
                            .symbol(symbol)
                            .and_then(|sym| literal(sym.initial_value.as_ref()));
                    }
                    reference.constant = Some(!self.is_time_varying(&name));
                    reference.id = Some(name);
                }
            },
        }
        reference
    }

    /// SBML's modifiers, which Antimony has no statement for: every species
    /// the kinetic law uses that is not a reactant or product must be one
    /// (rule 21121). An interaction naming the reaction is a modifier too,
    /// with its SBO term and label, and replaces the bare one for its species.
    fn modifiers(
        &self,
        reaction: &crate::model::Reaction,
    ) -> Vec<sbml_io::ModifierSpeciesReference> {
        let sides: Vec<&str> = reaction
            .reactants
            .iter()
            .chain(&reaction.products)
            .map(|p| p.species.as_str())
            .collect();
        let mut derived: BTreeSet<&str> = reaction
            .modifiers
            .iter()
            .map(|p| p.species.as_str())
            .collect();
        if let Some(law) = &reaction.kinetic_law {
            for name in law.symbols() {
                if self
                    .module
                    .symbol(name)
                    .is_some_and(|s| s.kind == SymbolKind::Species)
                {
                    derived.insert(name);
                }
            }
        }
        derived.retain(|name| !sides.contains(name));

        let mut out = Vec::new();
        for i in self
            .module
            .interactions
            .iter()
            .filter(|i| i.reaction == reaction.name)
        {
            derived.remove(i.species.as_str());
            out.push(sbml_io::ModifierSpeciesReference {
                name: i.name.clone(),
                species: i.species.clone(),
                meta: sbml_io::Meta {
                    sbo_term: Some(format!("SBO:{:07}", i.kind.sbo_term())),
                    ..sbml_io::Meta::default()
                },
                ..sbml_io::ModifierSpeciesReference::default()
            });
        }
        for name in derived {
            out.push(sbml_io::ModifierSpeciesReference {
                species: name.to_string(),
                ..sbml_io::ModifierSpeciesReference::default()
            });
        }
        out
    }

    /// The symbol carrying a reaction's or event's annotation: `J1 identity
    /// "…"` names the reaction, whose name is a symbol.
    fn annotation_meta(&self, id: &str) -> sbml_io::Meta {
        self.module.symbol(id).map(symbol_meta).unwrap_or_default()
    }

    fn reactions(&self) -> Vec<sbml_io::Reaction> {
        self.module
            .reactions
            .iter()
            .map(|r| sbml_io::Reaction {
                id: r.name.clone(),
                // `->` is reversible and `=>` is not (Delphi §9.1.1).
                reversible: Some(r.is_reversible()),
                reactants: r
                    .reactants
                    .iter()
                    .map(|p| self.species_reference(p))
                    .collect(),
                products: r
                    .products
                    .iter()
                    .map(|p| self.species_reference(p))
                    .collect(),
                modifiers: self.modifiers(r),
                kinetic_law: r.kinetic_law.as_ref().map(|law| sbml_io::KineticLaw {
                    math: Some(to_math(law)),
                    ..sbml_io::KineticLaw::default()
                }),
                meta: if r.name_generated {
                    sbml_io::Meta::default()
                } else {
                    self.annotation_meta(&r.name)
                },
                ..sbml_io::Reaction::default()
            })
            .collect()
    }

    fn events(&self) -> Vec<sbml_io::Event> {
        self.module
            .events
            .iter()
            .enumerate()
            .map(|(i, e)| {
                let id = e.name.clone().unwrap_or_else(|| format!("_E{}", i + 1));
                let math = |expr: &Option<Expr>| expr.as_ref().map(to_math);
                sbml_io::Event {
                    use_values_from_trigger_time: Some(event_flag(e.from_trigger.as_ref(), true)),
                    trigger: Some(sbml_io::Trigger {
                        // `t0` is SBML's initialValue.
                        initial_value: Some(event_flag(e.initial_value.as_ref(), true)),
                        persistent: Some(event_flag(e.persistent.as_ref(), true)),
                        math: Some(to_math(&e.trigger)),
                        ..sbml_io::Trigger::default()
                    }),
                    priority: math(&e.priority).map(|m| sbml_io::Priority {
                        math: Some(m),
                        ..sbml_io::Priority::default()
                    }),
                    delay: math(&e.delay).map(|m| sbml_io::Delay {
                        math: Some(m),
                        ..sbml_io::Delay::default()
                    }),
                    assignments: e
                        .assignments
                        .iter()
                        .map(|a| sbml_io::EventAssignment {
                            variable: a.target.clone(),
                            math: math(&a.expr),
                            ..sbml_io::EventAssignment::default()
                        })
                        .collect(),
                    meta: self.annotation_meta(&id),
                    id: Some(id),
                    ..sbml_io::Event::default()
                }
            })
            .collect()
    }
}

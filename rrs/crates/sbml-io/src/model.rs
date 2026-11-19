//! The SBML core: a document and its model, as plain data.
//!
//! Attributes are stored as written: `None` means the attribute was absent,
//! not that it has its default value. See the crate documentation.

use crate::layout::Layout;
use crate::math::Math;

/// An SBML document: the `<sbml>` element.
#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub level: u32,
    pub version: u32,
    /// Namespace declarations on `<sbml>`, in order: (prefix, URI), with no
    /// prefix for the default namespace. Kept so that preserved XML from other
    /// packages stays valid. May be left empty: the writer always declares the
    /// core namespace, and the Layout namespace when there are layouts.
    pub namespaces: Vec<(Option<String>, String)>,
    pub model: Option<Model>,
    pub meta: Meta,
}

impl Document {
    /// An empty document at the given level and version.
    pub fn new(level: u32, version: u32) -> Self {
        Document {
            level,
            version,
            namespaces: Vec::new(),
            model: None,
            meta: Meta::default(),
        }
    }
}

impl Default for Document {
    /// An empty Level 3 Version 2 document.
    fn default() -> Self {
        Document::new(3, 2)
    }
}

/// What any SBML element may carry, besides its own attributes and children.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Meta {
    pub metaid: Option<String>,
    /// As written, for example `SBO:0000176`.
    pub sbo_term: Option<String>,
    /// The content of `<notes>`, as raw XML.
    pub notes: Option<String>,
    /// The content of `<annotation>`, as raw XML.
    pub annotation: Option<String>,
    pub extra: Extra,
}

/// What `sbml-io` doesn't understand on an element, kept to be written back.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Extra {
    pub attributes: Vec<Attribute>,
    /// Child elements, each as raw XML, in document order. They are written
    /// back after the children `sbml-io` knows.
    pub elements: Vec<String>,
}

impl Extra {
    pub fn is_empty(&self) -> bool {
        self.attributes.is_empty() && self.elements.is_empty()
    }
}

/// An attribute kept as written.
#[derive(Debug, Clone, PartialEq)]
pub struct Attribute {
    /// The qualified name as written, for example `fbc:charge`.
    pub name: String,
    /// The namespace URI, if the name has a prefix.
    pub namespace: Option<String>,
    pub value: String,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Model {
    pub id: Option<String>,
    pub name: Option<String>,
    pub substance_units: Option<String>,
    pub time_units: Option<String>,
    pub volume_units: Option<String>,
    pub area_units: Option<String>,
    pub length_units: Option<String>,
    pub extent_units: Option<String>,
    pub conversion_factor: Option<String>,
    pub function_definitions: Vec<FunctionDefinition>,
    pub unit_definitions: Vec<UnitDefinition>,
    /// Level 2 Versions 2 to 4 only.
    pub compartment_types: Vec<CompartmentType>,
    /// Level 2 Versions 2 to 4 only.
    pub species_types: Vec<SpeciesType>,
    pub compartments: Vec<Compartment>,
    pub species: Vec<Species>,
    pub parameters: Vec<Parameter>,
    pub initial_assignments: Vec<InitialAssignment>,
    pub rules: Vec<Rule>,
    pub constraints: Vec<Constraint>,
    pub reactions: Vec<Reaction>,
    pub events: Vec<Event>,
    /// The Layout package. In Level 2 these live in the model's annotation.
    pub layouts: Vec<Layout>,
    /// `listOfLayouts` itself: its attributes, notes, annotation and other
    /// children, as written.
    pub layouts_meta: Meta,
    /// The Render package's global render information, which belongs to
    /// `listOfLayouts` (in Level 2, to its annotation).
    pub global_render: Option<crate::render::RenderList>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct FunctionDefinition {
    pub id: String,
    pub name: Option<String>,
    /// A `Math::Lambda`.
    pub math: Option<Math>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UnitDefinition {
    pub id: String,
    pub name: Option<String>,
    pub units: Vec<Unit>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Unit {
    /// A base unit, for example `mole` or `second`.
    pub kind: String,
    /// An integer in Level 2.
    pub exponent: Option<f64>,
    pub scale: Option<i32>,
    pub multiplier: Option<f64>,
    /// Level 2 Version 1 only.
    pub offset: Option<f64>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct CompartmentType {
    pub id: String,
    pub name: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpeciesType {
    pub id: String,
    pub name: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Compartment {
    pub id: String,
    pub name: Option<String>,
    /// An integer from 0 to 3 in Level 2.
    pub spatial_dimensions: Option<f64>,
    pub size: Option<f64>,
    pub units: Option<String>,
    /// Level 2 only.
    pub outside: Option<String>,
    pub constant: Option<bool>,
    /// Level 2 Versions 2 to 4 only.
    pub compartment_type: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Species {
    pub id: String,
    pub name: Option<String>,
    pub compartment: String,
    pub initial_amount: Option<f64>,
    pub initial_concentration: Option<f64>,
    pub substance_units: Option<String>,
    /// Level 2 Versions 1 and 2 only.
    pub spatial_size_units: Option<String>,
    pub has_only_substance_units: Option<bool>,
    pub boundary_condition: Option<bool>,
    /// Level 2 only.
    pub charge: Option<i32>,
    pub constant: Option<bool>,
    /// Level 3 only.
    pub conversion_factor: Option<String>,
    /// Level 2 Versions 2 to 4 only.
    pub species_type: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Parameter {
    pub id: String,
    pub name: Option<String>,
    pub value: Option<f64>,
    pub units: Option<String>,
    pub constant: Option<bool>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct InitialAssignment {
    pub symbol: String,
    pub math: Option<Math>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Rule {
    Assignment {
        variable: String,
        math: Option<Math>,
        meta: Meta,
    },
    Rate {
        variable: String,
        math: Option<Math>,
        meta: Meta,
    },
    Algebraic {
        math: Option<Math>,
        meta: Meta,
    },
}

impl Rule {
    /// The variable an assignment or rate rule sets.
    pub fn variable(&self) -> Option<&str> {
        match self {
            Rule::Assignment { variable, .. } | Rule::Rate { variable, .. } => Some(variable),
            Rule::Algebraic { .. } => None,
        }
    }

    pub fn math(&self) -> Option<&Math> {
        match self {
            Rule::Assignment { math, .. }
            | Rule::Rate { math, .. }
            | Rule::Algebraic { math, .. } => math.as_ref(),
        }
    }

    pub fn meta(&self) -> &Meta {
        match self {
            Rule::Assignment { meta, .. }
            | Rule::Rate { meta, .. }
            | Rule::Algebraic { meta, .. } => meta,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Constraint {
    pub math: Option<Math>,
    /// The content of `<message>`, as raw XML.
    pub message: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Reaction {
    pub id: String,
    pub name: Option<String>,
    pub reversible: Option<bool>,
    pub fast: Option<bool>,
    /// Level 3 only.
    pub compartment: Option<String>,
    pub reactants: Vec<SpeciesReference>,
    pub products: Vec<SpeciesReference>,
    pub modifiers: Vec<ModifierSpeciesReference>,
    pub kinetic_law: Option<KineticLaw>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpeciesReference {
    pub id: Option<String>,
    pub name: Option<String>,
    pub species: String,
    pub stoichiometry: Option<f64>,
    /// Level 2 only: stoichiometry given as MathML.
    pub stoichiometry_math: Option<Math>,
    /// Level 3 only.
    pub constant: Option<bool>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ModifierSpeciesReference {
    pub id: Option<String>,
    pub name: Option<String>,
    pub species: String,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct KineticLaw {
    pub math: Option<Math>,
    /// Level 3's `listOfLocalParameters`, or Level 2's `listOfParameters`.
    pub parameters: Vec<LocalParameter>,
    /// Level 2 Version 1 only.
    pub time_units: Option<String>,
    /// Level 2 Version 1 only.
    pub substance_units: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct LocalParameter {
    pub id: String,
    pub name: Option<String>,
    pub value: Option<f64>,
    pub units: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Event {
    pub id: Option<String>,
    pub name: Option<String>,
    pub use_values_from_trigger_time: Option<bool>,
    /// Level 2 Versions 1 and 2 only.
    pub time_units: Option<String>,
    pub trigger: Option<Trigger>,
    pub delay: Option<Delay>,
    /// Level 3 only.
    pub priority: Option<Priority>,
    pub assignments: Vec<EventAssignment>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Trigger {
    /// Level 3 only.
    pub initial_value: Option<bool>,
    /// Level 3 only.
    pub persistent: Option<bool>,
    pub math: Option<Math>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Delay {
    pub math: Option<Math>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Priority {
    pub math: Option<Math>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct EventAssignment {
    pub variable: String,
    pub math: Option<Math>,
    pub meta: Meta,
}

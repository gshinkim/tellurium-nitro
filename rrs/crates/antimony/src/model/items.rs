//! The model's items: symbols, reactions, events and the rest, and the
//! modules and document that hold them.

use super::{Constness, InteractionKind, Qualifier, ReactionKind, SymbolKind};
use crate::expr::Expr;
use crate::span::Span;

// --- annotations -------------------------------------------------------------

/// One CV term: `S1 identity "http://identifiers.org/chebi/CHEBI:4167"`.
#[derive(Clone, Debug, PartialEq)]
pub struct CvTerm {
    pub qualifier: Qualifier,
    pub resource: String,
    /// The word the source used (`part` or `hasPart`), so the file is written
    /// back as it was. The qualifier is the meaning; the spelling is the
    /// author's.
    pub spelling: String,
    pub span: Span,
}

/// One `creator1.*` group. The number is only a grouping key; 0 is the
/// unnumbered `creator.name` form.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Creator {
    pub index: u32,
    pub name: Option<String>,
    pub given_name: Option<String>,
    pub family_name: Option<String>,
    pub organization: Option<String>,
    pub email: Option<String>,
}

impl Creator {
    pub fn is_empty(&self) -> bool {
        self.name.is_none()
            && self.given_name.is_none()
            && self.family_name.is_none()
            && self.organization.is_none()
            && self.email.is_none()
    }
}

/// The Dublin Core side of an annotation: dates and people. Dates are kept
/// as written, not parsed: an annotation records what the file says.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Metadata {
    pub created: Option<String>,
    pub modified: Vec<String>,
    pub creators: Vec<Creator>,
    pub notes: Option<String>,
}

impl Metadata {
    pub fn is_empty(&self) -> bool {
        self.created.is_none()
            && self.notes.is_none()
            && self.modified.is_empty()
            && self.creators.is_empty()
    }

    /// The creator group with this number, made if there is none yet.
    pub fn creator_mut(&mut self, index: u32) -> &mut Creator {
        let i = match self.creators.iter().position(|c| c.index == index) {
            Some(i) => i,
            None => {
                self.creators.push(Creator {
                    index,
                    ..Creator::default()
                });
                self.creators.len() - 1
            }
        };
        &mut self.creators[i]
    }
}

// --- symbols -----------------------------------------------------------------

/// A rule on a symbol: `x := expr` or `x' = expr`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RuleKind {
    /// `:=`, recomputed continuously.
    Assignment,
    /// `'`, the derivative with respect to time.
    Rate,
}

impl RuleKind {
    pub fn name(self) -> &'static str {
        match self {
            RuleKind::Assignment => "assignment rule",
            RuleKind::Rate => "rate rule",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Rule {
    pub kind: RuleKind,
    /// `None` when the right-hand side was written empty: `p' = ;`.
    pub expr: Option<Expr>,
}

/// The unified symbol (Delphi §4.1): every named thing in a module, whose
/// kind is inferred and promoted as the file is read.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Symbol {
    pub name: String,
    pub kind: SymbolKind,
    pub constness: Constness,
    /// `$S1`, or a `const` species.
    pub is_boundary: bool,
    pub is_substance_only: bool,
    /// `S1 in cell`.
    pub compartment: Option<String>,
    /// `x = expr`: a value, or an initial assignment if not constant.
    pub initial_value: Option<Expr>,
    /// `p = ;` was written. Kept beside `initial_value`, which it never
    /// clears: an empty right-hand side must not destroy a declared value
    /// (Delphi §9.1.4 #3).
    pub has_empty_value: bool,
    pub rule: Option<Rule>,
    /// `S1 has substance_per_volume`.
    pub units: Option<String>,
    /// `S1 conversionFactor cf` (Delphi §5.8).
    pub conversion_factor: Option<String>,
    /// `S1 is "Glucose"`.
    pub display_name: Option<String>,
    pub notes: Option<String>,
    pub sbo_term: Option<u32>,
    pub cv_terms: Vec<CvTerm>,
    pub metadata: Metadata,
    /// The comment above the statement that introduced the symbol.
    pub comment: Option<String>,
    /// Never declared, only referred to: the typo case (Delphi §6.5). Set by
    /// the semantic pass.
    pub is_implicit: bool,
    /// The name before flattening (Delphi §2.2.2 seam 3). The same as `name`
    /// until modules are built.
    pub source_path: String,
    /// Where it was first mentioned.
    pub span: Span,
    /// Where it was first introduced; the same as `span` when parsed. The
    /// semantic pass reports conflicts against it.
    pub decl_span: Span,
}

impl Symbol {
    /// A new, unclassified symbol.
    pub fn new(name: impl Into<String>, span: Span) -> Symbol {
        let name = name.into();
        Symbol {
            source_path: name.clone(),
            name,
            span,
            decl_span: span,
            ..Symbol::default()
        }
    }

    /// Has a value, a rule, or both: Delphi §6.5's test for whether it is
    /// uninitialised.
    pub fn is_initialised(&self) -> bool {
        self.initial_value.is_some() || self.rule.is_some()
    }
}

// --- reactions ---------------------------------------------------------------

/// A species in a reaction, with its stoichiometry.
#[derive(Clone, Debug, PartialEq)]
pub struct Participant {
    pub species: String,
    /// An expression, not a number: `k S1 -> S2` is legal, and named
    /// stoichiometries (`S1 => S1_S2_stoichiometry S2`) are in the corpus.
    /// `None` means 1.
    pub stoichiometry: Option<Expr>,
    /// `$S1` written in the reaction.
    pub is_boundary: bool,
    pub span: Span,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Reaction {
    /// Always set: an unlabelled reaction is named `_J1`, `_J2`, … when it
    /// is read, by its position among all the reactions.
    pub name: String,
    /// The name was made up when the file was read, not written in it. The
    /// writer then leaves the label out, so the file comes back as written
    /// and reading it again makes the same name (SPEC.md R28).
    pub name_generated: bool,
    pub kind: ReactionKind,
    pub reactants: Vec<Participant>,
    pub products: Vec<Participant>,
    pub modifiers: Vec<Participant>,
    pub kinetic_law: Option<Expr>,
    /// `R16: J => M; ;`: a law written and empty, which is not the same as
    /// no law at all.
    pub has_empty_kinetic_law: bool,
    pub compartment: Option<String>,
    pub comment: Option<String>,
    pub span: Span,
}

impl Reaction {
    pub fn new(name: impl Into<String>, span: Span) -> Reaction {
        Reaction {
            name: name.into(),
            name_generated: false,
            kind: ReactionKind::Irreversible,
            reactants: Vec::new(),
            products: Vec::new(),
            modifiers: Vec::new(),
            kinetic_law: None,
            has_empty_kinetic_law: false,
            compartment: None,
            comment: None,
            span,
        }
    }

    /// `->` and `<->` are reversible; `=>` is not (Delphi §9.1.1).
    pub fn is_reversible(&self) -> bool {
        matches!(
            self.kind,
            ReactionKind::Reversible | ReactionKind::BiDirectional
        )
    }
}

/// `i1: G1 -| J1;`: a species affecting a reaction without being consumed
/// or produced by it. Held on the module, not the reaction, because it may
/// come before the reaction it names.
#[derive(Clone, Debug, PartialEq)]
pub struct Interaction {
    /// `None` when unlabelled. No name is invented: nothing refers to one.
    pub name: Option<String>,
    pub species: String,
    pub reaction: String,
    pub kind: InteractionKind,
    pub comment: Option<String>,
    pub span: Span,
}

// --- units -------------------------------------------------------------------

/// `unit substance = 1e-6 mole;`. Not a symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct UnitDef {
    pub name: String,
    /// The right-hand side as written. A unit expression juxtaposes
    /// (`mole litre / second`), which the expression grammar cannot say, so
    /// it is kept as text.
    pub definition: String,
    /// `time_unit is "time"`.
    pub display_name: Option<String>,
    pub cv_terms: Vec<CvTerm>,
    pub metadata: Metadata,
    pub comment: Option<String>,
    pub span: Span,
}

impl UnitDef {
    /// `unit time_unit = time_unit;`, a libantimony quirk. A unit resolver
    /// must not recurse on it.
    pub fn is_self_referential(&self) -> bool {
        self.definition.trim() == self.name
    }
}

// --- events ------------------------------------------------------------------

/// `S1 = 0` inside an event.
#[derive(Clone, Debug, PartialEq)]
pub struct EventAssignment {
    pub target: String,
    pub expr: Option<Expr>,
    pub span: Span,
}

/// `E1: at (x > 5), priority = 1: y = 3, x = r + 2`, or with a delay,
/// `at 1 after S1 < 0.1: S1 = 1`.
#[derive(Clone, Debug, PartialEq)]
pub struct Event {
    pub name: Option<String>,
    pub trigger: Expr,
    pub delay: Option<Expr>,
    pub priority: Option<Expr>,
    pub persistent: Option<Expr>,
    /// The `t0` modifier: whether a trigger already true at the start counts
    /// as having just become true.
    pub initial_value: Option<Expr>,
    pub from_trigger: Option<Expr>,
    /// May be empty: `E0: at time > 5.5: ;`.
    pub assignments: Vec<EventAssignment>,
    pub comment: Option<String>,
    pub span: Span,
}

// --- functions, rules, constraints, objectives --------------------------------

/// `function MM(S, Vm, Km) Vm*S/(Km + S) end`.
#[derive(Clone, Debug, PartialEq)]
pub struct Function {
    pub name: String,
    pub args: Vec<String>,
    pub body: Option<Expr>,
    pub comment: Option<String>,
    pub span: Span,
}

/// `_alg0: 0 = -1*k1 + S1 + S2;`. The name is a label, not a symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct AlgebraicRule {
    pub name: String,
    /// `None` when written empty: `_alg0: 0 = ;` (case 01244).
    pub expr: Option<Expr>,
    pub comment: Option<String>,
    pub span: Span,
}

/// What a constraint says.
#[derive(Clone, Debug, PartialEq)]
pub enum ConstraintKind {
    /// `constraint c0: R16 >= 0`.
    Expr(Expr),
    /// `-10 <= J2 <= 10`: a chained comparison is a range, here and only
    /// here (Delphi decision 18).
    Range {
        lower: Expr,
        subject: Expr,
        upper: Expr,
        lower_inclusive: bool,
        upper_inclusive: bool,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct Constraint {
    /// `None` for a bare flux bound such as `R19 <= 1000`.
    pub name: Option<String>,
    pub kind: ConstraintKind,
    pub message: Option<String>,
    pub comment: Option<String>,
    pub span: Span,
}

impl Constraint {
    pub fn is_range(&self) -> bool {
        matches!(self.kind, ConstraintKind::Range { .. })
    }
}

/// `OBJF: maximize R26 + R10;`, the fbc package's objective. The name is a
/// label, not a symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct Objective {
    pub name: Option<String>,
    pub maximise: bool,
    pub expr: Expr,
    pub comment: Option<String>,
    pub span: Span,
}

/// What a raw section is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RawKind {
    /// An annotation in a vocabulary this library does not know.
    Annotation,
}

/// One of Antimony's layout statements (Antimony 3, the reference's "Layout
/// and Render"): `model.layout = on`, `S1.pos = {28, 35}`,
/// `J0.S1.arc2.b1 = {55, 31}`, `J0.--.position = {520, 483}`,
/// `species.color = azure`. This library reads and keeps them; what they
/// draw is for the caller to decide.
#[derive(Clone, Debug, PartialEq)]
pub struct LayoutStatement {
    /// The dotted parts as written: `["J0", "S1", "arc2", "b1"]`; a null
    /// node's part is `"--"`.
    pub path: Vec<String>,
    pub value: LayoutValue,
    pub comment: Option<String>,
    pub span: Span,
}

/// The right-hand side of a [`LayoutStatement`].
#[derive(Clone, Debug, PartialEq)]
pub enum LayoutValue {
    Number(f64),
    /// `{x, y}`.
    Point(f64, f64),
    /// Any other braced list, such as `{S1, S2, S3}`, each item as written.
    List(Vec<String>),
    /// One or more words: `on`, `azure`, `bold_italic`, `blue ombre`.
    Word(String),
    /// A quoted string: `"#268F49CC"`.
    Text(String),
}

/// Text kept verbatim so that a load and save loses nothing.
#[derive(Clone, Debug, PartialEq)]
pub struct RawSection {
    pub kind: RawKind,
    pub text: String,
    pub span: Span,
}

// --- modules and the document ------------------------------------------------

/// A `model … end` block, or the statements of a file with none.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Module {
    /// `None` for the statements of a file without a `model` block.
    pub name: Option<String>,
    /// `case00001 is "case00001"`: the module is the subject, not a symbol.
    pub display_name: Option<String>,
    /// `model *m()`.
    pub is_main: bool,
    pub formal_args: Vec<String>,
    pub symbols: Vec<Symbol>,
    pub reactions: Vec<Reaction>,
    pub interactions: Vec<Interaction>,
    pub units: Vec<UnitDef>,
    pub objectives: Vec<Objective>,
    pub events: Vec<Event>,
    pub functions: Vec<Function>,
    pub algebraic_rules: Vec<AlgebraicRule>,
    pub constraints: Vec<Constraint>,
    pub sbo_term: Option<u32>,
    /// Unknown annotations, kept verbatim.
    pub unparsed_sections: Vec<RawSection>,
    /// Layout statements, in the order read.
    pub layout: Vec<LayoutStatement>,
    /// `model conversionFactor cf` (Delphi §5.8).
    pub conversion_factor: Option<String>,
    pub notes: Option<String>,
    pub cv_terms: Vec<CvTerm>,
    pub metadata: Metadata,
    pub span: Span,
}

impl Module {
    pub fn new(name: Option<String>, span: Span) -> Module {
        Module {
            name,
            span,
            ..Module::default()
        }
    }

    pub fn symbol(&self, name: &str) -> Option<&Symbol> {
        self.symbols.iter().find(|s| s.name == name)
    }

    pub fn symbol_mut(&mut self, name: &str) -> Option<&mut Symbol> {
        self.symbols.iter_mut().find(|s| s.name == name)
    }

    pub fn reaction(&self, name: &str) -> Option<&Reaction> {
        self.reactions.iter().find(|r| r.name == name)
    }

    pub fn unit(&self, name: &str) -> Option<&UnitDef> {
        self.units.iter().find(|u| u.name == name)
    }

    pub fn unit_mut(&mut self, name: &str) -> Option<&mut UnitDef> {
        self.units.iter_mut().find(|u| u.name == name)
    }

    pub fn function(&self, name: &str) -> Option<&Function> {
        self.functions.iter().find(|f| f.name == name)
    }

    /// The symbols of one kind, in order.
    pub fn symbols_of(&self, kind: SymbolKind) -> impl Iterator<Item = &Symbol> {
        self.symbols.iter().filter(move |s| s.kind == kind)
    }
}

/// A parsed file: its modules, one of which is the main one.
///
/// A list from the start, although it holds one module until modules are
/// built (Delphi §2.2.2 seam 1): nothing may assume there is only one.
#[derive(Clone, Debug, PartialEq)]
pub struct Document {
    pub modules: Vec<Module>,
    /// Index into `modules`: the `*` module, or else the first.
    pub main: usize,
}

impl Default for Document {
    /// One unnamed main module.
    fn default() -> Document {
        let mut module = Module::new(None, Span::default());
        module.is_main = true;
        Document {
            modules: vec![module],
            main: 0,
        }
    }
}

impl Document {
    pub fn main_module(&self) -> &Module {
        &self.modules[self.main]
    }

    pub fn main_module_mut(&mut self) -> &mut Module {
        &mut self.modules[self.main]
    }

    pub fn module(&self, name: &str) -> Option<&Module> {
        self.modules
            .iter()
            .find(|m| m.name.as_deref() == Some(name))
    }

    /// The index of the module with this name.
    pub fn module_index(&self, name: &str) -> Option<usize> {
        self.modules
            .iter()
            .position(|m| m.name.as_deref() == Some(name))
    }
}

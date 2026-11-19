//! The model (Delphi §4, `uAntModel`).
//!
//! Plain data: structs with public fields, `Vec` in source order, `Option`
//! for what may be absent (SPEC.md §3). Items refer to each other by name,
//! never by pointer or index, so editing the model can never leave a
//! reference dangling (SPEC.md §3.1); the `symbol`, `reaction`, `unit` and
//! `function` methods look a name up.
//!
//! Nothing here is resolved: a name is what the file wrote. Binding names
//! and checking them is the semantic pass's job.

mod items;
mod values;

pub use items::*;

/// What a symbol is. Antimony infers this as a file is read, and promotes it
/// as more is learned (Delphi §4.1, §6.2).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum SymbolKind {
    /// Seen, not yet classified.
    #[default]
    Unknown,
    Species,
    Parameter,
    Compartment,
    Reaction,
    Event,
    Function,
    Submodule,
    UnitDef,
    Constraint,
    Interaction,
    Deleted,
    /// A symbol standing as a stoichiometric coefficient:
    /// `reaction1: S1 => S1_S2_stoichiometry S2;` (Delphi decision 17). Not
    /// a parameter: libantimony does not treat it as one, and SBML needs it
    /// as the id of a `<speciesReference>`, not as a `<parameter>`.
    Stoichiometry,
}

impl SymbolKind {
    pub const ALL: [SymbolKind; 13] = [
        SymbolKind::Unknown,
        SymbolKind::Species,
        SymbolKind::Parameter,
        SymbolKind::Compartment,
        SymbolKind::Reaction,
        SymbolKind::Event,
        SymbolKind::Function,
        SymbolKind::Submodule,
        SymbolKind::UnitDef,
        SymbolKind::Constraint,
        SymbolKind::Interaction,
        SymbolKind::Deleted,
        SymbolKind::Stoichiometry,
    ];

    /// The word used in messages: `species`, `unit`, …
    pub fn name(self) -> &'static str {
        match self {
            SymbolKind::Unknown => "unknown",
            SymbolKind::Species => "species",
            SymbolKind::Parameter => "parameter",
            SymbolKind::Compartment => "compartment",
            SymbolKind::Reaction => "reaction",
            SymbolKind::Event => "event",
            SymbolKind::Function => "function",
            SymbolKind::Submodule => "submodule",
            SymbolKind::UnitDef => "unit",
            SymbolKind::Constraint => "constraint",
            SymbolKind::Interaction => "interaction",
            SymbolKind::Deleted => "deleted",
            SymbolKind::Stoichiometry => "stoichiometry",
        }
    }

    /// Whether a symbol of this kind may become `to` (Delphi §6.2).
    ///
    /// An unknown symbol may become anything. A parameter may become a
    /// species, a compartment or a stoichiometry, because `k1 = 0.1` looks
    /// the same as an initial value until `k1` turns up in a reaction or an
    /// `in` clause. Nothing else moves; staying the same is always allowed.
    pub fn can_promote_to(self, to: SymbolKind) -> bool {
        self == to
            || self == SymbolKind::Unknown
            || (self == SymbolKind::Parameter
                && matches!(
                    to,
                    SymbolKind::Species | SymbolKind::Compartment | SymbolKind::Stoichiometry
                ))
    }
}

/// A `const` or `var` qualifier, or neither.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Constness {
    #[default]
    Default,
    Const,
    Var,
}

impl Constness {
    pub fn name(self) -> &'static str {
        match self {
            Constness::Default => "default",
            Constness::Const => "const",
            Constness::Var => "var",
        }
    }
}

/// The arrow of a reaction.
///
/// **`->` is reversible and `=>` irreversible**, the reverse of what one
/// might guess. Verified against the corpus by checking each case's SBML
/// `reversible` attribute (Delphi §9.1.1).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ReactionKind {
    /// `->`
    #[default]
    Reversible,
    /// `=>`
    Irreversible,
    /// `<->`
    BiDirectional,
    /// `-|`
    Inhibition,
    /// `-o`
    Activation,
}

/// How an interaction's species affects its reaction (Delphi §9.2):
/// `i1: G1 -| J1;` inhibits, `G2 -o J1` activates, `G3 -( J1` modifies in an
/// unspecified way.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum InteractionKind {
    #[default]
    Unknown,
    Inhibits,
    Activates,
}

impl InteractionKind {
    /// The Antimony arrow: `-|`, `-o` or `-(`.
    pub fn arrow(self) -> &'static str {
        match self {
            InteractionKind::Inhibits => "-|",
            InteractionKind::Activates => "-o",
            InteractionKind::Unknown => "-(",
        }
    }

    /// The SBO term SBML gives the modifier, as libantimony writes it:
    /// 20 (inhibitor), 459 (stimulator), or 19 (modifier), which is also the
    /// honest term for "we were not told which".
    pub fn sbo_term(self) -> u32 {
        match self {
            InteractionKind::Inhibits => 20,
            InteractionKind::Activates => 459,
            InteractionKind::Unknown => 19,
        }
    }
}

/// An annotation qualifier (Delphi M9), named after the SBML qualifier.
///
/// Each has two accepted Antimony spellings (`part` and `hasPart` mean the
/// same), and the model records which one the source used. The last two are
/// model qualifiers, legal only on a module. `IsDescribedBy` is decided by
/// its subject: on a species it is `bqbiol:isDescribedBy`, on a model
/// `bqmodel:isDescribedBy`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Qualifier {
    Is,
    HasPart,
    IsPartOf,
    IsVersionOf,
    HasVersion,
    IsHomologTo,
    IsDescribedBy,
    IsEncodedBy,
    Encodes,
    OccursIn,
    HasProperty,
    IsPropertyOf,
    HasTaxon,
    ModelIs,
    IsDerivedFrom,
}

/// short spelling, long spelling, SBML element, model only
type QualifierRow = (&'static str, &'static str, &'static str, bool);

impl Qualifier {
    pub const ALL: [Qualifier; 15] = [
        Qualifier::Is,
        Qualifier::HasPart,
        Qualifier::IsPartOf,
        Qualifier::IsVersionOf,
        Qualifier::HasVersion,
        Qualifier::IsHomologTo,
        Qualifier::IsDescribedBy,
        Qualifier::IsEncodedBy,
        Qualifier::Encodes,
        Qualifier::OccursIn,
        Qualifier::HasProperty,
        Qualifier::IsPropertyOf,
        Qualifier::HasTaxon,
        Qualifier::ModelIs,
        Qualifier::IsDerivedFrom,
    ];

    /// `identity`'s long form is `biological_entity_is`, neither the SBML
    /// name nor a shortening; the vocabulary simply has three words for one
    /// idea. The element is what SBML gets.
    fn row(self) -> QualifierRow {
        match self {
            Qualifier::Is => ("identity", "biological_entity_is", "is", false),
            Qualifier::HasPart => ("part", "hasPart", "hasPart", false),
            Qualifier::IsPartOf => ("parthood", "isPartOf", "isPartOf", false),
            Qualifier::IsVersionOf => ("hypernym", "isVersionOf", "isVersionOf", false),
            Qualifier::HasVersion => ("version", "hasVersion", "hasVersion", false),
            Qualifier::IsHomologTo => ("homolog", "isHomologTo", "isHomologTo", false),
            Qualifier::IsDescribedBy => ("description", "isDescribedBy", "isDescribedBy", false),
            Qualifier::IsEncodedBy => ("encoder", "isEncodedBy", "isEncodedBy", false),
            Qualifier::Encodes => ("encodement", "encodes", "encodes", false),
            Qualifier::OccursIn => ("container", "occursIn", "occursIn", false),
            Qualifier::HasProperty => ("property", "hasProperty", "hasProperty", false),
            Qualifier::IsPropertyOf => ("propertyBearer", "isPropertyOf", "isPropertyOf", false),
            Qualifier::HasTaxon => ("taxon", "hasTaxon", "hasTaxon", false),
            Qualifier::ModelIs => ("model_entity_is", "model_entity_is", "is", true),
            Qualifier::IsDerivedFrom => ("origin", "isDerivedFrom", "isDerivedFrom", true),
        }
    }

    /// The qualifier either spelling names, case-sensitively. `None` means
    /// the word is not a qualifier, which is not an error: the parser keeps
    /// such an annotation as raw text, so a file using a vocabulary we do not
    /// know still loads and round-trips.
    pub fn from_word(word: &str) -> Option<Qualifier> {
        Qualifier::ALL.into_iter().find(|q| {
            let (short, long, _, _) = q.row();
            short == word || long == word
        })
    }

    /// The SBML local name: `is`, `hasPart`, …
    pub fn element(self) -> &'static str {
        self.row().2
    }

    /// Whether it may appear only on a module.
    pub fn is_model_only(self) -> bool {
        self.row().3
    }

    /// The spelling written when none was recorded: the shorter, more
    /// Antimony-like synonym, which is what the corpus uses.
    pub fn keyword(self) -> &'static str {
        self.row().0
    }
}

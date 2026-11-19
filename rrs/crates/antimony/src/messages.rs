//! The text of every message (Delphi `uAntStrings`).
//!
//! Kept in one place so the messages can be reviewed as a set. The wording
//! is the Delphi version's, so the two implementations can be compared on
//! their diagnostics too.

// --- lexical --------------------------------------------------------------

/// `Unexpected character '~'`.
pub fn unexpected_char(c: char) -> String {
    if c == '\'' {
        "Unexpected character ''''".into()
    } else {
        format!("Unexpected character '{c}'")
    }
}

pub const UNTERMINATED_STRING: &str = "Unterminated string literal";
pub const UNTERMINATED_COMMENT: &str = "Unterminated /* comment";

pub fn invalid_number(text: &str) -> String {
    format!("Invalid number literal \"{text}\"")
}

// --- syntax ---------------------------------------------------------------

pub fn syntax_error(detail: &str) -> String {
    format!("Syntax error: {detail}")
}

pub fn unexpected_token(found: &str, expected: &str) -> String {
    format!("Unexpected {found}, expected {expected}")
}

pub fn expected_identifier(found: &str) -> String {
    format!("Expected an identifier, found {found}")
}

pub fn expected_end(block: &str, starting_at: &str) -> String {
    format!("Missing \"end\" for the {block} starting at {starting_at}")
}

pub fn negative_stoichiometry(species: &str) -> String {
    format!("Negative stoichiometry on \"{species}\" is not supported")
}

pub fn empty_reaction(reaction: &str) -> String {
    format!("Reaction \"{reaction}\" has neither reactants nor products")
}

pub const DNA_NOT_SUPPORTED: &str =
    "DNA strand syntax (\"--\") is not supported; use libantimony for this file";

pub fn not_yet_supported(what: &str, milestone: &str) -> String {
    format!("{what} is not supported yet (planned for milestone {milestone})")
}

// --- semantic -------------------------------------------------------------

pub fn type_conflict(name: &str, wanted: &str, established: &str, at: &str) -> String {
    format!(
        "Cannot use \"{name}\" as a {wanted}; it was already established as a {established} at {at}"
    )
}

pub fn undefined_symbol(name: &str) -> String {
    format!("Undefined symbol \"{name}\"")
}

pub fn uninitialised_symbol(name: &str) -> String {
    format!("Symbol \"{name}\" has no value, assignment rule or rate rule")
}

/// The same condition as [`uninitialised_symbol`], when the caller asked for
/// missing values to be filled in. The other wording would be false of the
/// output: by the time anything is written, the symbol has a value.
pub fn initialised_to_default(name: &str, value: &str) -> String {
    format!("Symbol \"{name}\" had no value; initialised to {value}")
}

pub fn implicit_parameter(name: &str) -> String {
    format!(
        "\"{name}\" was never declared; it has been assumed to be a parameter. Check for a typo"
    )
}

pub fn duplicate_rule(name: &str, existing: &str, second: &str) -> String {
    format!("Symbol \"{name}\" already has {existing}; a second rule ({second}) is not allowed")
}

pub fn const_with_rule(name: &str, rule: &str) -> String {
    format!("Symbol \"{name}\" is const and cannot carry a {rule}")
}

pub fn not_a_compartment(name: &str, kind: &str) -> String {
    format!("\"{name}\" is used as a compartment but is a {kind}")
}

pub fn bad_conversion_factor(name: &str, kind: &str) -> String {
    format!("conversionFactor \"{name}\" must name a parameter, but it is a {kind}")
}

pub fn arity_mismatch(function: &str, expected: usize, supplied: usize) -> String {
    format!("Function \"{function}\" takes {expected} argument(s), but {supplied} were supplied")
}

pub fn recursive_function(function: &str) -> String {
    format!("Function \"{function}\" is defined recursively")
}

pub fn duplicate_formal_arg(name: &str) -> String {
    format!("Duplicate formal argument \"{name}\"")
}

pub fn circular_rule(chain: &str) -> String {
    format!("Circular assignment-rule dependency: {chain}")
}

pub fn undefined_module(name: &str) -> String {
    format!("Undefined module \"{name}\"")
}

pub fn circular_module(chain: &str) -> String {
    format!("Circular module instantiation: {chain}")
}

pub fn incompatible_mapping(from: &str, from_kind: &str, to: &str, to_kind: &str) -> String {
    format!("Cannot map \"{from}\" ({from_kind}) onto \"{to}\" ({to_kind})")
}

pub fn information_loss(what: &str) -> String {
    format!("Information would be lost: {what}")
}

pub fn model_only_qualifier(qualifier: &str) -> String {
    format!(
        "\"{qualifier}\" is a model qualifier and means nothing on an element; it is kept but will not be written to SBML"
    )
}

pub fn annotation_dropped(what: &str) -> String {
    format!("annotations on {what} are not written to SBML")
}

// --- infrastructure -------------------------------------------------------

pub fn file_not_found(name: &str) -> String {
    format!("File not found: {name}")
}

pub fn read_error(name: &str, reason: &str) -> String {
    format!("Could not read \"{name}\": {reason}")
}

/// Not in Delphi, which read a file whole and reported its decoder's
/// message with `read_error`.
pub const NOT_UTF8: &str = "The source is not valid UTF-8";

pub fn internal_error(what: &str) -> String {
    format!("Internal error: {what}")
}

// --- reporting ------------------------------------------------------------

pub fn error_summary(errors: usize, warnings: usize) -> String {
    format!("{errors} error(s), {warnings} warning(s)")
}

pub const NO_DIAGNOSTICS: &str = "No diagnostics";

//! Built-in function names (Delphi `uAntTypes`).
//!
//! A name on these lists, in call position, is a function the language
//! provides and is never turned into an implicit parameter (Delphi §6.5).
//! Getting the lists wrong is not cosmetic: a missing entry shows up as a
//! spurious `dcImplicitParameter` warning on a valid model.

/// The function list from the Antimony language documentation, verbatim and
/// in its own order. Delphi §5.3's list is an abridged paraphrase of it.
pub const DOCUMENTED_FUNCTIONS: [&str; 48] = [
    "abs",
    "and",
    "arccos",
    "arccosh",
    "arccot",
    "arccoth",
    "arccsc",
    "arccsch",
    "arcsec",
    "arcsech",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctanh",
    "ceiling",
    "cos",
    "cosh",
    "cot",
    "coth",
    "csc",
    "csch",
    "divide",
    "eq",
    "exp",
    "factorial",
    "floor",
    "geq",
    "gt",
    "leq",
    "ln",
    "log",
    "lt",
    "minus",
    "neq",
    "not",
    "or",
    "piecewise",
    "plus",
    "power",
    "root",
    "sec",
    "sech",
    "sin",
    "sinh",
    "tan",
    "tanh",
    "times",
    "xor",
];

/// Names outside the documented list that must nonetheless be accepted. The
/// numbers are occurrences in the 1788-case corpus.
pub const ADDITIONAL_FUNCTIONS: [&str; 16] = [
    // libantimony's own spellings of documented functions.
    "asin",  // 36, documented as arcsin
    "acos",  // 36, documented as arccos
    "atan",  // 36, documented as arctan
    "ceil",  // 57, documented as ceiling
    "sqrt",  // 34, documented only as root(2, x)
    "log10", // 34, documented only as log
    "pow",   // Delphi §5.3's spelling of power
    // csymbols: not functions, but they must never become parameters.
    "delay",    // 67
    "time",     // 732, always bare
    "avogadro", // 23
    "rateOf",   // 86
    // MathML operators the corpus uses that neither list names.
    "implies",  // 12
    "max",      // 6
    "min",      // 6
    "quotient", // 3
    "rem",      // 3
];

/// Whether `name` is a built-in function. Case-sensitive, as Antimony is.
pub fn is_builtin_function(name: &str) -> bool {
    DOCUMENTED_FUNCTIONS.contains(&name) || ADDITIONAL_FUNCTIONS.contains(&name)
}

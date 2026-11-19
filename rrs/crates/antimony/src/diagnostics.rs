//! Problems found in a model file (SPEC.md §5, Delphi `uAntDiagnostics`).
//!
//! Reading a file never fails and never panics on what is in it: the lexer,
//! the parser and the semantic passes add a [`Diagnostic`] and keep going, so
//! an editor can show every problem in a file at once.

use std::fmt;
use std::ops::Deref;

use crate::messages;
use crate::span::Span;

/// How serious a diagnostic is. Ordered: `Hint < Warning < Error`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Severity {
    Hint,
    Warning,
    Error,
}

impl Severity {
    /// `hint`, `warning` or `error`.
    pub fn name(self) -> &'static str {
        match self {
            Severity::Hint => "hint",
            Severity::Warning => "warning",
            Severity::Error => "error",
        }
    }
}

impl fmt::Display for Severity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Every diagnostic the library can produce. A host may match on these
/// rather than on message text: the text may change, the codes do not.
///
/// They are the Delphi `TAntDiagCode` values without the `dc` prefix;
/// [`Code::name`] gives the Delphi spelling, which the `tests/invalid`
/// fixtures use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Code {
    None,
    // Lexical.
    UnexpectedChar,
    UnterminatedString,
    UnterminatedComment,
    InvalidNumber,
    // Syntax.
    SyntaxError,
    UnexpectedToken,
    ExpectedIdentifier,
    ExpectedEnd,
    /// Delphi §5.7, decision 14.
    NegativeStoichiometry,
    EmptyReaction,
    /// DNA strand syntax, out of scope (Delphi N2).
    DnaNotSupported,
    /// A deferred feature; the message names its milestone.
    NotYetSupported,
    // Semantic.
    /// The promotion table of Delphi §6.2.
    TypeConflict,
    UndefinedSymbol,
    /// Declared or inferred, never given a value (Delphi §6.5).
    UninitialisedSymbol,
    /// Never declared at all: the typo case (Delphi §6.5).
    ImplicitParameter,
    DuplicateRule,
    ConstWithRule,
    NotACompartment,
    /// Delphi §5.8.
    BadConversionFactor,
    ArityMismatch,
    RecursiveFunction,
    DuplicateFormalArg,
    CircularRule,
    UndefinedModule,
    CircularModule,
    IncompatibleMapping,
    /// Delphi §9.1.4 #5: data is never dropped silently.
    InformationLoss,
    // Infrastructure.
    FileNotFound,
    ReadError,
    InternalError,
}

impl Code {
    pub const ALL: [Code; 32] = [
        Code::None,
        Code::UnexpectedChar,
        Code::UnterminatedString,
        Code::UnterminatedComment,
        Code::InvalidNumber,
        Code::SyntaxError,
        Code::UnexpectedToken,
        Code::ExpectedIdentifier,
        Code::ExpectedEnd,
        Code::NegativeStoichiometry,
        Code::EmptyReaction,
        Code::DnaNotSupported,
        Code::NotYetSupported,
        Code::TypeConflict,
        Code::UndefinedSymbol,
        Code::UninitialisedSymbol,
        Code::ImplicitParameter,
        Code::DuplicateRule,
        Code::ConstWithRule,
        Code::NotACompartment,
        Code::BadConversionFactor,
        Code::ArityMismatch,
        Code::RecursiveFunction,
        Code::DuplicateFormalArg,
        Code::CircularRule,
        Code::UndefinedModule,
        Code::CircularModule,
        Code::IncompatibleMapping,
        Code::InformationLoss,
        Code::FileNotFound,
        Code::ReadError,
        Code::InternalError,
    ];

    /// The Delphi name: `dcTypeConflict`.
    pub fn name(self) -> &'static str {
        match self {
            Code::None => "dcNone",
            Code::UnexpectedChar => "dcUnexpectedChar",
            Code::UnterminatedString => "dcUnterminatedString",
            Code::UnterminatedComment => "dcUnterminatedComment",
            Code::InvalidNumber => "dcInvalidNumber",
            Code::SyntaxError => "dcSyntaxError",
            Code::UnexpectedToken => "dcUnexpectedToken",
            Code::ExpectedIdentifier => "dcExpectedIdentifier",
            Code::ExpectedEnd => "dcExpectedEnd",
            Code::NegativeStoichiometry => "dcNegativeStoichiometry",
            Code::EmptyReaction => "dcEmptyReaction",
            Code::DnaNotSupported => "dcDNANotSupported",
            Code::NotYetSupported => "dcNotYetSupported",
            Code::TypeConflict => "dcTypeConflict",
            Code::UndefinedSymbol => "dcUndefinedSymbol",
            Code::UninitialisedSymbol => "dcUninitialisedSymbol",
            Code::ImplicitParameter => "dcImplicitParameter",
            Code::DuplicateRule => "dcDuplicateRule",
            Code::ConstWithRule => "dcConstWithRule",
            Code::NotACompartment => "dcNotACompartment",
            Code::BadConversionFactor => "dcBadConversionFactor",
            Code::ArityMismatch => "dcArityMismatch",
            Code::RecursiveFunction => "dcRecursiveFunction",
            Code::DuplicateFormalArg => "dcDuplicateFormalArg",
            Code::CircularRule => "dcCircularRule",
            Code::UndefinedModule => "dcUndefinedModule",
            Code::CircularModule => "dcCircularModule",
            Code::IncompatibleMapping => "dcIncompatibleMapping",
            Code::InformationLoss => "dcInformationLoss",
            Code::FileNotFound => "dcFileNotFound",
            Code::ReadError => "dcReadError",
            Code::InternalError => "dcInternalError",
        }
    }

    /// The code with this Delphi name.
    pub fn from_name(name: &str) -> Option<Code> {
        Code::ALL.into_iter().find(|code| code.name() == name)
    }

    /// The severity a diagnostic with this code has unless it is given one.
    /// The two uninitialised-symbol codes are warnings by design (Delphi
    /// §6.5): valid published models must load. Strict mode raises them.
    pub fn default_severity(self) -> Severity {
        match self {
            Code::None => Severity::Hint,
            Code::UninitialisedSymbol | Code::ImplicitParameter | Code::InformationLoss => {
                Severity::Warning
            }
            _ => Severity::Error,
        }
    }
}

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// One problem, where it is, and how serious.
#[derive(Clone, Debug, PartialEq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: Code,
    pub message: String,
    pub span: Span,
}

impl Diagnostic {
    /// A diagnostic with its code's default severity.
    pub fn new(code: Code, span: Span, message: impl Into<String>) -> Diagnostic {
        Diagnostic {
            severity: code.default_severity(),
            code,
            message: message.into(),
            span,
        }
    }

    /// The one-line form with a file name, as a console prints it and an
    /// editor parses it: `model.ant(3:17): error: message [dcSyntaxError]`.
    pub fn with_file(&self, file: &str) -> String {
        let place = match (self.span.is_unknown(), file.is_empty()) {
            (true, true) => String::new(),
            (true, false) => format!("{file}: "),
            (false, _) => format!("{file}({}): ", self.span.start),
        };
        format!("{place}{}: {} [{}]", self.severity, self.message, self.code)
    }
}

impl fmt::Display for Diagnostic {
    /// `(3:17): error: message [dcSyntaxError]`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.with_file(""))
    }
}

/// The diagnostics from reading one source, in the order they were found.
/// Reads as a slice of [`Diagnostic`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Diagnostics(Vec<Diagnostic>);

impl Diagnostics {
    pub fn new() -> Diagnostics {
        Diagnostics::default()
    }

    /// Adds a diagnostic with its code's default severity, and returns it.
    pub fn add(&mut self, code: Code, span: Span, message: impl Into<String>) -> &Diagnostic {
        self.push(Diagnostic::new(code, span, message))
    }

    /// Adds a diagnostic with the given severity, and returns it.
    pub fn add_with(
        &mut self,
        code: Code,
        severity: Severity,
        span: Span,
        message: impl Into<String>,
    ) -> &Diagnostic {
        self.push(Diagnostic {
            severity,
            ..Diagnostic::new(code, span, message)
        })
    }

    pub fn push(&mut self, diagnostic: Diagnostic) -> &Diagnostic {
        self.0.push(diagnostic);
        &self.0[self.0.len() - 1]
    }

    /// Moves every diagnostic in `other` to the end of this list.
    pub fn append(&mut self, other: &mut Diagnostics) {
        self.0.append(&mut other.0);
    }

    pub fn count_of(&self, severity: Severity) -> usize {
        self.0.iter().filter(|d| d.severity == severity).count()
    }

    pub fn has_errors(&self) -> bool {
        self.0.iter().any(|d| d.severity == Severity::Error)
    }

    pub fn has_any(&self, code: Code) -> bool {
        self.first_of(code).is_some()
    }

    pub fn first_of(&self, code: Code) -> Option<&Diagnostic> {
        self.0.iter().find(|d| d.code == code)
    }

    /// Puts the list in source order, so a problem list matches the editor.
    /// Diagnostics at the same place keep their order.
    pub fn sort_by_span(&mut self) {
        self.0
            .sort_by_key(|d| (d.span.start.line, d.span.start.column));
    }

    /// `2 error(s), 1 warning(s)`.
    pub fn summary(&self) -> String {
        messages::error_summary(
            self.count_of(Severity::Error),
            self.count_of(Severity::Warning),
        )
    }

    pub fn into_vec(self) -> Vec<Diagnostic> {
        self.0
    }
}

impl Deref for Diagnostics {
    type Target = [Diagnostic];

    fn deref(&self) -> &[Diagnostic] {
        &self.0
    }
}

impl<'a> IntoIterator for &'a Diagnostics {
    type Item = &'a Diagnostic;
    type IntoIter = std::slice::Iter<'a, Diagnostic>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl IntoIterator for Diagnostics {
    type Item = Diagnostic;
    type IntoIter = std::vec::IntoIter<Diagnostic>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<Vec<Diagnostic>> for Diagnostics {
    fn from(list: Vec<Diagnostic>) -> Diagnostics {
        Diagnostics(list)
    }
}

impl fmt::Display for Diagnostics {
    /// One line per diagnostic, or `No diagnostics`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return f.write_str(messages::NO_DIAGNOSTICS);
        }
        for diagnostic in &self.0 {
            writeln!(f, "{diagnostic}")?;
        }
        Ok(())
    }
}

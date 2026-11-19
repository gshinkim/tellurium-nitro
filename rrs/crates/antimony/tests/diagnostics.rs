//! Ported from the Delphi `Tests/uAntDiagnosticsTests.pas`.
//!
//! Changed for Rust: `FileNameIsStampedOnEachDiagnostic` becomes
//! `file_name_is_given_when_displaying`, because a diagnostic does not store
//! a file name (SPEC.md §9, R7). `OwnsAndFreesItsDiagnostics` is not ported:
//! ownership makes a leak impossible to write.

use antimony::messages;
use antimony::{Code, Diagnostics, Pos, Severity, Span};

fn span_at(line: u32, column: u32) -> Span {
    Span::from_pos(Pos::new(line, column, 0), 1)
}

// --- TDiagnosticsTests -----------------------------------------------------

#[test]
fn starts_empty() {
    let list = Diagnostics::new();
    assert_eq!(list.len(), 0);
    assert!(!list.has_errors());
    assert_eq!(list.to_string(), messages::NO_DIAGNOSTICS);
}

#[test]
fn add_uses_the_default_severity_for_the_code() {
    let mut list = Diagnostics::new();
    let d = list.add(Code::ImplicitParameter, span_at(1, 1), "k1");
    assert_eq!(d.severity, Severity::Warning);
    let d = list.add(Code::SyntaxError, span_at(1, 1), "boom");
    assert_eq!(d.severity, Severity::Error);
}

#[test]
fn explicit_severity_overrides_the_default() {
    // How strict mode raises Delphi §6.5's warnings.
    let mut list = Diagnostics::new();
    let d = list.add_with(
        Code::ImplicitParameter,
        Severity::Error,
        span_at(2, 3),
        "k1",
    );
    assert_eq!(d.severity, Severity::Error);
    assert!(list.has_errors());
}

#[test]
fn add_fmt_formats_the_template() {
    let mut list = Diagnostics::new();
    let d = list.add(
        Code::NegativeStoichiometry,
        span_at(4, 8),
        messages::negative_stoichiometry("A"),
    );
    assert!(d.message.contains("\"A\""), "{}", d.message);
    assert_eq!(d.code, Code::NegativeStoichiometry);
}

#[test]
fn counts_by_severity() {
    let mut list = Diagnostics::new();
    list.add(Code::SyntaxError, span_at(1, 1), "a");
    list.add(Code::SyntaxError, span_at(2, 1), "b");
    list.add(Code::ImplicitParameter, span_at(3, 1), "c");
    list.add_with(Code::None, Severity::Hint, span_at(4, 1), "d");

    assert_eq!(list.count_of(Severity::Error), 2);
    assert_eq!(list.count_of(Severity::Warning), 1);
    assert_eq!(list.count_of(Severity::Hint), 1);
    assert_eq!(list.len(), 4);
    assert_eq!(list.summary(), messages::error_summary(2, 1));
    assert_eq!(list.summary(), "2 error(s), 1 warning(s)");
}

#[test]
fn has_errors_ignores_warnings() {
    let mut list = Diagnostics::new();
    list.add(Code::UninitialisedSymbol, span_at(1, 1), "k1");
    list.add(Code::ImplicitParameter, span_at(2, 1), "k2");
    assert!(
        !list.has_errors(),
        "a model full of §6.5 warnings must still load"
    );
}

#[test]
fn finds_by_code() {
    let mut list = Diagnostics::new();
    list.add(Code::SyntaxError, span_at(1, 1), "a");
    list.add(Code::NegativeStoichiometry, span_at(2, 1), "b");

    assert!(list.has_any(Code::NegativeStoichiometry));
    assert!(!list.has_any(Code::TypeConflict));
    assert!(list.first_of(Code::SyntaxError).is_some());
    assert!(list.first_of(Code::TypeConflict).is_none());
}

#[test]
fn file_name_is_given_when_displaying() {
    let mut list = Diagnostics::new();
    let d = list.add(Code::SyntaxError, span_at(1, 1), "a");
    assert!(d.with_file("model.ant").starts_with("model.ant(1:1): "));
    let d = list.add(Code::SyntaxError, Span::default(), "b");
    assert_eq!(
        d.with_file("model.ant"),
        "model.ant: error: b [dcSyntaxError]"
    );
    assert_eq!(d.to_string(), "error: b [dcSyntaxError]");
}

#[test]
fn sorts_into_source_order() {
    let mut list = Diagnostics::new();
    list.add(Code::SyntaxError, span_at(9, 2), "third");
    list.add(Code::SyntaxError, span_at(1, 5), "second");
    list.add(Code::SyntaxError, span_at(1, 1), "first");

    list.sort_by_span();

    let messages: Vec<&str> = list.iter().map(|d| d.message.as_str()).collect();
    assert_eq!(messages, ["first", "second", "third"]);
}

#[test]
fn renders_file_position_severity_and_code() {
    let mut list = Diagnostics::new();
    let text = list
        .add(Code::SyntaxError, span_at(3, 17), "unexpected \"->\"")
        .with_file("model.ant");
    assert_eq!(
        text,
        "model.ant(3:17): error: unexpected \"->\" [dcSyntaxError]"
    );
}

// --- TDiagnosticPolicyTests ------------------------------------------------

#[test]
fn uninitialised_symbol_is_a_warning() {
    assert_eq!(
        Code::UninitialisedSymbol.default_severity(),
        Severity::Warning
    );
}

#[test]
fn implicit_parameter_is_a_warning() {
    assert_eq!(
        Code::ImplicitParameter.default_severity(),
        Severity::Warning
    );
}

#[test]
fn information_loss_is_a_warning() {
    // Delphi §9.1.4 #5: data is never dropped silently, but a loss is not a
    // load failure.
    assert_eq!(Code::InformationLoss.default_severity(), Severity::Warning);
}

#[test]
fn negative_stoichiometry_is_an_error() {
    assert_eq!(
        Code::NegativeStoichiometry.default_severity(),
        Severity::Error
    );
}

#[test]
fn syntax_problems_are_errors() {
    for code in [
        Code::SyntaxError,
        Code::UnterminatedString,
        Code::DnaNotSupported,
        Code::NotYetSupported,
    ] {
        assert_eq!(code.default_severity(), Severity::Error, "{code}");
    }
}

#[test]
fn every_code_has_a_name() {
    // The `tests/invalid` fixtures name codes, so the names are part of the
    // contract. They are all distinct and all read back.
    for code in Code::ALL {
        assert!(code.name().starts_with("dc"), "{code:?}");
        assert_eq!(Code::from_name(code.name()), Some(code));
    }
    assert_eq!(Code::DnaNotSupported.name(), "dcDNANotSupported");
}

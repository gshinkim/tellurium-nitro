//! Ported from the Delphi `Tests/uAntTypesTests.pas`.
//!
//! Not ported, because Rust cannot have the defect they guard against
//! (SPEC.md §9, R12): `SurvivesCommaDecimalLocaleOnInput`,
//! `SurvivesCommaDecimalLocaleOnOutput` and the locale setup of the number
//! tests. Rust's number conversions do not read the machine's locale, so the
//! same assertions without a locale are already in `parses_plain_decimal`
//! and `emits_dot_decimal_separator`. `RejectsGarbage`'s check that
//! `AntStrToFloat` raises has no counterpart: `number::parse` returns
//! `Option`.

use antimony::builtins::is_builtin_function;
use antimony::lexer::{Keyword, is_digit, is_exponent_marker, is_ident_char, is_ident_start};
use antimony::model::SymbolKind;
use antimony::number;
use antimony::{Pos, Span};

// --- TSourcePosTests -------------------------------------------------------

#[test]
fn unknown_is_unknown() {
    assert!(Pos::default().is_unknown());
    assert!(Span::default().is_unknown());
    assert_eq!(Pos::default().to_string(), "<unknown>");
}

#[test]
fn created_position_is_known() {
    let pos = Pos::new(3, 17, 42);
    assert!(!pos.is_unknown());
    assert_eq!((pos.line, pos.column, pos.offset), (3, 17, 42));
    assert_eq!(pos.to_string(), "3:17");
}

#[test]
fn span_length_is_offset_difference() {
    let span = Span::new(Pos::new(1, 1, 0), Pos::new(1, 6, 5));
    assert_eq!(span.len(), 5);
    assert_eq!(Span::default().len(), 0);
}

#[test]
fn span_from_pos_advances_column_and_offset() {
    let span = Span::from_pos(Pos::new(2, 4, 10), 3);
    assert_eq!(span.end.line, 2);
    assert_eq!(span.end.column, 7);
    assert_eq!(span.end.offset, 13);
    assert_eq!(span.len(), 3);
}

#[test]
fn union_covers_both_spans() {
    let left = Span::new(Pos::new(1, 1, 0), Pos::new(1, 4, 3));
    let right = Span::new(Pos::new(2, 1, 10), Pos::new(2, 5, 14));
    let union = left.union(right);
    assert_eq!((union.start.offset, union.end.offset), (0, 14));
    // Order must not matter.
    let union = right.union(left);
    assert_eq!((union.start.offset, union.end.offset), (0, 14));
}

#[test]
fn union_with_unknown_is_identity() {
    let span = Span::new(Pos::new(1, 1, 0), Pos::new(1, 4, 3));
    assert_eq!(span.union(Span::default()), span);
    assert_eq!(Span::default().union(span), span);
}

#[test]
fn slice_of_returns_the_original_text() {
    // Offsets 3..5 are the arrow.
    let span = Span::new(Pos::new(1, 4, 3), Pos::new(1, 6, 5));
    assert_eq!(span.slice_of("S1 -> S2; k1*S1"), "->");
}

#[test]
fn slice_of_unknown_span_is_empty() {
    assert_eq!(Span::default().slice_of("anything"), "");
}

// --- TKeywordTableTests ----------------------------------------------------

#[test]
fn known_keywords_resolve() {
    assert_eq!(Keyword::from_word("model"), Some(Keyword::Model));
    assert_eq!(Keyword::from_word("end"), Some(Keyword::End));
    assert_eq!(Keyword::from_word("species"), Some(Keyword::Species));
    assert_eq!(
        Keyword::from_word("substanceOnly"),
        Some(Keyword::SubstanceOnly)
    );
    assert_eq!(
        Keyword::from_word("conversionFactor"),
        Some(Keyword::ConversionFactor)
    );
    assert_eq!(Keyword::from_word("inf"), Some(Keyword::Inf));
}

#[test]
fn non_keywords_resolve_to_none() {
    assert_eq!(Keyword::from_word("S1"), None);
    assert_eq!(Keyword::from_word(""), None);
    assert_eq!(Keyword::from_word("modelling"), None);
}

#[test]
fn lookup_is_case_sensitive() {
    // Identifiers are case-sensitive, and so is the keyword table: "Model"
    // is an ordinary identifier.
    assert_eq!(Keyword::from_word("Model"), None);
    assert_eq!(Keyword::from_word("END"), None);
    assert_eq!(Keyword::from_word("substanceonly"), None);
}

#[test]
fn every_spelling_round_trips() {
    for keyword in Keyword::ALL {
        assert_eq!(
            Keyword::from_word(keyword.spelling()),
            Some(keyword),
            "keyword {} does not round-trip",
            keyword.spelling()
        );
    }
}

#[test]
fn builtin_functions_are_recognised() {
    for name in ["rateOf", "piecewise", "avogadro", "geq", "sqrt", "time"] {
        assert!(is_builtin_function(name), "{name}");
    }
}

#[test]
fn builtin_lookup_rejects_unknown_names() {
    assert!(!is_builtin_function("MM"));
    assert!(
        !is_builtin_function("Sqrt"),
        "lookup must be case-sensitive"
    );
    assert!(!is_builtin_function(""));
}

// --- TPromotionTests -------------------------------------------------------

#[test]
fn unknown_promotes_to_anything() {
    for kind in SymbolKind::ALL {
        assert!(SymbolKind::Unknown.can_promote_to(kind), "{}", kind.name());
    }
}

#[test]
fn parameter_promotes_to_species_and_compartment() {
    // `k1 = 0.1` looks like an initial concentration until k1 turns up in a
    // reaction (Delphi §6.2).
    assert!(SymbolKind::Parameter.can_promote_to(SymbolKind::Species));
    assert!(SymbolKind::Parameter.can_promote_to(SymbolKind::Compartment));
    assert!(!SymbolKind::Parameter.can_promote_to(SymbolKind::Reaction));
}

#[test]
fn species_does_not_demote() {
    assert!(!SymbolKind::Species.can_promote_to(SymbolKind::Parameter));
    assert!(!SymbolKind::Species.can_promote_to(SymbolKind::Compartment));
    assert!(!SymbolKind::Species.can_promote_to(SymbolKind::Reaction));
}

#[test]
fn compartment_is_terminal() {
    assert!(!SymbolKind::Compartment.can_promote_to(SymbolKind::Species));
    assert!(!SymbolKind::Compartment.can_promote_to(SymbolKind::Parameter));
    assert!(!SymbolKind::Compartment.can_promote_to(SymbolKind::Reaction));
}

#[test]
fn identity_is_always_legal() {
    for kind in SymbolKind::ALL {
        assert!(kind.can_promote_to(kind), "{}", kind.name());
    }
}

// --- TNumberConversionTests ------------------------------------------------

#[test]
fn parses_plain_decimal() {
    assert_eq!(number::parse("0.1"), Some(0.1));
    assert_eq!(number::parse("-2.5"), Some(-2.5));
    assert_eq!(number::parse("42"), Some(42.0));
}

#[test]
fn parses_scientific() {
    assert_eq!(number::parse("1.2e-3"), Some(0.0012));
    assert_eq!(number::parse("6.022E23"), Some(6.022e23));
}

#[test]
fn parses_inf_and_nan() {
    assert_eq!(number::parse("inf"), Some(f64::INFINITY));
    assert_eq!(number::parse("-inf"), Some(f64::NEG_INFINITY));
    assert!(number::parse("nan").is_some_and(f64::is_nan));
}

#[test]
fn rejects_garbage() {
    assert_eq!(number::parse(""), None);
    assert_eq!(number::parse("S1"), None);
    assert_eq!(number::parse("1.2.3"), None);
    assert_eq!(number::parse("not a number"), None);
}

#[test]
fn emits_dot_decimal_separator() {
    assert_eq!(number::format(0.1), "0.1");
    assert_eq!(number::format(-2.5), "-2.5");
    assert_eq!(number::format(42.0), "42");
}

#[test]
fn emits_inf_and_nan_spellings() {
    assert_eq!(number::format(f64::INFINITY), "inf");
    assert_eq!(number::format(f64::NEG_INFINITY), "-inf");
    assert_eq!(number::format(f64::NAN), "nan");
}

#[test]
fn round_trips_exactly() {
    for value in [
        0.1,
        1.0 / 3.0,
        6.02214076e23,
        1.2e-300,
        -7.5e18,
        0.0,
        1.0,
        std::f64::consts::E,
    ] {
        let text = number::format(value);
        assert_eq!(number::parse(&text), Some(value), "via {text}");
    }
}

#[test]
fn comma_is_not_a_decimal_separator() {
    // A comma separates declaration items in Antimony; it must never be
    // read as part of a number.
    assert_eq!(number::parse("0,1"), None);
}

// --- TCharClassTests -------------------------------------------------------

#[test]
fn identifier_start_accepts_letters_and_underscore() {
    assert!(is_ident_start('S'));
    assert!(is_ident_start('k'));
    assert!(is_ident_start('_'));
}

#[test]
fn identifier_start_rejects_digits() {
    assert!(!is_ident_start('1'));
    assert!(!is_ident_start('$'));
    assert!(!is_ident_start(' '));
}

#[test]
fn identifier_char_accepts_digits() {
    assert!(is_ident_char('1'));
    assert!(is_ident_char('_'));
    assert!(!is_ident_char('-'));
}

#[test]
fn non_ascii_is_not_an_identifier_char() {
    assert!(!is_ident_char('\u{03b1}')); // Greek small alpha
    assert!(!is_exponent_marker('\u{03b5}'));
    assert!(is_exponent_marker('e'));
    assert!(is_exponent_marker('E'));
    assert!(is_digit('7'));
}

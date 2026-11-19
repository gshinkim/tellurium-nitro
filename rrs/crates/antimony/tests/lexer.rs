//! Ported from the Delphi `Tests/uAntLexerTests.pas`.
//!
//! Changed for Rust: `NumbersAreLocaleIndependent` keeps only its
//! assertions, since Rust has no locale to set (SPEC.md §9, R12). The Delphi
//! `Kinds` helper compared token-kind names as a string; here the kinds are
//! compared as values. `EveryCorpusFileEndsWithEOF` and
//! `NoCorpusFileProducesAnUnknownToken` are checks inside the one corpus
//! sweep that the other corpus tests share, rather than sweeps of their own.
//! The tests after "Added in the Rust port" have no Delphi counterpart.
//!
//! The corpus tests need `ANTIMONY_CORPUS` and are skipped without it.

mod common;

use std::fmt::Write as _;
use std::time::Instant;

use antimony::Code;
use antimony::Diagnostics;
use antimony::lexer::{Keyword, Lexer, Token, TokenKind as K, tokenize};

fn lex(source: &str) -> (Vec<Token>, Diagnostics) {
    tokenize(source)
}

fn tokens(source: &str) -> Vec<Token> {
    tokenize(source).0
}

fn kinds(source: &str) -> Vec<K> {
    tokens(source).iter().map(|t| t.kind).collect()
}

fn has_keyword(tokens: &[Token], keyword: Keyword) -> bool {
    tokens.iter().any(|t| t.is_keyword(keyword))
}

// --- TLexerBasicsTests -----------------------------------------------------

#[test]
fn empty_source_yields_only_eof() {
    assert_eq!(kinds(""), [K::Eof]);
}

#[test]
fn whitespace_only_yields_only_eof() {
    assert_eq!(kinds("   \t  "), [K::Eof]);
    assert_eq!(kinds("\r\n\r\n"), [K::Eof]);
}

#[test]
fn eof_repeats_forever() {
    // The lexer parks on end of file instead of running past it.
    let mut lexer = Lexer::new("k1");
    lexer.next_token();
    for _ in 0..10 {
        assert_eq!(lexer.next_token().kind, K::Eof);
    }
}

#[test]
fn identifiers_are_case_sensitive() {
    let t = tokens("S1 s1 _x9");
    assert_eq!(
        (t[0].text.as_str(), t[1].text.as_str(), t[2].text.as_str()),
        ("S1", "s1", "_x9")
    );
}

#[test]
fn keywords_are_annotated_not_reserved() {
    // `conversionFactor` is a keyword in one statement position only, so the
    // lexer notes it and the parser decides. The kind stays Identifier.
    let t = tokens("conversionFactor Model species");
    assert_eq!(t[0].kind, K::Identifier);
    assert!(t[0].is_keyword(Keyword::ConversionFactor));
    assert!(
        !t[1].is_keyword(Keyword::Model),
        "\"Model\" is not \"model\""
    );
    assert!(t[2].is_keyword(Keyword::Species));
}

#[test]
fn peek_does_not_consume() {
    let mut lexer = Lexer::new("a b");
    assert_eq!(lexer.peek(0).text, "a");
    assert_eq!(lexer.peek(0).text, "a");
    assert_eq!(lexer.next_token().text, "a");
    assert_eq!(lexer.peek(0).text, "b");
}

#[test]
fn peek_reaches_arbitrary_depth() {
    // Telling a reaction from an assignment means looking ahead to the arrow.
    let mut lexer = Lexer::new("2 S1 + S2 -> S3");
    assert_eq!(lexer.peek(4).kind, K::ArrowRev);
    assert_eq!(lexer.peek(0).text, "2", "peeking deep must not consume");
}

#[test]
fn peek_past_end_returns_eof() {
    let mut lexer = Lexer::new("a");
    assert_eq!(lexer.peek(50).kind, K::Eof);
}

#[test]
fn try_consume_only_moves_on_match() {
    let mut lexer = Lexer::new("species S1");
    assert!(lexer.try_consume(K::Semicolon).is_none());
    assert!(lexer.try_consume_keyword(Keyword::Species));
    assert_eq!(lexer.peek(0).text, "S1");
    assert_eq!(lexer.current().text, "species");
}

// --- TLexerOperatorTests ---------------------------------------------------

#[test]
fn arrows_lex_as_single_tokens() {
    assert_eq!(
        kinds("-> => <-> -|"),
        [
            K::ArrowRev,
            K::ArrowIrrev,
            K::ArrowBiDir,
            K::Inhibits,
            K::Eof
        ]
    );
}

#[test]
fn reversible_arrow_is_the_single_dash() {
    // `->` is reversible and `=>` irreversible (Delphi §9.1.1). Getting this
    // backwards would invert every reaction in every model.
    assert_eq!(tokens("A -> B")[1].kind, K::ArrowRev);
    assert_eq!(tokens("A => B")[1].kind, K::ArrowIrrev);
}

#[test]
fn two_character_operators() {
    assert_eq!(
        kinds("&& || :="),
        [K::AndAnd, K::OrOr, K::ColonAssign, K::Eof]
    );
}

#[test]
fn comparison_operators() {
    assert_eq!(
        kinds("== != < <= > >="),
        [
            K::Equal,
            K::NotEqual,
            K::Less,
            K::LessEqual,
            K::Greater,
            K::GreaterEqual,
            K::Eof
        ]
    );
}

#[test]
fn assignment_versus_assignment_rule() {
    assert_eq!(
        kinds("k1 = 1"),
        [K::Identifier, K::Assign, K::Number, K::Eof]
    );
    assert_eq!(
        kinds("k1 := 1"),
        [K::Identifier, K::ColonAssign, K::Number, K::Eof]
    );
}

#[test]
fn activation_arrow_needs_a_standalone_o() {
    assert_eq!(
        kinds("S3 -o J1"),
        [K::Identifier, K::Activates, K::Identifier, K::Eof]
    );
    assert_eq!(
        kinds("x -old"),
        [K::Identifier, K::Minus, K::Identifier, K::Eof]
    );
}

#[test]
fn dash_dash_is_recognised_for_rejection() {
    // DNA strand syntax is not supported, but the lexer names it so the
    // parser can say so clearly (Delphi N2).
    assert_eq!(
        kinds("g1--g2"),
        [K::Identifier, K::DashDash, K::Identifier, K::Eof]
    );
}

#[test]
fn prime_and_dollar_and_at() {
    assert_eq!(
        kinds("S1' = 5"),
        [K::Identifier, K::Prime, K::Assign, K::Number, K::Eof]
    );
    assert_eq!(kinds("$Xo"), [K::Dollar, K::Identifier, K::Eof]);
    assert_eq!(kinds("@"), [K::At, K::Eof]);
}

#[test]
fn ellipsis_is_one_token() {
    assert_eq!(kinds("..."), [K::Ellipsis, K::Eof]);
    assert_eq!(kinds("A.B"), [K::Identifier, K::Dot, K::Identifier, K::Eof]);
    assert_eq!(
        kinds(".5"),
        [K::Number, K::Eof],
        "a leading-dot number still wins"
    );
}

#[test]
fn unknown_interaction_is_not_one_token() {
    // `-(` is the unknown-interaction arrow, but a minus before a bracket is
    // far more common (52 times in the corpus) and looks the same.
    assert_eq!(
        kinds("k*-(a)"),
        [
            K::Identifier,
            K::Star,
            K::Minus,
            K::LParen,
            K::Identifier,
            K::RParen,
            K::Eof
        ]
    );
}

#[test]
fn lone_ampersand_is_unknown_and_reported() {
    let (t, d) = lex("a & b");
    assert_eq!(t[1].kind, K::Unknown);
    assert!(d.has_any(Code::UnexpectedChar));
    assert_eq!(d[0].message, "Unexpected character '&'");
}

#[test]
fn unknown_character_does_not_stall_the_lexer() {
    assert_eq!(
        kinds("~~ok"),
        [K::Unknown, K::Unknown, K::Identifier, K::Eof]
    );
}

// --- TLexerNumberTests -----------------------------------------------------

#[test]
fn integer_and_decimal() {
    let t = tokens("42 0.5 3.");
    assert_eq!((t[0].number, t[1].number, t[2].number), (42.0, 0.5, 3.0));
}

#[test]
fn scientific() {
    let t = tokens("1.2e-3 6.022E23 1e5");
    assert_eq!(
        (t[0].number, t[1].number, t[2].number),
        (0.0012, 6.022e23, 100000.0)
    );
}

#[test]
fn leading_dot_number() {
    let t = tokens(".5");
    assert_eq!(t[0].kind, K::Number);
    assert_eq!(t[0].number, 0.5);
}

#[test]
fn minus_is_never_part_of_a_number() {
    // The most serious defect of the lexer before Delphi's: swallowing the
    // sign made `a-1` two operands with no operator between them.
    assert_eq!(kinds("a-1"), [K::Identifier, K::Minus, K::Number, K::Eof]);
    let t = tokens("-5");
    assert_eq!(t[0].kind, K::Minus);
    assert_eq!(t[1].number, 5.0);
}

#[test]
fn subtraction_of_two_numbers() {
    assert_eq!(kinds("2-3"), [K::Number, K::Minus, K::Number, K::Eof]);
}

#[test]
fn exponent_marker_without_digits_rolls_back() {
    let (t, d) = lex("2e");
    assert_eq!(
        t.iter().map(|t| t.kind).collect::<Vec<_>>(),
        [K::Number, K::Identifier, K::Eof]
    );
    assert_eq!(t[0].number, 2.0);
    assert_eq!(t[1].text, "e");
    assert!(!d.has_errors());
}

#[test]
fn inf_and_nan_are_number_literals() {
    // Four spellings, and only four: libantimony writes `INF` and `NaN`
    // (case 00950), while case 01810 uses `Inf` as a parameter name.
    let t = tokens("inf nan INF NaN");
    for (i, token) in t[..4].iter().enumerate() {
        assert_eq!(token.kind, K::Number, "{}", token.text);
        if i % 2 == 0 {
            assert_eq!(token.number, f64::INFINITY, "{}", token.text);
        } else {
            assert!(token.number.is_nan(), "{}", token.text);
        }
    }
}

#[test]
fn inf_lookalikes_remain_identifiers() {
    // Corpus case 01810 exists to test exactly this.
    assert_eq!(
        kinds("Inf INF_ inf_ NAN"),
        [
            K::Identifier,
            K::Identifier,
            K::Identifier,
            K::Identifier,
            K::Eof
        ]
    );
}

#[test]
fn numbers_are_locale_independent() {
    assert_eq!(tokens("0.1")[0].number, 0.1);
    assert_eq!(kinds("1,2"), [K::Number, K::Comma, K::Number, K::Eof]);
}

// --- TLexerStringTests -----------------------------------------------------

#[test]
fn simple_string() {
    let t = tokens("S1 is \"Glucose 6-phosphate\"");
    assert_eq!(t[2].kind, K::String);
    assert_eq!(t[2].text, "Glucose 6-phosphate");
}

#[test]
fn escapes_are_decoded() {
    assert_eq!(tokens(r#""a\"b\\c\nd""#)[0].text, "a\"b\\c\nd");
}

#[test]
fn unknown_escape_is_kept_verbatim() {
    // An unrecognised escape keeps both characters, so nothing is lost. `\t`
    // really is a tab, as the language says, even inside a Windows path.
    assert_eq!(tokens(r#""a\qb""#)[0].text, r"a\qb");
    assert_eq!(tokens(r#""C:\temp""#)[0].text, "C:\temp");
}

#[test]
fn unterminated_string_is_reported() {
    let (t, d) = lex("\"no end");
    assert!(d.has_any(Code::UnterminatedString));
    assert_eq!(t[0].text, "no end");
}

#[test]
fn string_does_not_span_lines() {
    let (t, d) = lex("\"abc\nk1 = 1");
    assert!(d.has_any(Code::UnterminatedString));
    // The next line is still read as ordinary statements.
    assert_eq!(t[2].text, "k1");
}

#[test]
fn span_covers_the_quotes() {
    let source = "x is \"hi\"";
    assert_eq!(
        tokens(source)[2].span.slice_of(source),
        "\"hi\"",
        "the span must recover the raw literal, quotes included"
    );
}

// --- TLexerTriviaTests -----------------------------------------------------

#[test]
fn line_comment_is_trivia_on_the_next_token() {
    let t = tokens("// Compartments and Species:\ncompartment c;");
    assert_eq!(t[0].text, "compartment");
    assert_eq!(t[0].trivia, "// Compartments and Species:");
}

#[test]
fn hash_starts_a_line_comment() {
    // The language reference gives `#` and `//` equally. The corpus has no
    // `#` at all, because libantimony always writes `//`.
    assert_eq!(
        kinds("S1 = 10  #The initial concentration of S1\nS2 = 2;"),
        [
            K::Identifier,
            K::Assign,
            K::Number,
            K::Newline,
            K::Identifier,
            K::Assign,
            K::Number,
            K::Semicolon,
            K::Eof
        ]
    );
    let t = tokens("# leading comment\nk = 1");
    assert_eq!(t[0].text, "k");
    assert_eq!(
        t[0].trivia, "# leading comment",
        "a # comment is trivia, exactly like //"
    );
}

#[test]
fn block_comment_is_trivia() {
    let t = tokens("/* two\n   lines */ k1 = 1");
    assert_eq!(t[0].text, "k1");
    assert!(t[0].has_trivia());
    assert!(t[0].trivia.contains("lines"));
}

#[test]
fn unterminated_block_comment_is_reported() {
    let (_, d) = lex("k1 = 1 /* and then nothing");
    assert!(d.has_any(Code::UnterminatedComment));
}

#[test]
fn trailing_comment_carries_to_the_next_statement() {
    // A comment belongs to the statement after it, so a comment at the end
    // of a line crosses the line break rather than being dropped.
    let t = tokens("S1 = 1; // about S1\nS2 = 2;");
    assert_eq!(t[5].text, "S2");
    assert_eq!(t[5].trivia, "// about S1");
}

#[test]
fn comment_at_end_of_file_is_not_lost() {
    let t = tokens("k1 = 1\n// the end");
    let last = t.last().unwrap();
    assert_eq!(last.kind, K::Eof);
    assert_eq!(
        last.trivia, "// the end",
        "the writer must be able to put it back"
    );
}

#[test]
fn blank_line_run_is_one_terminator() {
    assert_eq!(
        kinds("a\n\n\nb"),
        [K::Identifier, K::Newline, K::Identifier, K::Eof]
    );
}

#[test]
fn leading_newlines_produce_no_terminator() {
    assert_eq!(kinds("\n\na"), [K::Identifier, K::Eof]);
}

#[test]
fn line_continuation_suppresses_the_terminator() {
    // The corpus never uses this, so only a fixture covers it.
    assert_eq!(
        kinds("k = 1 + \\\n2"),
        [
            K::Identifier,
            K::Assign,
            K::Number,
            K::Plus,
            K::Number,
            K::Eof
        ]
    );
    assert_eq!(
        kinds("a\nb"),
        [K::Identifier, K::Newline, K::Identifier, K::Eof]
    );
}

// --- TLexerPositionTests ---------------------------------------------------

#[test]
fn line_and_column_track_across_lines() {
    let t = tokens("ab cd\n  ef");
    assert_eq!((t[0].span.start.line, t[0].span.start.column), (1, 1));
    assert_eq!((t[1].span.start.line, t[1].span.start.column), (1, 4));
    // t[2] is the newline terminator.
    assert_eq!((t[3].span.start.line, t[3].span.start.column), (2, 3));
    assert_eq!(t[3].text, "ef");
}

#[test]
fn span_slices_back_to_the_original_text() {
    // Every token can point at the text it came from.
    let source = "J1: 2 S1 + S2 -> S3; k1*S1^2*S2";
    for token in tokens(source) {
        if !matches!(token.kind, K::Eof | K::Newline) {
            assert_eq!(token.span.slice_of(source), token.text, "{token}");
            assert_eq!(
                &source[token.span.start.offset..token.span.end.offset],
                token.text
            );
        }
    }
    assert_eq!(tokens(source)[6].span.slice_of(source), "->");
}

#[test]
fn crlf_counts_as_one_line_break() {
    let t = tokens("a\r\nb");
    assert_eq!((t[2].span.start.line, t[2].span.start.column), (2, 1));
}

#[test]
fn block_comment_advances_the_line_number() {
    assert_eq!(
        tokens("/* one\ntwo */ k1")[0].span.start.line,
        2,
        "a token after a multi-line comment must not report line 1"
    );
}

// --- TLexerLayoutTests -----------------------------------------------------
//
// Changed for Rust: the Delphi lexer read `layout … end` and `layout { … }`
// as one token. Antimony has no such block; its layout is dotted statements,
// so `layout` is an ordinary word (docs/antimony-layout-design.md, D1).

#[test]
fn layout_is_an_ordinary_word() {
    assert_eq!(
        kinds("model.layout = on\nk1 = 1"),
        [
            K::Identifier,
            K::Dot,
            K::Identifier,
            K::Assign,
            K::Identifier,
            K::Newline,
            K::Identifier,
            K::Assign,
            K::Number,
            K::Eof
        ]
    );
}

#[test]
fn layout_followed_by_end_swallows_nothing() {
    assert_eq!(
        kinds("layout\nk1 = 1\nend"),
        [
            K::Identifier,
            K::Newline,
            K::Identifier,
            K::Assign,
            K::Number,
            K::Newline,
            K::Identifier,
            K::Eof
        ]
    );
}

// --- TLexerStatementShapeTests ---------------------------------------------
//
// The statement forms the corpus actually contains, pinned for the parser.

#[test]
fn reaction_with_stoichiometry() {
    assert_eq!(
        kinds("J1: 2 S1 + S2 -> S3; k1*S1"),
        [
            K::Identifier,
            K::Colon,
            K::Number,
            K::Identifier,
            K::Plus,
            K::Identifier,
            K::ArrowRev,
            K::Identifier,
            K::Semicolon,
            K::Identifier,
            K::Star,
            K::Identifier,
            K::Eof
        ]
    );
}

#[test]
fn main_module_marker() {
    assert_eq!(
        kinds("model *case00001()"),
        [
            K::Identifier,
            K::Star,
            K::Identifier,
            K::LParen,
            K::RParen,
            K::Eof
        ]
    );
}

#[test]
fn algebraic_rule() {
    assert_eq!(
        kinds("_alg0: 0 = -1*k1 + S1;"),
        [
            K::Identifier,
            K::Colon,
            K::Number,
            K::Assign,
            K::Minus,
            K::Number,
            K::Star,
            K::Identifier,
            K::Plus,
            K::Identifier,
            K::Semicolon,
            K::Eof
        ]
    );
}

#[test]
fn empty_kinetic_law() {
    assert_eq!(
        kinds("R16: J => M; ;"),
        [
            K::Identifier,
            K::Colon,
            K::Identifier,
            K::ArrowIrrev,
            K::Identifier,
            K::Semicolon,
            K::Semicolon,
            K::Eof
        ]
    );
}

#[test]
fn unit_attachment() {
    assert!(tokens("S1 has substance_per_volume;")[1].is_keyword(Keyword::Has));
}

#[test]
fn substance_only_species() {
    let t = tokens("substanceOnly species S1 in c;");
    assert!(t[0].is_keyword(Keyword::SubstanceOnly));
    assert!(t[1].is_keyword(Keyword::Species));
    assert!(t[3].is_keyword(Keyword::In));
}

#[test]
fn boundary_species_in_declaration_list() {
    assert_eq!(
        kinds("species G in Cell, $T in Cell;"),
        [
            K::Identifier,
            K::Identifier,
            K::Identifier,
            K::Identifier,
            K::Comma,
            K::Dollar,
            K::Identifier,
            K::Identifier,
            K::Identifier,
            K::Semicolon,
            K::Eof
        ]
    );
}

#[test]
fn event_with_after_and_modifiers() {
    let t = tokens("event1: at 1 after S1 < 0.1, priority=1: S1 = 1;");
    assert!(t[2].is_keyword(Keyword::At));
    assert!(t[4].is_keyword(Keyword::After));
    assert!(t[9].is_keyword(Keyword::Priority));
}

#[test]
fn event_modifiers_from_the_language_reference() {
    let t = tokens("E1: at (x>5), t0=false: y=3, x=r+2");
    assert!(has_keyword(&t, Keyword::At));
    assert!(has_keyword(&t, Keyword::T0), "t0 must be a keyword");

    let t = tokens("at (x>5), persistent=true, fromTrigger=false, priority=2: y=3");
    assert!(has_keyword(&t, Keyword::At));
    assert!(has_keyword(&t, Keyword::Persistent));
    assert!(has_keyword(&t, Keyword::FromTrigger));
    assert!(has_keyword(&t, Keyword::Priority));
    assert!(
        !has_keyword(&t, Keyword::Delay),
        "delay is not present here"
    );
}

#[test]
fn dna_is_a_keyword_so_it_can_be_rejected() {
    assert!(tokens("dna")[0].is_keyword(Keyword::Dna));
    assert_eq!(
        kinds("dna1: --P1--G1--"),
        [
            K::Identifier,
            K::Colon,
            K::DashDash,
            K::Identifier,
            K::DashDash,
            K::Identifier,
            K::DashDash,
            K::Eof
        ]
    );
}

#[test]
fn rate_rule() {
    assert_eq!(
        kinds("S1' = -k1*S1"),
        [
            K::Identifier,
            K::Prime,
            K::Assign,
            K::Minus,
            K::Identifier,
            K::Star,
            K::Identifier,
            K::Eof
        ]
    );
}

#[test]
fn dotted_name() {
    assert_eq!(
        kinds("A.S1"),
        [K::Identifier, K::Dot, K::Identifier, K::Eof]
    );
}

#[test]
fn conversion_factor_statement() {
    let t = tokens("S1 conversionFactor s1_cf;");
    assert!(t[1].is_keyword(Keyword::ConversionFactor));
    assert_eq!(t[2].text, "s1_cf");
}

#[test]
fn empty_rule_keeps_both_sides() {
    // `p = ;` is real corpus text (case 01234); the parser decides what it
    // means.
    assert_eq!(
        kinds("p = ;"),
        [K::Identifier, K::Assign, K::Semicolon, K::Eof]
    );
}

// --- Added in the Rust port ------------------------------------------------

#[test]
fn fenced_string_is_raw_and_keeps_line_numbers() {
    let source = "model notes ```\n <p>a \\ \"b\"</p>\r\n  x\n```\nk1";
    let (t, d) = lex(source);
    assert!(d.is_empty(), "{d}");
    assert_eq!(t[2].kind, K::String);
    assert_eq!(t[2].text, " <p>a \\ \"b\"</p>\n  x");
    // The Delphi version once lost a line for every fenced note.
    assert_eq!((t[4].text.as_str(), t[4].span.start.line), ("k1", 5));
}

#[test]
fn unterminated_fence_is_reported() {
    let (t, d) = lex("notes ```\nnever closed");
    assert!(d.has_any(Code::UnterminatedString));
    assert_eq!(t[1].text, "never closed");
}

#[test]
fn columns_count_characters_and_skip_a_byte_order_mark() {
    let t = tokens("\u{feff}x = \"\u{e9}t\u{e9}\" y");
    assert_eq!((t[0].span.start.line, t[0].span.start.column), (1, 1));
    assert_eq!(t[3].text, "y");
    assert_eq!(t[3].span.start.column, 11);
}

#[test]
fn multibyte_unknown_character_is_one_token() {
    let (t, d) = lex("a \u{2192} b");
    assert_eq!(t[1].kind, K::Unknown);
    assert_eq!(t[1].text, "\u{2192}");
    assert_eq!(t[2].text, "b");
    assert_eq!(d.len(), 1);
}

#[test]
fn backslash_before_a_line_break_in_a_string_keeps_the_string_open() {
    // The Delphi behaviour: an unknown escape keeps both characters, and a
    // line break is one of them. Positions after it are still right.
    let (t, d) = lex("x is \"a\\\nb\" y");
    assert!(d.is_empty(), "{d}");
    assert_eq!(t[2].text, "a\\\nb");
    assert_eq!((t[3].text.as_str(), t[3].span.start.line), ("y", 2));
}

#[test]
fn lone_carriage_return_is_a_line_break() {
    let t = tokens("a\rb");
    assert_eq!(t[1].kind, K::Newline);
    assert_eq!((t[2].span.start.line, t[2].span.start.column), (2, 1));
}

// --- TLexerCorpusTests -----------------------------------------------------

struct Sweep {
    files: usize,
    tokens: usize,
    arrows_rev: usize,
    arrows_irrev: usize,
    problems: Vec<String>,
}

fn sweep_corpus(test: &str) -> Option<Sweep> {
    let root = common::corpus_root(test)?;
    let mut sweep = Sweep {
        files: 0,
        tokens: 0,
        arrows_rev: 0,
        arrows_irrev: 0,
        problems: Vec::new(),
    };
    for (case, path) in common::corpus_files(&root) {
        let (tokens, diagnostics) = tokenize(&common::read(&path));
        assert!(tokens.len() > 1, "{case}");
        assert_eq!(tokens.last().map(|t| t.kind), Some(K::Eof), "{case}");
        for token in &tokens {
            match token.kind {
                K::Unknown => sweep
                    .problems
                    .push(format!("{case}: unknown token {token} at {}", token.span)),
                K::ArrowRev => sweep.arrows_rev += 1,
                K::ArrowIrrev => sweep.arrows_irrev += 1,
                _ => {}
            }
        }
        sweep.tokens += tokens.len();
        if let Some(first) = diagnostics.first() {
            sweep.problems.push(first.with_file(&case));
        }
        sweep.files += 1;
    }
    Some(sweep)
}

#[test]
fn every_corpus_file_lexes_without_diagnostics() {
    let Some(sweep) = sweep_corpus("every_corpus_file_lexes_without_diagnostics") else {
        return;
    };
    assert_eq!(sweep.files, common::CORPUS_SIZE);
    // The Delphi lexer makes exactly this many tokens from the corpus, so
    // the two agree on every file's token count, not just on the total
    // being plausible. (Regenerating the corpus would move it.)
    assert_eq!(sweep.tokens, 273_129, "tokens over the whole corpus");
    assert!(
        sweep.problems.is_empty(),
        "{} problem(s):\n{}",
        sweep.problems.len(),
        sweep.problems.join("\n")
    );
}

#[test]
fn arrow_counts_match_the_corpus() {
    // An independent count of the corpus text finds 969 `->` and 2165 `=>`.
    // A lexer that mis-split an arrow would move these.
    let Some(sweep) = sweep_corpus("arrow_counts_match_the_corpus") else {
        return;
    };
    assert_eq!(sweep.arrows_rev, 969, "reversible \"->\" arrows");
    assert_eq!(sweep.arrows_irrev, 2165, "irreversible \"=>\" arrows");
}

#[test]
fn lexes_five_thousand_lines_in_under_fifty_milliseconds() {
    // The shapes the corpus contains, repeated: the largest corpus file is
    // only 435 lines.
    let mut source = String::new();
    for i in 1..=500 {
        let j = i + 1;
        writeln!(source, "  // Reactions:").unwrap();
        writeln!(
            source,
            "  reaction{i}: S{i} => S{j}; compartment_*k{i}*S{i};"
        )
        .unwrap();
        writeln!(source, "  S{i} = 0.00015/compartment_;").unwrap();
        writeln!(source, "  S{i} has substance_per_volume;").unwrap();
        writeln!(source, "  k{i} = 1.2e-3;").unwrap();
        writeln!(source, "  p{i} := piecewise(2, S{i} != 3, x^2 - 1);").unwrap();
        writeln!(source, "  E{i}: at (S{i} > 5, priority=1): S{i} = 0;").unwrap();
        writeln!(source, "  species G{i} in Cell, $T{i} in Cell;").unwrap();
        writeln!(source, "  const compartment_, k{i};").unwrap();
        writeln!(source, "  S{i} is \"Species number {i}\";").unwrap();
    }
    assert!(source.lines().count() >= 5000);

    let start = Instant::now();
    let (tokens, _) = tokenize(&source);
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;

    assert!(tokens.len() > 40_000, "only {} tokens", tokens.len());
    // The budget is for optimised code. An unoptimised test build is several
    // times slower, so it gets a looser bound that still catches anything
    // quadratic.
    let budget = if cfg!(debug_assertions) { 500.0 } else { 50.0 };
    assert!(
        elapsed < budget,
        "lexing {} lines ({} tokens) took {elapsed:.1} ms; the budget is {budget} ms",
        source.lines().count(),
        tokens.len()
    );
    eprintln!("lexed {} tokens in {elapsed:.1} ms", tokens.len());
}

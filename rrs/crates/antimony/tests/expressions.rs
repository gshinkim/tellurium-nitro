//! Ported from the Delphi `Tests/uAntExprTests.pas`.
//!
//! Not ported, because ownership makes the defect impossible to write:
//! `FailureLeaksNothing`, `CloneSurvivesTheOriginalBeingFreed`,
//! `SettingAChildFreesTheOldOne` and `LeafNodesAllocateNoArgumentList`
//! (a leaf has no argument list in the `enum`). `EqualsTreeDistinguishesOperators`
//! loses its `nil` case: there is no null `Expr`.
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart. The
//! corpus tests need `ANTIMONY_CORPUS` and are skipped without it.

mod common;

use antimony::Code;
use antimony::Diagnostics;
use antimony::builtins::{DOCUMENTED_FUNCTIONS, is_builtin_function};
use antimony::expr::{
    BinaryOp, Constant, Expr, ExprKind, ExprParser, MAX_DEPTH, MAX_NESTING, UnaryOp,
    parse_expression,
};
use antimony::lexer::{Keyword, Lexer, TokenKind, tokenize};

/// Parses, panicking with the diagnostics if it fails.
fn p(text: &str) -> Expr {
    let (expr, diagnostics) = parse_expression(text);
    expr.unwrap_or_else(|| panic!("{text:?} did not parse:\n{diagnostics}"))
}

/// Parses and writes back: the shorthand the precedence tests use.
fn r(text: &str) -> String {
    p(text).to_string()
}

/// Parses something that must fail.
fn fails(text: &str) -> Diagnostics {
    let (expr, diagnostics) = parse_expression(text);
    assert!(expr.is_none(), "{text:?} should not parse");
    assert!(diagnostics.has_errors(), "{text:?} failed without an error");
    diagnostics
}

fn eval(text: &str) -> Option<f64> {
    p(text).evaluate()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

fn assert_folds(text: &str, expected: f64) {
    match eval(text) {
        Some(value) => assert!(close(value, expected), "{text} = {value}, not {expected}"),
        None => panic!("{text} did not fold"),
    }
}

// --- TPrecedenceTests -------------------------------------------------------

#[test]
fn multiplication_binds_tighter_than_addition() {
    assert_eq!(r("a + b * c"), "a + b * c");
    assert_eq!(r("(a + b) * c"), "(a + b) * c");
}

#[test]
fn addition_is_left_associative() {
    assert_eq!(r("a + b + c"), "a + b + c");
    assert_eq!(r("a + (b + c)"), "a + (b + c)");
}

#[test]
fn subtraction_is_left_associative() {
    // The bracket must survive: a - (b - c) is not a - b - c.
    assert_eq!(r("a - b - c"), "a - b - c");
    assert_eq!(r("a - (b - c)"), "a - (b - c)");
}

#[test]
fn division_is_left_associative() {
    assert_eq!(r("a / b / c"), "a / b / c");
    assert_eq!(r("a / (b / c)"), "a / (b / c)");
    assert_eq!(r("a / (b * c)"), "a / (b * c)");
}

#[test]
fn power_is_right_associative() {
    assert_eq!(r("a ^ b ^ c"), "a ^ b ^ c");
    assert_eq!(r("(a ^ b) ^ c"), "(a ^ b) ^ c");
}

#[test]
fn power_binds_tighter_than_multiplication() {
    assert_eq!(r("a * b ^ c"), "a * b ^ c");
    assert_eq!(r("k1*S1^2*S2"), "k1 * S1 ^ 2 * S2");
}

#[test]
fn unary_minus_binds_tighter_than_power() {
    // Delphi §5.3 puts unary +/- ABOVE ^, so -2^2 is (-2)^2 = 4. The text
    // stays `-2 ^ 2` because no bracket is needed to give this tree back, so
    // the claim is checked on the tree.
    let e = p("-2^2");
    let ExprKind::Binary {
        op: BinaryOp::Power,
        left,
        ..
    } = &e.kind
    else {
        panic!("not a power: {e:?}");
    };
    assert!(
        matches!(left.kind, ExprKind::Unary { .. }),
        "the unary minus must be INSIDE the power, on its left"
    );
    assert_eq!(e.to_string(), "-2 ^ 2");
}

#[test]
fn power_accepts_a_negative_exponent() {
    // `^-` occurs 42 times in the corpus.
    assert_eq!(r("x^-1"), "x ^ -1");
    assert_eq!(r("S1^-2.5"), "S1 ^ -2.5");
}

#[test]
fn unary_minus_binds_tighter_than_addition() {
    assert_eq!(r("-a + b"), "-a + b");
    assert_eq!(r("-(a + b)"), "-(a + b)");
    assert_eq!(r("-1*k1 + S1"), "-1 * k1 + S1");
}

#[test]
fn double_unary_minus() {
    // Must not come out as `--a`, which reads back as the DNA operator.
    // Delphi wrote `-(-a)`; the space does the same job without adding a
    // level of nesting (SPEC.md §5.3).
    assert_eq!(r("- -a"), "- -a");
    assert_eq!(r("-(-a)"), "- -a");
    assert!(p("- -a").same_tree(&p("-(-a)")));
}

#[test]
fn comparisons_bind_looser_than_arithmetic() {
    assert_eq!(r("a + b < c * d"), "a + b < c * d");
}

#[test]
fn equality_binds_looser_than_comparison() {
    assert_eq!(r("a < b == c > d"), "a < b == c > d");
}

#[test]
fn not_sits_between_and_and_equality() {
    // `!` takes the comparison to its right, but not the conjunction.
    let e = p("!a == b");
    let ExprKind::Unary {
        op: UnaryOp::Not,
        operand,
    } = &e.kind
    else {
        panic!("! must be the root: {e:?}");
    };
    assert!(
        matches!(
            operand.kind,
            ExprKind::Binary {
                op: BinaryOp::Equal,
                ..
            }
        ),
        "the == must be INSIDE the !"
    );
    assert_eq!(e.to_string(), "!a == b");

    assert!(
        matches!(
            p("!a && b").kind,
            ExprKind::Binary {
                op: BinaryOp::And,
                ..
            }
        ),
        "&& must be the root"
    );
    assert_eq!(r("!!a"), "! !a");
}

#[test]
fn and_binds_tighter_than_or() {
    assert_eq!(r("a || b && c"), "a || b && c");
    assert_eq!(r("(a || b) && c"), "(a || b) && c");
}

#[test]
fn word_operators_match_their_symbols() {
    assert_eq!(r("a or b and c"), r("a || b && c"));
    assert_eq!(r("not a"), "!a");
}

#[test]
fn conditional_is_lowest_and_right_associative() {
    assert_eq!(r("a > 1 ? b : c"), "a > 1 ? b : c");
    assert_eq!(r("a ? b : c ? d : e"), "a ? b : c ? d : e");
    assert_eq!(r("a ? b : (c ? d : e)"), "a ? b : c ? d : e");
    assert_eq!(r("(a ? b : c) ? d : e"), "(a ? b : c) ? d : e");
}

#[test]
fn modulo_binds_with_multiplication() {
    assert_eq!(r("a + b % c"), "a + b % c");
    assert_eq!(r("a % b * c"), "a % b * c");
}

#[test]
fn redundant_parentheses_are_discarded() {
    assert_eq!(r("(a + (b * c))"), "a + b * c");
    assert_eq!(r("((((a))))"), "a");
}

#[test]
fn necessary_parentheses_are_kept() {
    assert_eq!(r("(a + b) * (c - d)"), "(a + b) * (c - d)");
}

// --- TExprParseTests --------------------------------------------------------

#[test]
fn numbers_and_symbols() {
    assert_eq!(
        p("42").kind,
        ExprKind::Number {
            value: 42.0,
            units: None
        }
    );
    assert_eq!(p("S1").kind, ExprKind::Symbol("S1".into()));
}

#[test]
fn named_constants_are_not_symbols() {
    // So that inlining a function can substitute them (Delphi §9.1.4 #4).
    assert_eq!(p("pi").kind, ExprKind::Constant(Constant::Pi));
    assert_eq!(p("false").kind, ExprKind::Constant(Constant::False));
    assert!(p("pi + true").symbols().is_empty());
}

#[test]
fn function_call() {
    let e = p("piecewise(2, x > 1, 3)");
    let ExprKind::Call { name, args } = &e.kind else {
        panic!("not a call: {e:?}");
    };
    assert_eq!((name.as_str(), args.len()), ("piecewise", 3));
    assert_eq!(e.to_string(), "piecewise(2, x > 1, 3)");
}

#[test]
fn zero_argument_call() {
    let e = p("f()");
    assert_eq!(
        e.kind,
        ExprKind::Call {
            name: "f".into(),
            args: vec![]
        }
    );
    assert_eq!(e.to_string(), "f()");
}

#[test]
fn nested_calls() {
    assert_eq!(r("abs(sin(x)-cos(y))"), "abs(sin(x) - cos(y))");
}

#[test]
fn piecewise_is_an_ordinary_call() {
    assert_eq!(r("piecewise(2, P1!= pi, 3)"), "piecewise(2, P1 != pi, 3)");
}

#[test]
fn word_operator_in_call_position_is_a_call() {
    // The corpus has both `and(a, b)` and `a && b`; position decides.
    assert!(matches!(p("and(a, b)").kind, ExprKind::Call { ref name, .. } if name == "and"));
    assert!(matches!(
        p("a and b").kind,
        ExprKind::Binary {
            op: BinaryOp::And,
            ..
        }
    ));
    assert!(matches!(p("not(a)").kind, ExprKind::Call { .. }));
}

#[test]
fn dotted_name_is_one_symbol() {
    let e = p("A.S1 * 2");
    let ExprKind::Binary { left, .. } = &e.kind else {
        panic!("{e:?}");
    };
    assert_eq!(left.kind, ExprKind::Symbol("A.S1".into()));
    assert_eq!(e.to_string(), "A.S1 * 2");
}

#[test]
fn bare_time_is_a_symbol() {
    // `time` occurs 732 times in the corpus, never as `time()`.
    assert_eq!(p("time").kind, ExprKind::Symbol("time".into()));
}

#[test]
fn numeric_literal_can_carry_a_unit() {
    let e = p("1 mole");
    assert_eq!(
        e.kind,
        ExprKind::Number {
            value: 1.0,
            units: Some("mole".into())
        }
    );
    assert_eq!(e.to_string(), "1 mole", "must round-trip");

    // The unit is exactly one identifier: in the reference's own example the
    // `/` divides two numbers that each carry a unit.
    let e = p("40 foo/3 seconds");
    let ExprKind::Binary {
        op: BinaryOp::Divide,
        left,
        right,
    } = &e.kind
    else {
        panic!("{e:?}");
    };
    assert_eq!(
        left.kind,
        ExprKind::Number {
            value: 40.0,
            units: Some("foo".into())
        }
    );
    assert_eq!(
        right.kind,
        ExprKind::Number {
            value: 3.0,
            units: Some("seconds".into())
        }
    );
    assert_eq!(e.to_string(), "40 foo / 3 seconds");
}

#[test]
fn unit_suffix_stops_at_the_first_identifier() {
    // Corpus cases 00928 and 00929: `S1 = 0.00015 mole/C`, where C is the
    // compartment. Taking `/C` into the unit would invent a unit "mole/C"
    // and lose the division.
    let e = p("0.00015 mole/C");
    let ExprKind::Binary {
        op: BinaryOp::Divide,
        left,
        right,
    } = &e.kind
    else {
        panic!("the / must divide: {e:?}");
    };
    assert_eq!(
        left.kind,
        ExprKind::Number {
            value: 0.00015,
            units: Some("mole".into())
        },
        "the unit is \"mole\", not \"mole/C\""
    );
    assert_eq!(right.kind, ExprKind::Symbol("C".into()));
    assert_eq!(e.symbols(), ["C"], "C is a symbol to resolve");
}

#[test]
fn unit_suffix_does_not_eat_ordinary_arithmetic() {
    let e = p("2 * x");
    let ExprKind::Binary { left, .. } = &e.kind else {
        panic!("{e:?}");
    };
    assert!(matches!(left.kind, ExprKind::Number { units: None, .. }));
    assert_eq!(e.to_string(), "2 * x");
    assert_eq!(r("2 + x"), "2 + x");
    assert_eq!(r("k1*S1"), "k1 * S1");
}

#[test]
fn unit_suffix_can_be_turned_off_for_stoichiometry() {
    // `2 S1 + S2 -> S3`: in a participant list `2 S1` is a coefficient and a
    // species, so the statement parser turns units off there.
    let mut lexer = Lexer::new("2 S1");
    let mut parser = ExprParser::new(&mut lexer);
    parser.parse_unit_suffixes = false;
    let e = parser.parse().unwrap();
    assert_eq!(
        e.kind,
        ExprKind::Number {
            value: 2.0,
            units: None
        },
        "S1 is a species, not a unit"
    );
    assert_eq!(lexer.peek(0).text, "S1", "the species must still be there");
}

#[test]
fn kinetic_law_from_the_corpus() {
    assert_eq!(r("compartment_*k1*S1"), "compartment_ * k1 * S1");
    assert_eq!(r("0.00015/compartment_"), "0.00015 / compartment_");
}

#[test]
fn spans_cover_the_whole_expression() {
    let text = "k1 * S1";
    let e = p(text);
    assert_eq!(e.span.start.column, 1);
    assert_eq!(e.span.slice_of(text), text);
}

#[test]
fn empty_input_is_an_error_not_a_zero() {
    // `p = ;` is real corpus text, but deciding that a right-hand side is
    // empty is the statement parser's job; this parser refuses.
    fails("");
}

#[test]
fn unbalanced_paren_is_reported() {
    assert!(fails("(a + b").has_any(Code::UnexpectedToken));
}

#[test]
fn trailing_rubbish_is_reported() {
    fails("a + b )");
}

#[test]
fn dna_operator_is_rejected_by_name() {
    assert!(fails("g1--g2").has_any(Code::DnaNotSupported));
}

#[test]
fn missing_operand_is_reported() {
    fails("a *");
}

// --- TExprNodeTests ---------------------------------------------------------

#[test]
fn clone_is_a_deep_copy() {
    let e = p("k1 * S1 ^ 2 + f(x, y)");
    let copy = e.clone();
    assert!(e.same_tree(&copy));
    assert_eq!(e, copy);
    assert_eq!(e.to_string(), copy.to_string());
}

#[test]
fn unary_uses_operand_not_right() {
    let ExprKind::Unary { operand, .. } = p("-x").kind else {
        panic!("not unary");
    };
    assert_eq!(operand.kind, ExprKind::Symbol("x".into()));
}

#[test]
fn equals_tree_ignores_spans() {
    let a = p("a + b");
    let b = p("   a   +   b   ");
    assert!(a.same_tree(&b));
    assert_ne!(a, b, "== compares spans too");
}

#[test]
fn equals_tree_distinguishes_operators() {
    assert!(!p("a + b").same_tree(&p("a - b")));
}

#[test]
fn nan_literals_compare_equal() {
    // NaN is not equal to itself as a number, but two `nan` literals are the
    // same tree, or no model holding one could pass the round-trip test.
    assert!(p("nan").same_tree(&p("nan")));
}

#[test]
fn collect_symbols_dedupes_in_order() {
    assert_eq!(p("k1 * S1 + k1 * S2 - S1").symbols(), ["k1", "S1", "S2"]);
}

#[test]
fn collect_symbols_skips_constants_and_function_names() {
    // A function name taken for a symbol would become an invented parameter.
    assert_eq!(p("sin(x) + pi * abs(y)").symbols(), ["x", "y"]);
}

#[test]
fn node_count_and_depth() {
    assert_eq!((p("a + b").node_count(), p("a + b").depth()), (3, 2));
    assert_eq!((p("a").node_count(), p("a").depth()), (1, 1));
}

// --- TExprEvaluateTests -----------------------------------------------------

#[test]
fn arithmetic_folds() {
    assert_folds("2 + 3 * 4", 14.0);
}

#[test]
fn precedence_is_respected_when_folding() {
    assert_folds("(2 + 3) * 4", 20.0);
    assert_folds("2 ^ 3 ^ 2", 512.0); // right-associative: 2^(3^2)
}

#[test]
fn unary_minus_power_evaluates_as_specified() {
    // The arithmetic consequence of Delphi §5.3's ordering, pinned as a number.
    assert_folds("-2^2", 4.0);
}

#[test]
fn comparisons_yield_one_and_zero() {
    assert_folds("3 > 2", 1.0);
    assert_folds("3 < 2", 0.0);
}

#[test]
fn named_constants_evaluate() {
    assert_folds("pi", std::f64::consts::PI);
    assert_folds("true", 1.0);
}

#[test]
fn builtin_functions_evaluate() {
    assert_folds("sqrt(16)", 4.0);
    assert_folds("abs(0 - 3)", 3.0);
    assert_folds("pow(2, 10)", 1024.0);
}

#[test]
fn ceiling_does_not_overflow_on_large_values() {
    assert_folds("ceil(6.022e23)", 6.022e23);
    assert_folds("floor(2.7)", 2.0);
    assert_folds("ceil(2.1)", 3.0);
}

#[test]
fn unresolved_symbol_fails_rather_than_returning_zero() {
    // A typo that folded to zero would give a model that runs and is wrong.
    assert_eq!(eval("k1 * 2"), None);
}

#[test]
fn resolver_supplies_symbol_values() {
    let value = p("k1 * 2 + S1").evaluate_with(&|name| match name {
        "k1" => Some(5.0),
        "S1" => Some(1.0),
        _ => None,
    });
    assert_eq!(value, Some(11.0));
}

#[test]
fn division_by_zero_fails() {
    assert_eq!(eval("1 / 0"), None);
}

#[test]
fn conjunction_short_circuits() {
    assert_folds("0 && unknownSymbol", 0.0);
    assert_folds("1 || unknownSymbol", 1.0);
    assert_eq!(eval("1 && unknownSymbol"), None);
}

#[test]
fn piecewise_selects_lazily() {
    // A branch that is not chosen need not be computable.
    assert_folds("piecewise(1, 1 > 0, unknownSymbol)", 1.0);
    assert_folds("piecewise(unknownSymbol, 0 > 1, 9)", 9.0);
}

#[test]
fn is_constant_distinguishes_value_from_assignment() {
    // `x = 3` is a value and `x = 2*k` an initial assignment.
    assert!(p("2 * 3 + pi").is_constant());
    assert!(!p("2 * k").is_constant());
    assert!(!p("sin(2)").is_constant(), "a call may not be pure");
}

#[test]
fn unknown_function_is_not_evaluatable() {
    assert_eq!(eval("myUserFunction(1, 2)"), None);
}

#[test]
fn logical_operators_fold_in_both_spellings() {
    // and, or and not are both operators and functions, and the corpus uses
    // both. xor is a function only.
    assert_folds("1 && 1", 1.0);
    assert_folds("and(1, 1)", 1.0);
    assert_folds("and(1, 0)", 0.0);
    assert_folds("0 || 1", 1.0);
    assert_folds("or(0, 1)", 1.0);
    assert_folds("or(0, 0)", 0.0);
    assert_folds("!0", 1.0);
    assert_folds("not(0)", 1.0);
    assert_folds("xor(1, 0)", 1.0);
    assert_folds("xor(1, 1)", 0.0);
    assert_folds("implies(0, 0)", 1.0); // false implies anything
    assert_folds("implies(1, 0)", 0.0);
    // An argument that cannot be computed still defeats the call.
    assert_eq!(eval("and(1, unknownSymbol)"), None);
}

#[test]
fn every_documented_function_is_recognised() {
    let missing: Vec<&str> = DOCUMENTED_FUNCTIONS
        .iter()
        .copied()
        .filter(|f| !is_builtin_function(f))
        .collect();
    assert!(missing.is_empty(), "not recognised: {missing:?}");
    assert_eq!(DOCUMENTED_FUNCTIONS.len(), 48);
    // The six Delphi §5.3 left out.
    for f in ["ceiling", "coth", "csch", "sech", "divide", "power"] {
        assert!(is_builtin_function(f), "{f}");
    }
}

#[test]
fn every_documented_function_evaluates() {
    // One in-domain call per documented name. piecewise has its own test.
    let calls = [
        "abs(-2)",
        "and(1, 1)",
        "arccos(0.5)",
        "arccosh(2)",
        "arccot(2)",
        "arccoth(2)",
        "arccsc(2)",
        "arccsch(2)",
        "arcsec(2)",
        "arcsech(0.5)",
        "arcsin(0.5)",
        "arcsinh(2)",
        "arctan(0.5)",
        "arctanh(0.5)",
        "ceiling(2.1)",
        "cos(1)",
        "cosh(1)",
        "cot(1)",
        "coth(1)",
        "csc(1)",
        "csch(1)",
        "divide(6, 3)",
        "eq(1, 1)",
        "exp(1)",
        "factorial(5)",
        "floor(2.7)",
        "geq(2, 1)",
        "gt(2, 1)",
        "leq(1, 2)",
        "ln(1)",
        "log(100)",
        "lt(1, 2)",
        "minus(5, 3)",
        "neq(1, 2)",
        "not(0)",
        "or(0, 1)",
        "plus(1, 2)",
        "power(2, 3)",
        "root(2, 9)",
        "sec(1)",
        "sech(1)",
        "sin(1)",
        "sinh(1)",
        "tan(1)",
        "tanh(1)",
        "times(2, 3)",
        "xor(1, 0)",
    ];
    let unfolded: Vec<&str> = calls.into_iter().filter(|c| eval(c).is_none()).collect();
    assert!(unfolded.is_empty(), "did not fold: {unfolded:?}");

    assert_folds("factorial(5)", 120.0);
    assert_folds("divide(6, 3)", 2.0);
    assert_folds("power(2, 3)", 8.0);
    assert_folds("minus(5, 3)", 2.0);
    assert_folds("minus(5)", -5.0); // MathML minus is unary as well as binary
    assert_folds("ceiling(2.1)", 3.0);
    assert_folds("log(100)", 2.0); // a bare log is base 10, as libSBML reads it
    assert_folds("log(2, 8)", 3.0); // MathML log takes the base first
}

#[test]
fn reciprocal_functions_fold() {
    assert_folds("sec(0)", 1.0);
    assert_folds("coth(1)", 1.0 / 1f64.tanh());
    assert_eq!(eval("csc(0)"), None, "csc(0) is 1/0");
    assert_eq!(eval("csch(0)"), None);
}

#[test]
fn inverse_functions_fold() {
    // libantimony writes `asin` where the documentation says `arcsin`.
    assert_folds("asin(0.5)", 0.5f64.asin());
    assert_folds("arcsin(0.5)", 0.5f64.asin());
    assert_folds("arcsec(2)", 0.5f64.acos()); // arcsec(x) = arccos(1/x)
    assert_eq!(eval("arcsec(0)"), None);
}

#[test]
fn n_ary_math_ml_operators_fold() {
    assert_folds("plus(1, 2, 3, 4)", 10.0);
    assert_folds("times(2, 3, 4)", 24.0);
    assert_folds("xor(1, 1, 1)", 1.0); // n-ary xor is parity
    assert_folds("min(3, 1, 2)", 1.0);
    assert_folds("max(3, 1, 2)", 3.0);
}

// --- Added in the Rust port -------------------------------------------------

#[test]
fn bracketed_span_includes_both_brackets() {
    let text = "(a + b) * c";
    let ExprKind::Binary { left, .. } = p(text).kind else {
        panic!("not binary");
    };
    assert_eq!(left.span.slice_of(text), "(a + b)");
}

#[test]
fn unit_is_not_a_keyword_or_a_call() {
    // `at 5.7 after …`: the keyword is the event's, not a unit.
    let mut lexer = Lexer::new("5.7 after x");
    let e = ExprParser::new(&mut lexer).parse().unwrap();
    assert!(matches!(e.kind, ExprKind::Number { units: None, .. }));
    assert!(lexer.peek(0).is_keyword(Keyword::After));
    // `2 f(x)` is not 2 in units of f.
    let mut lexer = Lexer::new("2 f(x)");
    let e = ExprParser::new(&mut lexer).parse().unwrap();
    assert!(matches!(e.kind, ExprKind::Number { units: None, .. }));
}

#[test]
fn parser_stops_after_the_expression() {
    // The statement parser carries on from where the expression ends.
    let mut lexer = Lexer::new("k1*S1; S2 = 3");
    let e = ExprParser::new(&mut lexer).parse().unwrap();
    assert_eq!(e.to_string(), "k1 * S1");
    assert_eq!(lexer.peek(0).kind, TokenKind::Semicolon);
}

/// Expressions nested `n` deep in each way the limits count, as
/// (description, text, how many levels the limit sees).
fn deep_shapes(n: usize) -> Vec<(&'static str, String, usize)> {
    // The whole expression is one level of nesting, and each bracket, call,
    // unary operator, `^` or `?:` inside it is another. A tree's depth counts
    // its nodes from root to deepest leaf.
    vec![
        (
            "brackets",
            format!("{}a{}", "(".repeat(n - 1), ")".repeat(n - 1)),
            MAX_NESTING,
        ),
        (
            "calls",
            format!("{}a{}", "f(".repeat(n - 1), ")".repeat(n - 1)),
            MAX_NESTING,
        ),
        (
            "unary minus",
            format!("{}1", "- ".repeat(n - 1)),
            MAX_NESTING,
        ),
        ("not", format!("{}a", "!".repeat(n - 1)), MAX_NESTING),
        ("powers", format!("2{}", "^1".repeat(n - 1)), MAX_NESTING),
        (
            "conditionals",
            format!("{}1", "a ? 1 : ".repeat(n - 1)),
            MAX_NESTING,
        ),
        ("a sum", vec!["1"; n].join(" + "), MAX_DEPTH),
        ("a difference", vec!["1"; n].join(" - "), MAX_DEPTH),
    ]
}

#[test]
fn deep_nesting_is_an_error_not_a_crash() {
    for (shape, _, limit) in deep_shapes(1) {
        let at_limit = deep_shapes(limit)
            .into_iter()
            .find(|s| s.0 == shape)
            .unwrap()
            .1;
        let past = deep_shapes(limit + 1)
            .into_iter()
            .find(|s| s.0 == shape)
            .unwrap()
            .1;
        let (expr, d) = parse_expression(&at_limit);
        assert!(expr.is_some(), "{shape} at the limit should parse: {d}");
        let d = fails(&past);
        assert_eq!(
            d.len(),
            1,
            "{shape}: one diagnostic, not one per level: {d}"
        );
        assert!(d[0].message.contains("nested too deeply"), "{shape}: {d}");
    }
    // Far past the limit is no different.
    fails(&format!("{}a", "(".repeat(100_000)));
    fails(&vec!["1"; 100_000].join("+"));
}

#[test]
fn the_deepest_allowed_expressions_fit_in_a_small_stack() {
    // 1 MB is the smallest stack this library runs on (a Windows main
    // thread, WebAssembly). Everything that walks a tree recursively must
    // cope with the deepest trees the parser accepts, in an unoptimised
    // build, which uses the most stack.
    let worker = std::thread::Builder::new()
        .stack_size(1 << 20)
        .spawn(|| {
            let mut inputs: Vec<String> = deep_shapes(1)
                .into_iter()
                .map(|(shape, _, limit)| {
                    deep_shapes(limit)
                        .into_iter()
                        .find(|s| s.0 == shape)
                        .unwrap()
                        .1
                })
                .collect();
            // A long sum at the bottom of deep brackets.
            let n = MAX_NESTING - 1;
            inputs.push(format!(
                "{}{}{}",
                "(".repeat(n),
                vec!["1"; MAX_DEPTH - n].join(" * "),
                ")".repeat(n)
            ));
            for text in &inputs {
                let e = p(text);
                let again = p(&e.to_string());
                assert!(e.same_tree(&again), "{text}");
                assert_eq!(e, e.clone());
                let _ = (
                    e.evaluate(),
                    e.symbols(),
                    e.depth(),
                    e.node_count(),
                    e.is_constant(),
                );
            }
        })
        .unwrap();
    worker.join().expect("ran out of stack");
}

// --- TExprCorpusTests -------------------------------------------------------

/// The expressions in a source, found without a statement parser: the
/// right-hand side of each `=` and `:=`, up to a top-level `;`, `,`, line
/// break or unmatched `:`. Kinetic laws are not found this way. `unit`
/// lines are skipped: a unit definition is not an ordinary expression.
fn extract_expressions(source: &str) -> Vec<String> {
    let (tokens, _) = tokenize(source);
    let mut found = Vec::new();
    let mut unit_line = false;
    for (i, token) in tokens.iter().enumerate() {
        if token.kind == TokenKind::Newline {
            unit_line = false;
            continue;
        }
        if token.is_keyword(Keyword::Unit) {
            unit_line = true;
            continue;
        }
        if unit_line || !matches!(token.kind, TokenKind::Assign | TokenKind::ColonAssign) {
            continue;
        }
        let start = i + 1;
        let (mut depth, mut pending, mut j) = (0i32, 0, start);
        while j < tokens.len() {
            match tokens[j].kind {
                TokenKind::LParen => depth += 1,
                TokenKind::RParen => depth -= 1,
                TokenKind::Question => pending += 1,
                _ => {}
            }
            if depth < 0 {
                break;
            }
            if depth == 0
                && matches!(
                    tokens[j].kind,
                    TokenKind::Semicolon | TokenKind::Comma | TokenKind::Newline | TokenKind::Eof
                )
            {
                break;
            }
            // A top-level `:` ends an event modifier's value; one belonging
            // to `?:` does not.
            if depth == 0 && tokens[j].kind == TokenKind::Colon {
                if pending == 0 {
                    break;
                }
                pending -= 1;
            }
            j += 1;
        }
        if j > start {
            let text =
                source[tokens[start].span.start.offset..tokens[j - 1].span.end.offset].trim();
            if !text.is_empty() {
                found.push(text.to_string());
            }
        }
    }
    found
}

/// Every expression in the corpus, as (case, text).
fn corpus_expressions(test: &str) -> Option<Vec<(String, String)>> {
    let root = common::corpus_root(test)?;
    let mut all = Vec::new();
    for (case, path) in common::corpus_files(&root) {
        for text in extract_expressions(&common::read(&path)) {
            all.push((case.clone(), text));
        }
    }
    Some(all)
}

#[test]
fn every_corpus_expression_parses() {
    let Some(all) = corpus_expressions("every_corpus_expression_parses") else {
        return;
    };
    assert!(all.len() > 10_000, "only {} expressions found", all.len());
    let failed: Vec<String> = all
        .iter()
        .filter_map(|(case, text)| {
            let (expr, d) = parse_expression(text);
            (expr.is_none() || d.has_errors()).then(|| format!("{case}: {text:?}\n{d}"))
        })
        .collect();
    assert!(
        failed.is_empty(),
        "{} of {} failed. First: {}",
        failed.len(),
        all.len(),
        failed[0]
    );
}

#[test]
fn every_corpus_expression_round_trips() {
    // ToString(Parse(s)) must read back as the same tree.
    let Some(all) = corpus_expressions("every_corpus_expression_round_trips") else {
        return;
    };
    let mut failed = Vec::new();
    for (case, text) in &all {
        let Some(first) = parse_expression(text).0 else {
            continue;
        };
        let written = first.to_string();
        let ok = parse_expression(&written)
            .0
            .is_some_and(|second| first.same_tree(&second));
        if !ok {
            failed.push(format!("{case}: {text:?} written as {written:?}"));
        }
    }
    assert!(
        failed.is_empty(),
        "{} failed. First: {}",
        failed.len(),
        failed[0]
    );
}

#[test]
fn rendering_is_idempotent() {
    // Writing what was read, then reading and writing again, gives the same
    // text: brackets do not pile up.
    let Some(all) = corpus_expressions("rendering_is_idempotent") else {
        return;
    };
    let mut failed = Vec::new();
    for (case, text) in &all {
        let Some(first) = parse_expression(text).0 else {
            continue;
        };
        let once = first.to_string();
        let Some(second) = parse_expression(&once).0 else {
            continue;
        };
        let twice = second.to_string();
        if once != twice {
            failed.push(format!("{case}: {once:?} then {twice:?}"));
        }
    }
    assert!(
        failed.is_empty(),
        "{} grew. First: {}",
        failed.len(),
        failed[0]
    );
}

fn unknown_calls(e: &Expr, into: &mut Vec<String>) {
    if let ExprKind::Call { name, .. } = &e.kind {
        if !is_builtin_function(name) && !into.contains(name) {
            into.push(name.clone());
        }
    }
    for child in e.children() {
        unknown_calls(child, into);
    }
}

#[test]
fn no_corpus_call_is_an_unknown_function() {
    // A built-in missing from the list would become an invented parameter.
    let Some(all) = corpus_expressions("no_corpus_call_is_an_unknown_function") else {
        return;
    };
    let mut unknown = Vec::new();
    for (_, text) in &all {
        if let Some(e) = parse_expression(text).0 {
            unknown_calls(&e, &mut unknown);
        }
    }
    assert!(unknown.is_empty(), "not built-ins: {unknown:?}");
}

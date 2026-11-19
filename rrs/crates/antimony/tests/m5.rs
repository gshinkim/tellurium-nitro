//! Ported from the Delphi `Tests/uAntM5Tests.pas`: functions, events,
//! algebraic rules, constraints and annotations (Delphi M5 and M9).
//!
//! Changed for Rust: `EventReferencesResolve` checks that the assignment's
//! target name resolves to the symbol, since references are names (SPEC.md
//! §3.1). The constraint tests read `ConstraintKind` rather than Delphi's
//! nullable `Expr`/`Lower`/`Subject`/`Upper`.
//!
//! The corpus tests need `ANTIMONY_CORPUS` and `ANTIMONY_CORPUS_INDEX`.

mod common;

use antimony::expr::{BinaryOp, Expr, ExprKind, parse_expression};
use antimony::model::{Constraint, ConstraintKind, Event, Module, Qualifier, RawKind, SymbolKind};
use antimony::parser::parse_document;
use antimony::semantics::resolve;
use antimony::{Code, Diagnostics};

/// A parse, optionally analysed: the main module and the diagnostics.
struct Read {
    document: antimony::model::Document,
    diagnostics: Diagnostics,
}

impl Read {
    fn main(&self) -> &Module {
        self.document.main_module()
    }

    fn ev(&self, i: usize) -> &Event {
        self.main()
            .events
            .get(i)
            .unwrap_or_else(|| panic!("expected at least {} event(s)", i + 1))
    }

    fn constraint(&self, i: usize) -> &Constraint {
        &self.main().constraints[i]
    }

    fn has(&self, code: Code) -> bool {
        self.diagnostics.has_any(code)
    }

    fn assert_clean(&self) -> &Self {
        assert!(!self.diagnostics.has_errors(), "{}", self.diagnostics);
        self
    }
}

/// Parse only.
fn p(source: &str) -> Read {
    let (document, diagnostics) = parse_document(source);
    Read {
        document,
        diagnostics,
    }
}

/// Parse and analyse.
fn a(source: &str) -> Read {
    let loaded = antimony::parse(source);
    Read {
        document: loaded.document,
        diagnostics: loaded.diagnostics,
    }
}

fn text(e: &Option<Expr>) -> String {
    e.as_ref().expect("no expression").to_string()
}

// --- TEventTests ---------------------------------------------------------------

#[test]
fn simple_event() {
    let r = p("event1: at S1 < 0.1: S1 = 1;");
    r.assert_clean();
    let e = r.ev(0);
    assert_eq!(e.name.as_deref(), Some("event1"));
    assert_eq!(e.trigger.to_string(), "S1 < 0.1");
    assert_eq!(e.assignments.len(), 1);
    assert_eq!(e.assignments[0].target, "S1");
    assert_eq!(text(&e.assignments[0].expr), "1");
}

#[test]
fn unnamed_event() {
    let r = p("at (x > 5): y = 3, x = r + 2;");
    r.assert_clean();
    assert_eq!(r.ev(0).name, None);
    assert_eq!(r.ev(0).assignments.len(), 2);
}

#[test]
fn parenthesised_trigger() {
    // 200 corpus cases use the compound bracketed form.
    let r = p("E0: at (trig > 0) && (A0 < 1): B0 = 2;");
    r.assert_clean();
    assert_eq!(r.ev(0).trigger.to_string(), "trig > 0 && A0 < 1");
}

#[test]
fn delay_is_written_before_the_trigger() {
    // 328 corpus cases: the delay comes first, joined by `after`.
    let r = p("event1: at 1 after S1 < 0.1: S1 = 1;");
    r.assert_clean();
    assert_eq!(text(&r.ev(0).delay), "1");
    assert_eq!(r.ev(0).trigger.to_string(), "S1 < 0.1");
}

#[test]
fn delay_may_be_an_expression() {
    let r = p("E1: at 1*(5.3 + exp(exponentiale)) after trig > -5.2: P1 = 2;");
    r.assert_clean();
    assert_eq!(text(&r.ev(0).delay), "1 * (5.3 + exp(exponentiale))");
    assert_eq!(r.ev(0).trigger.to_string(), "trig > -5.2");
}

#[test]
fn all_five_modifiers() {
    // Every one occurs in the corpus.
    let r = p(
        "E0: at x > 1, priority = 5, t0 = false, persistent = false, fromTrigger = true, delay = 2: B0 = 2;",
    );
    r.assert_clean();
    let e = r.ev(0);
    assert_eq!(text(&e.priority), "5");
    assert_eq!(text(&e.initial_value), "false", "t0");
    assert_eq!(text(&e.persistent), "false");
    assert_eq!(text(&e.from_trigger), "true");
    assert_eq!(text(&e.delay), "2");
}

#[test]
fn multiple_assignments() {
    let r = p("event1: at 1 after S1 < 0.1: S2 = 0, S1 = 1;");
    r.assert_clean();
    let targets: Vec<&str> = r
        .ev(0)
        .assignments
        .iter()
        .map(|a| a.target.as_str())
        .collect();
    assert_eq!(targets, ["S2", "S1"]);
}

#[test]
fn empty_assignment_list_is_legal() {
    // A challenge case: libantimony writes it and cannot read it back.
    let r = p("E0: at time > 5.5: ;");
    r.assert_clean();
    assert!(r.ev(0).assignments.is_empty());
    assert_eq!(r.ev(0).trigger.to_string(), "time > 5.5");
}

#[test]
fn event_name_becomes_an_event_symbol() {
    let r = p("E1: at x > 1: y = 2;");
    r.assert_clean();
    assert_eq!(r.main().symbol("E1").unwrap().kind, SymbolKind::Event);
}

#[test]
fn event_references_resolve() {
    let r = a("x = 0; y = 0; r = 1;\nE1: at x > 5: y = r + 2;");
    r.assert_clean();
    assert!(!r.has(Code::ImplicitParameter), "{}", r.diagnostics);
    let target = &r.ev(0).assignments[0].target;
    assert_eq!(
        resolve(r.main(), target).map(|s| s.name.as_str()),
        Some("y")
    );
}

#[test]
fn unit_suffix_does_not_eat_the_after_keyword() {
    // Read as a unit, the delay would be "5.7 in units of after" and the
    // event would lose its trigger. A keyword is never a unit.
    let r = p("event1: at 5.7 after S1 < 0.75: S2 = 1;");
    r.assert_clean();
    let delay = r.ev(0).delay.as_ref().unwrap();
    assert_eq!(
        delay.kind,
        ExprKind::Number {
            value: 5.7,
            units: None
        }
    );
    assert_eq!(r.ev(0).trigger.to_string(), "S1 < 0.75");
}

// --- TFunctionTests --------------------------------------------------------------
//
// The corpus has no `function … end` blocks (libantimony inlines them), so
// these fixtures are the only coverage functions have.

#[test]
fn function_definition() {
    let r = p("function quadratic(x, a, b, c)\n  a*x^2 + b*x + c\nend");
    r.assert_clean();
    let functions = &r.main().functions;
    assert_eq!(functions.len(), 1);
    assert_eq!(functions[0].name, "quadratic");
    assert_eq!(functions[0].args.len(), 4);
    assert_eq!(text(&functions[0].body), "a * x ^ 2 + b * x + c");
}

#[test]
fn no_argument_function() {
    let r = p("function nothing()\n  42\nend");
    r.assert_clean();
    assert!(r.main().functions[0].args.is_empty());
}

#[test]
fn formal_arguments_are_local_to_the_body() {
    // The `bvar` leak (Delphi §9.1.4 #4): libantimony lets a function's
    // arguments escape as invented parameters.
    let r = a("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend");
    r.assert_clean();
    for name in ["S", "Vm", "Km"] {
        assert!(r.main().symbol(name).is_none(), "{name} must not leak");
    }
    assert!(!r.has(Code::ImplicitParameter));
}

#[test]
fn duplicate_formal_argument_is_reported() {
    assert!(p("function f(a, a)\n  a\nend").has(Code::DuplicateFormalArg));
}

#[test]
fn missing_end_is_reported() {
    assert!(p("function f(a)\n  a").has(Code::ExpectedEnd));
}

#[test]
fn call_with_correct_arity_is_silent() {
    let r = a("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend\nv = MM(1, 2, 3);");
    r.assert_clean();
    assert!(!r.has(Code::ArityMismatch));
}

#[test]
fn call_with_wrong_arity_is_reported() {
    let r = a("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend\nv = MM(1, 2);");
    assert!(r.has(Code::ArityMismatch), "{}", r.diagnostics);
    assert!(
        r.diagnostics
            .first_of(Code::ArityMismatch)
            .unwrap()
            .message
            .contains("MM")
    );
}

#[test]
fn builtins_are_exempt_from_arity_checking() {
    // piecewise is variadic.
    let r = a("v = piecewise(1, x > 0, 2, y > 0, 3);\nx = 1; y = 1;");
    assert!(!r.has(Code::ArityMismatch), "{}", r.diagnostics);
}

#[test]
fn recursive_function_is_reported() {
    // SBML function definitions cannot recurse at all.
    let r = a("function fact(n)\n  n*fact(n - 1)\nend");
    assert!(r.has(Code::RecursiveFunction), "{}", r.diagnostics);
}

#[test]
fn function_name_becomes_a_function_symbol() {
    let r = p("function f(a)\n  a\nend");
    r.assert_clean();
    assert_eq!(r.main().symbol("f").unwrap().kind, SymbolKind::Function);
}

// --- TAlgebraicRuleTests -----------------------------------------------------------

#[test]
fn named_algebraic_rule() {
    let r = p("_alg0: 0 = -1*k1 + S1 + S2;");
    r.assert_clean();
    let rules = &r.main().algebraic_rules;
    assert_eq!(rules.len(), 1);
    assert_eq!(rules[0].name, "_alg0");
    assert_eq!(text(&rules[0].expr), "-1 * k1 + S1 + S2");
}

#[test]
fn empty_algebraic_rule_is_legal() {
    // Case 01244 writes `_alg0: 0 = ;`.
    let r = p("_alg0: 0 = ;");
    r.assert_clean();
    assert_eq!(r.main().algebraic_rules.len(), 1);
    assert!(r.main().algebraic_rules[0].expr.is_none());
}

#[test]
fn rule_name_is_not_a_symbol() {
    let r = p("_alg0: 0 = X - p1;");
    r.assert_clean();
    assert!(
        r.main().symbol("_alg0").is_none(),
        "the label is not a symbol, as in libantimony"
    );
}

#[test]
fn rule_references_resolve() {
    let r = a("X = 1; p1 = 2;\n_alg0: 0 = X - p1;");
    r.assert_clean();
    assert!(!r.has(Code::ImplicitParameter), "{}", r.diagnostics);
}

// --- TConstraintTests ----------------------------------------------------------------

#[test]
fn named_constraint() {
    let r = p("constraint c0: R16 >= 0");
    r.assert_clean();
    assert_eq!(r.main().constraints.len(), 1);
    let c = r.constraint(0);
    assert_eq!(c.name.as_deref(), Some("c0"));
    assert!(matches!(&c.kind, ConstraintKind::Expr(e) if e.to_string() == "R16 >= 0"));
    assert!(!c.is_range());
}

#[test]
fn bare_flux_bound() {
    // The language reference's flux-balance constraints: no keyword.
    let r = p("R19 <= 1000");
    r.assert_clean();
    assert_eq!(r.main().constraints.len(), 1);
    assert!(
        matches!(&r.constraint(0).kind, ConstraintKind::Expr(e) if e.to_string() == "R19 <= 1000")
    );

    let r = p("0 <= J0");
    r.assert_clean();
    assert_eq!(r.main().constraints.len(), 1);
}

#[test]
fn chained_comparison_is_a_range() {
    // Decision 18: the expression parser builds ((-10 <= J2) <= 10); the
    // constraint takes it apart.
    let r = p("constraint b: -10 <= J2 <= 10");
    r.assert_clean();
    let ConstraintKind::Range {
        lower,
        subject,
        upper,
        ..
    } = &r.constraint(0).kind
    else {
        panic!("not a range: {:?}", r.constraint(0));
    };
    assert_eq!(lower.to_string(), "-10");
    assert_eq!(subject.to_string(), "J2");
    assert_eq!(upper.to_string(), "10");
}

#[test]
fn chained_comparison_keeps_strictness() {
    let r = p("constraint b: 0 < J2 <= 10");
    r.assert_clean();
    let ConstraintKind::Range {
        lower_inclusive,
        upper_inclusive,
        ..
    } = r.constraint(0).kind
    else {
        panic!("not a range");
    };
    assert!(!lower_inclusive);
    assert!(upper_inclusive);
}

#[test]
fn chained_comparison_elsewhere_stays_left_associative() {
    // The range reading is for constraints only.
    let e = parse_expression("a <= b <= c").0.unwrap();
    assert_eq!(e.to_string(), "a <= b <= c");
    let ExprKind::Binary { left, .. } = &e.kind else {
        panic!("{e:?}");
    };
    assert!(
        matches!(
            left.kind,
            ExprKind::Binary {
                op: BinaryOp::LessEqual,
                ..
            }
        ),
        "the left child is itself the comparison: ((a <= b) <= c)"
    );
}

#[test]
fn constraint_references_resolve() {
    let r = a("R16: A => B; k*A;\nA = 1; B = 0; k = 1;\nconstraint c0: R16 >= 0");
    r.assert_clean();
    assert!(!r.has(Code::ImplicitParameter), "{}", r.diagnostics);
}

// --- TAnnotationTests ------------------------------------------------------------------

#[test]
fn model_history_annotations_load() {
    let r = p(
        "model *m()\n  model created \"2012-12-21T11:05:21Z\"\n  model modified \"2012-12-21T10:05:21Z\"\nend",
    );
    r.assert_clean();
    let metadata = &r.main().metadata;
    assert_eq!(metadata.created.as_deref(), Some("2012-12-21T11:05:21Z"));
    assert_eq!(metadata.modified, ["2012-12-21T10:05:21Z"]);
    assert!(
        r.main().unparsed_sections.is_empty(),
        "a recognised annotation is parsed, not kept as text"
    );
}

#[test]
fn dotted_annotation_qualifier() {
    let r = p(
        "model *m()\n  model creator1.givenName \"Brett\"\n  model creator1.email \"bgoli@users.sourceforge.net\"\nend",
    );
    r.assert_clean();
    // Both fields belong to the same creator: the number groups them.
    let creators = &r.main().metadata.creators;
    assert_eq!(creators.len(), 1);
    assert_eq!(creators[0].index, 1);
    assert_eq!(creators[0].given_name.as_deref(), Some("Brett"));
    assert_eq!(
        creators[0].email.as_deref(),
        Some("bgoli@users.sourceforge.net")
    );
}

#[test]
fn symbol_annotation_loads() {
    // The corpus has none (libantimony drops CV terms), so only fixtures.
    let r = p(
        "S1 = 1;\nS1 identity \"http://identifiers.org/chebi/CHEBI:4167\"\nS1 part \"http://identifiers.org/chebi/CHEBI:17234\"",
    );
    r.assert_clean();
    let terms = &r.main().symbol("S1").unwrap().cv_terms;
    assert_eq!(terms.len(), 2);
    assert_eq!(terms[0].qualifier, Qualifier::Is);
    assert_eq!(terms[0].resource, "http://identifiers.org/chebi/CHEBI:4167");
    // The synonym the file used is kept, so a save gives it back as written.
    assert_eq!(terms[0].spelling, "identity");
    assert_eq!(terms[1].qualifier, Qualifier::HasPart);
    assert_eq!(terms[1].spelling, "part");
}

#[test]
fn annotation_text_is_preserved() {
    // An unknown qualifier is kept as raw text, written back verbatim, so a
    // file that loads today never stops loading.
    let r = p("model *m()\n  model somethingWeDoNotKnow \"a value\"\nend");
    r.assert_clean();
    let sections = &r.main().unparsed_sections;
    assert_eq!(sections.len(), 1);
    assert_eq!(sections[0].kind, RawKind::Annotation);
    assert!(sections[0].text.contains("a value"), "{}", sections[0].text);
}

#[test]
fn unit_display_name_does_not_invent_a_symbol() {
    // `time_unit is "time"` names a unit, not a symbol.
    let r = a("unit time_unit = second;\ntime_unit is \"time\";");
    r.assert_clean();
    assert!(r.main().symbol("time_unit").is_none());
    assert_eq!(
        r.main().unit("time_unit").unwrap().display_name.as_deref(),
        Some("time")
    );
    assert!(!r.has(Code::UninitialisedSymbol));
}

// --- TM5CorpusTests -------------------------------------------------------------------

#[test]
fn m5_subset_parses_and_classifies() {
    // Delphi M5's acceptance set: the 1629 cases with no SBML package.
    let Some(cases) = common::corpus_index("m5_subset_parses_and_classifies") else {
        return;
    };
    let m5: Vec<_> = cases.iter().filter(|c| c.is_m5()).collect();
    assert_eq!(m5.len(), 1629);
    let mut failed = Vec::new();
    for case in m5.iter().filter(|c| !c.expects_reject()) {
        let loaded = antimony::parse(&common::read(&case.path));
        if loaded.has_errors() {
            failed.push(format!("{}: {}", case.id, loaded.diagnostics));
        }
        for symbol in &loaded.document.main_module().symbols {
            assert_ne!(
                symbol.kind,
                SymbolKind::Unknown,
                "{}: {}",
                case.id,
                symbol.name
            );
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
fn no_corpus_case_is_unsupported() {
    // Every non-comp case reaches the semantic layer; only `comp` may still
    // report a construct as not supported.
    let Some(cases) = common::corpus_index("no_corpus_case_is_unsupported") else {
        return;
    };
    let mut unsupported = Vec::new();
    for case in cases
        .iter()
        .filter(|c| !c.has_package("comp") && !c.expects_reject())
    {
        let loaded = antimony::parse(&common::read(&case.path));
        if let Some(d) = loaded.diagnostics.first_of(Code::NotYetSupported) {
            unsupported.push(format!("{}: {}", case.id, d.message));
        }
    }
    assert!(
        unsupported.is_empty(),
        "{} unsupported. First: {}",
        unsupported.len(),
        unsupported[0]
    );
}

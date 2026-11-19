//! Ported from the Delphi `Tests/uAntFacadeTests.pas`.
//!
//! Delphi's `TAntimony` is an object that is loaded and then queried, and
//! its contract has three rules: a caller writes no tests for missing
//! values, a bad argument raises, a query before a load raises. Here
//! `antimony::parse` returns the model, and the queries are methods on it
//! (SPEC.md §6), so:
//!
//! - Rule 1 holds as it was: the value methods answer with each kind's
//!   default where the model is silent.
//! - Rule 2 becomes `None`: an id that is not of the kind asked about, or an
//!   index past the end (`get`), answers `None` rather than raising.
//! - Rule 3 cannot arise: there is no model to ask before `parse` returns
//!   one, and a failed load is `Err` from `Loaded::into_result`.
//!
//! Not ported, for those reasons: `AMissingFileIsADiagnosticNotAnException`
//! (the crate opens no files: reading one, and reporting a missing one, is
//! the caller's part; `parse_bytes` takes over the unreadable-file test),
//! `ReloadingDiscardsThePreviousModel` (nothing is reused between parses)
//! and `ARaisedMessageNamesTheMethod` (nothing is raised).
//! `EveryQueryRaisesBeforeAnyLoad` keeps its second half, a failed load.
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart.

use antimony::model::{Document, Module};
use antimony::{Code, Loaded, Options, flatten, parse, parse_bytes, parse_with, write_document};

/// Parses `source`, which must load without errors.
fn load(source: &str) -> Document {
    let loaded = parse(source);
    assert!(!loaded.has_errors(), "{}", loaded.diagnostics);
    loaded.document
}

fn module_of(document: &Document) -> &Module {
    document.main_module()
}

// --- loading ---------------------------------------------------------------------

#[test]
fn load_string_returns_true_for_a_good_model() {
    let loaded = parse("J1: S1 => S2; k1*S1;\nS1 = 1; S2 = 0; k1 = 0.1;");
    assert!(!loaded.has_errors());
    assert!(loaded.into_result().is_ok());
}

#[test]
fn a_model_error_is_returned_not_raised() {
    // An editor must be able to open a broken file and show what is wrong.
    let loaded = parse("J1: S1 => ~~~ ;");
    assert!(loaded.has_errors());
    assert!(!loaded.diagnostics.to_string().is_empty());
}

#[test]
fn warnings_do_not_fail_the_load() {
    let loaded = parse("J1: S1 => S2; k1*S1;");
    assert!(!loaded.has_errors());
    assert!(
        !loaded.diagnostics.is_empty(),
        "the warnings are still reported"
    );
}

#[test]
fn strict_mode_makes_warnings_fail() {
    let strict = Options {
        strict: true,
        ..Options::default()
    };
    assert!(parse_with("J1: S1 => S2; k1*S1;", &strict).has_errors());
}

#[test]
fn initialise_missing_reaches_the_semantic_pass() {
    let options = Options {
        initialise_missing: true,
        ..Options::default()
    };
    let loaded = parse_with("J1: S1 => S2; k1*S1;", &options);
    assert!(!loaded.has_errors());
    let text = loaded.diagnostics.to_string();
    assert!(text.contains("initialised"), "{text}");
}

#[test]
fn load_file_reads_from_disk() {
    // The caller reads the file; `parse_bytes` reads what it read.
    let path = std::env::temp_dir().join(format!("antimony_facade_{}.ant", std::process::id()));
    std::fs::write(&path, "S1 = 4;").unwrap();
    let bytes = std::fs::read(&path).unwrap();
    std::fs::remove_file(&path).ok();
    let loaded = parse_bytes(&bytes, &Options::default());
    assert!(!loaded.has_errors(), "{}", loaded.diagnostics);
    assert_eq!(module_of(&loaded.document).value("S1"), Some(4.0));
}

#[test]
fn an_unreadable_file_is_a_diagnostic_not_an_exception() {
    // 0xC0 and 0xC1 can never appear in UTF-8.
    let loaded = parse_bytes(
        &[0x53, 0x31, 0x20, 0x3D, 0x20, 0xC0, 0xC1, 0x3B],
        &Options::default(),
    );
    assert!(loaded.has_errors());
    assert!(loaded.diagnostics.has_any(Code::ReadError));
    assert!(!loaded.diagnostics.has_any(Code::FileNotFound));
    // Added in the port: it says where.
    let found = loaded.diagnostics.first_of(Code::ReadError).unwrap();
    assert_eq!((found.span.start.line, found.span.start.column), (1, 6));
}

#[test]
fn every_query_raises_before_any_load() {
    // A failed load gives no model to a caller who asks for one.
    let loaded = parse("J1: S1 => ~~~ ;");
    assert!(loaded.into_result().is_err());
}

// --- ownership -------------------------------------------------------------------

#[test]
fn flatten_returns_something_the_caller_owns() {
    let document = load("J1: S1 => S2; k1*S1;\nS1 = 1; k1 = 2;");
    let mut flat = flatten::flatten(&document);
    flat.reactions.clear();
    // Changing the copy leaves the document as it was.
    assert_eq!(module_of(&document).reactions.len(), 1);
    assert!(!antimony::write_document(&document).is_empty());
}

#[test]
fn an_unknown_module_id_raises() {
    let document = load("S1 = 1;");
    assert!(flatten::flatten_module(&document, "nosuchmodule").is_none());
    assert!(document.module("nosuchmodule").is_none());
}

// --- writing ---------------------------------------------------------------------

#[test]
fn sbml_and_antimony_round_trip() {
    let document = load("J1: S1 => S2; k1*S1;\nS1 = 1; S2 = 0; k1 = 0.1;");
    let text = write_document(&document);
    assert!(!text.is_empty());
    let again: Loaded = parse(&text);
    assert!(!again.has_errors(), "{}", again.diagnostics);
}

#[cfg(feature = "sbml")]
#[test]
fn sbml_defaults_to_the_main_module() {
    use antimony::sbml::{write_sbml, write_sbml_module};
    let document = load("model m1() S1 = 1; end\nmodel *m2() S2 = 2; end");
    // `*` names the main module, so that is what is written by default.
    assert!(write_sbml(&document).contains(r#"id="m2""#));
    let m1 = flatten::flatten_module(&document, "m1").unwrap();
    assert!(write_sbml_module(&m1).contains(r#"id="m1""#));
}

// --- the query surface -------------------------------------------------------------

/// Five reactions so "the 5th" is a real index, and the last one carries
/// stoichiometries worth asking about.
const QUERY_MODEL: &str = "compartment C = 1;
species A in C, B in C, P in C, Q in C, R in C, $X in C;
J1: A -> B; k*A;
J2: B -> A; k*B;
J3: A -> P; k*A;
J4: P -> Q; k*P;
J5: 2 A + 3 Q => R; k*A*Q;
A = 1; B = 0; P = 0; Q = 0; R = 0; X = 7;
k = 0.5;
const k;
";

#[test]
fn one_call_answers_a_stoichiometry_question() {
    // The stoichiometry of the 2nd reactant of the 5th reaction.
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.stoichiometry(&m.reactions[4].reactants[1]), 3.0);
}

#[test]
fn stoichiometry_knows_its_three_cases() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    // None written is 1.
    let first = &m.reactions[0].reactants[0];
    assert_eq!(m.stoichiometry(first), 1.0);
    assert_eq!(first.stoichiometry_text(), "1");
    // A number is itself.
    assert_eq!(m.stoichiometry(&m.reactions[4].reactants[0]), 2.0);

    // A named one is the symbol's value; its text is the name.
    let document = load("sr = 2;\nJ0: sr A -> B; k*A;\nA = 1; B = 0; k = 1;");
    let m = module_of(&document);
    let named = &m.reactions[0].reactants[0];
    assert_eq!(m.stoichiometry(named), 2.0);
    assert_eq!(named.stoichiometry_text(), "sr");
}

#[test]
fn reaction_queries_read_the_model() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.reactions.len(), 5);
    let j5 = &m.reactions[4];
    assert_eq!(j5.name, "J5");
    assert_eq!(j5.reactants.len(), 2, "A and Q");
    assert_eq!(j5.products.len(), 1, "R");
    assert_eq!(j5.reactants[1].species, "Q");
    assert_eq!(j5.products[0].species, "R");
    assert!(j5.kinetic_law.is_some());
    // `->` is reversible, `=>` is not.
    assert!(m.reactions[0].is_reversible(), "J1 uses ->");
    assert!(!j5.is_reversible(), "J5 uses =>");
}

#[test]
fn reactions_answer_by_id_as_well_as_index() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.reaction("J5"), Some(&m.reactions[4]));
}

#[test]
fn an_unnamed_reaction_still_has_an_id() {
    let document = load("S1 -> S2; k*S1;\nS1 = 1; S2 = 0; k = 1;");
    let m = module_of(&document);
    assert_eq!(m.reactions[0].name, "_J1");
    assert_eq!(m.reaction("_J1").unwrap().reactants[0].species, "S1");

    // Its position among all reactions, so naming one renumbers nothing.
    let document = load("Ja: S1 -> S2; k*S1;\nS2 -> S3; k*S2;\nS1 = 1; k = 1;");
    let names: Vec<_> = module_of(&document)
        .reactions
        .iter()
        .map(|r| r.name.as_str())
        .collect();
    assert_eq!(names, ["Ja", "_J2"]);

    // The author's own id wins; the made-up one steps past it.
    let document = load("_J1: S1 -> S2; k*S1;\nS2 -> S3; k*S2;\nS1 = 1; k = 1;");
    let names: Vec<_> = module_of(&document)
        .reactions
        .iter()
        .map(|r| r.name.as_str())
        .collect();
    assert_eq!(names, ["_J1", "_J2"]);
}

#[test]
fn quantity_kinds_have_their_own_lists() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.species().count(), 6);
    assert_eq!(m.compartments().count(), 1);
    assert!(m.parameters().count() >= 1, "k at least");
    assert_eq!(m.compartments().next().unwrap().name, "C");
    assert_eq!(m.symbol("A").unwrap().compartment.as_deref(), Some("C"));
    assert_eq!(m.compartment_size("C"), Some(1.0));
    assert_eq!(m.parameter_value("k"), Some(0.5));
    assert_eq!(m.concentration("X"), Some(7.0));
}

#[test]
fn boundary_is_not_constness() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.boundary_species().count(), 1);
    assert_eq!(m.floating_species().count(), 5);
    assert_eq!(m.boundary_species().next().unwrap().name, "X");
}

#[test]
fn boundary_ids_carry_no_dollar() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert!(m.species().all(|s| !s.name.contains('$')));
}

#[test]
fn absent_values_answer_with_their_default() {
    let loaded = parse("species S1 in C;\ncompartment C;\np = ;");
    assert!(!loaded.has_errors(), "{}", loaded.diagnostics);
    let m = module_of(&loaded.document);
    assert_eq!(
        m.concentration("S1"),
        Some(0.0),
        "a species with no value is 0"
    );
    // A size of 0 would divide every concentration by zero.
    assert_eq!(m.compartment_size("C"), Some(1.0));
    assert_eq!(
        m.parameter_value("p"),
        Some(0.0),
        "the warning is the signal"
    );
    assert!(
        loaded.diagnostics.has_any(Code::UninitialisedSymbol)
            || loaded.diagnostics.has_any(Code::ImplicitParameter)
    );

    // `a = b; b = a` is not a circular rule, so it loads; following it must
    // stop, and the answer is the default.
    let document = load("a = b; b = a;");
    assert_eq!(module_of(&document).parameter_value("a"), Some(0.0));
}

#[test]
fn concentration_and_amount_convert_through_the_volume() {
    let document = load("compartment C = 4;\nspecies S1 in C;\nS1 = 2;");
    let m = module_of(&document);
    // Declared as a concentration, so the amount is worked out.
    assert_eq!(m.concentration("S1"), Some(2.0));
    assert_eq!(m.amount("S1"), Some(8.0));
}

#[test]
fn the_three_expression_slots_are_separate() {
    let document = load("x := 3*y; y' = 2; y = 1; k = 2*j; j = 3;");
    let m = module_of(&document);
    let (x, y, k) = (
        m.symbol("x").unwrap(),
        m.symbol("y").unwrap(),
        m.symbol("k").unwrap(),
    );
    assert!(x.assignment_rule().is_some());
    assert!(x.rate_rule().is_none(), "x has an assignment, not a rate");
    assert!(y.rate_rule().is_some());
    assert!(y.assignment_rule().is_none());
    // As written, and as worked out.
    assert_eq!(k.initial_value.as_ref().unwrap().to_string(), "2 * j");
    assert_eq!(m.parameter_value("k"), Some(6.0));
}

// --- rule 2: a bad argument is None ------------------------------------------------

#[test]
fn an_out_of_range_index_raises() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert!(m.reactions.get(99).is_none());
    assert!(m.events.get(99).is_none());
    assert!(m.species().nth(99).is_none());
    assert!(m.parameters().nth(99).is_none());
    assert!(m.compartments().nth(99).is_none());
    // A real reaction, a participant past the end.
    assert!(m.reactions[4].reactants.get(99).is_none());
}

#[test]
fn an_unknown_id_raises() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.parameter_value("nosuch"), None);
    assert_eq!(m.concentration("nosuch"), None);
    assert_eq!(m.amount("nosuch"), None);
    assert_eq!(m.compartment_size("nosuch"), None);
    assert!(m.reaction("nosuch").is_none());
    assert!(m.symbol("nosuch").is_none());
    assert!(m.event("nosuch").is_none());
}

#[test]
fn the_wrong_kind_of_id_raises() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    assert_eq!(m.parameter_value("A"), None, "A is a species");
    assert_eq!(m.concentration("k"), None, "k is a parameter");
    assert_eq!(m.compartment_size("A"), None, "A is a species");
    // Floating and boundary are different lists.
    assert!(m.floating_species().all(|s| s.name != "X"));
    assert!(m.boundary_species().all(|s| s.name != "A"));
}

#[test]
fn counts_bound_a_loop_so_correct_code_never_raises() {
    let document = load(QUERY_MODEL);
    let m = module_of(&document);
    let total: f64 = m
        .reactions
        .iter()
        .flat_map(|r| &r.reactants)
        .map(|p| m.stoichiometry(p))
        .sum();
    // J1 to J4 have one reactant each; J5 is 2 A + 3 Q.
    assert_eq!(total, 9.0);
}

// --- Added in the Rust port --------------------------------------------------------

#[test]
fn a_byte_order_mark_is_skipped() {
    let loaded = parse_bytes(b"\xEF\xBB\xBFS1 = 4;", &Options::default());
    assert!(!loaded.has_errors(), "{}", loaded.diagnostics);
    assert_eq!(module_of(&loaded.document).value("S1"), Some(4.0));
}

#[test]
fn bytes_are_read_with_the_options() {
    let strict = Options {
        strict: true,
        ..Options::default()
    };
    assert!(parse_bytes(b"J1: S1 => S2; k1*S1;", &strict).has_errors());
}

#[test]
fn a_value_follows_an_assignment_rule_first() {
    // The rule holds at the start too, so it wins over an initial value.
    let document = load("x = 1; x := 2*k; k = 5;");
    assert_eq!(module_of(&document).value("x"), Some(10.0));
    // `time` has no value at the start, so neither does what uses it.
    let document = load("y := time + 1;");
    assert_eq!(module_of(&document).value("y"), None);
    assert_eq!(module_of(&document).parameter_value("y"), Some(0.0));
}

#[test]
fn an_amount_is_divided_by_the_size() {
    let document = load("compartment C = 4;\nsubstanceOnly species S1 in C;\nS1 = 2;");
    let m = module_of(&document);
    assert_eq!(m.amount("S1"), Some(2.0));
    assert_eq!(m.concentration("S1"), Some(0.5));
    // A size of 0 gives 0, not an infinity.
    let document = load("compartment C = 0;\nsubstanceOnly species S1 in C;\nS1 = 2;");
    assert_eq!(module_of(&document).concentration("S1"), Some(0.0));
}

#[test]
fn a_species_with_no_compartment_has_size_one() {
    let document = load("species S1;\nS1 = 3;");
    assert_eq!(module_of(&document).amount("S1"), Some(3.0));
}

#[test]
fn an_event_is_found_by_name() {
    let document = load("E1: at time > 5: S1 = 3;\nS1 = 0;");
    let m = module_of(&document);
    assert_eq!(m.event("E1").unwrap().trigger.to_string(), "time > 5");
}

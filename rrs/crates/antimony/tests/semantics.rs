//! Ported from the Delphi `Tests/uAntSemanticsTests.pas`.
//!
//! Changed for Rust: Delphi binds each reference to its symbol with a
//! pointer (`Binding`, `Compartment`, `Species`). Here references are names
//! (SPEC.md §3.1), so `ExpressionReferencesAreBound`,
//! `CompartmentReferenceIsBound`, `ParticipantSpeciesAreBound` and
//! `StoichiometrySymbolIsBound` check that the name resolves to that symbol
//! with `semantics::resolve`.
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart. The
//! corpus tests need `ANTIMONY_CORPUS` and `ANTIMONY_CORPUS_INDEX`.

mod common;

use antimony::expr::ExprKind;
use antimony::model::{Module, Symbol, SymbolKind};
use antimony::semantics::resolve;
use antimony::{Code, Loaded, Options, Severity};

struct Analysed(Loaded);

impl Analysed {
    fn main(&self) -> &Module {
        self.0.document.main_module()
    }

    fn sym(&self, name: &str) -> &Symbol {
        self.main()
            .symbol(name)
            .unwrap_or_else(|| panic!("no symbol named {name:?}"))
    }

    fn has(&self, code: Code) -> bool {
        self.0.diagnostics.has_any(code)
    }

    fn message(&self, code: Code) -> &str {
        &self.0.diagnostics.first_of(code).unwrap().message
    }

    fn assert_no_errors(&self) -> &Self {
        assert!(!self.0.has_errors(), "{}", self.0.diagnostics);
        self
    }

    fn value(&self, name: &str) -> f64 {
        self.sym(name)
            .initial_value
            .as_ref()
            .and_then(|e| e.evaluate())
            .unwrap_or_else(|| panic!("{name} has no value"))
    }
}

fn a(source: &str) -> Analysed {
    Analysed(antimony::parse(source))
}

fn a_with(source: &str, strict: bool, initialise_missing: bool) -> Analysed {
    Analysed(antimony::parse_with(
        source,
        &Options {
            strict,
            initialise_missing,
        },
    ))
}

// --- TResolutionTests ----------------------------------------------------------

#[test]
fn forward_reference_resolves() {
    // `S1 = k*2` may come before `k = 3`.
    let r = a("S1 = k*2;\nk = 3;");
    r.assert_no_errors();
    assert!(!r.sym("k").is_implicit, "k is declared, just later");
    assert!(!r.has(Code::ImplicitParameter));
}

#[test]
fn expression_references_are_bound() {
    let r = a("k = 3;\nS1 = k*2;");
    r.assert_no_errors();
    let value = r.sym("S1").initial_value.as_ref().unwrap();
    let ExprKind::Binary { left, .. } = &value.kind else {
        panic!("{value:?}");
    };
    let ExprKind::Symbol(name) = &left.kind else {
        panic!("{left:?}");
    };
    assert_eq!(resolve(r.main(), name), Some(r.sym("k")));
}

#[test]
fn builtin_names_are_not_symbols() {
    // `time` as a symbol would be invented as a parameter, 732 times.
    let r = a("S1 = time*2;");
    r.assert_no_errors();
    assert!(r.main().symbol("time").is_none());
    assert!(!r.has(Code::ImplicitParameter));
    assert!(resolve(r.main(), "time").is_none());
}

#[test]
fn call_names_are_not_resolved_as_symbols() {
    let r = a("S1 = sin(x);");
    r.assert_no_errors();
    assert!(
        r.main().symbol("sin").is_none(),
        "a function name is not a reference"
    );
    assert!(r.main().symbol("x").is_some(), "but its argument is");
}

#[test]
fn compartment_reference_is_bound() {
    let r = a("compartment c;\nspecies S1 in c;");
    r.assert_no_errors();
    let compartment = r.sym("S1").compartment.as_deref().unwrap();
    assert_eq!(resolve(r.main(), compartment), Some(r.sym("c")));
    assert_eq!(r.sym("c").kind, SymbolKind::Compartment);
}

#[test]
fn compartment_named_but_not_declared_is_inferred() {
    // The `in` clause says it is a compartment: inference, not invention.
    let r = a("species S1 in cell;");
    r.assert_no_errors();
    assert_eq!(r.sym("cell").kind, SymbolKind::Compartment);
    assert!(!r.has(Code::ImplicitParameter));
}

#[test]
fn participant_species_are_bound() {
    let r = a("J1: S1 => S2; k*S1;");
    r.assert_no_errors();
    let reaction = &r.main().reactions[0];
    assert_eq!(
        resolve(r.main(), &reaction.reactants[0].species),
        Some(r.sym("S1"))
    );
    assert_eq!(
        resolve(r.main(), &reaction.products[0].species),
        Some(r.sym("S2"))
    );
}

#[test]
fn stoichiometry_symbol_is_bound() {
    let r = a(
        "S1_S2_stoich := 2*p1;\np1 = 4;\nk = 1; S1 = 1; S2 = 0;\nreaction1: S1 => S1_S2_stoich S2; k*S1;",
    );
    r.assert_no_errors();
    assert!(!r.has(Code::ImplicitParameter), "{}", r.0.diagnostics);
    let stoichiometry = r.main().reactions[0].products[0].stoichiometry.as_ref();
    let Some(ExprKind::Symbol(name)) = stoichiometry.map(|e| &e.kind) else {
        panic!("{stoichiometry:?}");
    };
    assert_eq!(resolve(r.main(), name), Some(r.sym("S1_S2_stoich")));
}

// --- TTypeInferenceTests -------------------------------------------------------

#[test]
fn unclassified_becomes_parameter() {
    assert_eq!(a("k1 = 0.1;").sym("k1").kind, SymbolKind::Parameter);
}

#[test]
fn reaction_participant_becomes_species() {
    assert_eq!(
        a("k1 = 0.1;\nJ1: k1 => S2; 1;").sym("k1").kind,
        SymbolKind::Species,
        "appearing in a reaction promotes a parameter to a species"
    );
}

#[test]
fn parameter_promotes_to_species() {
    let r = a("k1 = 0.1;\nJ1: S1 => k1; 1;");
    r.assert_no_errors();
    assert_eq!(r.sym("k1").kind, SymbolKind::Species);
}

#[test]
fn parameter_promotes_to_compartment() {
    let r = a("v = 2;\nspecies S1 in v;");
    r.assert_no_errors();
    assert_eq!(r.sym("v").kind, SymbolKind::Compartment);
}

#[test]
fn compartment_does_not_become_species() {
    let r = a("compartment c;\nJ1: c => S2; 1;");
    assert!(r.has(Code::TypeConflict));
    assert_eq!(
        r.sym("c").kind,
        SymbolKind::Compartment,
        "an illegal promotion must not overwrite the established type"
    );
}

#[test]
fn nothing_is_left_unknown() {
    let r = a("compartment c;\nspecies S1 in c;\nJ1: S1 => S2; k*S1;\nk = 1;");
    for symbol in &r.main().symbols {
        assert_ne!(
            symbol.kind,
            SymbolKind::Unknown,
            "{} was left unclassified",
            symbol.name
        );
    }
}

#[test]
fn named_stoichiometry_is_not_a_parameter() {
    // The rule comes first, as in the corpus, so the symbol looks like an
    // ordinary formula until the participant list.
    let r = a(
        "S1_S2_stoichiometry := 2*p1;\np1 = 4; S1 = 1; S2 = 0; k = 1;\nreaction1: S1 => S1_S2_stoichiometry S2; k*S1;",
    );
    r.assert_no_errors();
    assert_eq!(r.sym("S1_S2_stoichiometry").kind, SymbolKind::Stoichiometry);
    assert_eq!(
        r.sym("p1").kind,
        SymbolKind::Parameter,
        "an ordinary parameter is unaffected"
    );
}

#[test]
fn module_name_is_not_a_symbol() {
    // Declaring one added a phantom parameter to 1542 corpus cases at once.
    let r = a("model *m()\n  k = 1;\nend\nm is \"A pretty name\"");
    r.assert_no_errors();
    assert!(
        r.main().symbol("m").is_none(),
        "the module name must not become a symbol"
    );
    assert_eq!(r.main().display_name.as_deref(), Some("A pretty name"));
}

// --- TUninitialisedSymbolTests -------------------------------------------------

#[test]
fn never_declared_is_an_implicit_parameter() {
    // The typo case, the more serious of the two.
    let r = a("A = 1; B = 0; k1 = 2;\nJ1: A => B; k1*typoo;");
    assert!(r.has(Code::ImplicitParameter), "{}", r.0.diagnostics);
    let message = r.message(Code::ImplicitParameter);
    assert!(
        message.contains("typoo"),
        "the message must name the symbol"
    );
    assert!(message.contains("typo"), "and point at the likely cause");
    assert!(r.sym("typoo").is_implicit);
}

#[test]
fn declared_but_valueless_is_uninitialised() {
    let r = a("species S1;");
    assert!(r.has(Code::UninitialisedSymbol));
    assert!(
        !r.has(Code::ImplicitParameter),
        "a declared symbol is not an invented one"
    );
}

#[test]
fn valued_symbols_are_silent() {
    let r = a("k1 = 1;\ncompartment c;\nc = 1;");
    assert!(!r.has(Code::UninitialisedSymbol), "{}", r.0.diagnostics);
    assert!(!r.has(Code::ImplicitParameter));
}

#[test]
fn a_symbol_with_only_a_rule_is_initialised() {
    let r = a("Vtot := Vmax*E;\nVmax = 1;\nE = 2;");
    assert!(
        !r.has(Code::UninitialisedSymbol),
        "a rule is a source of value: {}",
        r.0.diagnostics
    );
}

#[test]
fn warnings_do_not_fail_the_load() {
    let r = a("J1: A => B; k1*unknownThing;");
    assert!(r.has(Code::ImplicitParameter));
    assert!(!r.0.has_errors(), "a warning must not fail the load");
}

#[test]
fn strict_mode_promotes_them_to_errors() {
    let r = a_with("J1: A => B; k1*unknownThing;", true, false);
    assert!(r.0.has_errors(), "strict mode makes them errors");
    assert_eq!(
        r.0.diagnostics
            .first_of(Code::ImplicitParameter)
            .unwrap()
            .severity,
        Severity::Error
    );
}

#[test]
fn initialise_missing_is_off_by_default() {
    let r = a("species S1;");
    assert!(
        r.sym("S1").initial_value.is_none(),
        "nothing is invented by default"
    );
    assert!(r.has(Code::UninitialisedSymbol));
}

#[test]
fn initialise_missing_gives_a_species_zero() {
    assert_eq!(a_with("species S1;", false, true).value("S1"), 0.0);
}

#[test]
fn initialise_missing_gives_a_compartment_one() {
    // 0 would make every concentration a division by zero.
    let r = a_with("compartment C; species S1 in C;", false, true);
    assert_eq!(r.value("C"), 1.0, "a compartment gets 1, not 0");
    assert_eq!(r.value("S1"), 0.0, "the species still gets 0");
}

#[test]
fn initialise_missing_says_what_it_did() {
    let r = a_with("species S1;", false, true);
    assert!(
        r.message(Code::UninitialisedSymbol)
            .contains("initialised to 0")
    );
}

#[test]
fn initialise_missing_leaves_real_values_alone() {
    let r = a_with("species S1; S1 = 4; k1 = 0.1;", false, true);
    assert_eq!(r.value("S1"), 4.0, "an existing value is not overwritten");
    assert!(!r.has(Code::UninitialisedSymbol));
}

#[test]
fn an_undeclared_name_is_still_reported_as_a_typo() {
    let r = a_with("J1: A => B; k1*unknownThing;", false, true);
    assert!(
        r.has(Code::ImplicitParameter),
        "the typo warning survives filling in"
    );
    assert!(
        r.message(Code::ImplicitParameter)
            .contains("Check for a typo")
    );
}

// --- TValidationTests ------------------------------------------------------------

#[test]
fn const_with_a_rate_rule_is_reported() {
    let r = a("k1 = 1;\nconst k1;\nk1' = 2;");
    assert!(r.has(Code::ConstWithRule), "{}", r.0.diagnostics);
}

#[test]
fn const_with_an_assignment_rule_is_reported() {
    assert!(a("k1 = 1;\nconst k1;\nk1 := 2*x;").has(Code::ConstWithRule));
}

#[test]
fn species_in_a_non_compartment_is_reported() {
    let r = a("J1: v => X; 1;\nspecies S1 in v;");
    assert!(
        r.has(Code::NotACompartment) || r.has(Code::TypeConflict),
        "{}",
        r.0.diagnostics
    );
}

#[test]
fn conversion_factor_must_be_a_parameter() {
    let r = a("compartment c;\nspecies S1;\nS1 conversionFactor c;");
    assert!(r.has(Code::BadConversionFactor), "{}", r.0.diagnostics);
}

#[test]
fn direct_circular_rule_is_reported() {
    assert!(a("x := x + 1;").has(Code::CircularRule));
}

#[test]
fn indirect_circular_rule_names_the_path() {
    let r = a("a := b + 1;\nb := c + 1;\nc := a + 1;");
    assert_eq!(
        r.message(Code::CircularRule),
        "Circular assignment-rule dependency: a -> b -> c -> a"
    );
}

#[test]
fn rate_rule_may_mention_its_own_subject() {
    // `S1' = -k1*S1` is a differential equation, not a cycle.
    let r = a("k1 = 1;\nS1 = 10;\nS1' = -k1*S1;");
    assert!(!r.has(Code::CircularRule), "{}", r.0.diagnostics);
}

#[test]
fn acyclic_rule_chain_is_silent() {
    let r = a("a := b + 1;\nb := c + 1;\nc = 1;");
    assert!(!r.has(Code::CircularRule), "{}", r.0.diagnostics);
}

// --- Added in the Rust port ---------------------------------------------------------

#[test]
fn function_arguments_are_local() {
    // The `bvar` leak (Delphi §9.1.4 #4): a function's arguments are not
    // module parameters.
    let r =
        a("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend\nJ1: A -> B; MM(A, 1, 2)\nA = 1; B = 0");
    r.assert_no_errors();
    assert!(r.main().symbol("Vm").is_none());
    assert!(!r.has(Code::ImplicitParameter), "{}", r.0.diagnostics);
}

#[test]
fn a_declared_parameter_promotes_to_compartment() {
    // `parameter_promotes_to_compartment` above never reaches the promotion:
    // `v = 2` leaves v unclassified until the end of the pass. `formula`
    // makes it a parameter at once.
    let r = a("formula v = 2;\nspecies S1 in v;");
    r.assert_no_errors();
    assert_eq!(r.sym("v").kind, SymbolKind::Compartment);
}

#[test]
fn arity_and_recursion_are_checked() {
    let r = a("function f(x)\n  f(x) + 1\nend\nk = f(1, 2)");
    assert!(r.has(Code::RecursiveFunction));
    assert_eq!(
        r.message(Code::ArityMismatch),
        "Function \"f\" takes 1 argument(s), but 2 were supplied"
    );
}

#[test]
fn a_long_chain_of_rules_does_not_exhaust_the_stack() {
    // The cycle check is a loop, not recursion.
    let mut source = String::new();
    for i in 0..20_000 {
        source.push_str(&format!("x{i} := x{} + 1\n", i + 1));
    }
    source.push_str("x20000 := x0\n");
    let worker = std::thread::Builder::new()
        .stack_size(1 << 20)
        .spawn(move || a(&source).has(Code::CircularRule))
        .unwrap();
    assert!(worker.join().expect("ran out of stack"));
}

#[test]
fn diagnostics_come_parser_first_then_semantics() {
    let r = a("compartment c\nspecies c\nJ1: A -> B; typoo");
    let codes: Vec<Code> = r.0.diagnostics.iter().map(|d| d.code).collect();
    assert_eq!(codes[0], Code::TypeConflict);
    assert!(codes.contains(&Code::ImplicitParameter));
}

#[test]
fn into_result_fails_only_on_errors() {
    assert!(antimony::parse("J1: A -> B; k*typoo").into_result().is_ok());
    let errors = antimony::parse("k = (1").into_result().unwrap_err();
    assert!(errors.has_errors());
}

#[test]
fn flatten_is_a_copy_of_the_main_module() {
    let loaded = antimony::parse("model *m()\n  k = 1\nend");
    let flat = antimony::flatten::flatten(&loaded.document);
    assert_eq!(&flat, loaded.document.main_module());
    assert_eq!(
        antimony::flatten::flatten_module(&loaded.document, "m"),
        Some(flat)
    );
    assert!(antimony::flatten::flatten_module(&loaded.document, "x").is_none());
}

// --- TSemanticsCorpusTests -------------------------------------------------------------

/// The cases of a set that invent a parameter, with the names.
fn implicit_parameters(cases: &[&common::Case]) -> Vec<String> {
    let mut found = Vec::new();
    for case in cases.iter().filter(|c| !c.expects_reject()) {
        let loaded = antimony::parse(&common::read(&case.path));
        let names: Vec<&str> = loaded
            .diagnostics
            .iter()
            .filter(|d| d.code == Code::ImplicitParameter)
            .map(|d| d.message.as_str())
            .collect();
        if !names.is_empty() {
            found.push(format!("{}: {}", case.id, names.join("; ")));
        }
    }
    found
}

#[test]
fn core_subset_produces_no_implicit_parameters() {
    // These are valid models, so every name in them resolves; one warning
    // means a missing built-in or a construct the parser did not declare.
    let Some(cases) = common::corpus_index("core_subset_produces_no_implicit_parameters") else {
        return;
    };
    let core: Vec<_> = cases.iter().filter(|c| c.is_core()).collect();
    assert_eq!(core.iter().filter(|c| !c.expects_reject()).count(), 831);
    let found = implicit_parameters(&core);
    assert!(
        found.is_empty(),
        "{} cases. First: {}",
        found.len(),
        found[0]
    );
}

#[test]
fn whole_corpus_produces_no_implicit_parameters() {
    // `comp` is excluded: a dotted name into a module that is not built
    // cannot resolve.
    let Some(cases) = common::corpus_index("whole_corpus_produces_no_implicit_parameters") else {
        return;
    };
    let non_comp: Vec<_> = cases.iter().filter(|c| !c.has_package("comp")).collect();
    assert!(non_comp.len() > 1600);
    let found = implicit_parameters(&non_comp);
    assert!(
        found.is_empty(),
        "{} cases. First: {}",
        found.len(),
        found[0]
    );
}

#[test]
fn core_subset_classifies_every_symbol() {
    let Some(cases) = common::corpus_index("core_subset_classifies_every_symbol") else {
        return;
    };
    for case in cases.iter().filter(|c| c.is_core() && !c.expects_reject()) {
        let loaded = antimony::parse(&common::read(&case.path));
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
}

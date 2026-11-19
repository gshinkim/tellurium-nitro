//! Ported from the Delphi `Tests/uAntParserTests.pas`.
//!
//! Changed for Rust: `ForwardReferencesAreNotResolved` drops its check that
//! the compartment object is unbound, since references are names (SPEC.md
//! §3.1). The invalid-fixture tests run the whole pipeline, parser and
//! semantic pass, as Delphi's do (`antimony::parse`).
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart. The
//! corpus tests need `ANTIMONY_CORPUS` and `ANTIMONY_CORPUS_INDEX`.

mod common;

use antimony::model::{
    Constness, ConstraintKind, Document, Module, RawKind, Reaction, RuleKind, Symbol, SymbolKind,
};
use antimony::parser::parse_document;
use antimony::{Code, Diagnostics, Severity};

/// A parse: the document and its diagnostics.
struct Parsed {
    document: Document,
    diagnostics: Diagnostics,
}

impl Parsed {
    fn main(&self) -> &Module {
        self.document.main_module()
    }

    fn sym(&self, name: &str) -> &Symbol {
        self.main()
            .symbol(name)
            .unwrap_or_else(|| panic!("no symbol named {name:?}"))
    }

    fn rxn(&self, i: usize) -> &Reaction {
        self.main()
            .reactions
            .get(i)
            .unwrap_or_else(|| panic!("expected at least {} reaction(s)", i + 1))
    }

    fn assert_clean(&self) -> &Self {
        assert!(!self.diagnostics.has_errors(), "{}", self.diagnostics);
        self
    }
}

fn p(source: &str) -> Parsed {
    let (document, diagnostics) = parse_document(source);
    Parsed {
        document,
        diagnostics,
    }
}

/// Parses `source`, which must have no errors.
fn clean(source: &str) -> Parsed {
    let parsed = p(source);
    parsed.assert_clean();
    parsed
}

fn text(e: &Option<antimony::expr::Expr>) -> String {
    e.as_ref().expect("no expression").to_string()
}

// --- TModuleTests ------------------------------------------------------------

#[test]
fn empty_source_gives_an_empty_module() {
    let parsed = clean("");
    assert!(parsed.main().symbols.is_empty());
}

#[test]
fn model_block_names_the_module() {
    let parsed = clean("model cell()\n  k1 = 1;\nend");
    assert_eq!(parsed.main().name.as_deref(), Some("cell"));
}

#[test]
fn main_module_marker_is_recorded() {
    let parsed = clean("model *case00001()\nend");
    assert!(parsed.main().is_main);
    assert_eq!(parsed.main().name.as_deref(), Some("case00001"));
}

#[test]
fn formal_arguments_are_captured() {
    let parsed = clean("model side_reaction(S, k1)\nend");
    assert_eq!(parsed.main().formal_args, ["S", "k1"]);
}

#[test]
fn duplicate_formal_argument_is_reported() {
    assert!(
        p("model m(a, a)\nend")
            .diagnostics
            .has_any(Code::DuplicateFormalArg)
    );
}

#[test]
fn missing_end_is_reported() {
    assert!(
        p("model cell()\n  k1 = 1;")
            .diagnostics
            .has_any(Code::ExpectedEnd)
    );
}

#[test]
fn statements_may_follow_the_end() {
    // `case00001 is "case00001"` ends almost every corpus file, outside the
    // block and with no semicolon.
    let parsed = clean("model *m()\nend\nm is \"Case one\"");
    assert_eq!(parsed.main().display_name.as_deref(), Some("Case one"));
    assert!(parsed.main().symbols.is_empty(), "the module is no symbol");
}

#[test]
fn statements_may_precede_any_model() {
    assert_eq!(clean("k1 = 1;\nS1 = 2;").main().symbols.len(), 2);
}

#[test]
fn document_holds_a_list_of_modules() {
    let parsed = p("model *m()\nend");
    assert!(!parsed.document.modules.is_empty());
    assert_eq!(parsed.document.main, 0);
}

// --- TDeclarationTests -------------------------------------------------------

#[test]
fn compartment_declaration() {
    assert_eq!(
        clean("compartment compartment_;").sym("compartment_").kind,
        SymbolKind::Compartment
    );
}

#[test]
fn species_list_with_compartments() {
    let parsed = clean("species S1 in compartment_, S2 in compartment_;");
    assert_eq!(parsed.sym("S1").kind, SymbolKind::Species);
    assert_eq!(
        parsed.sym("S1").compartment.as_deref(),
        Some("compartment_")
    );
    assert_eq!(
        parsed.sym("S2").compartment.as_deref(),
        Some("compartment_")
    );
}

#[test]
fn boundary_marker_is_per_item() {
    let parsed = clean("species G in Cell, $T in Cell;");
    assert!(!parsed.sym("G").is_boundary);
    assert!(parsed.sym("T").is_boundary);
}

#[test]
fn substance_only_species() {
    let parsed = clean("substanceOnly species S1 in c;");
    assert!(parsed.sym("S1").is_substance_only);
    assert_eq!(parsed.sym("S1").kind, SymbolKind::Species);

    let parsed = clean("substanceOnly species $S4 in c;");
    assert!(parsed.sym("S4").is_substance_only);
    assert!(parsed.sym("S4").is_boundary);
}

#[test]
fn const_qualifier_list() {
    let parsed = clean("compartment c;\nk1 = 1;\nconst c, k1;");
    assert_eq!(parsed.sym("c").constness, Constness::Const);
    assert_eq!(parsed.sym("k1").constness, Constness::Const);
    assert_eq!(
        parsed.sym("c").kind,
        SymbolKind::Compartment,
        "the qualifier must not disturb the established type"
    );
}

#[test]
fn var_qualifier_list() {
    let parsed = clean("var x, y;");
    assert_eq!(parsed.sym("x").constness, Constness::Var);
    assert_eq!(parsed.sym("y").constness, Constness::Var);
}

#[test]
fn const_species_is_boundary() {
    assert!(clean("species S1;\nconst S1;").sym("S1").is_boundary);
}

#[test]
fn declaration_with_initial_value() {
    let parsed = clean("species S1 in c = 5;");
    assert_eq!(parsed.sym("S1").compartment.as_deref(), Some("c"));
    assert_eq!(text(&parsed.sym("S1").initial_value), "5");
}

#[test]
fn in_placement_on_its_own_statement() {
    let parsed = clean("S1 in cell = 5;");
    assert_eq!(parsed.sym("S1").compartment.as_deref(), Some("cell"));
    assert_eq!(text(&parsed.sym("S1").initial_value), "5");
}

#[test]
fn unit_attachment() {
    assert_eq!(
        clean("S1 has substance_per_volume;")
            .sym("S1")
            .units
            .as_deref(),
        Some("substance_per_volume")
    );
}

#[test]
fn display_name() {
    assert_eq!(
        clean("S1 is \"Glucose 6-phosphate\";")
            .sym("S1")
            .display_name
            .as_deref(),
        Some("Glucose 6-phosphate")
    );
}

#[test]
fn notes_statement() {
    assert_eq!(
        clean("m notes \"Originally created by me\";")
            .sym("m")
            .notes
            .as_deref(),
        Some("Originally created by me")
    );
}

#[test]
fn conversion_factor_statement() {
    assert_eq!(
        clean("S1 conversionFactor s1_cf;")
            .sym("S1")
            .conversion_factor
            .as_deref(),
        Some("s1_cf")
    );
}

// --- TAssignmentAndRuleTests -------------------------------------------------

#[test]
fn initial_assignment_stores_an_expression() {
    let parsed = clean("S1 = 0.00015/compartment_;");
    assert_eq!(
        text(&parsed.sym("S1").initial_value),
        "0.00015 / compartment_"
    );
    assert!(parsed.sym("S1").rule.is_none());
}

#[test]
fn assignment_rule_uses_colon_equals() {
    let parsed = clean("Vtot := Vmax*E;");
    let rule = parsed.sym("Vtot").rule.as_ref().unwrap();
    assert_eq!(rule.kind, RuleKind::Assignment);
    assert_eq!(text(&rule.expr), "Vmax * E");
    assert!(parsed.sym("Vtot").initial_value.is_none());
}

#[test]
fn rate_rule_uses_a_prime() {
    let parsed = clean("S1' = -k1*S1;");
    let rule = parsed.sym("S1").rule.as_ref().unwrap();
    assert_eq!(rule.kind, RuleKind::Rate);
    assert_eq!(text(&rule.expr), "-k1 * S1");
}

#[test]
fn empty_right_hand_side_keeps_the_declared_value() {
    // libantimony loses the 3 here (case 01234); we must not.
    let parsed = clean("p = 3;\np = ;");
    assert_eq!(
        text(&parsed.sym("p").initial_value),
        "3",
        "the 3 must survive"
    );
    assert!(parsed.sym("p").has_empty_value);
}

#[test]
fn empty_right_hand_side_is_recorded() {
    // 31 core cases write `compartment_ = ;`.
    let parsed = clean("compartment_ = ;");
    assert!(parsed.sym("compartment_").initial_value.is_none());
    assert!(parsed.sym("compartment_").has_empty_value);
}

#[test]
fn forward_references_are_not_resolved() {
    // The parser binds nothing, which is what lets `S1 = k*2` come first.
    let parsed = clean("S1 = k*2;\nk = 3;");
    assert!(parsed.main().symbol("k").is_some());
    assert_eq!(parsed.main().symbols[0].name, "S1");
}

// --- TReactionTests -----------------------------------------------------------

#[test]
fn named_reaction_with_kinetic_law() {
    let parsed = clean("reaction1: S1 => S2; compartment_*k1*S1;");
    let r = parsed.rxn(0);
    assert_eq!(r.name, "reaction1");
    assert_eq!(r.reactants.len(), 1);
    assert_eq!(r.reactants[0].species, "S1");
    assert_eq!(r.products[0].species, "S2");
    assert_eq!(text(&r.kinetic_law), "compartment_ * k1 * S1");
}

#[test]
fn unnamed_reaction() {
    // Named `_J1` when read, so a simulator and every writer agree on it.
    let parsed = clean("S1 -> S2; k1*S1;");
    assert_eq!(parsed.rxn(0).name, "_J1");
    assert_eq!(parsed.rxn(0).reactants.len(), 1);
}

#[test]
fn arrow_direction_follows_the_corpus() {
    // `->` is reversible and `=>` irreversible (Delphi §9.1.1).
    assert!(
        clean("J1: S1 -> S2;").rxn(0).is_reversible(),
        "-> is reversible"
    );
    assert!(
        !clean("J1: S1 => S2;").rxn(0).is_reversible(),
        "=> is irreversible"
    );
}

#[test]
fn numeric_stoichiometry() {
    let parsed = clean("reaction1: S1 => 2 S2; k*S1;");
    let r = parsed.rxn(0);
    assert_eq!(text(&r.products[0].stoichiometry), "2");
    assert_eq!(r.products[0].species, "S2");
    assert!(r.reactants[0].stoichiometry.is_none(), "absent means one");
}

#[test]
fn named_stoichiometry_is_an_expression() {
    let parsed = clean("reaction1: S1 => S1_S2_stoichiometry S2; k*S1;");
    let product = &parsed.rxn(0).products[0];
    assert_eq!(text(&product.stoichiometry), "S1_S2_stoichiometry");
    assert_eq!(product.species, "S2");
    assert_eq!(
        parsed.sym("S1_S2_stoichiometry").kind,
        SymbolKind::Stoichiometry
    );
}

#[test]
fn boundary_species_in_a_reaction() {
    let parsed = clean("$Xo -> S1; k*Xo;");
    assert!(parsed.rxn(0).reactants[0].is_boundary);
    assert!(parsed.sym("Xo").is_boundary);
}

#[test]
fn empty_reactant_side() {
    let parsed = clean("J0: -> Xref X; k1;");
    let r = parsed.rxn(0);
    assert!(r.reactants.is_empty());
    assert_eq!(r.products.len(), 1);
    assert_eq!(text(&r.products[0].stoichiometry), "Xref");
    assert_eq!(r.products[0].species, "X");
}

#[test]
fn empty_product_side() {
    let parsed = clean("pp9_v2: pp9_mrna => ; rate*pp9_mrna;");
    let r = parsed.rxn(0);
    assert_eq!((r.reactants.len(), r.products.len()), (1, 0));
    assert!(r.kinetic_law.is_some());
}

#[test]
fn reaction_with_no_participants_at_all() {
    // A challenge case: libantimony writes it and cannot read it back. We
    // accept it, with a hint.
    let parsed = clean("J0:  -> ; 2;");
    assert!(parsed.rxn(0).reactants.is_empty() && parsed.rxn(0).products.is_empty());
    let hint = parsed.diagnostics.first_of(Code::EmptyReaction).unwrap();
    assert_eq!(hint.severity, Severity::Hint);
}

#[test]
fn empty_kinetic_law_is_distinct_from_absent() {
    // 34 fbc cases write `R16: J => M; ;`.
    let parsed = clean("R16: J => M; ;");
    assert!(parsed.rxn(0).kinetic_law.is_none());
    assert!(parsed.rxn(0).has_empty_kinetic_law);

    let parsed = clean("R17: J => M;");
    assert!(parsed.rxn(0).kinetic_law.is_none());
    assert!(!parsed.rxn(0).has_empty_kinetic_law);
}

#[test]
fn multiple_participants_on_both_sides() {
    let parsed = clean("reaction1: S1 + S2 => S3 + S4; k1*S1*S2;");
    assert_eq!(parsed.rxn(0).reactants.len(), 2);
    assert_eq!(parsed.rxn(0).products.len(), 2);
    assert_eq!(parsed.rxn(0).products[1].species, "S4");
}

#[test]
fn reaction_becomes_a_reaction_symbol() {
    let parsed = clean("reaction1: S1 => S2; k*S1;");
    assert_eq!(parsed.sym("reaction1").kind, SymbolKind::Reaction);
    assert_eq!(parsed.sym("S1").kind, SymbolKind::Species);
}

// --- TUnitDefinitionTests -----------------------------------------------------

#[test]
fn simple_unit_definition() {
    let parsed = clean("unit volume = litre;");
    assert_eq!(parsed.main().units.len(), 1);
    assert_eq!(parsed.main().units[0].name, "volume");
    assert_eq!(parsed.main().units[0].definition, "litre");
}

#[test]
fn compound_unit_definition_is_kept_verbatim() {
    assert_eq!(
        clean("unit substance_per_volume = mole / litre;")
            .main()
            .units[0]
            .definition,
        "mole / litre"
    );
    assert_eq!(
        clean("unit area = metre^2;").main().units[0].definition,
        "metre^2"
    );
}

#[test]
fn juxtaposed_unit_definition_does_not_fail() {
    // A unit expression juxtaposes, which the expression grammar cannot say.
    assert_eq!(
        clean("unit x = mole litre / (second metre metre);")
            .main()
            .units[0]
            .definition,
        "mole litre / (second metre metre)"
    );
}

#[test]
fn self_referential_unit_is_flagged() {
    assert!(clean("unit time_unit = time_unit;").main().units[0].is_self_referential());
}

// --- TFluxBalanceTests ---------------------------------------------------------

#[test]
fn named_objective() {
    let parsed = clean("OBJF: maximize R26 + R10;");
    let objectives = &parsed.main().objectives;
    assert_eq!(objectives.len(), 1);
    assert_eq!(objectives[0].name.as_deref(), Some("OBJF"));
    assert!(objectives[0].maximise);
    assert_eq!(objectives[0].expr.to_string(), "R26 + R10");
}

#[test]
fn unnamed_objective() {
    let parsed = clean("maximize J1;");
    assert_eq!(parsed.main().objectives.len(), 1);
    assert_eq!(parsed.main().objectives[0].name, None);
}

#[test]
fn minimize_and_maximize_are_distinct() {
    let parsed = clean("obj1: minimize J2;");
    assert!(!parsed.main().objectives[0].maximise);
    assert_eq!(parsed.main().objectives[0].name.as_deref(), Some("obj1"));
}

#[test]
fn objective_name_is_not_a_symbol() {
    // libantimony reports no `OBJF` symbol for case 01612.
    let parsed = clean("R26: A => B; k*A;\nOBJF: maximize R26;");
    assert!(parsed.main().symbol("OBJF").is_none());
    assert!(parsed.main().symbol("R26").is_some());
}

#[test]
fn objective_takes_an_expression() {
    assert_eq!(
        clean("obj3: maximize J1 + 3*J2 + 4*J3^2 + 5*J4*J5;")
            .main()
            .objectives[0]
            .expr
            .to_string(),
        "J1 + 3 * J2 + 4 * J3 ^ 2 + 5 * J4 * J5"
    );
}

#[test]
fn dotted_sbo_term_does_not_invent_a_symbol() {
    let parsed = clean("fb_0 = 0;\nfb_0.sboTerm = 625");
    assert!(parsed.main().symbol("fb_0.sboTerm").is_none());
    assert_eq!(parsed.sym("fb_0").sbo_term, Some(625));
    assert_eq!(parsed.main().symbols.len(), 1);
}

#[test]
fn space_separated_sbo_term() {
    assert_eq!(clean("A sboTerm = 236").sym("A").sbo_term, Some(236));
}

// --- TRecoveryTests --------------------------------------------------------------

#[test]
fn one_bad_statement_costs_one_diagnostic() {
    let parsed = p("k1 = 1;\n@@@ ;\nk2 = 2;");
    assert!(parsed.diagnostics.has_errors());
    assert!(
        parsed.diagnostics.count_of(Severity::Error) <= 3,
        "{}",
        parsed.diagnostics
    );
}

#[test]
fn parsing_continues_after_an_error() {
    let parsed = p("k1 = 1;\n@@@ ;\nk2 = 2;");
    assert!(parsed.main().symbol("k1").is_some());
    assert!(
        parsed.main().symbol("k2").is_some(),
        "the statement after the error must still be parsed"
    );
}

#[test]
fn deferred_constructs_name_their_milestone() {
    for source in ["import \"other.txt\";", "A: subcell();"] {
        let parsed = p(source);
        let d = parsed
            .diagnostics
            .first_of(Code::NotYetSupported)
            .unwrap_or_else(|| panic!("{source}: {}", parsed.diagnostics));
        assert!(d.message.contains("M10"), "{source}: {}", d.message);
    }
}

#[test]
fn dna_syntax_is_rejected_by_name() {
    let parsed = p("dna1: --P1--G1--");
    assert!(
        parsed.diagnostics.has_any(Code::DnaNotSupported),
        "{}",
        parsed.diagnostics
    );
}

#[test]
fn layout_statement_is_read_and_kept() {
    // Antimony 3's layout (tests/layout.rs has the rest).
    let parsed = clean("k1 = 1;\nmodel.layout = on\nk2 = 2;");
    assert_eq!(parsed.main().layout.len(), 1);
    assert!(parsed.main().unparsed_sections.is_empty());
    assert!(parsed.main().symbol("k2").is_some());
}

// --- Added in the Rust port ---------------------------------------------------------

#[test]
fn a_qualifier_may_carry_a_value() {
    // libantimony accepts these; Delphi did not (SPEC.md R36).
    let parsed = clean("const k1 = 0.5;\nvar k = 3, j = 2;\nconst S1 in C = 3;");
    assert_eq!(parsed.sym("k1").constness, Constness::Const);
    assert_eq!(text(&parsed.sym("k1").initial_value), "0.5");
    assert_eq!(parsed.sym("k").constness, Constness::Var);
    assert_eq!(text(&parsed.sym("j").initial_value), "2");
    assert_eq!(parsed.sym("S1").compartment.as_deref(), Some("C"));
    assert_eq!(text(&parsed.sym("S1").initial_value), "3");
}

#[test]
fn a_dollar_assignment_is_a_const_one() {
    // libantimony reads `$B = 4` as `const B; B = 4` (SPEC.md R36).
    let parsed = clean("$B = 4;");
    let b = parsed.sym("B");
    assert_eq!((b.constness, b.is_boundary), (Constness::Const, true));
    assert_eq!(text(&b.initial_value), "4");
    // Used in a reaction, it is a boundary species.
    let parsed = clean("$B = 4;\nJ1: $B -> C; k*B;");
    assert_eq!(parsed.sym("B").kind, SymbolKind::Species);
    assert!(parsed.sym("B").is_boundary);
    // An unnamed reaction may still begin with `$`.
    assert_eq!(clean("$A -> B; k;").rxn(0).reactants[0].species, "A");
}

#[test]
fn constraints_and_chained_ranges() {
    let parsed = clean("constraint c0: R16 >= 0\n-10 <= J2 < 10\nJ1 <= 1000");
    let constraints = &parsed.main().constraints;
    assert_eq!(constraints.len(), 3);
    assert_eq!(constraints[0].name.as_deref(), Some("c0"));
    assert!(matches!(&constraints[0].kind, ConstraintKind::Expr(e) if e.to_string() == "R16 >= 0"));
    let ConstraintKind::Range {
        lower,
        subject,
        upper,
        lower_inclusive,
        upper_inclusive,
    } = &constraints[1].kind
    else {
        panic!("not a range: {:?}", constraints[1]);
    };
    assert_eq!(
        (lower.to_string(), subject.to_string(), upper.to_string()),
        ("-10".into(), "J2".into(), "10".into())
    );
    assert!(*lower_inclusive && !*upper_inclusive);
    assert!(!constraints[2].is_range());
}

#[test]
fn events_with_delay_and_modifiers() {
    let parsed = clean(
        "E1: at 2 after time > 5, priority = 1, t0 = false, persistent = true: S1 = 0, k = k + 1",
    );
    let e = &parsed.main().events[0];
    assert_eq!(e.name.as_deref(), Some("E1"));
    assert_eq!(e.trigger.to_string(), "time > 5");
    assert_eq!(text(&e.delay), "2");
    assert_eq!(text(&e.priority), "1");
    assert_eq!(text(&e.initial_value), "false");
    assert_eq!(text(&e.persistent), "true");
    assert_eq!(e.assignments.len(), 2);
    assert_eq!(e.assignments[1].target, "k");
    assert_eq!(parsed.sym("E1").kind, SymbolKind::Event);

    // An empty assignment list is a challenge case.
    assert!(
        clean("E0: at time > 5.5: ;").main().events[0]
            .assignments
            .is_empty()
    );
}

#[test]
fn functions_before_the_model_belong_to_it() {
    // BIOMD0000000064: functions written before the `model` block it uses
    // them in must not be left in a module of their own.
    let parsed = clean("function f(x, y)\n  x*y\nend\nmodel *m()\n  J1: A -> B; f(k, A)\nend");
    assert_eq!(parsed.document.modules.len(), 1);
    assert_eq!(parsed.main().functions[0].args, ["x", "y"]);
    assert_eq!(text(&parsed.main().functions[0].body), "x * y");
}

#[test]
fn interactions_are_not_reactions() {
    let parsed = clean("J1: A -> B; k*A\ni1: G1 -| J1\nG2 -o J1\nG3 -( J1");
    let interactions = &parsed.main().interactions;
    assert_eq!(interactions.len(), 3);
    assert_eq!(interactions[0].name.as_deref(), Some("i1"));
    assert_eq!(interactions[2].reaction, "J1");
    assert_eq!(parsed.main().reactions.len(), 1);
    assert_eq!(parsed.sym("G3").kind, SymbolKind::Species);
}

#[test]
fn annotations_are_filed_on_their_subjects() {
    let source = "model *m()\n  S1 = 1\n  S1 identity \"http://a\",\n     \"http://b\"\n  unit u = mole\n  u part \"http://c\"\nend\nm created \"2012-12-21T11:05:21Z\"\nm creator1.familyName \"Smith\"\nS1 unknownQualifier \"x\"";
    let parsed = clean(source);
    let s1 = parsed.sym("S1");
    assert_eq!(s1.cv_terms.len(), 2, "one term per resource");
    assert_eq!(s1.cv_terms[1].resource, "http://b");
    assert_eq!(parsed.main().units[0].cv_terms[0].spelling, "part");
    assert!(
        parsed.main().symbol("u").is_none(),
        "a unit is not a symbol"
    );
    let metadata = &parsed.main().metadata;
    assert_eq!(metadata.created.as_deref(), Some("2012-12-21T11:05:21Z"));
    assert_eq!(metadata.creators[0].family_name.as_deref(), Some("Smith"));
    let raw = &parsed.main().unparsed_sections;
    assert_eq!(raw.len(), 1);
    assert_eq!(raw[0].kind, RawKind::Annotation);
    assert_eq!(raw[0].text, "S1 unknownQualifier \"x\"");
}

#[test]
fn generated_reaction_ids_skip_taken_names() {
    // The second reaction would be `_J2`, but a symbol has that name; the
    // third would be `_J3`, but the second took it.
    let parsed = clean("_J2 = 1\nA -> B; 1\nC -> D; 1\nE -> F; 1");
    let names: Vec<&str> = parsed
        .main()
        .reactions
        .iter()
        .map(|r| r.name.as_str())
        .collect();
    assert_eq!(names, ["_J1", "_J3", "_J4"]);
}

#[test]
fn comments_attach_to_what_they_precede() {
    let parsed = clean("// the rate\nk1 = 1\n// the reaction\nJ1: A -> B; k1*A");
    assert_eq!(parsed.sym("k1").comment.as_deref(), Some("// the rate"));
    assert_eq!(parsed.rxn(0).comment.as_deref(), Some("// the reaction"));
}

#[test]
fn type_conflict_names_both_places() {
    let parsed = p("compartment c\nspecies c");
    let d = parsed.diagnostics.first_of(Code::TypeConflict).unwrap();
    assert_eq!(d.span.start.line, 2);
    assert!(d.message.contains("1:13-14"), "{}", d.message);
}

/// Mangled corpus files must give diagnostics, never a panic.
#[test]
fn mangled_input_never_panics() {
    let Some(root) = common::corpus_root("mangled_input_never_panics") else {
        return;
    };
    // A small deterministic generator, so a failure can be reproduced.
    let mut seed: u64 = 0x2545_f491_4f6c_dd1d;
    let mut random = |n: usize| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % n.max(1) as u64) as usize
    };
    let junk = [
        "(", ")", ":", ";", "->", "=", "\"", "```", "/*", "end", "model", "at", ",", "-", "$", "'",
        "\n",
    ];
    for (case, path) in common::corpus_files(&root).into_iter().step_by(7) {
        let source = common::read(&path);
        for _ in 0..8 {
            let mut bytes: Vec<char> = source.chars().collect();
            for _ in 0..1 + random(6) {
                let at = random(bytes.len());
                match random(3) {
                    0 => {
                        bytes.truncate(at);
                    }
                    1 if at < bytes.len() => {
                        bytes.remove(at);
                    }
                    _ => {
                        let insert = junk[random(junk.len())];
                        for (i, c) in insert.chars().enumerate() {
                            bytes.insert((at + i).min(bytes.len()), c);
                        }
                    }
                }
            }
            let mangled: String = bytes.into_iter().collect();
            let result = std::panic::catch_unwind(|| parse_document(&mangled));
            assert!(result.is_ok(), "{case} panicked on:\n{mangled}");
        }
    }
}

// --- TParserCorpusTests ---------------------------------------------------------------

#[test]
fn all_core_subset_cases_parse_without_errors() {
    let Some(cases) = common::corpus_index("all_core_subset_cases_parse_without_errors") else {
        return;
    };
    let core: Vec<_> = cases.iter().filter(|c| c.is_core()).collect();
    assert_eq!(core.len(), 836);
    let (mut rejected, mut failed) = (0, Vec::new());
    for case in &core {
        let (_, diagnostics) = parse_document(&common::read(&case.path));
        if case.expects_reject() {
            // The five negative-stoichiometry cases must fail, by name.
            rejected += 1;
            assert!(
                diagnostics.has_any(Code::NegativeStoichiometry),
                "{} should be rejected with dcNegativeStoichiometry: {diagnostics}",
                case.id
            );
        } else if diagnostics.has_errors() {
            failed.push(format!("{}: {diagnostics}", case.id));
        }
    }
    assert_eq!(rejected, 5, "the expect_reject cases in the core subset");
    assert!(
        failed.is_empty(),
        "{} core cases failed. First: {}",
        failed.len(),
        failed[0]
    );
}

#[test]
fn challenge_cases_parse() {
    // Text libantimony writes and cannot read back: we must read it.
    let Some(cases) = common::corpus_index("challenge_cases_parse") else {
        return;
    };
    let challenges: Vec<_> = cases.iter().filter(|c| c.is_challenge()).collect();
    assert_eq!(challenges.len(), 10);
    for case in challenges {
        let (_, diagnostics) = parse_document(&common::read(&case.path));
        assert!(!diagnostics.has_errors(), "{}: {diagnostics}", case.id);
    }
}

#[test]
fn core_subset_symbol_counts_are_sane() {
    // Case 00001 is one compartment, two species and one reaction.
    let Some(root) = common::corpus_root("core_subset_symbol_counts_are_sane") else {
        return;
    };
    let parsed = clean(&common::read(&root.join("00001/00001-antimony.ant")));
    let main = parsed.main();
    assert_eq!(main.name.as_deref(), Some("case00001"));
    assert!(main.is_main);
    assert_eq!(main.symbols_of(SymbolKind::Compartment).count(), 1);
    assert_eq!(main.symbols_of(SymbolKind::Species).count(), 2);
    assert_eq!(main.reactions.len(), 1);
    assert_eq!(main.units.len(), 4);
    assert_eq!(
        parsed.sym("compartment_").display_name.as_deref(),
        Some("compartment")
    );
}

#[test]
fn every_non_comp_case_parses() {
    // Beyond the Delphi parser tests: Delphi's figure of 1657 parseable
    // cases outside the `comp` package. Of the 1663, the five flagged
    // `expect_reject` and case 01622 have a negative stoichiometry, which is
    // rejected by name (01622, an fbc case, is not flagged in the index, but
    // the Delphi converter rejects it the same way). The `comp` cases use
    // module syntax that is not built (`sub1.p80 * conv`, `sub1.E0: …`).
    let Some(cases) = common::corpus_index("every_non_comp_case_parses") else {
        return;
    };
    let mut failed = Vec::new();
    let mut parsed = 0;
    for case in cases.iter().filter(|c| !c.has_package("comp")) {
        let (_, diagnostics) = parse_document(&common::read(&case.path));
        if case.expects_reject() || case.id == "01622" {
            assert!(
                diagnostics.has_any(Code::NegativeStoichiometry),
                "{}: {diagnostics}",
                case.id
            );
        } else if diagnostics.has_errors() {
            failed.push(format!("{}: {diagnostics}", case.id));
        } else {
            parsed += 1;
        }
    }
    assert_eq!(parsed + failed.len(), 1657);
    assert!(
        failed.is_empty(),
        "{} failed. First: {}",
        failed.len(),
        failed[0]
    );
}

// --- TInvalidCorpusTests ------------------------------------------------------------------

/// Every fixture. They cover semantic diagnostics too (arity, recursion,
/// type conflicts), so each is read by the whole pipeline.
fn fixtures() -> Vec<std::path::PathBuf> {
    common::files_with(&common::invalid_dir(), "ant")
}

fn diagnose(file: &std::path::Path) -> Diagnostics {
    antimony::parse(&common::read(file)).diagnostics
}

#[test]
fn every_fixture_produces_its_expected_diagnostics() {
    let mut missing = String::new();
    for file in fixtures() {
        let expected = common::read_expected(&file.with_extension("expected")).unwrap();
        let diagnostics = diagnose(&file);
        for e in &expected {
            let found = diagnostics
                .iter()
                .any(|d| d.code.name() == e.code && d.span.start.line == e.line);
            if !found {
                missing += &format!(
                    "{}: expected {} at line {}; got:\n{diagnostics}\n",
                    file.display(),
                    e.code,
                    e.line
                );
            }
        }
    }
    assert!(missing.is_empty(), "{missing}");
}

#[test]
fn no_fixture_produces_an_unexpected_code() {
    // An error the fixture did not predict means it is not testing what it
    // says. Hints and warnings are commentary, and allowed.
    let mut surprises = String::new();
    for file in fixtures() {
        let expected = common::read_expected(&file.with_extension("expected")).unwrap();
        let diagnostics = diagnose(&file);
        for d in diagnostics.iter().filter(|d| d.severity == Severity::Error) {
            if !expected.iter().any(|e| e.code == d.code.name()) {
                surprises += &format!("{}: unexpected {d}\n", file.display());
            }
        }
    }
    assert!(surprises.is_empty(), "{surprises}");
}

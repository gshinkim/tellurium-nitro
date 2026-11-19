//! Ported from the Delphi `Tests/uAntWriterTests.pas`.
//!
//! Changed for Rust: `NumbersAreLocaleIndependent` keeps only its assertion,
//! since Rust has no locale to set (SPEC.md R12). The round trips use the
//! comparer in `common/compare.rs`, which checks more than Delphi's did.
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart. The
//! corpus tests need `ANTIMONY_CORPUS` and `ANTIMONY_CORPUS_INDEX`.

mod common;

use antimony::model::Document;
use antimony::parser::parse_document;
use antimony::{WriteOptions, write_document};

/// Parses `source` (which must parse) and writes it back.
fn w(source: &str) -> String {
    write_document(&parsed(source))
}

fn parsed(source: &str) -> Document {
    let (document, diagnostics) = parse_document(source);
    assert!(
        !diagnostics.has_errors(),
        "the source did not parse: {diagnostics}"
    );
    document
}

/// The differences between a document and the one its written text reads
/// back as, or why the text did not read back.
fn round_trip(document: &Document) -> Result<(), String> {
    let text = write_document(document);
    let (second, diagnostics) = parse_document(&text);
    if diagnostics.has_errors() {
        return Err(format!(
            "the written text did not parse back:\n{text}\n{diagnostics}"
        ));
    }
    let differences = common::compare::differences(document.main_module(), second.main_module());
    if differences.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "the model changed:\n{}\n--- written text ---\n{text}",
            differences.join("\n")
        ))
    }
}

fn assert_round_trips(source: &str) {
    if let Err(problem) = round_trip(&parsed(source)) {
        panic!("{problem}");
    }
}

// --- TWriterOutputTests ----------------------------------------------------------

#[test]
fn module_header_and_end() {
    // The `*` may appear where the source had none: a document with one
    // module has one main module, and the model is what round-trips.
    let text = w("model cell()\n  k1 = 1;\nend");
    assert!(text.contains("cell()"), "{text}");
    assert!(text.contains("end"), "{text}");
}

#[test]
fn main_module_marker() {
    assert!(w("model *m()\nend").contains("model *m()"));
}

#[test]
fn species_and_compartment_declarations() {
    let text = w("compartment c;\nspecies S1 in c;");
    assert!(text.contains("compartment c;"), "{text}");
    assert!(text.contains("species S1 in c;"), "{text}");
}

#[test]
fn boundary_and_substance_only() {
    let text = w("substanceOnly species $S1 in c;");
    assert!(text.contains("substanceOnly species $S1"), "{text}");
}

#[test]
fn reaction_arrows_match_the_corpus_convention() {
    // `->` reversible, `=>` irreversible.
    assert!(w("J1: A -> B; k*A;").contains("A -> B"));
    assert!(w("J1: A => B; k*A;").contains("A => B"));
}

#[test]
fn stoichiometry_is_emitted() {
    assert!(w("J1: A => 2 B; k*A;").contains("=> 2 B"));
}

#[test]
fn qualifiers_are_grouped() {
    let text = w("a = 1; b = 2;\nconst a, b;");
    assert!(text.contains("const a, b;"), "{text}");
}

#[test]
fn unit_definition_is_verbatim() {
    let text = w("unit x = mole litre / (second metre metre);");
    assert!(
        text.contains("unit x = mole litre / (second metre metre);"),
        "{text}"
    );
}

#[test]
fn numbers_are_locale_independent() {
    let text = w("k = 0.1;");
    assert!(text.contains("0.1") && !text.contains("0,1"), "{text}");
}

// --- TRoundTripTests -------------------------------------------------------------
//
// Round-trip closure (Delphi §9.1.4 #1), on the degenerate forms that break a
// writer which only handles the typical case.

#[test]
fn simple_model() {
    assert_round_trips(
        "model *m()\n  compartment c;\n  species S1 in c, S2 in c;\n  J1: S1 => S2; c*k1*S1;\n  S1 = 10; S2 = 0; c = 1; k1 = 0.1;\n  const c, k1;\nend",
    );
}

#[test]
fn empty_right_hand_side_survives() {
    // libantimony loses the value along with the empty rule.
    assert_round_trips("p = 3;\np = ;");
    assert_round_trips("compartment_ = ;");
}

#[test]
fn empty_kinetic_law_survives() {
    assert_round_trips("R16: J => M; ;");
    assert_round_trips("R17: J => M;");
}

#[test]
fn reaction_with_no_participants() {
    assert_round_trips("J0:  -> ; 2;");
}

#[test]
fn event_with_no_assignments() {
    assert_round_trips("E0: at time > 5.5: ;");
}

#[test]
fn negative_stoichiometry_is_never_emitted_bare() {
    // A bare `-` in a participant is rejected, so it is bracketed.
    let text = w("J1: A => (0 - 2) B; 1;");
    assert!(!text.contains("=> -"), "a bare minus: {text}");
    assert_round_trips("J1: A => (0 - 2) B; 1;");
}

#[test]
fn double_unary_minus_does_not_become_the_dna_operator() {
    assert_round_trips("k = - -5;");
    assert_round_trips("k = 1 - -5;");
}

#[test]
fn named_stoichiometry() {
    assert_round_trips(
        "S1_S2_stoichiometry := 2*p1;\np1 = 4;\nreaction1: S1 => S1_S2_stoichiometry S2; k*S1;",
    );
}

#[test]
fn event_with_delay_and_modifiers() {
    assert_round_trips(
        "E0: at 5.7 after S1 < 0.75, priority = 5, t0 = false, persistent = false, fromTrigger = true: S2 = 1, S3 = 2;",
    );
}

#[test]
fn function_definition() {
    assert_round_trips("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend");
}

#[test]
fn algebraic_rule_including_empty() {
    assert_round_trips("_alg0: 0 = -1*k1 + S1 + S2;");
    assert_round_trips("_alg0: 0 = ;");
}

#[test]
fn constraints_and_ranges() {
    assert_round_trips("constraint c0: R16 >= 0");
    assert_round_trips("constraint b: -10 <= J2 <= 10");
    assert_round_trips("constraint b: 0 < J2 <= 10");
}

#[test]
fn objectives_and_sbo_terms() {
    assert_round_trips("OBJF: maximize R26 + R10;");
    assert_round_trips("obj1: minimize J2;");
    assert_round_trips("fb_0 = 0;\nfb_0.sboTerm = 625");
}

#[test]
fn units_and_display_names() {
    assert_round_trips(
        "unit volume = litre;\nunit substance_per_volume = mole / litre;\nunit time_unit = time_unit;\nS1 has substance_per_volume;\nS1 is \"Glucose 6-phosphate\";\ntime_unit is \"time\";",
    );
}

#[test]
fn conversion_factor() {
    assert_round_trips("S1 conversionFactor s1_cf;\ns1_cf = 2;");
}

#[test]
fn infinity_and_nan() {
    assert_round_trips("P = INF; Q = -INF; R = NaN;");
    assert_round_trips("P = inf; R = nan;");
}

#[test]
fn numeric_literal_with_units() {
    assert_round_trips("S1 = 0.00015 mole/C;\ncompartment C;");
    assert_round_trips("x = 40 foo/3 seconds;");
}

// --- TCommentTests --------------------------------------------------------------

#[test]
fn comment_above_a_statement_is_preserved() {
    let text = w("// Reactions:\nJ1: A => B; k*A;");
    assert!(text.contains("// Reactions:"), "{text}");
}

#[test]
fn no_section_headers_are_invented() {
    // A file keeps the headings it came with rather than gaining a set on
    // every save.
    let text = w("J1: A => B; k*A;");
    assert!(
        !text.contains("// Reactions") && !text.contains("// Species"),
        "{text}"
    );
}

// Changed for Rust: Delphi re-emitted a `layout { … }` block, which Antimony
// does not have; its layout statements are re-emitted instead (D1 of
// docs/antimony-layout-design.md; tests/layout.rs has the rest).
#[test]
fn layout_statements_are_re_emitted() {
    let text = w("k1 = 1;\nmodel.layout = on\nS1.pos = {1, 2}");
    assert!(
        text.contains("model.layout = on\nS1.pos = {1, 2}"),
        "{text}"
    );
}

#[test]
fn annotation_is_re_emitted_verbatim() {
    let text = w("model *m()\n  model created \"2012-12-21T11:05:21Z\"\nend");
    assert!(text.contains("2012-12-21T11:05:21Z"), "{text}");
}

// --- Added in the Rust port -------------------------------------------------------

#[test]
fn qualifiers_with_values_round_trip() {
    // `$B = 4` is written `B = 4; const $B;`, and `k in C` on something
    // that is not a species keeps its compartment (SPEC.md R36).
    let text = w("$B = 4;\nconst k1 = 0.5;\nvar S1 in C = 3;");
    assert!(text.contains("const $B, k1;"), "{text}");
    assert!(text.contains("S1 in C;"), "{text}");
    assert_round_trips("$B = 4;\nconst k1 = 0.5;\nvar S1 in C = 3;");
    assert_round_trips("$B = 4;\nJ1: $B -> D; k1*B;");
}

#[test]
fn an_unlabelled_reaction_stays_unlabelled() {
    // Delphi wrote the generated `_J1:` label, and reading that back added a
    // symbol `_J1` the model did not have (SPEC.md R28).
    let text = w("S1 = 1; S2 = 0; k1 = 1\nS1 -> S2; k1*S1");
    assert!(!text.contains("_J1"), "{text}");
    assert_round_trips("S1 = 1; S2 = 0; k1 = 1\nS1 -> S2; k1*S1\nJ2: S2 -> S1; k1*S2");
}

#[test]
fn awkward_strings_round_trip() {
    // Quotes, backslashes and line breaks in a value, and a note with all of
    // them: the cases that broke Delphi's writer on real BioModels.
    assert_round_trips(
        "S1 is \"a \\\"quoted\\\" name\\nover two lines with a \\\\ backslash\";\nS1 = 1",
    );
    assert_round_trips(
        "model *m()\n  model notes ```\n  <p>a \"note\" with \\ and   trailing spaces   </p>\n```\n  S1 notes \"plain\"\nend",
    );
    assert_round_trips(
        "model *m()\n  model creator1.organization \"Institute of Parasitology and Biomedicine \\\"López-Neyra\\\", CSIC.\"\n  ATP identity \"http://a\",\n    \"http://b\"\n  ATP part \"http://c\"\nend",
    );
}

#[test]
fn interactions_round_trip() {
    assert_round_trips("J1: A -> B; k*A\ni1: G1 -| J1\nG2 -o J1\nG3 -( J1");
}

#[test]
fn every_option_reads_back_the_same() {
    let document = parsed(
        "model *m()\n  compartment c\n  species S1 in c\n  J1: S1 -> S2; k*S1\n  function f(x)\n    x\n  end\nend",
    );
    for options in [
        WriteOptions {
            indent: "\t".into(),
            explicit_declarations: true,
            blank_line_between_groups: false,
        },
        WriteOptions {
            indent: String::new(),
            explicit_declarations: false,
            blank_line_between_groups: true,
        },
    ] {
        let text = antimony::write::write_document_with(&document, &options);
        let (second, d) = parse_document(&text);
        assert!(!d.has_errors(), "{d}\n{text}");
        // Without explicit declarations, S1 is still a species (it is in a
        // reaction) and c a compartment (S1 is in it).
        let changed = common::compare::differences(document.main_module(), second.main_module());
        if options.explicit_declarations {
            assert!(changed.is_empty(), "{changed:?}\n{text}");
        }
    }
}

#[test]
fn comments_come_back_on_what_they_described() {
    let source = "// the compartment\ncompartment c\n// a rate\nk1 = 1\n// the reaction\nJ1: A -> B; k1*A\n// an event\nE1: at time > 1: k1 = 2";
    let first = parsed(source);
    let second = parsed(&write_document(&first));
    let comments = |d: &Document| {
        let m = d.main_module();
        (
            m.symbol("c").unwrap().comment.clone(),
            m.symbol("k1").unwrap().comment.clone(),
            m.reactions[0].comment.clone(),
            m.events[0].comment.clone(),
        )
    };
    assert_eq!(comments(&first), comments(&second));
    assert_eq!(comments(&second).2.as_deref(), Some("// the reaction"));
}

#[test]
fn an_analysed_model_round_trips() {
    // The writer on a model after the semantic pass, with implicit
    // parameters and inferred compartments.
    let loaded = antimony::parse("species S1 in cell\nJ1: S1 -> S2; k*S1*typo\nS1 = 1");
    let text = write_document(&loaded.document);
    let again = antimony::parse(&text);
    let changed =
        common::compare::differences(loaded.document.main_module(), again.document.main_module());
    assert!(changed.is_empty(), "{changed:?}\n{text}");
}

// --- TWriterCorpusTests ------------------------------------------------------------

/// Round-trips each case; returns (cases tried, failures).
fn round_trip_cases(cases: &[&common::Case], analyse: bool) -> (usize, Vec<String>) {
    let mut tried = 0;
    let mut failed = Vec::new();
    for case in cases {
        let source = common::read(&case.path);
        let first = if analyse {
            let loaded = antimony::parse(&source);
            if loaded.has_errors() {
                continue;
            }
            loaded.document
        } else {
            let (document, d) = parse_document(&source);
            if d.has_errors() {
                continue;
            }
            document
        };
        tried += 1;
        let text = write_document(&first);
        let second = if analyse {
            let loaded = antimony::parse(&text);
            (loaded.document, loaded.diagnostics)
        } else {
            parse_document(&text)
        };
        if second.1.has_errors() {
            failed.push(format!("{} did not read back: {}", case.id, second.1));
            continue;
        }
        let changed = common::compare::differences(first.main_module(), second.0.main_module());
        if !changed.is_empty() {
            failed.push(format!("{} changed: {}", case.id, changed.join("; ")));
        }
    }
    (tried, failed)
}

fn non_comp(cases: &[common::Case]) -> Vec<&common::Case> {
    cases
        .iter()
        .filter(|c| !c.expects_reject() && !c.has_package("comp"))
        .collect()
}

#[test]
fn every_corpus_file_round_trips() {
    // Delphi M6's acceptance criterion: parse, write, parse gives the same
    // model, over all 1657 parseable non-comp cases.
    let Some(cases) = common::corpus_index("every_corpus_file_round_trips") else {
        return;
    };
    let (tried, failed) = round_trip_cases(&non_comp(&cases), false);
    assert_eq!(tried, 1657);
    assert!(
        failed.is_empty(),
        "{} of {tried} failed. First: {}",
        failed.len(),
        failed[0]
    );
}

#[test]
fn every_analysed_corpus_file_round_trips() {
    // Beyond Delphi: the same after the semantic pass, as `antimony::parse`
    // gives a model to a caller.
    let Some(cases) = common::corpus_index("every_analysed_corpus_file_round_trips") else {
        return;
    };
    let (tried, failed) = round_trip_cases(&non_comp(&cases), true);
    assert_eq!(tried, 1657);
    assert!(
        failed.is_empty(),
        "{} of {tried} failed. First: {}",
        failed.len(),
        failed[0]
    );
}

#[test]
fn challenge_cases_round_trip() {
    let Some(cases) = common::corpus_index("challenge_cases_round_trip") else {
        return;
    };
    let challenges: Vec<_> = cases.iter().filter(|c| c.is_challenge()).collect();
    let (tried, failed) = round_trip_cases(&challenges, false);
    assert_eq!(tried, 10);
    assert!(failed.is_empty(), "{}", failed.join("\n"));
}

#[test]
fn writing_is_idempotent() {
    // The fixed point comes after one save: the first write groups the
    // statements, which reorders the symbols; after that nothing moves. So
    // the check is that the second and third writes are the same text.
    let Some(cases) = common::corpus_index("writing_is_idempotent") else {
        return;
    };
    let mut unsettled = Vec::new();
    for case in non_comp(&cases) {
        let (first, d) = parse_document(&common::read(&case.path));
        if d.has_errors() {
            continue;
        }
        let once = write_document(&first);
        let twice = write_document(&parse_document(&once).0);
        let thrice = write_document(&parse_document(&twice).0);
        if twice != thrice {
            unsettled.push(case.id.clone());
        }
    }
    assert!(
        unsettled.is_empty(),
        "{} had not settled by the second write. First: {}",
        unsettled.len(),
        unsettled[0]
    );
}

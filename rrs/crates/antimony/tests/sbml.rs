//! Ported from the Delphi `Tests/uAntSBMLTests.pas`: flattening, the SBML
//! structure, and MathML. Only with the `sbml` feature.
//!
//! Changed for Rust:
//! - Every document a test writes must also pass `sbml_io::validate` and
//!   read back with `sbml_io::read`; the helper checks both.
//! - The flatten tests that checked Delphi's pointers (a deep copy, references
//!   and bindings into the copy) check the same promise with names, since
//!   references are names here (SPEC.md R3).
//! - Where the Delphi tests matched XML text with attributes in the Delphi
//!   writer's order, these match `sbml-io`'s order (`metaid` first); the
//!   attributes are the same.
//! - `NumbersAreLocaleIndependent` keeps only its assertion: Rust has no
//!   locale to set (SPEC.md R12).
//!
//! The tests after "Added in the Rust port" have no Delphi counterpart. The
//! corpus test needs `ANTIMONY_CORPUS` and `ANTIMONY_CORPUS_INDEX`; with
//! `ANTIMONY_SBML_OUT` set to a folder it also saves each case's SBML there,
//! as `<case>.xml`, for checking with libSBML.

#![cfg(feature = "sbml")]

mod common;

use antimony::flatten::{flatten, flatten_module};
use antimony::model::{Document, InteractionKind, Module, SymbolKind};
use antimony::sbml::{to_sbml, write_sbml, xml_escape};
use antimony::{Code, write_document};

/// Reads `source`, which must load without errors.
fn m(source: &str) -> Document {
    let loaded = antimony::parse(source);
    assert!(
        !loaded.has_errors(),
        "the source did not load: {}",
        loaded.diagnostics
    );
    loaded.document
}

/// The SBML for `source`, checked to be valid and to read back.
fn x(source: &str) -> String {
    sbml_of(&m(source))
}

fn sbml_of(document: &Document) -> String {
    let sbml = to_sbml(&flatten(document));
    let issues = sbml_io::validate(&sbml);
    assert!(issues.is_empty(), "invalid SBML: {issues:?}");
    let text = write_sbml(document);
    if let Err(problems) = sbml_io::read(&text) {
        panic!("the SBML does not read back: {problems:?}\n{text}");
    }
    text
}

/// The SBML for `source` as `sbml-io` data.
fn model(source: &str) -> sbml_io::Model {
    to_sbml(&flatten(&m(source))).model.expect("a model")
}

/// How many times `sub` occurs in `text`. "One hasPart element" and "two"
/// both contain the string; the bag tests need a count.
fn occurrences(text: &str, sub: &str) -> usize {
    text.matches(sub).count()
}

#[track_caller]
fn has(xml: &str, sub: &str) {
    assert!(xml.contains(sub), "expected {sub:?} in:\n{xml}");
}

#[track_caller]
fn lacks(xml: &str, sub: &str) {
    assert!(!xml.contains(sub), "did not expect {sub:?} in:\n{xml}");
}

// --- TFlattenTests -----------------------------------------------------------

#[test]
fn flatten_returns_an_equal_module() {
    let document =
        m("compartment c; species S1 in c;\nJ1: S1 => S2; k*S1;\nk = 1; c = 1; S1 = 5; S2 = 0;");
    let flat = flatten(&document);
    assert_eq!(flat.symbols.len(), document.main_module().symbols.len());
    assert_eq!(flat.reactions.len(), document.main_module().reactions.len());
    assert_eq!(flat.symbol("S1").unwrap().compartment.as_deref(), Some("c"));
}

#[test]
fn flatten_is_a_deep_copy() {
    // Delphi compared object identity. Ownership makes the copy a copy; what
    // is left to check is that changing it leaves the original alone.
    let document = m("k = 1;");
    let mut flat = flatten(&document);
    flat.symbol_mut("k").unwrap().initial_value = None;
    assert!(
        document
            .main_module()
            .symbol("k")
            .unwrap()
            .initial_value
            .is_some()
    );
}

#[test]
fn the_original_is_untouched() {
    let document = m("k = 1;");
    let mut flat = flatten(&document);
    flat.symbol_mut("k").unwrap().name = "renamed".into();
    assert!(document.main_module().symbol("k").is_some());
}

#[test]
fn references_point_into_the_copy() {
    // A reference is a name, so it resolves in whichever module holds it.
    let flat = flatten(&m("compartment c; species S1 in c;\nc = 1; S1 = 0;"));
    let compartment = flat.symbol("S1").unwrap().compartment.clone().unwrap();
    assert_eq!(
        flat.symbol(&compartment).unwrap().kind,
        SymbolKind::Compartment
    );
}

#[test]
fn expression_bindings_point_into_the_copy() {
    let flat = flatten(&m("k = 2;\nx = k*3;"));
    let value = flat.symbol("x").unwrap().initial_value.as_ref().unwrap();
    assert_eq!(value.symbols(), ["k"]);
    assert!(flat.symbol("k").is_some());
}

#[test]
fn flatten_by_name_and_by_default() {
    let document = m("model *m()\n  k = 1;\nend");
    assert_eq!(
        flatten_module(&document, "m").unwrap().name.as_deref(),
        Some("m")
    );
    assert_eq!(flatten(&document).name.as_deref(), Some("m"));
    assert!(flatten_module(&document, "nosuchmodule").is_none());
}

// --- TSBMLStructureTests -------------------------------------------------------

#[test]
fn document_header_is_l3v2() {
    let xml = x("model *m()\nend");
    has(
        &xml,
        r#"xmlns="http://www.sbml.org/sbml/level3/version2/core""#,
    );
    has(&xml, r#"level="3" version="2""#);
    has(&xml, r#"<model id="m""#);
}

#[test]
fn compartments_species_and_parameters() {
    let xml = x("compartment c; species S1 in c;\nc = 1; S1 = 5; k = 0.1;");
    has(&xml, r#"<compartment id="c""#);
    has(&xml, r#"spatialDimensions="3""#);
    has(&xml, r#"<species id="S1""#);
    has(&xml, r#"compartment="c""#);
    has(&xml, r#"<parameter id="k""#);
}

#[test]
fn reversibility_follows_the_arrow() {
    has(&x("J1: A -> B; k*A;"), r#"reversible="true""#);
    has(&x("J1: A => B; k*A;"), r#"reversible="false""#);
}

#[test]
fn stoichiometry_is_an_attribute() {
    let xml = x("J1: A => 2 B; k*A;");
    has(
        &xml,
        r#"<speciesReference species="B" stoichiometry="2" constant="true"/>"#,
    );
    // An absent stoichiometry is 1.
    has(
        &xml,
        r#"<speciesReference species="A" stoichiometry="1" constant="true"/>"#,
    );
}

#[test]
fn variable_stoichiometry_gets_an_id() {
    // L3V2 removed StoichiometryMath: an id on the reference plus a matching
    // assignment rule, which the named stoichiometry already carries.
    let xml = x("S1_S2_stoich := 2*p1;\np1 = 4;\nreaction1: S1 => S1_S2_stoich S2; k*S1;");
    has(
        &xml,
        r#"<speciesReference id="S1_S2_stoich" species="S2" constant="false"/>"#,
    );
    has(&xml, r#"<assignmentRule variable="S1_S2_stoich">"#);
}

#[test]
fn constant_value_is_an_attribute_not_an_initial_assignment() {
    let xml = x("k = 0.1;");
    has(&xml, r#"value="0.1""#);
    lacks(&xml, "initialAssignment");
}

#[test]
fn non_constant_value_becomes_an_initial_assignment() {
    let xml = x("c = 2;\nS1 = 0.00015/c;");
    has(&xml, r#"<initialAssignment symbol="S1">"#);
    has(&xml, "<divide/>");
}

#[test]
fn assignment_and_rate_rules() {
    let xml = x("Vtot := Vmax*E;\nS1' = -k1*S1;");
    has(&xml, r#"<assignmentRule variable="Vtot">"#);
    has(&xml, r#"<rateRule variable="S1">"#);
}

#[test]
fn events_and_triggers() {
    let xml = x("E1: at 2 after S1 > 5, priority = 3: S2 = 0;");
    has(&xml, r#"<event id="E1""#);
    has(&xml, "<trigger");
    has(&xml, "<delay>");
    has(&xml, "<priority>");
    has(&xml, r#"<eventAssignment variable="S2">"#);
}

#[test]
fn function_definition_becomes_a_lambda() {
    let xml = x("function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend");
    has(&xml, r#"<functionDefinition id="MM">"#);
    has(&xml, "<lambda>");
    has(&xml, "<bvar><ci> S </ci></bvar>");
}

#[test]
fn conversion_factor_is_standard_l3v2() {
    let xml = x("species S1;\nS1 conversionFactor cf;\ncf = 2;");
    has(&xml, r#"conversionFactor="cf""#);
}

#[test]
fn boundary_species_and_substance_only() {
    let xml = x("substanceOnly species $S1 in c;\nc = 1; S1 = 3;");
    has(&xml, r#"hasOnlySubstanceUnits="true""#);
    has(&xml, r#"boundaryCondition="true""#);
    // substanceOnly means the value is an amount.
    has(&xml, r#"initialAmount="3""#);
}

#[test]
fn xml_is_escaped() {
    assert_eq!(
        xml_escape(r#"a & b <c> "d""#),
        "a &amp; b &lt;c&gt; &quot;d&quot;"
    );
}

#[test]
fn numbers_are_locale_independent() {
    let xml = x("k = 0.1;");
    has(&xml, r#"value="0.1""#);
    lacks(&xml, "0,1");
}

#[test]
fn identifiers_are_never_rewritten() {
    // A converter must not rename anything the author wrote, even an id
    // libantimony renamed on import to avoid a keyword.
    let xml = x(
        "compartment compartment_;\nspecies S1 in compartment_;\ncompartment_ = 1; S1 = 1; k = 1;\ncompartment_ is \"compartment\";\nJ1: S1 => S2; compartment_*k*S1;",
    );
    has(&xml, r#"<compartment id="compartment_""#);
    has(&xml, r#"compartment="compartment_""#);
    has(&xml, "<ci> compartment_ </ci>");
    // The original SBML id survives where it always was: as the name.
    has(&xml, r#"name="compartment""#);
}

#[test]
fn display_name_stays_a_name() {
    let xml = x("species S34;\nS34 = 1;\nS34 is \"Ethyl Alcohol\";");
    has(&xml, r#"id="S34""#);
    has(&xml, r#"name="Ethyl Alcohol""#);
}

// What the Delphi SBML comparison tier found.

#[test]
fn negative_literal_keeps_its_sign() {
    // `-0.1` is a minus over a number; reading the root as a number wrote
    // `value="0"` (corpus case 00124).
    let xml = x("k3 = -0.1;");
    has(&xml, r#"value="-0.1""#);
    lacks(&xml, r#"value="0""#);
}

#[test]
fn constant_expression_stays_an_expression() {
    // `3/2` could be folded, but that rounds what the author wrote (00912).
    let xml = x("c = 3/2;");
    has(&xml, r#"<initialAssignment symbol="c">"#);
    has(&xml, "<divide/>");
    lacks(&xml, r#"size="1.5""#);
}

#[test]
fn species_with_no_compartment_gets_a_default_one() {
    let xml = x("species X;\nX = 0;");
    has(&xml, r#"id="default_compartment""#);
    has(&xml, r#"sboTerm="SBO:0000410""#);
    has(&xml, r#"compartment="default_compartment""#);
}

#[test]
fn named_stoichiometry_carries_its_value() {
    // The id on a species reference is the quantity, so its value goes on
    // the reference (case 01434 had an id and no number anywhere).
    let xml = x("A1_sr = 2;\nJ0: A1_sr A -> ; -1;");
    has(&xml, r#"id="A1_sr""#);
    has(&xml, r#"stoichiometry="2""#);
}

#[test]
fn named_stoichiometry_is_constant_unless_a_rule_moves_it() {
    let xml = x("A1_sr = 2;\nJ0: A1_sr A -> ; -1;");
    has(
        &xml,
        r#"<speciesReference id="A1_sr" species="A" stoichiometry="2" constant="true"/>"#,
    );

    let xml = x("Xref := 3;\nJ0: -> Xref X; k1;");
    has(&xml, r#"constant="false""#);
}

#[test]
fn a_rule_with_no_expression_is_not_emitted() {
    // `p' = ;` is kept so the Antimony round trip closes, but a rule
    // without math is not SBML (cases 01236 and 01416).
    let xml = x("var p;\np = 3;\np' = ;");
    lacks(&xml, "<rateRule");
    lacks(&xml, "<listOfRules>");
}

// What the Delphi Tellurium tier found.

#[test]
fn events_carry_their_required_attributes() {
    // Required in L3 with no defaults (libSBML 21225/21226).
    let xml = x("E1: at time > 5: S1 = 3;");
    has(&xml, r#"useValuesFromTriggerTime="true""#);
    has(&xml, r#"initialValue="true""#);
    has(&xml, r#"persistent="true""#);

    // What the source does say still wins.
    let xml = x("E1: at time > 5, fromTrigger=false: S1 = 3;");
    has(&xml, r#"useValuesFromTriggerTime="false""#);
}

#[test]
fn modifiers_are_inferred_from_the_kinetic_law() {
    // S3 sets the rate but is neither consumed nor produced (rule 21121,
    // case 00063). It must be declared a species: a bare name in a law is a
    // parameter, and a parameter is not a modifier.
    let xml =
        x("species S1, S2, S3;\nreaction1: S1 => S2; k1*S1*S3;\nS1 = 1; S2 = 0; S3 = 2; k1 = 1.5;");
    has(&xml, "<listOfModifiers>");
    has(&xml, r#"<modifierSpeciesReference species="S3"/>"#);
    // A reactant is not also a modifier.
    lacks(&xml, r#"<modifierSpeciesReference species="S1"/>"#);
}

#[test]
fn a_unit_exponent_is_not_part_of_the_kind() {
    // Rule 20410 restricts kind to the base units (81 cases were invalid).
    let xml = x("unit area = metre^2;");
    has(&xml, r#"kind="metre" exponent="2""#);
    lacks(&xml, "metre^2");

    // The denominator's exponent is negated, not dropped.
    let xml = x("unit substance_per_area = mole / metre^2;");
    has(&xml, r#"kind="mole" exponent="1""#);
    has(&xml, r#"kind="metre" exponent="-2""#);
}

#[test]
fn an_unrepresentable_unit_is_empty_not_wrong() {
    let xml = x("unit odd = 1e-3 mole / (gram * 0.00027777 second);");
    has(&xml, r#"<unitDefinition id="odd"/>"#);

    // A reference to another unit definition is not a base unit either.
    let xml = x("unit volume = litre;\nunit v2 = volume;");
    has(&xml, r#"<unitDefinition id="v2"/>"#);
}

#[test]
fn the_sbml_prefix_is_declared() {
    // `<math>` makes MathML the default namespace, so a unit on a number
    // needs the declared `sbml:` prefix (case 00928).
    let xml = x("E1: at time > 5: S1 = 0.00015 mole;");
    has(&xml, "xmlns:sbml=");
    has(&xml, r#"sbml:units="mole""#);
}

#[test]
fn a_flux_bound_is_not_a_core_constraint() {
    // A bound names a reaction, which core SBML does not allow in a
    // constraint (case 01188).
    let xml = x("J0: A -> B; k1*A;\nk1 = 1;\nJ0 <= 1000;");
    lacks(&xml, "<listOfConstraints>");

    // Delphi passed whether or not an ordinary constraint went out; here
    // it does.
    let xml = x("S1 = 1;\nS1 < 10;");
    has(&xml, "<constraint>");
}

// Annotations and CV terms (Delphi M9).

#[test]
fn model_history_becomes_dcterms() {
    let xml = x(r#"model *m()
  S1 = 1;
  model created "2012-12-21T11:05:21Z"
  model modified "2013-01-15T12:25:55Z"
  model creator1.givenName "Brett"
  model creator1.familyName "Olivier"
  model creator1.email "bgoli@example.org"
  model creator1.organization "VU university Amsterdam"
end"#);
    has(
        &xml,
        "<dcterms:W3CDTF>2012-12-21T11:05:21Z</dcterms:W3CDTF>",
    );
    has(&xml, "<dcterms:created");
    has(&xml, "<dcterms:modified");
    has(&xml, "<vCard:Family>Olivier</vCard:Family>");
    has(&xml, "<vCard:Given>Brett</vCard:Given>");
    has(&xml, "<vCard:EMAIL>bgoli@example.org</vCard:EMAIL>");
    has(
        &xml,
        "<vCard:Orgname>VU university Amsterdam</vCard:Orgname>",
    );
}

#[test]
fn cv_terms_carry_their_qualifier() {
    let xml = x("S1 = 1;\nS1 identity \"http://identifiers.org/chebi/CHEBI:4167\"");
    has(&xml, "<bqbiol:is>");
    has(
        &xml,
        r#"<rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:4167"/>"#,
    );
    has(&xml, r##"rdf:about="#meta_S1""##);
    has(&xml, r#"metaid="meta_S1""#);
}

#[test]
fn terms_sharing_a_qualifier_share_a_bag() {
    // `part` and `hasPart` are two spellings of one qualifier.
    let xml = x(
        "S1 = 1;\nS1 part \"http://identifiers.org/chebi/CHEBI:1\"\nS1 hasPart \"http://identifiers.org/chebi/CHEBI:2\"",
    );
    assert_eq!(occurrences(&xml, "<bqbiol:hasPart>"), 1, "{xml}");
    assert_eq!(occurrences(&xml, "<rdf:li rdf:resource="), 2, "{xml}");
}

#[test]
fn model_qualifiers_use_the_model_namespace() {
    // By subject: `description` on a model is bqmodel, on a species bqbiol.
    // `taxon` stays bqbiol even on a model.
    let xml = x(r#"model *m()
  S1 = 1;
  model model_entity_is "http://identifiers.org/biomodels.db/BIOMD1"
  model description "http://identifiers.org/pubmed/1833774"
  model taxon "http://identifiers.org/taxonomy/8292"
  S1 description "http://identifiers.org/pubmed/9999999"
end"#);
    has(&xml, "<bqmodel:is>");
    has(&xml, "<bqmodel:isDescribedBy>");
    has(&xml, "<bqbiol:hasTaxon>");
    has(&xml, "<bqbiol:isDescribedBy>");
}

#[test]
fn notes_become_xhtml() {
    let xml =
        x("model *m()\n  S1 = 1;\n  model notes \"Maximize single objective function.\"\nend");
    has(&xml, "<notes>");
    has(
        &xml,
        r#"<p xmlns="http://www.w3.org/1999/xhtml">Maximize single objective function.</p>"#,
    );
}

#[test]
fn a_metaid_is_emitted_only_when_it_is_needed() {
    let xml = x("S1 = 1; S2 = 2;\nS1 identity \"http://identifiers.org/chebi/CHEBI:4167\"");
    has(&xml, r#"metaid="meta_S1""#);
    lacks(&xml, r#"metaid="meta_S2""#);

    let xml = x("S1 = 1;");
    lacks(&xml, "metaid");
    lacks(&xml, "<annotation>");
    lacks(&xml, "<notes>");
}

#[test]
fn annotation_text_is_escaped() {
    // Prose is escaped; `"` is left alone in element content.
    let xml = x(r#"model *m()
  S1 = 1;
  model notes "a < b & c > d, say \"hello\""
end"#);
    has(&xml, "a &lt; b &amp; c &gt; d");
    lacks(&xml, "&quot;");
}

#[test]
fn xhtml_notes_are_emitted_as_markup() {
    let xml = x(
        "model *m()\n  S1 = 1;\n  model notes ```\n<p>See <a href=\"http://example.org/x\">the paper</a>.</p>\n```\nend",
    );
    has(&xml, r#"<body xmlns="http://www.w3.org/1999/xhtml">"#);
    has(
        &xml,
        r#"<p>See <a href="http://example.org/x">the paper</a>.</p>"#,
    );
    lacks(&xml, "&lt;");
    lacks(&xml, "&quot;");
}

#[test]
fn reactions_events_and_units_carry_annotations_too() {
    let xml = x(r#"model *m()
  unit u = second;
  S1 = 1; S2 = 0; k = 1;
  J1: S1 -> S2; k*S1;
  E1: at (S1 > 5): S2 = 3;
  J1 identity "http://identifiers.org/ec-code/1.1.1.1"
  E1 description "http://identifiers.org/pubmed/12345"
  u identity "http://identifiers.org/uo/UO:0000010"
end"#);
    has(&xml, r#"metaid="meta_J1""#);
    has(
        &xml,
        r#"<rdf:li rdf:resource="http://identifiers.org/ec-code/1.1.1.1"/>"#,
    );
    has(&xml, r#"metaid="meta_E1""#);
    has(
        &xml,
        r#"<rdf:li rdf:resource="http://identifiers.org/pubmed/12345"/>"#,
    );
    has(&xml, r#"metaid="meta_u""#);
    has(
        &xml,
        r#"<rdf:li rdf:resource="http://identifiers.org/uo/UO:0000010"/>"#,
    );
}

#[test]
fn annotating_a_unit_does_not_invent_a_parameter() {
    let xml = x("unit u = second;\nS1 = 1;\nu identity \"http://identifiers.org/uo/UO:0000010\"");
    lacks(&xml, r#"<parameter id="u""#);
    has(&xml, r#"<unitDefinition metaid="meta_u" id="u""#);
}

#[test]
fn one_qualifier_may_carry_several_resources() {
    let xml = x(
        "ATP = 1;\nATP identity \"http://identifiers.org/chebi/CHEBI:15422\",\n             \"http://identifiers.org/kegg.compound/C00002\"",
    );
    assert_eq!(occurrences(&xml, "<bqbiol:is>"), 1, "{xml}");
    assert_eq!(occurrences(&xml, "<rdf:li rdf:resource="), 2, "{xml}");
    has(&xml, "CHEBI:15422");
    has(&xml, "C00002");
}

const MULTI_LINE_NOTES: &str = "model *m()\n  S1 = 1;\n  model notes ```\n<p>First \"quoted\" line.</p>\n<p>Second line.</p>\n```\nend";

#[test]
fn multi_line_notes_survive_an_antimony_round_trip() {
    let first = x(MULTI_LINE_NOTES);
    let antimony = write_document(&m(MULTI_LINE_NOTES));
    has(&antimony, "notes ```");
    assert_eq!(
        x(&antimony),
        first,
        "the note did not survive our own writer"
    );
}

#[test]
fn a_fenced_string_carries_multi_line_notes() {
    let xml = x(
        "model *m()\n  S1 = 1;\n  model notes ```\n<p>Converted from BIOMD0000000042.</p>\n<p>Second paragraph.</p>\n```\nend",
    );
    has(&xml, "<notes>");
    has(&xml, "<p>Converted from BIOMD0000000042.</p>");
    has(&xml, "Second paragraph");
}

// What the Delphi curated BioModels tier found.

#[test]
fn modifier_derivation_is_case_sensitive() {
    // BIOMD0000000557 holds both `duspn` and `DUSPn`.
    let xml = x(
        "compartment C = 1;\nspecies duspn in C, DUSPn in C;\nduspn = 1; DUSPn = 0; k = 1;\nJ1:  => DUSPn; k*duspn*C;",
    );
    has(&xml, r#"<modifierSpeciesReference species="duspn"/>"#);
    lacks(&xml, r#"<modifierSpeciesReference species="DUSPn"/>"#);
}

#[test]
fn a_declared_default_compartment_is_not_duplicated() {
    // Two compartments with one id is libSBML 10301.
    let xml = x(
        "compartment default_compartment;\nspecies As, Bs;\nAs = 1; Bs = 0; k = 1;\nR1: As -> Bs; k*As;",
    );
    assert_eq!(occurrences(&xml, "<compartment "), 1, "{xml}");
    has(&xml, r#"compartment="default_compartment""#);
}

#[test]
fn file_level_functions_reach_the_model() {
    // BIOMD0000000064 called functions that were nowhere in the document.
    let xml =
        x("function twice(x)\n  x*2\nend\n\nmodel *m()\n  S1 = 1;\n  J1: S1 => ; twice(S1);\nend");
    has(&xml, "<listOfFunctionDefinitions>");
    has(&xml, r#"<functionDefinition id="twice""#);
}

#[test]
fn the_model_itself_can_carry_an_sbo_term() {
    let document = m("model *m()\n  S1 = 1;\n  model.sboTerm = 167\nend");
    let module = document.main_module();
    assert_eq!(module.sbo_term, Some(167));
    // And no phantom parameter called "model".
    assert!(module.symbol("model").is_none());
}

#[test]
fn an_annotation_by_the_models_own_name_makes_no_parameter() {
    let source = "model *m()\n  S1 = 1;\n  m identity \"http://identifiers.org/biomodels.db/MODEL0001\"\nend";
    let document = m(source);
    let module = document.main_module();
    assert!(
        module.symbol("m").is_none(),
        "no phantom parameter called m"
    );
    assert_eq!(module.cv_terms.len(), 1, "the model's own annotation");
    let xml = x(source);
    assert!(!xml.contains("<parameter id=\"m\""), "{xml}");
    assert!(xml.contains("biomodels.db/MODEL0001"), "{xml}");
}

#[test]
fn a_functions_annotation_reaches_the_sbml() {
    let source = "function mm(x)\n  x/(1+x)\nend\n\nmodel *m()\n  S1 = 1;\n  \
         mm identity \"http://identifiers.org/sbo/SBO:0000028\"\n  \
         mm notes \"Michaelis-Menten\"\nend";
    let xml = x(source);
    let function = xml
        .split("<functionDefinition")
        .nth(1)
        .and_then(|f| f.split("</functionDefinition>").next())
        .expect("a function definition");
    has(function, "http://identifiers.org/sbo/SBO:0000028");
    has(function, "Michaelis-Menten");
    // Annotating the function adds no parameter for it.
    assert!(!xml.contains("<parameter id=\"mm\""), "{xml}");
}

#[test]
fn sbo_terms_reach_the_sbml() {
    // Every `X.sboTerm = n` must come out as an sboTerm attribute, with or
    // without other annotations on the same element.
    let source = "function mm(x)\n  x/(1+x)\nend\n\nmodel *m()\n  \
         compartment cyt = 1;\n  species Glc in cyt = 1, G6P in cyt = 0;\n  \
         J0: Glc -> G6P; k1*mm(Glc);\n  k1 = 0.1;\n  E1: at time > 1: k1 = 0.2;\n  \
         model.sboTerm = 4\n  cyt.sboTerm = 290\n  Glc.sboTerm = 247\n  \
         J0.sboTerm = 176\n  k1.sboTerm = 9\n  E1.sboTerm = 231\n  mm.sboTerm = 28\n  \
         Glc identity \"http://identifiers.org/chebi/CHEBI:4167\"\nend";
    // Valid, and reads back.
    let xml = x(source);
    assert_eq!(occurrences(&xml, "sboTerm="), 7, "{xml}");

    let sbml = to_sbml(&flatten(&m(source)));
    let sbo = |meta: &sbml_io::Meta| meta.sbo_term.clone();
    let model = sbml.model.expect("a model");
    assert_eq!(sbo(&model.meta).as_deref(), Some("SBO:0000004"));
    assert_eq!(
        sbo(&model.compartments[0].meta).as_deref(),
        Some("SBO:0000290")
    );
    let species = |id: &str| model.species.iter().find(|s| s.id == id).unwrap();
    // With a CV term beside it, and without.
    assert_eq!(sbo(&species("Glc").meta).as_deref(), Some("SBO:0000247"));
    assert_eq!(sbo(&species("G6P").meta), None);
    assert_eq!(
        sbo(&model.reactions[0].meta).as_deref(),
        Some("SBO:0000176")
    );
    assert_eq!(
        sbo(&model.parameters[0].meta).as_deref(),
        Some("SBO:0000009")
    );
    assert_eq!(sbo(&model.events[0].meta).as_deref(), Some("SBO:0000231"));
    assert_eq!(
        sbo(&model.function_definitions[0].meta).as_deref(),
        Some("SBO:0000028")
    );
}

#[test]
fn notes_keep_their_trailing_whitespace() {
    let source = "model *m()\n  S1 = 1;\n  model notes ```\na line with trailing space   \nand another\n```\nend";
    assert_eq!(
        x(source),
        x(&write_document(&m(source))),
        "the note changed on the way through our own writer"
    );
}

#[test]
fn line_numbers_survive_a_fenced_note() {
    // The note spans lines 2-5; `k` is undeclared and on line 7.
    let loaded =
        antimony::parse("model *m()\n  model notes ```\none\ntwo\n```\n  S1 = 1;\n  x := k;\nend");
    let found = loaded
        .diagnostics
        .first_of(Code::ImplicitParameter)
        .unwrap_or_else(|| panic!("{}", loaded.diagnostics));
    assert_eq!(found.span.start.line, 7, "{}", loaded.diagnostics);
}

// Interactions.

/// The species appear in the law: a modifier the law never mentions would
/// mean nothing in SBML.
const INTERACTION_MODEL: &str = "model *m()
  S1 = 1; S2 = 0; k = 1; G1 = 1; G2 = 1; G3 = 1;
  J1: S1 -> S2; k*S1*G1*G2*G3;
  i1: G1 -| J1;
  i2: G2 -o J1;
  G3 -( J1;
end";

#[test]
fn interactions_become_modifiers_with_sbo_terms() {
    let xml = x(INTERACTION_MODEL);
    has(
        &xml,
        r#"<modifierSpeciesReference sboTerm="SBO:0000020" name="i1" species="G1"/>"#,
    );
    has(
        &xml,
        r#"<modifierSpeciesReference sboTerm="SBO:0000459" name="i2" species="G2"/>"#,
    );
    // Unlabelled: no name, and none invented.
    has(
        &xml,
        r#"<modifierSpeciesReference sboTerm="SBO:0000019" species="G3"/>"#,
    );
    // Each species once, though each is in the law as well.
    assert_eq!(occurrences(&xml, "<modifierSpeciesReference"), 3, "{xml}");
}

#[test]
fn an_interaction_is_not_a_reaction() {
    let loaded = antimony::parse(INTERACTION_MODEL);
    let module = loaded.document.main_module();
    assert_eq!(module.reactions.len(), 1, "J1 only");
    assert_eq!(module.interactions.len(), 3);
    assert_eq!(module.interactions[0].reaction, "J1");
    assert_eq!(module.interactions[0].species, "G1");
    assert_eq!(module.interactions[0].kind, InteractionKind::Inhibits);
    assert_eq!(module.interactions[1].kind, InteractionKind::Activates);
    assert_eq!(module.interactions[2].kind, InteractionKind::Unknown);
    assert_eq!(module.symbol("J1").unwrap().kind, SymbolKind::Reaction);
    assert_eq!(module.symbol("i1").unwrap().kind, SymbolKind::Interaction);
    // A label is a name, not a quantity: no warning that it has no value.
    assert!(
        !loaded.diagnostics.has_any(Code::UninitialisedSymbol),
        "{}",
        loaded.diagnostics
    );
}

#[test]
fn interactions_round_trip() {
    let antimony = write_document(&m(INTERACTION_MODEL));
    has(&antimony, "i1: G1 -| J1;");
    has(&antimony, "i2: G2 -o J1;");
    has(&antimony, "G3 -( J1;");
    assert_eq!(
        x(&antimony),
        x(INTERACTION_MODEL),
        "the interactions did not survive our own writer"
    );
}

// --- TMathMLTests ----------------------------------------------------------------

#[test]
fn arithmetic_becomes_apply() {
    let xml = x("J1: A => B; k*A + 2;");
    has(&xml, "<plus/>");
    has(&xml, "<times/>");
    has(&xml, "<ci> k </ci>");
    has(&xml, "<cn> 2 </cn>");
}

#[test]
fn time_is_a_csymbol() {
    let xml = x("x := time*2;");
    has(
        &xml,
        r#"definitionURL="http://www.sbml.org/sbml/symbols/time""#,
    );
    lacks(&xml, "<ci> time </ci>");
}

#[test]
fn avogadro_is_a_csymbol() {
    has(
        &x("x := avogadro;"),
        r#"definitionURL="http://www.sbml.org/sbml/symbols/avogadro""#,
    );
}

#[test]
fn delay_and_rate_of_are_csymbols() {
    has(
        &x("x := delay(y, 1);"),
        r#"definitionURL="http://www.sbml.org/sbml/symbols/delay""#,
    );
    has(
        &x("p2 = rateOf(p1);"),
        r#"definitionURL="http://www.sbml.org/sbml/symbols/rateOf""#,
    );
}

#[test]
fn named_constants() {
    let xml = x("x := pi + exponentiale;");
    has(&xml, "<pi/>");
    has(&xml, "<exponentiale/>");
}

#[test]
fn infinity_and_nan() {
    let xml = x("P := INF; R := NaN; Q := -INF;");
    has(&xml, "<infinity/>");
    has(&xml, "<notanumber/>");
}

#[test]
fn sqrt_becomes_root_degree_two() {
    let xml = x("x := sqrt(y);");
    has(&xml, "<root/>");
    has(&xml, "<degree>");
    lacks(&xml, "<sqrt/>");
}

#[test]
fn log10_becomes_log_with_base() {
    let xml = x("x := log10(y);");
    has(&xml, "<log/>");
    has(&xml, "<logbase>");
}

#[test]
fn piecewise_becomes_pieces() {
    let xml = x("x := piecewise(2, y > 1, 3);");
    has(&xml, "<piecewise>");
    has(&xml, "<piece>");
    has(&xml, "<otherwise>");
}

#[test]
fn conditional_becomes_piecewise() {
    let xml = x("x := y > 1 ? 2 : 3;");
    has(&xml, "<piecewise>");
    has(&xml, "<otherwise>");
}

#[test]
fn user_function_call_uses_ci() {
    let xml = x("function f(a)\n  a*2\nend\nx := f(3);");
    has(&xml, "<ci> f </ci>");
}

#[test]
fn unit_bearing_literal_carries_its_unit() {
    has(&x("x := 5 mole;"), r#"sbml:units="mole""#);
}

// --- Added in the Rust port ------------------------------------------------------

#[test]
fn a_rule_or_event_makes_a_parameter_variable() {
    // SBML rules 20903 and 21213; the Delphi writer read the declaration
    // alone and wrote `constant="true"` (SPEC.md R30).
    let model = model(
        "x := 2*k; k = 1;\ny = 1; y' = 1;\nz = 0;\nE1: at time > 1: z = 1;\ncompartment c = 1; c := 2;\nw = 3;",
    );
    let constant = |id: &str| {
        model
            .parameters
            .iter()
            .find(|p| p.id == id)
            .and_then(|p| p.constant)
            .or_else(|| {
                model
                    .compartments
                    .iter()
                    .find(|c| c.id == id)
                    .and_then(|c| c.constant)
            })
    };
    assert_eq!(constant("x"), Some(false));
    assert_eq!(constant("y"), Some(false));
    assert_eq!(constant("z"), Some(false));
    assert_eq!(constant("c"), Some(false));
    assert_eq!(constant("w"), Some(true));
    assert_eq!(constant("k"), Some(true));
}

#[test]
fn root_and_log_take_their_qualifier_first() {
    // `root(n, x)` and `log(b, x)`: the first argument is the degree or
    // base, as a qualifier. Delphi wrote both as plain arguments, which is
    // not valid MathML (SPEC.md R31).
    use sbml_io::math::Math;
    let model = model("a := root(3, y);\nb := log(2, y);\nc := log(y);\nd := root(y);\ny = 8;");
    let math = |id: &str| {
        model
            .rules
            .iter()
            .find_map(|r| match r {
                sbml_io::Rule::Assignment { variable, math, .. } if variable == id => math.clone(),
                _ => None,
            })
            .unwrap()
    };
    let ident = |name: &str| Box::new(Math::Ident(name.into()));
    let three = Box::new(Math::Number {
        value: sbml_io::math::Number::Real(3.0),
        units: None,
    });
    let two = Box::new(Math::Number {
        value: sbml_io::math::Number::Real(2.0),
        units: None,
    });
    assert_eq!(
        math("a"),
        Math::Root {
            degree: Some(three),
            arg: ident("y")
        }
    );
    assert_eq!(
        math("b"),
        Math::Log {
            base: Some(two),
            arg: ident("y")
        }
    );
    assert_eq!(
        math("c"),
        Math::Log {
            base: None,
            arg: ident("y")
        }
    );
    assert_eq!(
        math("d"),
        Math::Root {
            degree: None,
            arg: ident("y")
        }
    );
}

#[test]
fn an_unnamed_event_gets_an_id() {
    // SBML events need no id, but the Delphi writer gave each one `_E1`,
    // `_E2`, … and the differential tier compares by id.
    let model = model("at time > 1: x = 1;\nat time > 2: x = 2;\nx = 0;");
    let ids: Vec<_> = model.events.iter().map(|e| e.id.clone().unwrap()).collect();
    assert_eq!(ids, ["_E1", "_E2"]);
}

#[test]
fn modifiers_are_listed_in_order_of_name() {
    let model = model("species A, B, Z, M;\nJ1: A => B; Z*M*A;\nA = 1; B = 0; Z = 1; M = 1;");
    let names: Vec<_> = model.reactions[0]
        .modifiers
        .iter()
        .map(|m| m.species.as_str())
        .collect();
    assert_eq!(names, ["M", "Z"]);
}

#[test]
fn an_empty_module_writes_an_empty_model() {
    let module = Module::default();
    let document = to_sbml(&module);
    assert_eq!((document.level, document.version), (3, 2));
    assert!(sbml_io::validate(&document).is_empty());
    sbml_io::read(&sbml_io::write(&document)).unwrap();
}

#[test]
fn every_corpus_file_converts_to_valid_sbml() {
    // Every case that loads, `comp` aside (flattening does not yet expand
    // submodels), gives a document that passes `sbml_io::validate` and
    // reads back. libSBML's checks are run on the saved files (SPEC.md §8).
    let Some(cases) = common::corpus_index("every_corpus_file_converts_to_valid_sbml") else {
        return;
    };
    let out = std::env::var_os("ANTIMONY_SBML_OUT").map(std::path::PathBuf::from);
    if let Some(dir) = &out {
        std::fs::create_dir_all(dir).expect("cannot create ANTIMONY_SBML_OUT");
    }
    let mut tried = 0;
    let mut failed = Vec::new();
    for case in cases.iter().filter(|c| !c.has_package("comp")) {
        let loaded = antimony::parse(&common::read(&case.path));
        if loaded.has_errors() {
            continue;
        }
        tried += 1;
        let sbml = to_sbml(&flatten(&loaded.document));
        let issues = sbml_io::validate(&sbml);
        if !issues.is_empty() {
            failed.push(format!("{}: {}", case.id, issues[0]));
            continue;
        }
        let text = sbml_io::write(&sbml);
        if let Err(problems) = sbml_io::read(&text) {
            failed.push(format!("{}: does not read back: {problems:?}", case.id));
            continue;
        }
        if let Some(dir) = &out {
            std::fs::write(dir.join(format!("{}.xml", case.id)), &text).expect("cannot save");
        }
    }
    assert_eq!(tried, 1657);
    assert!(
        failed.is_empty(),
        "{} of {tried} failed:\n{}",
        failed.len(),
        failed.join("\n")
    );
}

fn losses(source: &str) -> Vec<String> {
    let flat = flatten(&m(source));
    antimony::sbml::information_loss(&flat)
        .iter()
        .map(|d| {
            assert_eq!(d.code, antimony::Code::InformationLoss);
            assert_eq!(d.severity, antimony::Severity::Warning);
            d.message.clone()
        })
        .collect()
}

#[test]
fn what_the_sbml_cannot_carry_is_reported() {
    // Every kind of annotation written: nothing to report.
    let written = "function mm(x)\n  x/(1+x)\nend\n\nmodel *m()\n  \
         compartment cyt = 1;\n  species Glc in cyt = 1;\n  J0: Glc -> ; k1*mm(Glc);\n  \
         k1 = 0.1;\n  E1: at time > 1: k1 = 0.2;\n  model.sboTerm = 4\n  \
         cyt.sboTerm = 290\n  J0.sboTerm = 176\n  E1.sboTerm = 231\n  \
         Glc identity \"http://identifiers.org/chebi/CHEBI:4167\"\n  \
         mm identity \"http://identifiers.org/sbo/SBO:0000028\"\n  \
         k1 notes \"a rate constant\"\n  \
         m identity \"http://identifiers.org/biomodels.db/MODEL0001\"\nend";
    assert_eq!(losses(written), Vec::<String>::new());

    // A vocabulary the library does not know, and an interaction, which
    // becomes a modifier with no annotation of its own: each named.
    let lost = losses(
        "S1 = 1\nS1 frobnicates \"http://example.org/x\"\nJ0: -> S1; 1\n\
         J1: S1 -| J0\nJ1 identity \"http://example.org/y\"",
    );
    assert_eq!(
        lost,
        [
            "annotations on J1 are not written to SBML",
            "annotations on S1 are not written to SBML"
        ]
    );
}

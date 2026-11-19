//! Antimony 3's layout statements (docs/antimony-layout-design.md, step 1):
//! read, kept in order, and written back. The examples are the reference's
//! ("Layout and Render", https://tellurium.readthedocs.io/en/latest/antimony.html)
//! and libAntimony 3.2.1's own output.

use antimony::model::{LayoutStatement, LayoutValue, Module};
use antimony::{Loaded, write_document};

fn load(source: &str) -> Loaded {
    antimony::parse(source)
}

/// A load that must have no errors.
fn clean(source: &str) -> Loaded {
    let loaded = load(source);
    assert!(!loaded.has_errors(), "{}", loaded.diagnostics);
    loaded
}

fn main(loaded: &Loaded) -> &Module {
    loaded.document.main_module()
}

fn value(module: &Module, name: &str) -> String {
    module
        .symbol(name)
        .and_then(|s| s.initial_value.as_ref())
        .map(|e| e.to_string())
        .unwrap_or_else(|| panic!("{name} has no value"))
}

/// The statements as `path = value` text, as the writer would put them.
fn statements(module: &Module) -> Vec<(String, LayoutValue)> {
    module
        .layout
        .iter()
        .map(|s: &LayoutStatement| (s.path.join("."), s.value.clone()))
        .collect()
}

fn point(x: f64, y: f64) -> LayoutValue {
    LayoutValue::Point(x, y)
}

fn word(w: &str) -> LayoutValue {
    LayoutValue::Word(w.into())
}

// --- the fault this step fixes --------------------------------------------------

#[test]
fn values_after_model_layout_are_kept() {
    // The old `layout … end` rule swallowed everything after `layout`.
    let loaded =
        clean("J0: S1 -> S2; k1*S1\nk1 = 0.1\nmodel.layout = on\nS1.pos = {28, 35}\nS2 = 3\n");
    let m = main(&loaded);
    assert_eq!(value(m, "S2"), "3");
    assert_eq!(value(m, "k1"), "0.1");
    assert_eq!(
        statements(m),
        [
            ("model.layout".into(), word("on")),
            ("S1.pos".into(), point(28.0, 35.0)),
        ]
    );
}

#[test]
fn layout_inside_a_model_block_leaves_its_end_alone() {
    let loaded = clean(
        "model m\n  J0: S1 -> S2; k1*S1\n  model.layout = on\n  S1.pos = {28, 35}\n  k1 = 0.1; S2 = 3\nend\n",
    );
    let m = loaded.document.module("m").expect("model m");
    assert_eq!(value(m, "S2"), "3");
    assert_eq!(value(m, "k1"), "0.1");
    assert_eq!(m.layout.len(), 2);
}

#[test]
fn layout_statements_declare_no_symbols() {
    let loaded = clean("J0: S1 -> S2; k1*S1\nS1.pos = {1, 2}\nJ0.S1.b1 = {3, 4}\nk1 = 1\n");
    let m = main(&loaded);
    for name in ["S1.pos", "pos", "J0.S1.b1", "J0.S1", "b1"] {
        assert!(m.symbol(name).is_none(), "{name} became a symbol");
    }
}

// --- every form in the reference --------------------------------------------------

#[test]
fn every_statement_form_reads() {
    let loaded = clean(concat!(
        "J0: S1 -> S2; k1*S1\n",
        "k1 = 1\n",
        "model.layout = on\n",
        "model.autolayout.maxNumConnectedEdges = 5\n",
        "model.layout.style = calm\n",
        "model.layout.align_top = {S1, S2}\n",
        "model.layout.size = {505, 404}\n",
        "model.layout.width = 505\n",
        "model.layout.background = fuchsia\n",
        "model.layout.fontstyle = italic\n",
        "species.size = {68, 101}\n",
        "species.color = blue\n",
        "reaction.size = {3, 3}\n",
        "compartment.linewidth = 2\n",
        "S1.pos = {28, 35}\n",
        "S1.x = 28\n",
        "S1.y = -35.5\n",
        "S1.color = \"#268F49CC\"\n",
        "S1.lineColor = black\n",
        "S1.fontstyle = bold_italic\n",
        "S1.font = sans_serif\n",
        "S1.shape = ellipse\n",
        "J0.position = {80.54, 153.07}\n",
        "J0.S1.species_end = {50.89, 74.02}\n",
        "J0.S1.rxn_end = {72.09, 165.65}\n",
        "J0.S1.b1 = {71.04, 129.94}\n",
        "J0.S1.arc2.b2 = {177.02, 245.98}\n",
        "J0.S1.seg2.species_end = {685, 1008}\n",
        "J0.S1.arc2.seg3.species_end = {740, 992}\n",
        "J0.--.position = {520, 483}\n",
        "S1.position.J0 = {349.69, 1024.82}\n",
        "S1.size.J0 = {80, 50}\n",
    ));
    let got = statements(main(&loaded));
    let want: Vec<(String, LayoutValue)> = vec![
        ("model.layout".into(), word("on")),
        (
            "model.autolayout.maxNumConnectedEdges".into(),
            LayoutValue::Number(5.0),
        ),
        ("model.layout.style".into(), word("calm")),
        (
            "model.layout.align_top".into(),
            LayoutValue::List(vec!["S1".into(), "S2".into()]),
        ),
        ("model.layout.size".into(), point(505.0, 404.0)),
        ("model.layout.width".into(), LayoutValue::Number(505.0)),
        ("model.layout.background".into(), word("fuchsia")),
        ("model.layout.fontstyle".into(), word("italic")),
        ("species.size".into(), point(68.0, 101.0)),
        ("species.color".into(), word("blue")),
        ("reaction.size".into(), point(3.0, 3.0)),
        ("compartment.linewidth".into(), LayoutValue::Number(2.0)),
        ("S1.pos".into(), point(28.0, 35.0)),
        ("S1.x".into(), LayoutValue::Number(28.0)),
        ("S1.y".into(), LayoutValue::Number(-35.5)),
        ("S1.color".into(), LayoutValue::Text("#268F49CC".into())),
        ("S1.lineColor".into(), word("black")),
        ("S1.fontstyle".into(), word("bold_italic")),
        ("S1.font".into(), word("sans_serif")),
        ("S1.shape".into(), word("ellipse")),
        ("J0.position".into(), point(80.54, 153.07)),
        ("J0.S1.species_end".into(), point(50.89, 74.02)),
        ("J0.S1.rxn_end".into(), point(72.09, 165.65)),
        ("J0.S1.b1".into(), point(71.04, 129.94)),
        ("J0.S1.arc2.b2".into(), point(177.02, 245.98)),
        ("J0.S1.seg2.species_end".into(), point(685.0, 1008.0)),
        ("J0.S1.arc2.seg3.species_end".into(), point(740.0, 992.0)),
        ("J0.--.position".into(), point(520.0, 483.0)),
        ("S1.position.J0".into(), point(349.69, 1024.82)),
        ("S1.size.J0".into(), point(80.0, 50.0)),
    ];
    assert_eq!(got, want);
}

#[test]
fn several_statements_on_one_line() {
    let loaded = clean("S1 = 1\nS1.x = 28; S1.y = 35; S1 = 2\n");
    let m = main(&loaded);
    assert_eq!(m.layout.len(), 2);
    assert_eq!(value(m, "S1"), "2");
}

#[test]
fn a_style_name_of_several_words() {
    let loaded = clean("model.layout.style = blue ombre\n");
    assert_eq!(statements(main(&loaded))[0].1, word("blue ombre"));
}

#[test]
fn a_braced_list_may_break_lines() {
    let loaded = clean("model.layout.align_right = {S3,\n  S4, S5}\nk1 = 1\n");
    let m = main(&loaded);
    assert_eq!(
        statements(m)[0].1,
        LayoutValue::List(vec!["S3".into(), "S4".into(), "S5".into()])
    );
    assert_eq!(value(m, "k1"), "1");
}

#[test]
fn comments_stay_with_their_statement() {
    let loaded = clean("# General layout options\nmodel.layout = on\n");
    let text = write_document(&loaded.document);
    assert!(
        text.contains("General layout options\nmodel.layout = on"),
        "{text}"
    );
}

// --- what is not a layout statement -----------------------------------------------

#[test]
fn sbo_terms_are_not_layout() {
    let loaded = clean("S1 = 1\nS1.sboTerm = 247\nmodel.sboTerm = 4\n");
    let m = main(&loaded);
    assert!(m.layout.is_empty());
    assert_eq!(m.symbol("S1").unwrap().sbo_term, Some(247));
    assert_eq!(m.sbo_term, Some(4));
}

#[test]
fn a_dotted_name_without_a_layout_word_is_not_layout() {
    // As before this step: a dotted name is one symbol.
    let loaded = load("A.k1 = 3\n");
    let m = main(&loaded);
    assert!(m.layout.is_empty());
    assert!(m.symbol("A.k1").is_some());
}

#[test]
fn a_symbol_called_layout_is_still_a_symbol() {
    let loaded = clean("layout = 3\nk1 = layout * 2\n");
    let m = main(&loaded);
    assert_eq!(value(m, "layout"), "3");
    assert!(m.layout.is_empty());
}

// --- errors -----------------------------------------------------------------------

#[test]
fn a_statement_with_no_value_is_one_error_and_the_rest_reads() {
    let loaded = load("S1.pos = \nk1 = 1\n");
    assert!(loaded.has_errors());
    assert_eq!(
        loaded
            .diagnostics
            .iter()
            .filter(|d| d.severity == antimony::Severity::Error)
            .count(),
        1
    );
    assert_eq!(value(main(&loaded), "k1"), "1");
    assert!(main(&loaded).layout.is_empty());
}

#[test]
fn an_unclosed_brace_is_an_error() {
    let loaded = load("S1.pos = {28, 35\nk1 = 1\n");
    assert!(loaded.has_errors());
}

// --- writing back -------------------------------------------------------------------

/// libAntimony 3.2.1's output for a model with aliases, a null node and a
/// compartment (2026-10-22).
const LIBANTIMONY: &str = "\
J0: -> S1; k0
J1: S1 -> S2; k1*S1
J2: S1 -> S3; k2*S1
k0 = 1; k1 = 1; k2 = 1; S1 = 0; S2 = 0; S3 = 0
compartment C = 1; S1 in C; S2 in C; S3 in C

# General layout options
model.layout = on
model.layout.size = {838.83, 410}

// Individual element layout information
C.position = {10, 10}
C.size = {500, 400}
S1.position = {100, 100}
S1.position.J2 = {100, 300}
J0.position = {72.5, 66.5}
J0.S1.species_end = {116.93, 106.58}
J0.S1.b1 = {112.62, 102.82}
J0.S1.b2 = {112.62, 102.82}
J0.--.species_end = {57.53, 54.68}
J0.--.b1 = {52.38, 50.18}
J0.--.position = {20, 20}
J0.--.color = white
J0.--.lineColor = black
J0.--.lineWidth = 2
";

#[test]
fn libantimony_output_reads_cleanly() {
    let loaded = clean(LIBANTIMONY);
    assert_eq!(main(&loaded).layout.len(), 16);
}

#[test]
fn written_statements_read_back_the_same() {
    let first = clean(LIBANTIMONY);
    let text = write_document(&first.document);
    let second = clean(&text);
    assert_eq!(statements(main(&first)), statements(main(&second)));
    // The first write reorders a file (write.rs); after that, writing again
    // gives the same text.
    let again = write_document(&second.document);
    assert_eq!(write_document(&clean(&again).document), again);
}

#[test]
fn written_form() {
    let loaded = clean(concat!(
        "S1 = 1\n",
        "S1.pos = {28, -3.5}\n",
        "S1.color = \"#268F49CC\"\n",
        "model.layout.align_top = {S1}\n",
        "model.layout.style = blue ombre\n",
        "J0.--.position = {1, 2}\n",
    ));
    let text = write_document(&loaded.document);
    for line in [
        "S1.pos = {28, -3.5}",
        "S1.color = \"#268F49CC\"",
        "model.layout.align_top = {S1}",
        "model.layout.style = blue ombre",
        "J0.--.position = {1, 2}",
    ] {
        assert!(
            text.lines().any(|l| l.trim() == line),
            "{line} missing from:\n{text}"
        );
    }
}

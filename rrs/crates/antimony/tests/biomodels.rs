//! The BioModels tier (Delphi `Tools/biomodels_tier.py`; SPEC.md §8).
//!
//! The corpus is synthetic: the SBML Test Suite exercises SBML one feature
//! at a time, and its Antimony was written by libantimony. Curated BioModels
//! are real published models, with prose notes, author lists, CV terms and
//! vocabulary the corpus has no reason to contain; the Delphi tier found
//! five defects there that every other tier had missed.
//!
//! `ANTIMONY_BIOMODELS` is a folder of `<id>.ant` files, made from the
//! curated SBML by `tools/biomodels_to_antimony.py` (libantimony converts
//! 708 of the 717; the 9 it cannot are skipped, as in Delphi). Unset, the
//! test is skipped with a note. For each model:
//!
//! - **loads**: no error diagnostics;
//! - **round-trip**: the Antimony this crate writes reads back as the same
//!   model (the deep comparison of `common/compare.rs`), and the second
//!   write equals the third;
//! - with the `sbml` feature, **sbml**: the SBML passes `sbml_io::validate`
//!   and reads back, and is the same SBML after a trip through this crate's
//!   own Antimony writer (Delphi's check). With `ANTIMONY_SBML_OUT` set the
//!   SBML is saved there as `<id>.xml`, for `tools/libsbml_check.py`.

mod common;

use std::path::PathBuf;

use antimony::{parse, write_document};

/// The number of curated models libantimony converts.
const CONVERTED: usize = 708;

/// What is wrong with one model, or `None`.
fn check(source: &str) -> Option<String> {
    let loaded = parse(source);
    if loaded.has_errors() {
        let first = loaded
            .diagnostics
            .iter()
            .find(|d| d.severity == antimony::Severity::Error);
        return Some(format!("loads: {}", first.unwrap()));
    }
    let once = write_document(&loaded.document);
    let again = parse(&once);
    if again.has_errors() {
        return Some(format!(
            "round-trip: our Antimony does not read back: {}",
            again.diagnostics
        ));
    }
    let differences =
        common::compare::differences(loaded.document.main_module(), again.document.main_module());
    if !differences.is_empty() {
        return Some(format!("round-trip: {}", differences[0]));
    }
    let twice = write_document(&again.document);
    if write_document(&parse(&twice).document) != twice {
        return Some("round-trip: the second write differs from the third".into());
    }
    sbml(&loaded.document, &again.document)
}

#[cfg(feature = "sbml")]
fn sbml(first: &antimony::model::Document, again: &antimony::model::Document) -> Option<String> {
    use antimony::flatten::flatten;
    use antimony::sbml::{to_sbml, write_sbml};
    let issues = sbml_io::validate(&to_sbml(&flatten(first)));
    if !issues.is_empty() {
        return Some(format!("sbml: {}", issues[0]));
    }
    let text = write_sbml(first);
    if let Err(problems) = sbml_io::read(&text) {
        return Some(format!("sbml: does not read back: {problems:?}"));
    }
    if write_sbml(again) != text {
        return Some("sbml: differs after a trip through our Antimony writer".into());
    }
    None
}

#[cfg(not(feature = "sbml"))]
fn sbml(_: &antimony::model::Document, _: &antimony::model::Document) -> Option<String> {
    None
}

/// Saves a model's SBML to `ANTIMONY_SBML_OUT`, if it is set.
#[cfg(feature = "sbml")]
fn save(id: &str, source: &str) {
    if let Some(dir) = std::env::var_os("ANTIMONY_SBML_OUT") {
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).expect("cannot create ANTIMONY_SBML_OUT");
        let text = antimony::sbml::write_sbml(&parse(source).document);
        std::fs::write(dir.join(format!("{id}.xml")), text).expect("cannot save");
    }
}

#[cfg(not(feature = "sbml"))]
fn save(_: &str, _: &str) {}

#[test]
fn every_curated_biomodel_loads_round_trips_and_converts() {
    let Some(dir) = std::env::var_os("ANTIMONY_BIOMODELS") else {
        eprintln!("note: ANTIMONY_BIOMODELS is not set; the BioModels tier skipped");
        return;
    };
    let files = common::files_with(&PathBuf::from(dir), "ant");
    let mut failed = Vec::new();
    for path in &files {
        let id = path.file_stem().unwrap().to_string_lossy().into_owned();
        let source = common::read(path);
        match check(&source) {
            Some(problem) => failed.push(format!("{id}: {problem}")),
            None => save(&id, &source),
        }
    }
    assert_eq!(files.len(), CONVERTED, "not the converted set");
    assert!(
        failed.is_empty(),
        "{} of {} failed:\n{}",
        failed.len(),
        files.len(),
        failed.join("\n")
    );
}

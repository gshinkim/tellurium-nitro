//! The differential tier (SPEC.md §8): this crate against the Delphi
//! implementation, case by case over the corpus.
//!
//! ```text
//! cargo run --release -p antimony-rs --features sbml --example differential
//! cargo run --release -p antimony-rs --features sbml --example differential -- 00001 01234
//! ```
//!
//! Needs `ANTIMONY_CORPUS`, `ANTIMONY_CORPUS_INDEX` and `ANTIMONY_ORACLE`
//! (the Delphi `ant2sbml.exe`, so in practice Windows). The oracle's output
//! is cached under `target/antimony-oracle/`; delete that folder after
//! rebuilding the oracle.
//!
//! Two things are compared, on every case of the corpus:
//!
//! - **Diagnostics**, every one, by line and code: the oracle's on stderr
//!   against `antimony::parse`'s. This covers the `comp` cases and the ones
//!   the oracle rejects as well.
//! - **The SBML**, where the oracle wrote some and the case is outside the
//!   `comp` package: the oracle's against this crate's `write_sbml` of the
//!   same file. Both are read with `sbml-io`, so they are compared as
//!   models, never as text, and then element by element (the model's own
//!   attributes and annotation, each function, unit, compartment, species,
//!   parameter, initial assignment, rule, constraint, reaction and event),
//!   every attribute and every MathML node. The only allowance is whitespace
//!   in raw notes and annotations. A divergence shows where the two first
//!   differ.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::{env, fs};

/// One case to compare.
struct Case {
    id: String,
    path: PathBuf,
    comp: bool,
}

/// What the oracle made of a case.
struct Oracle {
    /// `None` when it rejected the case.
    sbml: Option<String>,
    /// Its diagnostics, as "line code" (line 0 for none), sorted.
    diagnostics: Vec<String>,
}

fn var(name: &str) -> String {
    env::var(name).unwrap_or_else(|_| {
        eprintln!("{name} is not set; see the documentation at the top of this file");
        std::process::exit(2);
    })
}

/// Every case of the index, or those named on the command line.
fn cases(root: &Path, index: &Path, only: &[String]) -> Vec<Case> {
    let text = fs::read_to_string(index).expect("cannot read the index");
    let json: serde_json::Value = serde_json::from_str(&text).expect("index is not JSON");
    json.as_array()
        .expect("index is not a list")
        .iter()
        .map(|item| {
            let id = item["case"].as_str().unwrap_or_default().to_string();
            let packages = item["packages"].as_array();
            Case {
                path: root.join(&id).join(format!("{id}-antimony.ant")),
                comp: packages.is_some_and(|p| p.iter().any(|v| v == "comp")),
                id,
            }
        })
        .filter(|case| only.is_empty() || only.contains(&case.id))
        .collect()
}

/// A diagnostic line of the oracle's, `file(3:17): error: message [dcCode]`,
/// as "line code".
fn oracle_diagnostic(line: &str) -> Option<String> {
    let line = line.trim_end();
    let code = line.strip_suffix(']')?.rsplit_once('[')?.1;
    let at = line
        .find(".ant(")
        .and_then(|i| line[i + 5..].split_once(':'))
        .and_then(|(n, _)| n.parse::<u32>().ok())
        .unwrap_or(0);
    Some(format!("{at} {code}"))
}

/// The oracle's output for a case, from the cache or by running it.
fn run_oracle(oracle: &str, cache: &Path, case: &Case) -> Oracle {
    let sbml_file = cache.join(format!("{}.xml", case.id));
    let errors_file = cache.join(format!("{}.err", case.id));
    let (sbml, errors) = match fs::read_to_string(&errors_file) {
        Ok(errors) => (fs::read_to_string(&sbml_file).ok(), errors),
        Err(_) => {
            let output = Command::new(oracle)
                .arg(&case.path)
                .arg("--quiet")
                .output()
                .expect("cannot run the oracle");
            let errors = String::from_utf8_lossy(&output.stderr).into_owned();
            let sbml = output
                .status
                .success()
                .then(|| String::from_utf8(output.stdout).expect("the oracle wrote invalid UTF-8"));
            match &sbml {
                Some(text) => fs::write(&sbml_file, text).ok(),
                None => fs::remove_file(&sbml_file).ok(),
            };
            fs::write(&errors_file, &errors).ok();
            (sbml, errors)
        }
    };
    let mut diagnostics: Vec<String> = errors.lines().filter_map(oracle_diagnostic).collect();
    diagnostics.sort();
    Oracle { sbml, diagnostics }
}

/// Removes what may differ in form but not in meaning: whitespace in raw
/// notes and annotations, collapsed and dropped between tags. Delphi
/// indents both and ends their lines with CRLF.
fn normalise(model: &mut sbml_io::Model) {
    fn text(raw: &mut Option<String>) {
        if let Some(t) = raw {
            *t = t
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
                .replace("> <", "><");
        }
    }
    fn meta(m: &mut sbml_io::Meta) {
        text(&mut m.notes);
        text(&mut m.annotation);
    }
    meta(&mut model.meta);
    model
        .unit_definitions
        .iter_mut()
        .for_each(|u| meta(&mut u.meta));
    model
        .compartments
        .iter_mut()
        .for_each(|c| meta(&mut c.meta));
    model.species.iter_mut().for_each(|s| meta(&mut s.meta));
    model.parameters.iter_mut().for_each(|p| meta(&mut p.meta));
    model.events.iter_mut().for_each(|e| meta(&mut e.meta));
    model.reactions.iter_mut().for_each(|r| meta(&mut r.meta));
}

/// A model as comparable lines, one per element, keyed by what identifies
/// it: its id, a rule's variable, or its position.
fn elements(model: &sbml_io::Model) -> BTreeMap<String, String> {
    let mut facts = BTreeMap::new();
    let mut add = |key: String, value: String| {
        facts.insert(key, value);
    };
    let header = sbml_io::Model {
        meta: model.meta.clone(),
        id: model.id.clone(),
        name: model.name.clone(),
        conversion_factor: model.conversion_factor.clone(),
        ..sbml_io::Model::default()
    };
    add("model".into(), format!("{header:?}"));
    for f in &model.function_definitions {
        add(format!("function {}", f.id), format!("{f:?}"));
    }
    for u in &model.unit_definitions {
        add(format!("unit {}", u.id), format!("{u:?}"));
    }
    for c in &model.compartments {
        add(format!("compartment {}", c.id), format!("{c:?}"));
    }
    for s in &model.species {
        add(format!("species {}", s.id), format!("{s:?}"));
    }
    for p in &model.parameters {
        add(format!("parameter {}", p.id), format!("{p:?}"));
    }
    for a in &model.initial_assignments {
        add(format!("initial assignment {}", a.symbol), format!("{a:?}"));
    }
    let mut algebraic = 0;
    for r in &model.rules {
        let key = match r.variable() {
            Some(v) => format!("rule {v}"),
            None => {
                algebraic += 1;
                format!("algebraic rule {algebraic}")
            }
        };
        add(key, format!("{r:?}"));
    }
    for (i, c) in model.constraints.iter().enumerate() {
        add(format!("constraint {}", i + 1), format!("{c:?}"));
    }
    for r in &model.reactions {
        add(format!("reaction {}", r.id), format!("{r:?}"));
    }
    for (i, e) in model.events.iter().enumerate() {
        let key = e.id.clone().unwrap_or_else(|| format!("#{}", i + 1));
        add(format!("event {key}"), format!("{e:?}"));
    }
    facts
}

/// Where two descriptions of an element first differ, with some context.
fn first_difference(a: &str, b: &str) -> (String, String) {
    let at = a
        .char_indices()
        .zip(b.chars())
        .find(|((_, x), y)| x != y)
        .map(|((i, _), _)| i)
        .unwrap_or(a.len().min(b.len()));
    let window = |s: &str| {
        let mut start = at.saturating_sub(60);
        while !s.is_char_boundary(start) {
            start -= 1;
        }
        let mut end = (at + 100).min(s.len());
        while !s.is_char_boundary(end) {
            end += 1;
        }
        format!("…{}…", &s[start..end])
    };
    (window(a), window(b))
}

/// A model read from SBML text, normalised.
fn read_model(xml: &str, whose: &str, case: &str) -> sbml_io::Model {
    let mut model = match sbml_io::read(xml) {
        Ok(read) => read.document.model.unwrap_or_default(),
        Err(e) => panic!("{case}: sbml-io cannot read {whose} SBML: {e:?}"),
    };
    normalise(&mut model);
    model
}

fn main() {
    let root = PathBuf::from(var("ANTIMONY_CORPUS"));
    let index = PathBuf::from(var("ANTIMONY_CORPUS_INDEX"));
    let oracle = var("ANTIMONY_ORACLE");
    let only: Vec<String> = env::args().skip(1).collect();
    let cache = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/antimony-oracle");
    fs::create_dir_all(&cache).expect("cannot make the cache folder");

    let cases = cases(&root, &index, &only);
    let (mut agree, mut models, mut facts, mut messages) = (0, 0, 0, 0);
    let mut divergent = Vec::new();
    for case in &cases {
        let oracle = run_oracle(&oracle, &cache, case);
        let source = fs::read_to_string(&case.path).expect("cannot read the case");
        let loaded = antimony::parse(&source);
        let mut differences = Vec::new();

        let mut ours_diagnostics: Vec<String> = loaded
            .diagnostics
            .iter()
            .map(|d| format!("{} {}", d.span.start.line, d.code.name()))
            .collect();
        ours_diagnostics.sort();
        messages += ours_diagnostics.len().max(oracle.diagnostics.len());
        if ours_diagnostics != oracle.diagnostics {
            differences.push((
                "diagnostics".to_string(),
                Some(ours_diagnostics.join(", ")),
                Some(oracle.diagnostics.join(", ")),
            ));
        }

        if let (Some(xml), false) = (&oracle.sbml, case.comp) {
            models += 1;
            let ours_xml = antimony::sbml::write_sbml(&loaded.document);
            let ours = elements(&read_model(&ours_xml, "our", &case.id));
            let theirs = elements(&read_model(xml, "the oracle's", &case.id));
            facts += ours.len().max(theirs.len());
            for key in ours.keys().chain(theirs.keys()) {
                let (a, b) = (ours.get(key), theirs.get(key));
                if a != b && !differences.iter().any(|(k, _, _)| k == key) {
                    let (a, b) = match (a, b) {
                        (Some(a), Some(b)) => {
                            let (a, b) = first_difference(a, b);
                            (Some(a), Some(b))
                        }
                        _ => (a.cloned(), b.cloned()),
                    };
                    differences.push((key.clone(), a, b));
                }
            }
        }

        if differences.is_empty() {
            agree += 1;
        } else {
            divergent.push((case.id.clone(), differences));
        }
    }

    for (id, differences) in &divergent {
        println!("{id}:");
        for (key, ours, theirs) in differences.iter().take(5) {
            println!("  {key}\n    ours:   {ours:?}\n    Delphi: {theirs:?}");
        }
    }
    println!(
        "{} cases: {agree} agree, {} diverge ({messages} diagnostics and, over {models} models, {facts} facts compared)",
        cases.len(),
        divergent.len()
    );
    if !divergent.is_empty() {
        std::process::exit(1);
    }
}

//! Shared by the integration tests: finding the external corpus and its
//! index, and the hand-written invalid fixtures.

#![allow(dead_code)]

pub mod compare;

use std::fs;
use std::path::{Path, PathBuf};

/// The number of `NNNNN-antimony.ant` files the corpus generator wrote.
pub const CORPUS_SIZE: usize = 1788;

/// The corpus folder from `ANTIMONY_CORPUS`, or `None` (with a note on
/// stderr) when it is not set, so the calling test can skip.
pub fn corpus_root(test: &str) -> Option<PathBuf> {
    match std::env::var_os("ANTIMONY_CORPUS") {
        Some(root) => Some(PathBuf::from(root)),
        None => {
            eprintln!("note: ANTIMONY_CORPUS is not set; {test} skipped");
            None
        }
    }
}

/// Every case's Antimony file, as (case id, path), sorted by case id.
pub fn corpus_files(root: &PathBuf) -> Vec<(String, PathBuf)> {
    let mut files = Vec::new();
    for entry in fs::read_dir(root).unwrap_or_else(|e| panic!("{}: {e}", root.display())) {
        let case = entry.unwrap().path();
        let Some(id) = case.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        let file = case.join(format!("{id}-antimony.ant"));
        if file.is_file() {
            files.push((id.to_string(), file));
        }
    }
    files.sort();
    files
}

/// Reads a file as UTF-8.
pub fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

// --- the corpus index ---------------------------------------------------------

/// One case of the corpus index (`index.json`).
#[derive(Debug)]
pub struct Case {
    pub id: String,
    /// The `.ant` file.
    pub path: PathBuf,
    pub components: Vec<String>,
    pub packages: Vec<String>,
    /// Why libantimony lost data writing this case (Delphi §9.1.3).
    pub lossy: Vec<String>,
    /// Non-empty: libantimony cannot read its own text here; we must.
    pub challenge: String,
    /// Non-empty: the case must be rejected (Delphi §5.7).
    pub expect_reject: String,
    /// `conversionFactor` statements the corpus generator added to repair
    /// what libantimony drops (Delphi §5.8).
    pub conversion_factors: u64,
}

/// Delphi §9.1's core component set: the M3/M4 acceptance cases use these
/// components and no package.
const CORE_COMPONENTS: [&str; 9] = [
    "Species",
    "Reaction",
    "Parameter",
    "Compartment",
    "InitialAssignment",
    "AssignmentRule",
    "RateRule",
    "CSymbolTime",
    "StoichiometryMath",
];

impl Case {
    /// The M3/M4 set: core components only, no package (836 cases).
    pub fn is_core(&self) -> bool {
        self.packages.is_empty()
            && self
                .components
                .iter()
                .all(|c| CORE_COMPONENTS.contains(&c.as_str()))
    }

    /// The M5 set: the core plus functions, events, csymbols, algebraic
    /// rules and constraints, still with no package (1629 cases).
    pub fn is_m5(&self) -> bool {
        self.packages.is_empty()
            && self.components.iter().all(|c| {
                CORE_COMPONENTS.contains(&c.as_str())
                    || ["FunctionDefinition", "AlgebraicRule", "Constraint"].contains(&c.as_str())
                    || c.starts_with("Event")
                    || c.starts_with("CSymbol")
            })
    }

    pub fn has_package(&self, package: &str) -> bool {
        self.packages.iter().any(|p| p == package)
    }

    pub fn is_challenge(&self) -> bool {
        !self.challenge.is_empty()
    }

    pub fn expects_reject(&self) -> bool {
        !self.expect_reject.is_empty()
    }

    pub fn is_lossy(&self) -> bool {
        !self.lossy.is_empty()
    }
}

/// The corpus index, from `ANTIMONY_CORPUS_INDEX`, with each case's file
/// under `ANTIMONY_CORPUS`. `None` (with a note) when either is not set.
pub fn corpus_index(test: &str) -> Option<Vec<Case>> {
    let root = corpus_root(test)?;
    let Some(index) = std::env::var_os("ANTIMONY_CORPUS_INDEX") else {
        eprintln!("note: ANTIMONY_CORPUS_INDEX is not set; {test} skipped");
        return None;
    };
    let text = read(Path::new(&index));
    let json: serde_json::Value = serde_json::from_str(&text).expect("index.json is not JSON");
    let strings = |item: &serde_json::Value, key: &str| -> Vec<String> {
        item[key]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|v| v.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default()
    };
    let text_of =
        |item: &serde_json::Value, key: &str| item[key].as_str().unwrap_or_default().to_string();
    let cases = json
        .as_array()
        .expect("index.json is not a list")
        .iter()
        .map(|item| {
            let id = text_of(item, "case");
            Case {
                path: root.join(&id).join(format!("{id}-antimony.ant")),
                id,
                components: strings(item, "components"),
                packages: strings(item, "packages"),
                lossy: strings(item, "lossy"),
                challenge: text_of(item, "challenge"),
                expect_reject: text_of(item, "expect_reject"),
                conversion_factors: item["conversion_factors"].as_u64().unwrap_or(0),
            }
        })
        .collect();
    Some(cases)
}

// --- the invalid fixtures -------------------------------------------------------

/// One expected diagnostic: its Delphi code name (`dcTypeConflict`) and the
/// 1-based line it is reported on.
#[derive(Debug, PartialEq)]
pub struct Expected {
    pub code: String,
    pub line: u32,
}

/// Reads a `.expected` file: `#` lines are comments, every other non-blank
/// line is `<code> <line>`.
pub fn read_expected(path: &Path) -> Result<Vec<Expected>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut expected = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let bad = || {
            format!(
                "{}:{}: expected `<code> <line>`, got {line:?}",
                path.display(),
                i + 1
            )
        };
        let mut parts = line.split_whitespace();
        let (Some(code), Some(n), None) = (parts.next(), parts.next(), parts.next()) else {
            return Err(bad());
        };
        if !code.starts_with("dc") {
            return Err(bad());
        }
        let line = n.parse::<u32>().map_err(|_| bad())?;
        expected.push(Expected {
            code: code.to_string(),
            line,
        });
    }
    Ok(expected)
}

/// The folder of invalid fixtures.
pub fn invalid_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/invalid")
}

/// The files in `dir` with extension `ext`, sorted by name.
pub fn files_with(dir: &Path, ext: &str) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == ext))
        .collect();
    files.sort();
    files
}

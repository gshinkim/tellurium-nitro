//! The test data is in place and readable.
//!
//! `tests/invalid/` holds one Antimony file per diagnostic, each with a
//! `.expected` file naming the codes and lines it must produce. This file
//! only checks the set is well formed and names real codes; `parser.rs`
//! checks the parser against it.
//!
//! The corpus lives outside the repository and is found through
//! `ANTIMONY_CORPUS` (SPEC.md §8). Without it the corpus tests are skipped
//! with a note, so a fresh clone still passes.

mod common;

use std::fs;

use antimony::Code;
use common::{Expected, files_with, invalid_dir, read_expected};

#[test]
fn invalid_fixtures_are_well_formed() {
    let dir = invalid_dir();
    let sources = files_with(&dir, "ant");
    let expectations = files_with(&dir, "expected");
    assert_eq!(sources.len(), 12, "the Delphi set has 12 fixtures");
    assert_eq!(
        sources
            .iter()
            .map(|p| p.with_extension(""))
            .collect::<Vec<_>>(),
        expectations
            .iter()
            .map(|p| p.with_extension(""))
            .collect::<Vec<_>>(),
        "every .ant needs a .expected, and the other way round"
    );
    for source in &sources {
        let lines = fs::read_to_string(source).unwrap().lines().count() as u32;
        let expected = read_expected(&source.with_extension("expected")).unwrap();
        assert!(
            !expected.is_empty(),
            "{}: names no diagnostic",
            source.display()
        );
        // One past the last line is end of file, where an unterminated block
        // is reported (expected-end.ant).
        for e in &expected {
            assert!(
                Code::from_name(&e.code).is_some(),
                "{}: {} is not a diagnostic code",
                source.display(),
                e.code
            );
            assert!(
                (1..=lines + 1).contains(&e.line),
                "{}: line {} is outside the file ({lines} lines)",
                source.display(),
                e.line
            );
        }
    }
}

#[test]
fn expected_file_format() {
    let dir = std::env::temp_dir().join(format!("antimony-fixtures-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let good = dir.join("good.expected");
    fs::write(&good, "# a comment\n\ndcTypeConflict 3\ndcExpectedEnd 4\n").unwrap();
    assert_eq!(
        read_expected(&good).unwrap(),
        vec![
            Expected {
                code: "dcTypeConflict".into(),
                line: 3
            },
            Expected {
                code: "dcExpectedEnd".into(),
                line: 4
            },
        ]
    );
    for bad in [
        "dcTypeConflict",
        "dcTypeConflict three",
        "TypeConflict 3",
        "dcA 1 2",
    ] {
        let path = dir.join("bad.expected");
        fs::write(&path, bad).unwrap();
        assert!(read_expected(&path).is_err(), "{bad:?} should be rejected");
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn corpus_is_reachable() {
    let Some(root) = common::corpus_root("corpus_is_reachable") else {
        return;
    };
    assert_eq!(
        common::corpus_files(&root).len(),
        common::CORPUS_SIZE,
        "Antimony files under {}",
        root.display()
    );
}

#[test]
fn corpus_index_is_readable() {
    let Some(cases) = common::corpus_index("corpus_index_is_readable") else {
        return;
    };
    assert_eq!(cases.len(), common::CORPUS_SIZE);
    // Delphi §9.1's milestone sets.
    assert_eq!(cases.iter().filter(|c| c.is_core()).count(), 836);
    assert_eq!(cases.iter().filter(|c| c.is_m5()).count(), 1629);
    assert_eq!(cases.iter().filter(|c| c.has_package("comp")).count(), 125);
    assert_eq!(cases.iter().filter(|c| c.is_challenge()).count(), 10);
    assert_eq!(cases.iter().filter(|c| c.expects_reject()).count(), 5);
    assert_eq!(cases.iter().filter(|c| c.is_lossy()).count(), 9);
    assert!(
        cases.iter().all(|c| c.path.is_file()),
        "a case without its file"
    );
}

/// Ported from the Delphi `Tests/uAntCorpusTests.pas`: the rest of what
/// `corpus_index_is_readable` does not already check.
#[test]
fn corpus_sets_are_consistent() {
    let Some(cases) = common::corpus_index("corpus_sets_are_consistent") else {
        return;
    };
    // The M5 set is exactly the cases with no SBML package, and contains the
    // core set.
    for case in &cases {
        assert_eq!(case.is_m5(), case.packages.is_empty(), "{}", case.id);
        assert!(!case.is_core() || case.is_m5(), "{}", case.id);
    }
    // The rejected cases left the challenge set; if a case were in both, one
    // acceptance criterion could not be met.
    for case in cases.iter().filter(|c| c.expects_reject()) {
        assert_eq!(case.expect_reject, "negative-stoichiometry", "{}", case.id);
        assert!(!case.is_challenge(), "{}", case.id);
    }
    // The conversionFactor repair of Delphi §5.8: 117 statements in 78 cases.
    // If the corpus is regenerated without it, this notices.
    let repaired: Vec<_> = cases.iter().filter(|c| c.conversion_factors > 0).collect();
    assert_eq!(repaired.len(), 78);
    assert_eq!(
        repaired.iter().map(|c| c.conversion_factors).sum::<u64>(),
        117
    );
}

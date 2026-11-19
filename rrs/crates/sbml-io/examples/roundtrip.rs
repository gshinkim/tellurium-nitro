//! Reads every `.xml` file under a folder, writes each one back into a
//! second folder, and checks that writing is stable: reading what was
//! written and writing it again gives the same text. Also reports what
//! `validate` finds. Exits with 1 if any file can't be read, isn't stable,
//! has warnings or validation issues.
//!
//! ```text
//! cargo run -p sbml-io --release --example roundtrip -- <in> <out>
//! ```

use std::path::{Path, PathBuf};

fn files(dir: &Path, out: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(dir).unwrap().flatten() {
        let p = e.path();
        if p.is_dir() {
            files(&p, out);
        } else if p.extension().is_some_and(|x| x == "xml") {
            out.push(p);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (input, output) = (Path::new(&args[1]), Path::new(&args[2]));
    let mut all = Vec::new();
    files(input, &mut all);
    all.sort();
    let (mut ok, mut unread, mut unstable, mut warned, mut invalid) = (0, 0, 0, 0, 0);
    for p in &all {
        let rel = p.strip_prefix(input).unwrap();
        let xml = std::fs::read_to_string(p).unwrap();
        let read = match sbml_io::read(&xml) {
            Ok(r) => r,
            Err(e) => {
                unread += 1;
                println!("UNREAD {}: {}", rel.display(), e[0]);
                continue;
            }
        };
        let issues = sbml_io::validate(&read.document);
        if !issues.is_empty() {
            invalid += 1;
            println!("INVALID {}: {}", rel.display(), issues[0]);
        }
        if !read.warnings.is_empty() {
            warned += 1;
            println!("WARN {}: {}", rel.display(), read.warnings[0]);
        }
        let first = sbml_io::write(&read.document);
        let second = match sbml_io::read(&first) {
            Ok(r) => sbml_io::write(&r.document),
            Err(e) => format!("unreadable: {}", e[0]),
        };
        let dest = output.join(rel);
        std::fs::create_dir_all(dest.parent().unwrap()).unwrap();
        std::fs::write(&dest, &first).unwrap();
        if first == second {
            ok += 1;
        } else {
            unstable += 1;
            println!("UNSTABLE {}", rel.display());
        }
    }
    println!(
        "{} files: {ok} stable, {unstable} unstable, {unread} not read, {warned} with warnings, {invalid} with validation issues",
        all.len()
    );
    if unread + unstable + warned + invalid > 0 {
        std::process::exit(1);
    }
}

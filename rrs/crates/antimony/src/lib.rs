//! Read and write [Antimony](https://tellurium.readthedocs.io/en/latest/antimony.html),
//! the human-readable language for biochemical models, and convert it to
//! SBML. Pure Rust with no dependencies, so it also works in WebAssembly.
//!
//! It is a port of a complete Antimony implementation written in Delphi,
//! and is checked against it case by case: over a corpus of 1788 models and
//! 708 curated BioModels it reports the same problems on the same lines and
//! writes Antimony that reads back as the same model, and over the corpus
//! it writes the same SBML element by element. `SPEC.md` records the design
//! and every place it departs from the Delphi version.
//!
//! # In five lines
//!
//! ```
//! let loaded = antimony::parse("J1: S1 -> S2; k1*S1\nS1 = 10; S2 = 0; k1 = 0.1");
//! assert!(!loaded.has_errors());
//! let model = loaded.document.main_module();
//! assert_eq!(model.concentration("S1"), Some(10.0));
//! assert_eq!(model.reactions[0].kinetic_law.as_ref().unwrap().to_string(), "k1 * S1");
//! ```
//!
//! # Reading, and what can go wrong
//!
//! [`parse`] never fails and never panics. A malformed file is not an error
//! for the caller to unwrap — files in the wild are malformed in ordinary
//! ways, and an editor must be able to open one and show what is wrong — so
//! it returns a [`Loaded`]: the document, as much of it as could be read,
//! and the [`Diagnostics`] found on the way, each with a [`Code`], a
//! [`Severity`], a message and a [`Span`] (line and column from 1, byte
//! offset from 0).
//!
//! ```
//! let loaded = antimony::parse("J1: A -> B; k1*typoo\nA = 1; B = 0; k1 = 2");
//! // A typo is a warning, not a failure: the model still loads...
//! assert!(!loaded.has_errors());
//! // ...and the warning says where. Switch on the code, not the text.
//! let warning = loaded.diagnostics.first_of(antimony::Code::ImplicitParameter).unwrap();
//! assert_eq!(warning.span.start.line, 1);
//! assert_eq!(warning.to_string(),
//!     "(1:16): warning: \"typoo\" was never declared; it has been assumed to be a parameter. \
//!      Check for a typo [dcImplicitParameter]");
//!
//! // A caller who wants no model unless it is error-free:
//! let document = antimony::parse("J1: S1 => ~~~ ;").into_result();
//! assert!(document.is_err());
//! ```
//!
//! Two diagnostics mean different things: `dcUninitialisedSymbol` is a
//! symbol declared or used but never given a value; `dcImplicitParameter` is
//! a name never declared at all, invented from a use — usually a typo, and
//! the more serious of the two. [`Options`] (with [`parse_with`]) makes both
//! errors (`strict`), or gives each such symbol a value, 0 or 1 for a
//! compartment (`initialise_missing`; the warning is still given).
//!
//! The crate never opens a file: it also runs in a browser, where there is
//! no file system. Read the file yourself and pass its text to [`parse`], or
//! its bytes to [`parse_bytes`], which reports bytes that are not UTF-8 as a
//! diagnostic at the first bad one.
//!
//! # The model
//!
//! A [`Document`] holds the [`Module`](model::Module)s of a
//! file; `main_module()` is the `model *name()` one, or the only one. A
//! module is plain data with public fields — `symbols`, `reactions`,
//! `events`, `functions`, `units`, `constraints`, … — in the order the file
//! gave them, with `Option` wherever the file may have said nothing. Items
//! refer to each other by name (a species' `compartment` is the
//! compartment's id), looked up with `symbol`, `reaction`, `event`,
//! `function` and `unit`. Every [`Symbol`](model::Symbol) has a
//! [`SymbolKind`](model::SymbolKind) that the file implies: `J1: S1 -> S2`
//! makes `S1` and `S2` species before anything declares them.
//!
//! The fields say what the file wrote. The value methods answer what it is
//! worth at the start, worked out through the other symbols, with these
//! defaults where the file says nothing, so a caller writes no test for a
//! missing value:
//!
//! | Method | When the model is silent |
//! |---|---|
//! | `concentration`, `amount` | 0 |
//! | `compartment_size` | 1 — a size of 0 would divide every concentration by zero |
//! | `parameter_value` | 0, and the semantic pass warned about it |
//! | `stoichiometry` | 1 |
//!
//! Each answers `None` for an id that is not of its kind, which is a mistake
//! in the question rather than a fact about the model. `value` is the raw
//! answer, `None` when there is nothing to compute.
//!
//! ```
//! let loaded = antimony::parse("
//!     compartment C = 4;
//!     species A in C, $X in C;
//!     J1: 2 A => X; k*A;
//!     A = 2; X = 0; k = 2*j; j = 3;
//!     sr = 3; J2: sr A -> ; k;
//! ");
//! let m = loaded.document.main_module();
//! assert_eq!(m.amount("A"), Some(8.0));             // a concentration of 2 in a size of 4
//! assert_eq!(m.parameter_value("k"), Some(6.0));    // 2*j, with j = 3
//! assert_eq!(m.parameter_value("A"), None);         // A is a species
//! assert_eq!(m.symbol("k").unwrap().initial_value.as_ref().unwrap().to_string(), "2 * j");
//!
//! // `$X` is a boundary species. The `$` is not part of its id.
//! let floating: Vec<_> = m.floating_species().map(|s| s.name.as_str()).collect();
//! let boundary: Vec<_> = m.boundary_species().map(|s| s.name.as_str()).collect();
//! assert_eq!((floating, boundary), (vec!["A"], vec!["X"]));
//!
//! // A participant's stoichiometry, as a value and as written.
//! let j1 = m.reaction("J1").unwrap();
//! assert_eq!(m.stoichiometry(&j1.reactants[0]), 2.0);
//! assert!(!j1.is_reversible());                    // `=>`; `->` is reversible
//! let j2 = m.reaction("J2").unwrap();
//! assert_eq!(m.stoichiometry(&j2.reactants[0]), 3.0);
//! assert_eq!(j2.reactants[0].stoichiometry_text(), "sr");
//! ```
//!
//! A symbol's value comes from one of three places, which mean different
//! things: `k = 2*j` is its initial value (`initial_value`), `k := 2*j` an
//! assignment rule that holds at all times (`assignment_rule()`), and
//! `k' = 2*j` a rate rule, its derivative (`rate_rule()`).
//!
//! An unlabelled reaction is named `_J1`, `_J2`, … by its position among all
//! the reactions, so every reaction has an id, as a simulator needs; an id
//! the author wrote is never changed, and the made-up names step past it.
//! The Antimony writer leaves a made-up name out.
//!
//! # Writing
//!
//! [`write_document`] writes a document back as Antimony (and
//! [`write_module`] one module; [`WriteOptions`] for the layout). What it
//! writes reads back as the same model. The first write groups the
//! statements by kind, so it reorders a file; after that the text does not
//! change.
//!
//! ```
//! let loaded = antimony::parse("J1: S1 -> S2; k1*S1; S1 = 10; k1 = 0.1");
//! let text = antimony::write_document(&loaded.document);
//! assert!(text.contains("J1: S1 -> S2; k1 * S1;"));
//! assert_eq!(antimony::parse(&text).document.main_module().reactions.len(), 1);
//! ```
//!
//! With the `sbml` feature, `sbml::write_sbml` writes the main module as SBML
//! Level 3 Version 2, using the `sbml-io` crate, and `sbml::to_sbml` gives it
//! as data to inspect or change first:
//!
//! ```
//! # #[cfg(feature = "sbml")]
//! # {
//! let loaded = antimony::parse("J1: S1 -> S2; k1*S1; S1 = 10; S2 = 0; k1 = 0.1");
//! let xml = antimony::sbml::write_sbml(&loaded.document);
//! assert!(xml.contains(r#"<species id="S1" compartment="default_compartment""#));
//! # }
//! ```
//!
//! # Lower levels
//!
//! The layers are public, for an editor that wants them: the lexer's tokens
//! (for highlighting), the expression parser, the parser before the semantic
//! pass, and the semantic pass itself.
//!
//! ```
//! use antimony::lexer::{tokenize, TokenKind};
//!
//! let (tokens, diagnostics) = tokenize("J1: S1 -> S2; k1*S1");
//! assert!(diagnostics.is_empty());
//! assert_eq!(tokens[3].kind, TokenKind::ArrowRev);
//! assert_eq!(tokens[3].span.to_string(), "1:8-10");
//!
//! let (expr, _) = antimony::expr::parse_expression("(k1*S1^2)/(1 + S1)");
//! let expr = expr.unwrap();
//! assert_eq!(expr.to_string(), "k1 * S1 ^ 2 / (1 + S1)");
//! assert_eq!(expr.symbols(), ["k1", "S1"]);
//! assert_eq!(expr.evaluate_with(&|name| if name == "k1" { Some(2.0) } else { Some(3.0) }), Some(4.5));
//! ```
//!
//! # Not built
//!
//! - **Reading SBML.** This crate writes SBML; it never reads it.
//! - **Submodules, `import` and the `comp` package.** A file with several
//!   `model` blocks is read, and each can be written, but submodules are not
//!   expanded, and `import` is reported as not supported.
//! - **The `fbc` package.** Objectives and flux bounds are read and kept, and
//!   written back as Antimony, but not written to SBML.
//! - **DNA strand syntax.** Recognised, so it is reported by name.
//! - **Unit algebra.** A unit definition is kept as text. SBML output
//!   understands base units, `^` exponents and a leading multiplier, and
//!   writes an empty definition for anything more.
//! - **Layout and render statements** (Antimony 3: `model.layout = on`,
//!   `S1.pos = {28, 35}`). Read into [`model::Module::layout`], kept and
//!   written back, but not turned into SBML Layout or Render.
//!
//! # Conventions
//!
//! - `->` is reversible and `=>` is not.
//! - There are no reserved words: `compartment`, `at` or `is` can be names.
//!   Where a word stands decides what it means.
//! - An id is never rewritten, in Antimony or SBML.
//! - An index into `floating_species()` is not an index into `species()`;
//!   use ids across lists.

#![forbid(unsafe_code)]

pub mod builtins;
pub mod diagnostics;
pub mod expr;
pub mod flatten;
pub mod lexer;
pub mod messages;
pub mod model;
pub mod number;
pub mod parser;
#[cfg(feature = "sbml")]
pub mod sbml;
pub mod semantics;
pub mod span;
pub mod write;

pub use diagnostics::{Code, Diagnostic, Diagnostics, Severity};
pub use semantics::Options;
pub use span::{Pos, Span};
pub use write::{WriteOptions, write_document, write_module};

use model::Document;

// The README's examples are run as doc tests (the last needs `sbml`).
#[cfg(all(doctest, feature = "sbml"))]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

/// This crate's version.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A source read into a model, and what was found on the way.
///
/// A malformed file is not an error to unwrap: the document holds as much
/// as could be read, and the problems are in `diagnostics`. Use
/// [`Loaded::into_result`] to treat any error as failure.
#[derive(Clone, Debug)]
pub struct Loaded {
    pub document: Document,
    /// In the order found: the parser's, then the semantic pass's.
    pub diagnostics: Diagnostics,
}

impl Loaded {
    pub fn has_errors(&self) -> bool {
        self.diagnostics.has_errors()
    }

    /// The document if nothing was an error, else the diagnostics.
    pub fn into_result(self) -> Result<Document, Diagnostics> {
        if self.has_errors() {
            Err(self.diagnostics)
        } else {
            Ok(self.document)
        }
    }
}

/// Reads Antimony: parses it, then resolves and checks it. Never fails and
/// never panics; see [`Loaded`].
pub fn parse(source: &str) -> Loaded {
    parse_with(source, &Options::default())
}

/// As [`parse`], for a file's bytes, which are read as UTF-8 (a byte order
/// mark is skipped). Bytes that are not UTF-8 are not parsed: the document
/// is empty, and a `dcReadError` points at the first bad byte. The crate
/// never opens files itself (it also runs in a browser), so reading them,
/// and reporting a missing one, is the caller's part.
pub fn parse_bytes(bytes: &[u8], options: &Options) -> Loaded {
    match std::str::from_utf8(bytes) {
        Ok(source) => parse_with(source, options),
        Err(e) => {
            let good = e.valid_up_to();
            // The prefix is valid UTF-8 by definition.
            let prefix = std::str::from_utf8(&bytes[..good]).unwrap_or_default();
            let at = span::LineIndex::new(prefix).pos(prefix, good);
            let mut diagnostics = Diagnostics::default();
            diagnostics.add(Code::ReadError, Span::from_pos(at, 1), messages::NOT_UTF8);
            Loaded {
                document: Document::default(),
                diagnostics,
            }
        }
    }
}

/// As [`parse`], with [`Options`]: strict mode, or filling in missing values.
pub fn parse_with(source: &str, options: &Options) -> Loaded {
    let (mut document, mut diagnostics) = parser::parse_document(source);
    let mut found = semantics::analyse(&mut document, options);
    diagnostics.append(&mut found);
    Loaded {
        document,
        diagnostics,
    }
}

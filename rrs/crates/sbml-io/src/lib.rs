//! Read and write [SBML](https://sbml.org) as plain Rust data: Level 3
//! (Versions 1 and 2) and Level 2 (Versions 1 to 5), all of core, and the
//! Layout and Render packages. Pure Rust, so it also works in WebAssembly.
//!
//! A document is a tree of ordinary structs with public fields: `Vec` for
//! lists, `Option` for attributes that may be absent, enums where SBML has
//! alternatives. There are three entry points: [`read`], [`write`] and
//! [`validate`].
//!
//! ```
//! let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
//! <sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
//!   <model id="m">
//!     <listOfCompartments>
//!       <compartment id="cell" size="1" constant="true"/>
//!     </listOfCompartments>
//!   </model>
//! </sbml>"#;
//! let mut doc = sbml_io::read(xml).unwrap().document;
//! let model = doc.model.as_mut().unwrap();
//! model.species.push(sbml_io::Species {
//!     id: "S1".into(),
//!     compartment: "cell".into(),
//!     initial_concentration: Some(10.0),
//!     ..Default::default()
//! });
//! assert!(sbml_io::validate(&doc).is_empty());
//! assert!(sbml_io::write(&doc).contains(r#"<species id="S1""#));
//! ```
//!
//! **Attributes are stored as written.** `None` means the attribute was
//! absent, not that it has its default value; applications that need the
//! effective value apply the defaults of the document's level.
//!
//! **Nothing is lost.** Notes and annotations are kept as raw XML, and so is
//! anything `sbml-io` doesn't understand: other packages' elements and
//! attributes, unknown MathML. All of it is written back. Kept child elements
//! are written after the children `sbml-io` knows.

pub mod layout;
pub mod math;
pub mod miriam;
pub mod model;
pub mod ns;
mod read;
pub mod render;
mod validate;
mod write;

pub use layout::*;
pub use math::{Constant, Math, Number, Op, Symbol};
pub use model::*;
pub use read::{Diagnostic, Read, read};
pub use validate::{Issue, validate};
pub use write::{fmt_f64, write};

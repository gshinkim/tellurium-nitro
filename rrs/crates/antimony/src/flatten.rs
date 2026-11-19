//! Flattening (Delphi `uAntFlatten`, Delphi §2.2.2 seam 2).
//!
//! For now this is a copy of the module. That is the point: the SBML writer
//! is written against what flattening returns, never against the parse
//! tree, so expanding submodules later would change this function and no
//! consumer of it. The original is untouched, so an editor can keep both.

use crate::model::{Document, Module};

/// The main module, flattened.
pub fn flatten(document: &Document) -> Module {
    document.main_module().clone()
}

/// The named module, flattened, or `None` if there is no such module.
pub fn flatten_module(document: &Document, name: &str) -> Option<Module> {
    document.module(name).cloned()
}

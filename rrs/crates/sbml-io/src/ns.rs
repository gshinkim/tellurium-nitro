//! Namespace URIs.

pub const MATHML: &str = "http://www.w3.org/1998/Math/MathML";
pub const XHTML: &str = "http://www.w3.org/1999/xhtml";
pub const XSI: &str = "http://www.w3.org/2001/XMLSchema-instance";
/// The Layout package, Level 3 Version 1.
pub const LAYOUT: &str = "http://www.sbml.org/sbml/level3/version1/layout/version1";
/// Layout in a Level 2 annotation.
pub const LAYOUT_L2: &str = "http://projects.eml.org/bcb/sbml/level2";
/// The Render package, Level 3 Version 1.
pub const RENDER: &str = "http://www.sbml.org/sbml/level3/version1/render/version1";
/// Render in Level 2 annotations.
pub const RENDER_L2: &str = "http://projects.eml.org/bcb/sbml/render/level2";

/// The core namespace of a level and version, if `sbml-io` supports it.
pub fn core(level: u32, version: u32) -> Option<&'static str> {
    Some(match (level, version) {
        (2, 1) => "http://www.sbml.org/sbml/level2",
        (2, 2) => "http://www.sbml.org/sbml/level2/version2",
        (2, 3) => "http://www.sbml.org/sbml/level2/version3",
        (2, 4) => "http://www.sbml.org/sbml/level2/version4",
        (2, 5) => "http://www.sbml.org/sbml/level2/version5",
        (3, 1) => "http://www.sbml.org/sbml/level3/version1/core",
        (3, 2) => "http://www.sbml.org/sbml/level3/version2/core",
        _ => return None,
    })
}

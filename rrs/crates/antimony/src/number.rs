//! Reading and writing numbers (Delphi `AntTryStrToFloat`, `AntFloatToStr`).
//!
//! Rust's own conversions never depend on the machine's locale, so the
//! Delphi version's comma-decimal hazard cannot arise here (SPEC.md §9, R12).
//! Everything in this crate that turns text into a number or back goes
//! through these two functions.

/// Reads a number: surrounding whitespace, an optional sign, digits with an
/// optional `.` and exponent, or `inf`, `infinity` and `nan` in any case.
/// A comma is never a decimal separator.
pub fn parse(text: &str) -> Option<f64> {
    let text = text.trim();
    if text.is_empty() {
        return None;
    }
    text.parse::<f64>().ok()
}

/// The shortest text that [`parse`] reads back as exactly `value`.
///
/// Uses Antimony's spellings `inf`, `-inf` and `nan` for the non-finite
/// values, and exponent form (`6.02214076e23`) outside 1e-5 to 1e15, the
/// same thresholds `sbml-io` uses.
pub fn format(value: f64) -> String {
    if value.is_nan() {
        "nan".into()
    } else if value.is_infinite() {
        if value > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        }
    } else if value != 0.0 && (value.abs() >= 1e15 || value.abs() < 1e-5) {
        format!("{value:e}")
    } else {
        format!("{value}")
    }
}

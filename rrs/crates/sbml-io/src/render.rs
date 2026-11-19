//! The Render package (Level 3 Render Version 1), which in Level 2 files
//! lives in annotations: the global render information in the annotation of
//! `listOfLayouts`, a layout's own in the layout's annotation.
//!
//! The model is the file's, without interpretation: which style applies to a
//! glyph, and what its look is, is for the application to work out. Values
//! that may be relative (`50%`) or both (`5 + 50%`) are [`RelAbs`], kept as
//! written; colours are kept as written (`#RRGGBB`, `#RRGGBBAA`, a colour or
//! gradient id, or `none`).

use crate::layout::BoundingBox;
use crate::model::Meta;

/// `listOfGlobalRenderInformation` or `listOfRenderInformation`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RenderList {
    pub items: Vec<RenderInformation>,
    /// The list's own attributes (`versionMajor`, …) and children, such as
    /// `defaultValues`, kept as written.
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct RenderInformation {
    pub id: String,
    pub name: Option<String>,
    pub program_name: Option<String>,
    pub program_version: Option<String>,
    pub reference_render_information: Option<String>,
    pub background_color: Option<String>,
    pub color_definitions: Vec<ColorDefinition>,
    pub gradients: Vec<Gradient>,
    pub line_endings: Vec<LineEnding>,
    pub styles: Vec<Style>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ColorDefinition {
    pub id: String,
    pub name: Option<String>,
    pub value: String,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Gradient {
    pub id: String,
    pub name: Option<String>,
    pub spread_method: Option<String>,
    pub kind: GradientKind,
    pub stops: Vec<GradientStop>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq)]
pub enum GradientKind {
    Linear {
        x1: Option<RelAbs>,
        y1: Option<RelAbs>,
        z1: Option<RelAbs>,
        x2: Option<RelAbs>,
        y2: Option<RelAbs>,
        z2: Option<RelAbs>,
    },
    Radial {
        cx: Option<RelAbs>,
        cy: Option<RelAbs>,
        cz: Option<RelAbs>,
        r: Option<RelAbs>,
        fx: Option<RelAbs>,
        fy: Option<RelAbs>,
        fz: Option<RelAbs>,
    },
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct GradientStop {
    pub offset: Option<RelAbs>,
    pub stop_color: Option<String>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct LineEnding {
    pub id: String,
    pub enable_rotational_mapping: Option<bool>,
    pub bounding_box: Option<BoundingBox>,
    pub group: Option<Group>,
    pub meta: Meta,
}

/// A style: the look given to the glyphs it names by id, role or type.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Style {
    pub id: Option<String>,
    pub name: Option<String>,
    pub role_list: Vec<String>,
    pub type_list: Vec<String>,
    /// Local styles only.
    pub id_list: Vec<String>,
    pub group: Option<Group>,
    pub meta: Meta,
}

/// The attributes that set how something is drawn. Every element that draws
/// may carry them; a group passes them on to what it holds. All are as
/// written.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Presentation {
    pub stroke: Option<String>,
    pub stroke_width: Option<f64>,
    /// As written, for example `5,4`.
    pub stroke_dasharray: Option<String>,
    pub fill: Option<String>,
    pub fill_rule: Option<String>,
    pub font_family: Option<String>,
    pub font_size: Option<RelAbs>,
    pub font_weight: Option<String>,
    pub font_style: Option<String>,
    pub text_anchor: Option<String>,
    pub vtext_anchor: Option<String>,
    pub start_head: Option<String>,
    pub end_head: Option<String>,
    /// As written: six or twelve numbers.
    pub transform: Option<String>,
}

/// `g`: a group of primitives and further groups.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Group {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub children: Vec<Primitive>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Primitive {
    Group(Group),
    Rectangle(Rectangle),
    Ellipse(Ellipse),
    Polygon(Polygon),
    Curve(RenderCurve),
    Text(Text),
    Image(Image),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Rectangle {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub x: Option<RelAbs>,
    pub y: Option<RelAbs>,
    pub z: Option<RelAbs>,
    pub width: Option<RelAbs>,
    pub height: Option<RelAbs>,
    pub rx: Option<RelAbs>,
    pub ry: Option<RelAbs>,
    pub ratio: Option<f64>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Ellipse {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub cx: Option<RelAbs>,
    pub cy: Option<RelAbs>,
    pub cz: Option<RelAbs>,
    pub rx: Option<RelAbs>,
    pub ry: Option<RelAbs>,
    pub ratio: Option<f64>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polygon {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub elements: Vec<RenderPoint>,
    pub meta: Meta,
}

/// `curve` in a style: an open line through points, with line endings.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RenderCurve {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub elements: Vec<RenderPoint>,
    pub meta: Meta,
}

/// A point of a polygon or curve; with `bezier`, the end of a cubic Bézier
/// from the point before it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RenderPoint {
    pub x: Option<RelAbs>,
    pub y: Option<RelAbs>,
    pub z: Option<RelAbs>,
    pub bezier: Option<BezierBase>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BezierBase {
    pub x1: Option<RelAbs>,
    pub y1: Option<RelAbs>,
    pub z1: Option<RelAbs>,
    pub x2: Option<RelAbs>,
    pub y2: Option<RelAbs>,
    pub z2: Option<RelAbs>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Text {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub x: Option<RelAbs>,
    pub y: Option<RelAbs>,
    pub z: Option<RelAbs>,
    pub text: String,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Image {
    pub id: Option<String>,
    pub name: Option<String>,
    pub presentation: Presentation,
    pub x: Option<RelAbs>,
    pub y: Option<RelAbs>,
    pub z: Option<RelAbs>,
    pub width: Option<RelAbs>,
    pub height: Option<RelAbs>,
    pub href: Option<String>,
    pub meta: Meta,
}

/// A coordinate or length that may be absolute, relative to the box it is
/// drawn in, or both: `10`, `50%`, `5 + 50%`. Kept as written.
#[derive(Debug, Clone, PartialEq)]
pub struct RelAbs(pub String);

impl RelAbs {
    /// An absolute value.
    pub fn abs(x: f64) -> Self {
        RelAbs(crate::fmt_f64(x))
    }

    /// A percentage of the box.
    pub fn rel(percent: f64) -> Self {
        RelAbs(format!("{}%", crate::fmt_f64(percent)))
    }

    /// The absolute part and the percentage, or `None` if it doesn't read.
    pub fn parse(&self) -> Option<(f64, f64)> {
        let s: String = self.0.chars().filter(|c| !c.is_whitespace()).collect();
        let num = |t: &str| crate::read::parse_double(t);
        let Some(body) = s.strip_suffix('%') else {
            return Some((num(&s)?, 0.0));
        };
        // The last sign that starts a term, not one in an exponent.
        let b = body.as_bytes();
        let split = (1..b.len()).rev().find(|&i| {
            matches!(b[i], b'+' | b'-') && !matches!(b[i - 1], b'e' | b'E' | b'+' | b'-')
        });
        match split {
            None => Some((0.0, num(body)?)),
            Some(i) => {
                let rel = if b[i] == b'+' {
                    &body[i + 1..]
                } else {
                    &body[i..]
                };
                Some((num(&body[..i])?, num(rel)?))
            }
        }
    }

    /// Its value in a box of length `size`.
    pub fn resolve(&self, size: f64) -> Option<f64> {
        self.parse().map(|(a, r)| a + r / 100.0 * size)
    }
}

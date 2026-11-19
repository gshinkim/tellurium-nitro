//! The Layout package (Level 3 Layout Version 1), which in Level 2 files
//! lives in the model's annotation.

use crate::model::Meta;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Layout {
    pub id: String,
    pub name: Option<String>,
    pub dimensions: Dimensions,
    pub compartment_glyphs: Vec<CompartmentGlyph>,
    pub species_glyphs: Vec<SpeciesGlyph>,
    pub reaction_glyphs: Vec<ReactionGlyph>,
    pub text_glyphs: Vec<TextGlyph>,
    /// `listOfAdditionalGraphicalObjects`.
    pub additional: Vec<AdditionalGlyph>,
    /// The layout's own render information (the Render package): in Level 3
    /// `listOfRenderInformation`, in Level 2 in the layout's annotation.
    pub render: Option<crate::render::RenderList>,
    pub meta: Meta,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: Option<f64>,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Point { x, y, z: None }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Dimensions {
    pub width: f64,
    pub height: f64,
    pub depth: Option<f64>,
}

impl Dimensions {
    pub fn new(width: f64, height: f64) -> Self {
        Dimensions {
            width,
            height,
            depth: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BoundingBox {
    pub id: Option<String>,
    pub position: Point,
    pub dimensions: Dimensions,
}

impl BoundingBox {
    pub fn new(x: f64, y: f64, width: f64, height: f64) -> Self {
        BoundingBox {
            id: None,
            position: Point::new(x, y),
            dimensions: Dimensions::new(width, height),
        }
    }
}

/// A curve: a chain of line segments and cubic Béziers.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Curve {
    pub segments: Vec<Segment>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Segment {
    Line {
        start: Point,
        end: Point,
    },
    Bezier {
        start: Point,
        end: Point,
        base1: Point,
        base2: Point,
    },
}

impl Segment {
    pub fn start(&self) -> Point {
        match self {
            Segment::Line { start, .. } | Segment::Bezier { start, .. } => *start,
        }
    }

    pub fn end(&self) -> Point {
        match self {
            Segment::Line { end, .. } | Segment::Bezier { end, .. } => *end,
        }
    }
}

/// A plain `graphicalObject`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GraphicalObject {
    pub id: String,
    pub bounding_box: Option<BoundingBox>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct CompartmentGlyph {
    pub id: String,
    pub compartment: Option<String>,
    pub order: Option<f64>,
    pub bounding_box: Option<BoundingBox>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpeciesGlyph {
    pub id: String,
    pub species: Option<String>,
    pub bounding_box: Option<BoundingBox>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReactionGlyph {
    pub id: String,
    pub reaction: Option<String>,
    pub bounding_box: Option<BoundingBox>,
    pub curve: Option<Curve>,
    pub species_reference_glyphs: Vec<SpeciesReferenceGlyph>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpeciesReferenceGlyph {
    pub id: String,
    pub species_glyph: String,
    pub species_reference: Option<String>,
    pub role: Option<Role>,
    pub bounding_box: Option<BoundingBox>,
    pub curve: Option<Curve>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct TextGlyph {
    pub id: String,
    pub text: Option<String>,
    pub graphical_object: Option<String>,
    pub origin_of_text: Option<String>,
    pub bounding_box: Option<BoundingBox>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct GeneralGlyph {
    pub id: String,
    pub reference: Option<String>,
    pub bounding_box: Option<BoundingBox>,
    pub curve: Option<Curve>,
    pub reference_glyphs: Vec<ReferenceGlyph>,
    pub sub_glyphs: Vec<AdditionalGlyph>,
    pub meta: Meta,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReferenceGlyph {
    pub id: String,
    pub glyph: String,
    pub reference: Option<String>,
    pub role: Option<String>,
    pub bounding_box: Option<BoundingBox>,
    pub curve: Option<Curve>,
    pub meta: Meta,
}

/// What `listOfAdditionalGraphicalObjects` and `listOfSubGlyphs` may hold.
#[derive(Debug, Clone, PartialEq)]
pub enum AdditionalGlyph {
    Graphical(GraphicalObject),
    General(GeneralGlyph),
    Compartment(CompartmentGlyph),
    Species(SpeciesGlyph),
    Reaction(ReactionGlyph),
    Text(TextGlyph),
}

impl Layout {
    /// The size that holds everything drawn, from the origin: the furthest
    /// right and lowest of the glyphs' boxes and curves. For a layout read
    /// without its `<dimensions>`.
    pub fn extent(&self) -> Dimensions {
        let mut far = (0.0f64, 0.0f64);
        let add_box = |b: &Option<BoundingBox>, far: &mut (f64, f64)| {
            if let Some(b) = b {
                far.0 = far.0.max(b.position.x + b.dimensions.width);
                far.1 = far.1.max(b.position.y + b.dimensions.height);
            }
        };
        let add_curve = |c: &Option<Curve>, far: &mut (f64, f64)| {
            for s in c.iter().flat_map(|c| &c.segments) {
                for p in [s.start(), s.end()] {
                    far.0 = far.0.max(p.x);
                    far.1 = far.1.max(p.y);
                }
            }
        };
        for g in &self.compartment_glyphs {
            add_box(&g.bounding_box, &mut far);
        }
        for g in &self.species_glyphs {
            add_box(&g.bounding_box, &mut far);
        }
        for g in &self.text_glyphs {
            add_box(&g.bounding_box, &mut far);
        }
        for g in &self.reaction_glyphs {
            add_box(&g.bounding_box, &mut far);
            add_curve(&g.curve, &mut far);
            for r in &g.species_reference_glyphs {
                add_box(&r.bounding_box, &mut far);
                add_curve(&r.curve, &mut far);
            }
        }
        for g in &self.additional {
            let b = match g {
                AdditionalGlyph::Graphical(g) => &g.bounding_box,
                AdditionalGlyph::General(g) => &g.bounding_box,
                AdditionalGlyph::Compartment(g) => &g.bounding_box,
                AdditionalGlyph::Species(g) => &g.bounding_box,
                AdditionalGlyph::Reaction(g) => &g.bounding_box,
                AdditionalGlyph::Text(g) => &g.bounding_box,
            };
            add_box(b, &mut far);
        }
        Dimensions::new(far.0, far.1)
    }
}

impl AdditionalGlyph {
    pub fn id(&self) -> &str {
        match self {
            AdditionalGlyph::Graphical(g) => &g.id,
            AdditionalGlyph::General(g) => &g.id,
            AdditionalGlyph::Compartment(g) => &g.id,
            AdditionalGlyph::Species(g) => &g.id,
            AdditionalGlyph::Reaction(g) => &g.id,
            AdditionalGlyph::Text(g) => &g.id,
        }
    }
}

/// A species reference glyph's role.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Role {
    Substrate,
    Product,
    SideSubstrate,
    SideProduct,
    Modifier,
    Activator,
    Inhibitor,
    Undefined,
    /// A role the package doesn't define, kept as written.
    Other(String),
}

impl Role {
    const NAMED: [(Role, &str); 8] = [
        (Role::Substrate, "substrate"),
        (Role::Product, "product"),
        (Role::SideSubstrate, "sidesubstrate"),
        (Role::SideProduct, "sideproduct"),
        (Role::Modifier, "modifier"),
        (Role::Activator, "activator"),
        (Role::Inhibitor, "inhibitor"),
        (Role::Undefined, "undefined"),
    ];

    pub fn as_str(&self) -> &str {
        match self {
            Role::Other(s) => s,
            r => Self::NAMED.iter().find(|(n, _)| n == r).unwrap().1,
        }
    }

    pub fn parse(s: &str) -> Role {
        Self::NAMED
            .iter()
            .find(|(_, n)| *n == s)
            .map(|(r, _)| r.clone())
            .unwrap_or_else(|| Role::Other(s.to_string()))
    }
}

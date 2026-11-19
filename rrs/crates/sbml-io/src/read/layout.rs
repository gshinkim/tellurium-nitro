//! Reading the Layout package. The same code reads Level 3 layouts and the
//! Level 2 ones kept in annotations; `pkg` is the namespace in each case.

use roxmltree::Node;

use super::{E, Reader};
use crate::layout::*;
use crate::model::Meta;
use crate::ns;
use crate::render::RenderList;

type P<'p> = Option<&'p str>;

impl<'a, 'i: 'a> Reader<'a, 'i> {
    /// `listOfLayouts`: its layouts, its own content, and the global render
    /// information (in Level 2, out of its annotation).
    pub(super) fn layout_list(
        &mut self,
        n: Node<'a, 'i>,
        lns: &str,
    ) -> (Vec<Layout>, Meta, Option<RenderList>) {
        let l2 = lns == ns::LAYOUT_L2;
        let mut layouts = Vec::new();
        let mut global = None;
        let mut rest = None;
        let mut meta = self.meta(E::new(n), Some(lns), |r, c| {
            if r.name_in(c, lns) == Some("layout") {
                let mut l = r.layout(c, Some(lns));
                if l.id.is_empty() {
                    // Reactome leaves it out; SBML Layout requires it.
                    l.id = r.unused_id("layout");
                    let message = format!(
                        "a layout has no id, which SBML Layout requires; it is called '{}' here",
                        l.id
                    );
                    r.warn(c, &message);
                }
                layouts.push(l);
                return true;
            }
            if !l2 && r.name_in(c, ns::RENDER) == Some("listOfGlobalRenderInformation") {
                global = Some(r.render_list(c, ns::RENDER));
                return true;
            }
            if l2 && rest.is_none() && r.is_annotation(c, lns) {
                rest = Some(r.annotation_without(c, "listOfGlobalRenderInformation", &mut global));
                return true;
            }
            false
        });
        if let Some(rest) = rest {
            meta.annotation = rest;
        }
        (layouts, meta, global)
    }

    /// An `<annotation>`: SBML's, or written in the layout namespace.
    fn is_annotation(&self, c: Node<'a, 'i>, lns: &str) -> bool {
        (self.core_name(c) == Some("annotation")) || self.name_in(c, lns) == Some("annotation")
    }

    /// A Level 2 annotation with a render list taken out of it into `slot`;
    /// the rest is kept. `None` if nothing is left.
    fn annotation_without(
        &mut self,
        n: Node<'a, 'i>,
        list: &str,
        slot: &mut Option<RenderList>,
    ) -> Option<String> {
        let is_list = |c: Option<Node>| {
            c.is_some_and(|c| {
                c.is_element()
                    && c.tag_name().namespace() == Some(ns::RENDER_L2)
                    && c.tag_name().name() == list
            })
        };
        let mut rest = String::new();
        let mut took = false;
        for c in n.children() {
            // The whitespace before the list goes with it.
            if c.is_text()
                && c.text().is_some_and(|t| t.trim().is_empty())
                && is_list(c.next_sibling())
            {
                continue;
            }
            if is_list(Some(c)) && slot.is_none() {
                *slot = Some(self.render_list(c, ns::RENDER_L2));
                took = true;
            } else {
                rest.push_str(&self.raw(c));
            }
        }
        if took && rest.trim().is_empty() {
            None
        } else {
            Some(rest)
        }
    }

    pub(super) fn layout(&mut self, n: Node<'a, 'i>, pkg: P) -> Layout {
        let lns = pkg.unwrap_or(ns::LAYOUT);
        let l2 = lns == ns::LAYOUT_L2;
        let mut rest = None;
        let mut render = None;
        let mut e = E::new(n);
        let mut l = Layout {
            // Checked by `layout_list`, which can name one that is missing.
            id: e
                .take_pkg(pkg, "id")
                .map(str::to_string)
                .unwrap_or_default(),
            name: e.take_pkg(pkg, "name").map(str::to_string),
            ..Layout::default()
        };
        let mut dims = None;
        l.meta = self.meta(e, pkg, |r, c| {
            match r.name_in(c, lns) {
                Some("dimensions") => dims = Some(r.dimensions(c, pkg)),
                Some("listOfCompartmentGlyphs") => {
                    l.compartment_glyphs = r.list(c, lns, &["compartmentGlyph"], |r, c| {
                        r.compartment_glyph(c, pkg)
                    })
                }
                Some("listOfSpeciesGlyphs") => {
                    l.species_glyphs =
                        r.list(c, lns, &["speciesGlyph"], |r, c| r.species_glyph(c, pkg))
                }
                Some("listOfReactionGlyphs") => {
                    l.reaction_glyphs =
                        r.list(c, lns, &["reactionGlyph"], |r, c| r.reaction_glyph(c, pkg))
                }
                Some("listOfTextGlyphs") => {
                    l.text_glyphs = r.list(c, lns, &["textGlyph"], |r, c| r.text_glyph(c, pkg))
                }
                Some("listOfAdditionalGraphicalObjects") => {
                    l.additional = r.list(c, lns, ADDITIONAL, |r, c| r.additional(c, pkg))
                }
                _ if !l2 && r.name_in(c, ns::RENDER) == Some("listOfRenderInformation") => {
                    render = Some(r.render_list(c, ns::RENDER))
                }
                _ if l2 && rest.is_none() && r.is_annotation(c, lns) => {
                    rest = Some(r.annotation_without(c, "listOfRenderInformation", &mut render))
                }
                _ => return false,
            }
            true
        });
        if let Some(rest) = rest {
            l.meta.annotation = rest;
        }
        l.render = render;
        match dims {
            Some(d) => l.dimensions = d,
            None => {
                // As Reactome and MINERVA write: the size is what holds
                // everything drawn.
                l.dimensions = l.extent();
                let named = if l.id.is_empty() {
                    String::new()
                } else {
                    format!(" '{}'", l.id)
                };
                let message = format!(
                    "the layout{named} has no <dimensions> (its width and height), which SBML Layout requires; its size is taken from what it draws ({} by {})",
                    crate::fmt_f64(l.dimensions.width),
                    crate::fmt_f64(l.dimensions.height)
                );
                self.warn(n, &message)
            }
        }
        l
    }

    fn additional(&mut self, n: Node<'a, 'i>, pkg: P) -> AdditionalGlyph {
        match n.tag_name().name() {
            "generalGlyph" => AdditionalGlyph::General(self.general_glyph(n, pkg)),
            "compartmentGlyph" => AdditionalGlyph::Compartment(self.compartment_glyph(n, pkg)),
            "speciesGlyph" => AdditionalGlyph::Species(self.species_glyph(n, pkg)),
            "reactionGlyph" => AdditionalGlyph::Reaction(self.reaction_glyph(n, pkg)),
            "textGlyph" => AdditionalGlyph::Text(self.text_glyph(n, pkg)),
            _ => {
                let mut e = E::new(n);
                let id = self.req_pkg(&mut e, pkg, "id");
                let mut bounding_box = None;
                let meta = self.meta(e, pkg, |r, c| {
                    r.bounding_box_child(c, pkg, &mut bounding_box)
                });
                AdditionalGlyph::Graphical(GraphicalObject {
                    id,
                    bounding_box,
                    meta,
                })
            }
        }
    }

    fn compartment_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> CompartmentGlyph {
        let mut e = E::new(n);
        let mut g = CompartmentGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            compartment: e.take_pkg(pkg, "compartment").map(str::to_string),
            order: self.double_pkg(&mut e, pkg, "order"),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            r.bounding_box_child(c, pkg, &mut g.bounding_box)
        });
        g
    }

    fn species_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> SpeciesGlyph {
        let mut e = E::new(n);
        let mut g = SpeciesGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            species: e.take_pkg(pkg, "species").map(str::to_string),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            r.bounding_box_child(c, pkg, &mut g.bounding_box)
        });
        g
    }

    fn reaction_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> ReactionGlyph {
        let lns = pkg.unwrap_or(ns::LAYOUT);
        let mut e = E::new(n);
        let mut g = ReactionGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            reaction: e.take_pkg(pkg, "reaction").map(str::to_string),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            if r.name_in(c, lns) == Some("listOfSpeciesReferenceGlyphs") {
                g.species_reference_glyphs = r.list(c, lns, &["speciesReferenceGlyph"], |r, c| {
                    r.species_reference_glyph(c, pkg)
                });
                return true;
            }
            r.bounding_box_child(c, pkg, &mut g.bounding_box) || r.curve_child(c, pkg, &mut g.curve)
        });
        g
    }

    fn species_reference_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> SpeciesReferenceGlyph {
        let mut e = E::new(n);
        let mut g = SpeciesReferenceGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            species_glyph: self.req_pkg(&mut e, pkg, "speciesGlyph"),
            species_reference: e.take_pkg(pkg, "speciesReference").map(str::to_string),
            role: e.take_pkg(pkg, "role").map(Role::parse),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            r.bounding_box_child(c, pkg, &mut g.bounding_box) || r.curve_child(c, pkg, &mut g.curve)
        });
        g
    }

    fn text_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> TextGlyph {
        let mut e = E::new(n);
        let mut g = TextGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            text: e.take_pkg(pkg, "text").map(str::to_string),
            graphical_object: e.take_pkg(pkg, "graphicalObject").map(str::to_string),
            origin_of_text: e.take_pkg(pkg, "originOfText").map(str::to_string),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            r.bounding_box_child(c, pkg, &mut g.bounding_box)
        });
        g
    }

    fn general_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> GeneralGlyph {
        let lns = pkg.unwrap_or(ns::LAYOUT);
        let mut e = E::new(n);
        let mut g = GeneralGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            reference: e.take_pkg(pkg, "reference").map(str::to_string),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            match r.name_in(c, lns) {
                Some("listOfReferenceGlyphs") => {
                    g.reference_glyphs = r.list(c, lns, &["referenceGlyph"], |r, c| {
                        r.reference_glyph(c, pkg)
                    })
                }
                Some("listOfSubGlyphs") => {
                    g.sub_glyphs = r.list(c, lns, ADDITIONAL, |r, c| r.additional(c, pkg))
                }
                _ => {
                    return r.bounding_box_child(c, pkg, &mut g.bounding_box)
                        || r.curve_child(c, pkg, &mut g.curve);
                }
            }
            true
        });
        g
    }

    fn reference_glyph(&mut self, n: Node<'a, 'i>, pkg: P) -> ReferenceGlyph {
        let mut e = E::new(n);
        let mut g = ReferenceGlyph {
            id: self.req_pkg(&mut e, pkg, "id"),
            glyph: self.req_pkg(&mut e, pkg, "glyph"),
            reference: e.take_pkg(pkg, "reference").map(str::to_string),
            role: e.take_pkg(pkg, "role").map(str::to_string),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| {
            r.bounding_box_child(c, pkg, &mut g.bounding_box) || r.curve_child(c, pkg, &mut g.curve)
        });
        g
    }

    pub(super) fn bounding_box_child(
        &mut self,
        c: Node<'a, 'i>,
        pkg: P,
        slot: &mut Option<BoundingBox>,
    ) -> bool {
        if self.name_in(c, pkg.unwrap_or(ns::LAYOUT)) != Some("boundingBox") {
            return false;
        }
        let lns = pkg.unwrap_or(ns::LAYOUT);
        let mut e = E::new(c);
        let mut bb = BoundingBox {
            id: e.take_pkg(pkg, "id").map(str::to_string),
            ..Default::default()
        };
        let (mut pos, mut dims) = (None, None);
        for k in c.children().filter(|k| k.is_element()) {
            match self.name_in(k, lns) {
                Some("position") => pos = Some(self.point(k, pkg)),
                Some("dimensions") => dims = Some(self.dimensions(k, pkg)),
                _ => self.warn(k, "dropped from <boundingBox>"),
            }
        }
        match (pos, dims) {
            (Some(p), Some(d)) => (bb.position, bb.dimensions) = (p, d),
            _ => {
                // The glyph it belongs to, by its id in any namespace.
                let glyph = c
                    .parent_element()
                    .and_then(|g| g.attributes().find(|a| a.name() == "id"))
                    .map_or_else(|| "a glyph".to_owned(), |a| format!("'{}'", a.value()));
                let message = format!(
                    "the bounding box of {glyph} lacks its <position> or <dimensions>, which SBML Layout requires"
                );
                self.error(c, &message)
            }
        }
        *slot = Some(bb);
        true
    }

    fn curve_child(&mut self, c: Node<'a, 'i>, pkg: P, slot: &mut Option<Curve>) -> bool {
        let lns = pkg.unwrap_or(ns::LAYOUT);
        if self.name_in(c, lns) != Some("curve") {
            return false;
        }
        let mut curve = Curve::default();
        for k in c.children().filter(|k| k.is_element()) {
            if self.name_in(k, lns) != Some("listOfCurveSegments") {
                self.warn(k, "dropped from <curve>");
                continue;
            }
            curve.segments = self.list(k, lns, &["curveSegment"], |r, s| r.segment(s, pkg));
        }
        *slot = Some(curve);
        true
    }

    fn segment(&mut self, n: Node<'a, 'i>, pkg: P) -> Segment {
        let lns = pkg.unwrap_or(ns::LAYOUT);
        let kind = n
            .attribute((ns::XSI, "type"))
            .map(|t| t.rsplit(':').next().unwrap_or(t));
        let mut points = [None; 4];
        for k in n.children().filter(|k| k.is_element()) {
            let i = match self.name_in(k, lns) {
                Some("start") => 0,
                Some("end") => 1,
                Some("basePoint1") => 2,
                Some("basePoint2") => 3,
                _ => {
                    self.warn(k, "dropped from <curveSegment>");
                    continue;
                }
            };
            points[i] = Some(self.point(k, pkg));
        }
        let [start, end, b1, b2] = points;
        let (start, end) = match (start, end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                self.error(
                    n,
                    "a curve segment lacks its <start> or <end>, which SBML Layout requires",
                );
                (Point::default(), Point::default())
            }
        };
        let bezier = match kind {
            Some("CubicBezier") => true,
            Some("LineSegment") => false,
            _ => b1.is_some() || b2.is_some(),
        };
        if bezier {
            // A Bézier missing a base point degenerates towards its ends.
            Segment::Bezier {
                start,
                end,
                base1: b1.unwrap_or(start),
                base2: b2.unwrap_or(end),
            }
        } else {
            Segment::Line { start, end }
        }
    }

    fn point(&mut self, n: Node<'a, 'i>, pkg: P) -> Point {
        let mut e = E::new(n);
        let x = self.double_pkg(&mut e, pkg, "x");
        let y = self.double_pkg(&mut e, pkg, "y");
        if x.is_none() || y.is_none() {
            self.error(n, &format!("<{}> needs x and y", n.tag_name().name()));
        }
        Point {
            x: x.unwrap_or(0.0),
            y: y.unwrap_or(0.0),
            z: self.double_pkg(&mut e, pkg, "z"),
        }
    }

    fn dimensions(&mut self, n: Node<'a, 'i>, pkg: P) -> Dimensions {
        let mut e = E::new(n);
        let w = self.double_pkg(&mut e, pkg, "width");
        let h = self.double_pkg(&mut e, pkg, "height");
        if w.is_none() || h.is_none() {
            self.error(n, "<dimensions> needs width and height");
        }
        Dimensions {
            width: w.unwrap_or(0.0),
            height: h.unwrap_or(0.0),
            depth: self.double_pkg(&mut e, pkg, "depth"),
        }
    }
}

const ADDITIONAL: &[&str] = &[
    "graphicalObject",
    "generalGlyph",
    "compartmentGlyph",
    "speciesGlyph",
    "reactionGlyph",
    "textGlyph",
];

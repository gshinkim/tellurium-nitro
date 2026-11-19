//! Writing the Layout package: in Level 3 as package elements with the
//! layout prefix, in Level 2 into the model's annotation.

use super::{Attrs, W, attr, fmt_f64};
use crate::layout::*;
use crate::model::{Meta, Model};
use crate::ns;

impl W {
    /// The `layout:` prefix in Level 3; nothing in Level 2.
    fn lp(&self) -> String {
        if self.level >= 3 {
            format!("{}:", self.layout_prefix)
        } else {
            String::new()
        }
    }

    pub(crate) fn l3_layouts(&mut self, m: &Model) {
        if !super::has_layout_list(m) {
            return;
        }
        let p = self.lp();
        let mut a = Vec::new();
        if !self.root_declares(ns::XSI) {
            a.push(("xmlns:xsi".to_string(), ns::XSI.to_string()));
        }
        self.annotated(
            &format!("{p}listOfLayouts"),
            &mut a,
            &m.layouts_meta,
            None,
            |w| {
                for l in &m.layouts {
                    w.layout(l);
                }
                if let Some(g) = &m.global_render {
                    w.render_list("listOfGlobalRenderInformation", g);
                }
            },
        );
    }

    /// Level 2 layouts, as raw XML for the model's annotation.
    pub(crate) fn l2_layouts(&mut self, m: &Model) -> String {
        let saved = std::mem::take(&mut self.out);
        // `annotated` writes the annotation one level in from the model.
        self.out.push('\n');
        self.push_depth();
        self.push_depth();
        let mut a = vec![
            ("xmlns".to_string(), ns::LAYOUT_L2.to_string()),
            ("xmlns:xsi".to_string(), ns::XSI.to_string()),
        ];
        let head = m.global_render.as_ref().map(|g| {
            self.l2_render_head(
                "listOfGlobalRenderInformation",
                g,
                m.layouts_meta.annotation.is_some(),
            )
        });
        self.annotated("listOfLayouts", &mut a, &m.layouts_meta, head, |w| {
            for l in &m.layouts {
                w.layout(l);
            }
        });
        self.pop_depth();
        self.pop_depth();
        // No newline after it: what follows in the annotation keeps its own
        // whitespace, so the annotation reads back unchanged.
        let mut xml = std::mem::replace(&mut self.out, saved);
        xml.truncate(xml.trim_end().len());
        xml
    }

    /// An attribute of a layout element.
    fn la<T: super::ToAttr + ?Sized>(&self, a: &mut Attrs, name: &str, value: &T) {
        attr(a, &format!("{}{name}", self.lp()), value);
    }

    fn glyph(&mut self, tag: &str, mut a: Attrs, meta: &Meta, body: impl FnOnce(&mut Self)) {
        let tag = format!("{}{tag}", self.lp());
        self.annotated(&tag, &mut a, meta, None, body);
    }

    fn layout(&mut self, l: &Layout) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &l.id);
        self.la(&mut a, "name", &l.name);
        let head = match &l.render {
            Some(r) if self.level < 3 => {
                Some(self.l2_render_head("listOfRenderInformation", r, l.meta.annotation.is_some()))
            }
            _ => None,
        };
        let tag = format!("{}layout", self.lp());
        self.annotated(&tag, &mut a, &l.meta, head, |w| {
            w.dimensions("dimensions", &l.dimensions);
            let p = w.lp();
            w.list(
                &format!("{p}listOfCompartmentGlyphs"),
                &l.compartment_glyphs,
                |w, g| w.compartment_glyph(g),
            );
            w.list(
                &format!("{p}listOfSpeciesGlyphs"),
                &l.species_glyphs,
                |w, g| w.species_glyph(g),
            );
            w.list(
                &format!("{p}listOfReactionGlyphs"),
                &l.reaction_glyphs,
                |w, g| w.reaction_glyph(g),
            );
            w.list(&format!("{p}listOfTextGlyphs"), &l.text_glyphs, |w, g| {
                w.text_glyph(g)
            });
            w.list(
                &format!("{p}listOfAdditionalGraphicalObjects"),
                &l.additional,
                |w, g| w.additional(g),
            );
            if let (Some(r), true) = (&l.render, w.level >= 3) {
                w.render_list("listOfRenderInformation", r);
            }
        });
    }

    fn additional(&mut self, g: &AdditionalGlyph) {
        match g {
            AdditionalGlyph::Graphical(g) => {
                let mut a = Vec::new();
                self.la(&mut a, "id", &g.id);
                self.glyph("graphicalObject", a, &g.meta, |w| {
                    w.bounding_box(&g.bounding_box)
                });
            }
            AdditionalGlyph::General(g) => self.general_glyph(g),
            AdditionalGlyph::Compartment(g) => self.compartment_glyph(g),
            AdditionalGlyph::Species(g) => self.species_glyph(g),
            AdditionalGlyph::Reaction(g) => self.reaction_glyph(g),
            AdditionalGlyph::Text(g) => self.text_glyph(g),
        }
    }

    fn compartment_glyph(&mut self, g: &CompartmentGlyph) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &g.id);
        self.la(&mut a, "compartment", &g.compartment);
        self.la(&mut a, "order", &g.order);
        self.glyph("compartmentGlyph", a, &g.meta, |w| {
            w.bounding_box(&g.bounding_box)
        });
    }

    fn species_glyph(&mut self, g: &SpeciesGlyph) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &g.id);
        self.la(&mut a, "species", &g.species);
        self.glyph("speciesGlyph", a, &g.meta, |w| {
            w.bounding_box(&g.bounding_box)
        });
    }

    fn reaction_glyph(&mut self, g: &ReactionGlyph) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &g.id);
        self.la(&mut a, "reaction", &g.reaction);
        self.glyph("reactionGlyph", a, &g.meta, |w| {
            w.bounding_box(&g.bounding_box);
            w.curve(&g.curve);
            let p = w.lp();
            w.list(
                &format!("{p}listOfSpeciesReferenceGlyphs"),
                &g.species_reference_glyphs,
                |w, s| {
                    let mut a = Vec::new();
                    w.la(&mut a, "id", &s.id);
                    w.la(&mut a, "speciesGlyph", &s.species_glyph);
                    w.la(&mut a, "speciesReference", &s.species_reference);
                    w.la(
                        &mut a,
                        "role",
                        &s.role.as_ref().map(|r| r.as_str().to_string()),
                    );
                    w.glyph("speciesReferenceGlyph", a, &s.meta, |w| {
                        w.bounding_box(&s.bounding_box);
                        w.curve(&s.curve);
                    });
                },
            );
        });
    }

    fn text_glyph(&mut self, g: &TextGlyph) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &g.id);
        self.la(&mut a, "text", &g.text);
        self.la(&mut a, "graphicalObject", &g.graphical_object);
        self.la(&mut a, "originOfText", &g.origin_of_text);
        self.glyph("textGlyph", a, &g.meta, |w| w.bounding_box(&g.bounding_box));
    }

    fn general_glyph(&mut self, g: &GeneralGlyph) {
        let mut a = Vec::new();
        self.la(&mut a, "id", &g.id);
        self.la(&mut a, "reference", &g.reference);
        self.glyph("generalGlyph", a, &g.meta, |w| {
            w.bounding_box(&g.bounding_box);
            w.curve(&g.curve);
            let p = w.lp();
            w.list(
                &format!("{p}listOfReferenceGlyphs"),
                &g.reference_glyphs,
                |w, r| {
                    let mut a = Vec::new();
                    w.la(&mut a, "id", &r.id);
                    w.la(&mut a, "glyph", &r.glyph);
                    w.la(&mut a, "reference", &r.reference);
                    w.la(&mut a, "role", &r.role);
                    w.glyph("referenceGlyph", a, &r.meta, |w| {
                        w.bounding_box(&r.bounding_box);
                        w.curve(&r.curve);
                    });
                },
            );
            w.list(&format!("{p}listOfSubGlyphs"), &g.sub_glyphs, |w, s| {
                w.additional(s)
            });
        });
    }

    pub(super) fn bounding_box(&mut self, bb: &Option<BoundingBox>) {
        let Some(bb) = bb else { return };
        let p = self.lp();
        let id = match &bb.id {
            Some(id) => format!(" {p}id=\"{}\"", super::escape(id, true)),
            None => String::new(),
        };
        self.line(&format!("<{p}boundingBox{id}>"));
        self.push_depth();
        self.point("position", &bb.position);
        self.dimensions("dimensions", &bb.dimensions);
        self.pop_depth();
        self.line(&format!("</{p}boundingBox>"));
    }

    fn point(&mut self, tag: &str, pt: &Point) {
        let p = self.lp();
        let z =
            pt.z.map(|z| format!(" {p}z=\"{}\"", fmt_f64(z)))
                .unwrap_or_default();
        self.line(&format!(
            "<{p}{tag} {p}x=\"{}\" {p}y=\"{}\"{z}/>",
            fmt_f64(pt.x),
            fmt_f64(pt.y)
        ));
    }

    fn dimensions(&mut self, tag: &str, d: &Dimensions) {
        let p = self.lp();
        let depth = d
            .depth
            .map(|z| format!(" {p}depth=\"{}\"", fmt_f64(z)))
            .unwrap_or_default();
        self.line(&format!(
            "<{p}{tag} {p}width=\"{}\" {p}height=\"{}\"{depth}/>",
            fmt_f64(d.width),
            fmt_f64(d.height)
        ));
    }

    fn curve(&mut self, c: &Option<Curve>) {
        let Some(c) = c else { return };
        let p = self.lp();
        self.line(&format!("<{p}curve>"));
        self.push_depth();
        self.list(&format!("{p}listOfCurveSegments"), &c.segments, |w, s| {
            let kind = match s {
                Segment::Line { .. } => "LineSegment",
                Segment::Bezier { .. } => "CubicBezier",
            };
            w.line(&format!("<{p}curveSegment xsi:type=\"{kind}\">"));
            w.push_depth();
            w.point("start", &s.start());
            w.point("end", &s.end());
            if let Segment::Bezier { base1, base2, .. } = s {
                w.point("basePoint1", base1);
                w.point("basePoint2", base2);
            }
            w.pop_depth();
            w.line(&format!("</{p}curveSegment>"));
        });
        self.pop_depth();
        self.line(&format!("</{p}curve>"));
    }
}

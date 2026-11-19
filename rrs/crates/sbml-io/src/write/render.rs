//! Writing the Render package: in Level 3 as package elements with the
//! render prefix declared on `<sbml>`; in Level 2 into annotations, with the
//! Level 2 Render namespace declared on the list.

use super::{Attrs, W, attr, escape};
use crate::ns;
use crate::render::*;

impl W {
    /// The `render:` prefix.
    fn rp(&self) -> String {
        format!("{}:", self.render_prefix)
    }

    /// An attribute of a render element.
    fn ra<T: super::ToAttr + ?Sized>(&self, a: &mut Attrs, name: &str, value: &T) {
        attr(a, &format!("{}{name}", self.rp()), value);
    }

    fn ra_rel(&self, a: &mut Attrs, name: &str, value: &Option<RelAbs>) {
        self.ra(a, name, &value.as_ref().map(|v| v.0.clone()));
    }

    /// `listOfGlobalRenderInformation` or `listOfRenderInformation`.
    pub(crate) fn render_list(&mut self, tag: &str, list: &RenderList) {
        let mut a = Vec::new();
        if self.level < 3 {
            a.push((
                format!("xmlns:{}", self.render_prefix),
                ns::RENDER_L2.to_string(),
            ));
        }
        let tag = format!("{}{tag}", self.rp());
        self.annotated(&tag, &mut a, &list.meta, None, |w| {
            for ri in &list.items {
                w.render_information(ri);
            }
        });
    }

    /// A Level 2 render list, as raw XML to go first in the annotation of an
    /// element written at the current depth.
    pub(crate) fn l2_render_head(
        &mut self,
        tag: &str,
        list: &RenderList,
        has_rest: bool,
    ) -> String {
        let saved = std::mem::take(&mut self.out);
        self.out.push('\n');
        self.push_depth();
        self.push_depth();
        self.render_list(tag, list);
        self.pop_depth();
        self.pop_depth();
        let mut xml = std::mem::replace(&mut self.out, saved);
        xml.truncate(xml.trim_end().len());
        if !has_rest {
            xml.push('\n');
            xml.push_str(&"  ".repeat(self.depth + 1));
        }
        xml
    }

    fn render_information(&mut self, ri: &RenderInformation) {
        let p = self.rp();
        let mut a = Vec::new();
        self.ra(&mut a, "id", &ri.id);
        self.ra(&mut a, "name", &ri.name);
        self.ra(&mut a, "programName", &ri.program_name);
        self.ra(&mut a, "programVersion", &ri.program_version);
        self.ra(
            &mut a,
            "referenceRenderInformation",
            &ri.reference_render_information,
        );
        self.ra(&mut a, "backgroundColor", &ri.background_color);
        self.annotated(
            &format!("{p}renderInformation"),
            &mut a,
            &ri.meta,
            None,
            |w| {
                w.list(
                    &format!("{p}listOfColorDefinitions"),
                    &ri.color_definitions,
                    |w, c| {
                        let mut a = Vec::new();
                        w.ra(&mut a, "id", &c.id);
                        w.ra(&mut a, "name", &c.name);
                        w.ra(&mut a, "value", &c.value);
                        w.element(&format!("{p}colorDefinition"), a, &c.meta, |_| {});
                    },
                );
                w.list(
                    &format!("{p}listOfGradientDefinitions"),
                    &ri.gradients,
                    |w, g| w.gradient(g),
                );
                w.list(
                    &format!("{p}listOfLineEndings"),
                    &ri.line_endings,
                    |w, l| w.line_ending(l),
                );
                w.list(&format!("{p}listOfStyles"), &ri.styles, |w, s| w.style(s));
            },
        );
    }

    fn gradient(&mut self, g: &Gradient) {
        let p = self.rp();
        let mut a = Vec::new();
        self.ra(&mut a, "id", &g.id);
        self.ra(&mut a, "name", &g.name);
        self.ra(&mut a, "spreadMethod", &g.spread_method);
        let tag = match &g.kind {
            GradientKind::Linear {
                x1,
                y1,
                z1,
                x2,
                y2,
                z2,
            } => {
                for (n, v) in [
                    ("x1", x1),
                    ("y1", y1),
                    ("z1", z1),
                    ("x2", x2),
                    ("y2", y2),
                    ("z2", z2),
                ] {
                    self.ra_rel(&mut a, n, v);
                }
                "linearGradient"
            }
            GradientKind::Radial {
                cx,
                cy,
                cz,
                r,
                fx,
                fy,
                fz,
            } => {
                for (n, v) in [
                    ("cx", cx),
                    ("cy", cy),
                    ("cz", cz),
                    ("r", r),
                    ("fx", fx),
                    ("fy", fy),
                    ("fz", fz),
                ] {
                    self.ra_rel(&mut a, n, v);
                }
                "radialGradient"
            }
        };
        self.annotated(&format!("{p}{tag}"), &mut a, &g.meta, None, |w| {
            for s in &g.stops {
                let mut a = Vec::new();
                w.ra_rel(&mut a, "offset", &s.offset);
                w.ra(&mut a, "stop-color", &s.stop_color);
                w.element(&format!("{p}stop"), a, &s.meta, |_| {});
            }
        });
    }

    fn line_ending(&mut self, l: &LineEnding) {
        let p = self.rp();
        let mut a = Vec::new();
        self.ra(&mut a, "id", &l.id);
        self.ra(
            &mut a,
            "enableRotationalMapping",
            &l.enable_rotational_mapping,
        );
        self.annotated(&format!("{p}lineEnding"), &mut a, &l.meta, None, |w| {
            w.bounding_box(&l.bounding_box);
            if let Some(g) = &l.group {
                w.group(g);
            }
        });
    }

    fn style(&mut self, s: &Style) {
        let p = self.rp();
        let mut a = Vec::new();
        self.ra(&mut a, "id", &s.id);
        self.ra(&mut a, "name", &s.name);
        for (name, list) in [
            ("roleList", &s.role_list),
            ("typeList", &s.type_list),
            ("idList", &s.id_list),
        ] {
            if !list.is_empty() {
                self.ra(&mut a, name, &list.join(" "));
            }
        }
        self.annotated(&format!("{p}style"), &mut a, &s.meta, None, |w| {
            if let Some(g) = &s.group {
                w.group(g);
            }
        });
    }

    fn presentation(&self, a: &mut Attrs, p: &Presentation) {
        self.ra(a, "stroke", &p.stroke);
        self.ra(a, "stroke-width", &p.stroke_width);
        self.ra(a, "stroke-dasharray", &p.stroke_dasharray);
        self.ra(a, "fill", &p.fill);
        self.ra(a, "fill-rule", &p.fill_rule);
        self.ra(a, "font-family", &p.font_family);
        self.ra_rel(a, "font-size", &p.font_size);
        self.ra(a, "font-weight", &p.font_weight);
        self.ra(a, "font-style", &p.font_style);
        self.ra(a, "text-anchor", &p.text_anchor);
        self.ra(a, "vtext-anchor", &p.vtext_anchor);
        self.ra(a, "startHead", &p.start_head);
        self.ra(a, "endHead", &p.end_head);
        self.ra(a, "transform", &p.transform);
    }

    /// `id`, `name` and the presentation attributes, which every drawing
    /// element has.
    fn drawn(&self, id: &Option<String>, name: &Option<String>, p: &Presentation) -> Attrs {
        let mut a = Vec::new();
        self.ra(&mut a, "id", id);
        self.ra(&mut a, "name", name);
        self.presentation(&mut a, p);
        a
    }

    fn group(&mut self, g: &Group) {
        let mut a = self.drawn(&g.id, &g.name, &g.presentation);
        let tag = format!("{}g", self.rp());
        self.annotated(&tag, &mut a, &g.meta, None, |w| {
            for c in &g.children {
                w.primitive(c);
            }
        });
    }

    fn primitive(&mut self, prim: &Primitive) {
        let p = self.rp();
        match prim {
            Primitive::Group(g) => self.group(g),
            Primitive::Rectangle(r) => {
                let mut a = self.drawn(&r.id, &r.name, &r.presentation);
                for (n, v) in [
                    ("x", &r.x),
                    ("y", &r.y),
                    ("z", &r.z),
                    ("width", &r.width),
                    ("height", &r.height),
                    ("rx", &r.rx),
                    ("ry", &r.ry),
                ] {
                    self.ra_rel(&mut a, n, v);
                }
                self.ra(&mut a, "ratio", &r.ratio);
                self.element(&format!("{p}rectangle"), a, &r.meta, |_| {});
            }
            Primitive::Ellipse(e) => {
                let mut a = self.drawn(&e.id, &e.name, &e.presentation);
                for (n, v) in [
                    ("cx", &e.cx),
                    ("cy", &e.cy),
                    ("cz", &e.cz),
                    ("rx", &e.rx),
                    ("ry", &e.ry),
                ] {
                    self.ra_rel(&mut a, n, v);
                }
                self.ra(&mut a, "ratio", &e.ratio);
                self.element(&format!("{p}ellipse"), a, &e.meta, |_| {});
            }
            Primitive::Polygon(g) => {
                let a = self.drawn(&g.id, &g.name, &g.presentation);
                self.element(&format!("{p}polygon"), a, &g.meta, |w| {
                    w.points(&g.elements)
                });
            }
            Primitive::Curve(c) => {
                let a = self.drawn(&c.id, &c.name, &c.presentation);
                self.element(&format!("{p}curve"), a, &c.meta, |w| w.points(&c.elements));
            }
            Primitive::Image(i) => {
                let mut a = self.drawn(&i.id, &i.name, &i.presentation);
                for (n, v) in [
                    ("x", &i.x),
                    ("y", &i.y),
                    ("z", &i.z),
                    ("width", &i.width),
                    ("height", &i.height),
                ] {
                    self.ra_rel(&mut a, n, v);
                }
                self.ra(&mut a, "href", &i.href);
                self.element(&format!("{p}image"), a, &i.meta, |_| {});
            }
            Primitive::Text(t) => {
                let mut a = self.drawn(&t.id, &t.name, &t.presentation);
                for (n, v) in [("x", &t.x), ("y", &t.y), ("z", &t.z)] {
                    self.ra_rel(&mut a, n, v);
                }
                self.text_element(&format!("{p}text"), a, &t.meta, &t.text);
            }
        }
    }

    fn points(&mut self, points: &[RenderPoint]) {
        let p = self.rp();
        self.list(&format!("{p}listOfElements"), points, |w, pt| {
            let mut a: Attrs = vec![(
                "xsi:type".to_string(),
                if pt.bezier.is_some() {
                    "RenderCubicBezier"
                } else {
                    "RenderPoint"
                }
                .to_string(),
            )];
            for (n, v) in [("x", &pt.x), ("y", &pt.y), ("z", &pt.z)] {
                w.ra_rel(&mut a, n, v);
            }
            if let Some(b) = &pt.bezier {
                for (n, v) in [
                    ("basePoint1_x", &b.x1),
                    ("basePoint1_y", &b.y1),
                    ("basePoint1_z", &b.z1),
                    ("basePoint2_x", &b.x2),
                    ("basePoint2_y", &b.y2),
                    ("basePoint2_z", &b.z2),
                ] {
                    w.ra_rel(&mut a, n, v);
                }
            }
            w.element(&format!("{p}element"), a, &pt.meta, |_| {});
        });
    }

    /// An element whose content is text: written on one line, so the text is
    /// kept exactly. Its notes, annotation and kept elements, which Render
    /// doesn't give a text, are left out.
    fn text_element(&mut self, tag: &str, mut a: Attrs, meta: &crate::model::Meta, text: &str) {
        attr(&mut a, "metaid", &meta.metaid);
        attr(&mut a, "sboTerm", &meta.sbo_term);
        for x in &meta.extra.attributes {
            a.push((x.name.clone(), x.value.clone()));
        }
        let mut line = format!("<{tag}");
        for (k, v) in &a {
            line.push_str(&format!(" {k}=\"{}\"", escape(v, true)));
        }
        line.push_str(&format!(">{}</{tag}>", escape(text, false)));
        self.line(&line);
    }
}

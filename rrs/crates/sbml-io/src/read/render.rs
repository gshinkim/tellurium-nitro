//! Reading the Render package. The same code reads Level 3 render
//! information and the Level 2 kind kept in annotations; `pkg` is the Render
//! namespace in each case. Level 2 files from older tools write polygons and
//! curves as `listOfCurveSegments`; those are read into the same points.

use roxmltree::Node;

use super::{E, Reader};
use crate::layout::BoundingBox;
use crate::ns;
use crate::render::*;

type P<'p> = Option<&'p str>;

impl<'a, 'i: 'a> Reader<'a, 'i> {
    /// The local name of an element of Render: in its namespace, or, as some
    /// tools write the drawing inside a style, with no prefix in the core
    /// namespace. libSBML reads both.
    fn rname<'n>(&self, n: Node<'n, 'i>, rns: &str) -> Option<&'n str> {
        self.name_in(n, rns).or_else(|| self.core_name(n))
    }

    /// `list`, with items named as [`rname`](Self::rname) accepts them.
    fn rlist<T>(
        &mut self,
        n: Node<'a, 'i>,
        rns: &str,
        names: &[&str],
        mut item: impl FnMut(&mut Self, Node<'a, 'i>) -> T,
    ) -> Vec<T> {
        let mut out = Vec::new();
        let list_name = n.tag_name().name();
        for c in n.children().filter(|c| c.is_element()) {
            if self.rname(c, rns).is_some_and(|name| names.contains(&name)) {
                out.push(item(self, c));
            } else {
                let msg = format!("<{}> in <{list_name}> is dropped", c.tag_name().name());
                self.warn(c, &msg);
            }
        }
        out
    }

    /// `listOfGlobalRenderInformation` or `listOfRenderInformation`.
    pub(super) fn render_list(&mut self, n: Node<'a, 'i>, pkg: &str) -> RenderList {
        let mut items = Vec::new();
        let meta = self.meta(E::new(n), Some(pkg), |r, c| {
            if r.name_in(c, pkg) == Some("renderInformation") {
                items.push(r.render_information(c, Some(pkg)));
                true
            } else {
                false
            }
        });
        RenderList { items, meta }
    }

    fn render_information(&mut self, n: Node<'a, 'i>, pkg: P) -> RenderInformation {
        let rns = pkg.unwrap_or(ns::RENDER);
        let mut e = E::new(n);
        let mut ri = RenderInformation {
            id: self.req_pkg(&mut e, pkg, "id"),
            name: s(e.take_pkg(pkg, "name")),
            program_name: s(e.take_pkg(pkg, "programName")),
            program_version: s(e.take_pkg(pkg, "programVersion")),
            reference_render_information: s(e.take_pkg(pkg, "referenceRenderInformation")),
            background_color: s(e.take_pkg(pkg, "backgroundColor")),
            ..Default::default()
        };
        ri.meta = self.meta(e, pkg, |r, c| {
            match r.rname(c, rns) {
                Some("listOfColorDefinitions") => {
                    ri.color_definitions = r.rlist(c, rns, &["colorDefinition"], |r, c| {
                        let mut e = E::new(c);
                        ColorDefinition {
                            id: r.req_pkg(&mut e, pkg, "id"),
                            name: s(e.take_pkg(pkg, "name")),
                            value: r.req_pkg(&mut e, pkg, "value"),
                            meta: r.meta(e, pkg, |_, _| false),
                        }
                    })
                }
                Some("listOfGradientDefinitions") => {
                    ri.gradients = r.rlist(c, rns, &["linearGradient", "radialGradient"], |r, c| {
                        r.gradient(c, pkg)
                    })
                }
                Some("listOfLineEndings") => {
                    ri.line_endings = r.rlist(c, rns, &["lineEnding"], |r, c| r.line_ending(c, pkg))
                }
                Some("listOfStyles") => {
                    ri.styles = r.rlist(c, rns, &["style"], |r, c| r.style(c, pkg))
                }
                _ => return false,
            }
            true
        });
        ri
    }

    fn gradient(&mut self, n: Node<'a, 'i>, pkg: P) -> Gradient {
        let rns = pkg.unwrap_or(ns::RENDER);
        let mut e = E::new(n);
        let id = self.req_pkg(&mut e, pkg, "id");
        let name = s(e.take_pkg(pkg, "name"));
        let spread_method = s(e.take_pkg(pkg, "spreadMethod"));
        let mut ra = |name| rel(e.take_pkg(pkg, name));
        let kind = if n.tag_name().name() == "radialGradient" {
            GradientKind::Radial {
                cx: ra("cx"),
                cy: ra("cy"),
                cz: ra("cz"),
                r: ra("r"),
                fx: ra("fx"),
                fy: ra("fy"),
                fz: ra("fz"),
            }
        } else {
            GradientKind::Linear {
                x1: ra("x1"),
                y1: ra("y1"),
                z1: ra("z1"),
                x2: ra("x2"),
                y2: ra("y2"),
                z2: ra("z2"),
            }
        };
        let mut stops = Vec::new();
        let meta = self.meta(e, pkg, |r, c| {
            if r.rname(c, rns) != Some("stop") {
                return false;
            }
            let mut e = E::new(c);
            stops.push(GradientStop {
                offset: rel(e.take_pkg(pkg, "offset")),
                stop_color: s(e.take_pkg(pkg, "stop-color")),
                meta: r.meta(e, pkg, |_, _| false),
            });
            true
        });
        Gradient {
            id,
            name,
            spread_method,
            kind,
            stops,
            meta,
        }
    }

    fn line_ending(&mut self, n: Node<'a, 'i>, pkg: P) -> LineEnding {
        let rns = pkg.unwrap_or(ns::RENDER);
        // The bounding box is the Layout package's: in Level 2 the Level 2
        // layout namespace, or, in older files, Render's own.
        let lns = if rns == ns::RENDER_L2 {
            ns::LAYOUT_L2
        } else {
            ns::LAYOUT
        };
        let mut e = E::new(n);
        let mut le = LineEnding {
            id: self.req_pkg(&mut e, pkg, "id"),
            enable_rotational_mapping: self.parsed(
                &mut e,
                pkg,
                "enableRotationalMapping",
                "true or false",
                |v| match v {
                    "true" | "1" => Some(true),
                    "false" | "0" => Some(false),
                    _ => None,
                },
            ),
            ..Default::default()
        };
        let mut bb: Option<BoundingBox> = None;
        le.meta = self.meta(e, pkg, |r, c| {
            if r.rname(c, rns) == Some("g") {
                le.group = Some(r.group(c, pkg));
                return true;
            }
            r.bounding_box_child(c, Some(lns), &mut bb)
                || r.bounding_box_child(c, Some(rns), &mut bb)
        });
        le.bounding_box = bb;
        le
    }

    fn style(&mut self, n: Node<'a, 'i>, pkg: P) -> Style {
        let rns = pkg.unwrap_or(ns::RENDER);
        let mut e = E::new(n);
        let list = |v: Option<&str>| -> Vec<String> {
            v.map(|v| v.split_whitespace().map(str::to_string).collect())
                .unwrap_or_default()
        };
        let mut st = Style {
            id: s(e.take_pkg(pkg, "id")),
            name: s(e.take_pkg(pkg, "name")),
            role_list: list(e.take_pkg(pkg, "roleList")),
            type_list: list(e.take_pkg(pkg, "typeList")),
            id_list: list(e.take_pkg(pkg, "idList")),
            ..Default::default()
        };
        st.meta = self.meta(e, pkg, |r, c| {
            if r.rname(c, rns) == Some("g") && st.group.is_none() {
                st.group = Some(r.group(c, pkg));
                true
            } else {
                false
            }
        });
        st
    }

    fn presentation(&mut self, e: &mut E<'a, 'i>, pkg: P) -> Presentation {
        Presentation {
            stroke: s(e.take_pkg(pkg, "stroke")),
            stroke_width: self.double_pkg(e, pkg, "stroke-width"),
            stroke_dasharray: s(e.take_pkg(pkg, "stroke-dasharray")),
            fill: s(e.take_pkg(pkg, "fill")),
            fill_rule: s(e.take_pkg(pkg, "fill-rule")),
            font_family: s(e.take_pkg(pkg, "font-family")),
            font_size: rel(e.take_pkg(pkg, "font-size")),
            font_weight: s(e.take_pkg(pkg, "font-weight")),
            font_style: s(e.take_pkg(pkg, "font-style")),
            text_anchor: s(e.take_pkg(pkg, "text-anchor")),
            vtext_anchor: s(e.take_pkg(pkg, "vtext-anchor")),
            start_head: s(e.take_pkg(pkg, "startHead")),
            end_head: s(e.take_pkg(pkg, "endHead")),
            transform: s(e.take_pkg(pkg, "transform")),
        }
    }

    fn group(&mut self, n: Node<'a, 'i>, pkg: P) -> Group {
        let mut e = E::new(n);
        let mut g = Group {
            id: s(e.take_pkg(pkg, "id")),
            name: s(e.take_pkg(pkg, "name")),
            presentation: self.presentation(&mut e, pkg),
            ..Default::default()
        };
        g.meta = self.meta(e, pkg, |r, c| match r.primitive(c, pkg) {
            Some(p) => {
                g.children.push(p);
                true
            }
            None => false,
        });
        g
    }

    fn primitive(&mut self, n: Node<'a, 'i>, pkg: P) -> Option<Primitive> {
        let rns = pkg.unwrap_or(ns::RENDER);
        let name = self.rname(n, rns)?;
        if name == "g" {
            return Some(Primitive::Group(self.group(n, pkg)));
        }
        if !matches!(
            name,
            "rectangle" | "ellipse" | "polygon" | "curve" | "text" | "image"
        ) {
            return None;
        }
        let mut e = E::new(n);
        let id = s(e.take_pkg(pkg, "id"));
        let nm = s(e.take_pkg(pkg, "name"));
        let presentation = self.presentation(&mut e, pkg);
        let mut ra = |name| rel(e.take_pkg(pkg, name));
        Some(match name {
            "rectangle" => {
                let mut r = Rectangle {
                    id,
                    name: nm,
                    presentation,
                    x: ra("x"),
                    y: ra("y"),
                    z: ra("z"),
                    width: ra("width"),
                    height: ra("height"),
                    rx: ra("rx"),
                    ry: ra("ry"),
                    ..Default::default()
                };
                r.ratio = self.double_pkg(&mut e, pkg, "ratio");
                r.meta = self.meta(e, pkg, |_, _| false);
                Primitive::Rectangle(r)
            }
            "ellipse" => {
                let mut el = Ellipse {
                    id,
                    name: nm,
                    presentation,
                    cx: ra("cx"),
                    cy: ra("cy"),
                    cz: ra("cz"),
                    rx: ra("rx"),
                    ry: ra("ry"),
                    ..Default::default()
                };
                el.ratio = self.double_pkg(&mut e, pkg, "ratio");
                el.meta = self.meta(e, pkg, |_, _| false);
                Primitive::Ellipse(el)
            }
            "text" => {
                let (x, y, z) = (ra("x"), ra("y"), ra("z"));
                let text: String = n
                    .children()
                    .filter(|c| c.is_text())
                    .filter_map(|c| c.text())
                    .collect();
                Primitive::Text(Text {
                    id,
                    name: nm,
                    presentation,
                    x,
                    y,
                    z,
                    text,
                    meta: self.meta(e, pkg, |_, _| false),
                })
            }
            "image" => {
                let (x, y, z, width, height) =
                    (ra("x"), ra("y"), ra("z"), ra("width"), ra("height"));
                let href = s(e.take_pkg(pkg, "href"));
                Primitive::Image(Image {
                    id,
                    name: nm,
                    presentation,
                    x,
                    y,
                    z,
                    width,
                    height,
                    href,
                    meta: self.meta(e, pkg, |_, _| false),
                })
            }
            _ => {
                let mut elements = Vec::new();
                let meta = self.meta(e, pkg, |r, c| r.points_child(c, pkg, &mut elements));
                if name == "polygon" {
                    Primitive::Polygon(Polygon {
                        id,
                        name: nm,
                        presentation,
                        elements,
                        meta,
                    })
                } else {
                    Primitive::Curve(RenderCurve {
                        id,
                        name: nm,
                        presentation,
                        elements,
                        meta,
                    })
                }
            }
        })
    }

    /// `listOfElements`, or the older `listOfCurveSegments`.
    fn points_child(&mut self, c: Node<'a, 'i>, pkg: P, out: &mut Vec<RenderPoint>) -> bool {
        let rns = pkg.unwrap_or(ns::RENDER);
        match self.rname(c, rns) {
            Some("listOfElements") => {
                *out = self.rlist(c, rns, &["element"], |r, k| r.render_point(k, pkg));
                true
            }
            Some("listOfCurveSegments") => {
                for k in c.children().filter(|k| k.is_element()) {
                    if self.rname(k, rns) != Some("curveSegment") {
                        self.warn(k, "dropped from <listOfCurveSegments>");
                        continue;
                    }
                    self.curve_segment(k, rns, out);
                }
                true
            }
            _ => false,
        }
    }

    fn render_point(&mut self, n: Node<'a, 'i>, pkg: P) -> RenderPoint {
        let kind = n
            .attribute((ns::XSI, "type"))
            .map(|t| t.rsplit(':').next().unwrap_or(t));
        let mut e = E::new(n);
        let mut ra = |name| rel(e.take_pkg(pkg, name));
        let (x, y, z) = (ra("x"), ra("y"), ra("z"));
        let base = BezierBase {
            x1: ra("basePoint1_x"),
            y1: ra("basePoint1_y"),
            z1: ra("basePoint1_z"),
            x2: ra("basePoint2_x"),
            y2: ra("basePoint2_y"),
            z2: ra("basePoint2_z"),
        };
        let bezier = kind == Some("RenderCubicBezier")
            || (kind.is_none() && (base.x1.is_some() || base.x2.is_some()));
        // xsi:type is the XSI namespace's, so it isn't taken above; it is
        // written from `bezier`, not kept.
        let used_type = n.attributes().position(|a| a.namespace() == Some(ns::XSI));
        if let Some(i) = used_type {
            e.used[i] = true;
        }
        RenderPoint {
            x,
            y,
            z,
            bezier: bezier.then_some(base),
            meta: self.meta(e, pkg, |_, _| false),
        }
    }

    /// An older `curveSegment`: its start, if it is the first, then its end.
    fn curve_segment(&mut self, n: Node<'a, 'i>, rns: &str, out: &mut Vec<RenderPoint>) {
        let kind = n
            .attribute((ns::XSI, "type"))
            .map(|t| t.rsplit(':').next().unwrap_or(t));
        let mut pts: [Option<RenderPoint>; 4] = Default::default();
        for k in n.children().filter(|k| k.is_element()) {
            let i = match self.rname(k, rns) {
                Some("start") => 0,
                Some("end") => 1,
                Some("basePoint1") => 2,
                Some("basePoint2") => 3,
                _ => continue,
            };
            let mut e = E::new(k);
            pts[i] = Some(RenderPoint {
                x: rel(e.take_pkg(Some(rns), "x")),
                y: rel(e.take_pkg(Some(rns), "y")),
                z: rel(e.take_pkg(Some(rns), "z")),
                ..Default::default()
            });
        }
        let [start, end, b1, b2] = pts;
        if out.is_empty() {
            out.push(start.unwrap_or_default());
        }
        let mut end = end.unwrap_or_default();
        if kind == Some("CubicBezier") || b1.is_some() || b2.is_some() {
            let (b1, b2) = (b1.unwrap_or_default(), b2.unwrap_or_default());
            end.bezier = Some(BezierBase {
                x1: b1.x,
                y1: b1.y,
                z1: b1.z,
                x2: b2.x,
                y2: b2.y,
                z2: b2.z,
            });
        }
        out.push(end);
    }
}

fn s(v: Option<&str>) -> Option<String> {
    v.map(str::to_string)
}

fn rel(v: Option<&str>) -> Option<RelAbs> {
    v.map(|v| RelAbs(v.to_string()))
}

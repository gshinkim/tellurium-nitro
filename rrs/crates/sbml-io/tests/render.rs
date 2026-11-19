//! The Render package, read from the files in tests/data. Those were written
//! by python-libsbml (tests/scripts/make_render_fixtures.py): the same model,
//! layout and render information (a global one and the layout's own, between
//! them using every part of Render) in Level 3 Version 1, and in Level 2
//! Version 4, where Render lives in annotations.

use sbml_io::render::*;
use sbml_io::*;

fn load(name: &str) -> Document {
    let path = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    let read = read(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert!(read.warnings.is_empty(), "{:?}", read.warnings);
    read.document
}

fn rel(s: &str) -> Option<RelAbs> {
    Some(RelAbs(s.to_string()))
}

fn owned(s: &str) -> Option<String> {
    Some(s.to_string())
}

fn group(style: &Style) -> &Group {
    style.group.as_ref().unwrap()
}

fn style<'a>(ri: &'a RenderInformation, id: &str) -> &'a Style {
    ri.styles
        .iter()
        .find(|s| s.id.as_deref() == Some(id))
        .unwrap()
}

#[test]
fn level_3_render_is_read_in_full() {
    let doc = load("render-l3v1.xml");
    let m = doc.model.as_ref().unwrap();
    // Nothing of Render is left as raw XML.
    assert!(m.meta.extra.is_empty(), "{:?}", m.meta.extra);
    assert_eq!(m.layouts_meta, Meta::default());
    let l = &m.layouts[0];
    assert!(l.meta.extra.is_empty(), "{:?}", l.meta.extra);

    // The global render information.
    let global = m.global_render.as_ref().unwrap();
    assert_eq!(global.items.len(), 1);
    let g = &global.items[0];
    assert_eq!(g.id, "global1");
    assert_eq!(g.name.as_deref(), Some("Global style"));
    assert_eq!(g.program_name.as_deref(), Some("make_render_fixtures.py"));
    assert_eq!(g.background_color.as_deref(), Some("#F8F9FAFF"));
    let colours: Vec<_> = g
        .color_definitions
        .iter()
        .map(|c| (c.id.as_str(), c.value.as_str()))
        .collect();
    assert_eq!(
        colours,
        [
            ("black", "#000000"),
            ("fill", "#eef6ff"),
            ("teal", "#007777cc"),
            ("none_like", "#ffffff00")
        ]
    );

    let lin = &g.gradients[0];
    assert_eq!(lin.id, "lin");
    assert_eq!(lin.spread_method.as_deref(), Some("reflect"));
    assert!(matches!(lin.kind, GradientKind::Linear { .. }));
    let stops: Vec<_> = lin
        .stops
        .iter()
        .map(|s| (s.offset.clone().unwrap().0, s.stop_color.clone().unwrap()))
        .collect();
    assert_eq!(
        stops,
        [
            ("0".to_string(), "#FFFFFF".to_string()),
            ("100%".to_string(), "fill".to_string())
        ]
    );
    let GradientKind::Radial { fx, fy, .. } = &g.gradients[1].kind else {
        panic!("not radial")
    };
    assert_eq!((fx.clone(), fy.clone()), (rel("40%"), rel("40%")));

    let arrow = &g.line_endings[0];
    assert_eq!(arrow.id, "arrow");
    assert_eq!(arrow.enable_rotational_mapping, None);
    assert_eq!(
        arrow.bounding_box.as_ref().unwrap().position,
        Point::new(-12.0, -5.0)
    );
    let Primitive::Polygon(tri) = &arrow.group.as_ref().unwrap().children[0] else {
        panic!("no polygon")
    };
    assert_eq!(tri.presentation.fill.as_deref(), Some("black"));
    let pts: Vec<_> = tri
        .elements
        .iter()
        .map(|p| (p.x.clone().unwrap().0, p.y.clone().unwrap().0))
        .collect();
    assert_eq!(
        pts,
        [
            ("0".to_string(), "0".to_string()),
            ("100%".to_string(), "50%".to_string()),
            ("0".to_string(), "100%".to_string())
        ]
    );
    assert_eq!(g.line_endings[1].enable_rotational_mapping, Some(false));

    // Styles by type, with every presentation attribute a group may have.
    let sp = style(g, "species");
    assert_eq!(sp.type_list, ["SPECIESGLYPH"]);
    let p = &group(sp).presentation;
    assert_eq!(
        *p,
        Presentation {
            stroke: owned("black"),
            stroke_width: Some(1.5),
            fill: owned("fill"),
            font_family: owned("sans-serif"),
            font_size: rel("12"),
            font_weight: owned("bold"),
            font_style: owned("italic"),
            text_anchor: owned("middle"),
            vtext_anchor: owned("middle"),
            ..Default::default()
        }
    );
    let Primitive::Rectangle(r) = &group(sp).children[0] else {
        panic!("no rectangle")
    };
    assert_eq!((r.width.clone(), r.rx.clone()), (rel("100%"), rel("8")));
    assert_eq!(r.presentation.stroke_dasharray.as_deref(), Some("5 , 4"));

    // Styles by role.
    let lines = style(g, "reaction_lines");
    assert_eq!(lines.role_list, ["product", "substrate"]);
    assert_eq!(group(lines).presentation.end_head.as_deref(), Some("arrow"));
    let modifiers = &group(style(g, "modifier_lines")).presentation;
    assert_eq!(modifiers.stroke_dasharray.as_deref(), Some("6 , 4"));
    let Primitive::Ellipse(e) = &group(style(g, "compartments")).children[0] else {
        panic!("no ellipse")
    };
    assert_eq!(e.ry, rel("-2+50%"));

    // The layout's own, by id, with the other primitives.
    let local = &l.render.as_ref().unwrap().items[0];
    assert_eq!(
        local.reference_render_information.as_deref(),
        Some("global1")
    );
    let just_p = style(local, "just_P");
    assert_eq!(just_p.id_list, ["sg_P"]);
    let Primitive::Polygon(hex) = &group(just_p).children[0] else {
        panic!("no polygon")
    };
    assert_eq!(hex.elements.len(), 7);
    assert!(hex.elements[..6].iter().all(|p| p.bezier.is_none()));
    let b = hex.elements[6].bezier.as_ref().unwrap();
    assert_eq!((b.x1.clone(), b.y2.clone()), (rel("10%"), rel("70%")));

    let just_e = style(local, "just_E");
    assert_eq!(just_e.id_list, ["sg_E", "tg_E"]);
    let kids = &group(just_e).children;
    let Primitive::Group(inner) = &kids[0] else {
        panic!("no inner group")
    };
    assert_eq!(inner.presentation.stroke_width, Some(3.0));
    let Primitive::Curve(c) = &inner.children[0] else {
        panic!("no curve")
    };
    assert_eq!(c.presentation.end_head.as_deref(), Some("arrow"));
    assert_eq!(c.elements.len(), 2);
    let Primitive::Text(t) = &kids[1] else {
        panic!("no text")
    };
    assert_eq!((t.text.as_str(), t.x.clone()), ("E", rel("50%")));
    let Primitive::Image(i) = &kids[2] else {
        panic!("no image")
    };
    assert_eq!(i.href.as_deref(), Some("enzyme.png"));
    assert_eq!(
        style(local, "note")
            .group
            .as_ref()
            .unwrap()
            .presentation
            .font_family
            .as_deref(),
        Some("serif")
    );
}

#[test]
fn level_2_render_is_read_as_level_3() {
    // The same render information, from annotations.
    let l3 = load("render-l3v1.xml");
    let l2 = load("render-l2v4.xml");
    let (m3, m2) = (l3.model.as_ref().unwrap(), l2.model.as_ref().unwrap());
    assert_eq!(m2.global_render, m3.global_render);
    assert_eq!(m2.layouts[0].render, m3.layouts[0].render);
    // Taken out of the annotations, which are then empty.
    assert_eq!(m2.layouts_meta.annotation, None);
    assert_eq!(m2.layouts[0].meta.annotation, None);
    assert_eq!(m2.meta.annotation, None);
}

#[test]
fn written_and_read_back_unchanged() {
    for name in ["render-l3v1.xml", "render-l2v4.xml"] {
        let doc = load(name);
        let text = write(&doc);
        let again = read(&text).unwrap();
        assert!(again.warnings.is_empty(), "{name}: {:?}", again.warnings);
        assert_eq!(again.document, doc, "{name}");
        assert_eq!(write(&again.document), text, "{name}: not stable");
    }
}

#[test]
fn render_made_in_code_declares_its_namespace() {
    let mut doc = Document::new(3, 1);
    let mut m = Model::default();
    let mut l = Layout {
        id: "layout1".into(),
        dimensions: Dimensions::new(100.0, 100.0),
        ..Default::default()
    };
    let ri = RenderInformation {
        id: "pd_render".into(),
        color_definitions: vec![ColorDefinition {
            id: "c1".into(),
            value: "#FF0000".into(),
            ..Default::default()
        }],
        styles: vec![Style {
            id_list: vec!["sg_A".into()],
            group: Some(Group {
                presentation: Presentation {
                    fill: owned("c1"),
                    ..Default::default()
                },
                children: vec![Primitive::Ellipse(Ellipse {
                    cx: Some(RelAbs::rel(50.0)),
                    cy: Some(RelAbs::rel(50.0)),
                    rx: Some(RelAbs::rel(50.0)),
                    ..Default::default()
                })],
                ..Default::default()
            }),
            ..Default::default()
        }],
        ..Default::default()
    };
    l.render = Some(RenderList {
        items: vec![ri],
        ..Default::default()
    });
    m.layouts.push(l);
    doc.model = Some(m);
    let text = write(&doc);
    assert!(
        text.contains(&format!("xmlns:render=\"{}\"", ns::RENDER)),
        "{text}"
    );
    assert!(text.contains("render:required=\"false\""), "{text}");
    assert!(text.contains("<render:ellipse render:cx=\"50%\""), "{text}");
    let back = read(&text).unwrap().document;
    assert_eq!(
        back.model.unwrap().layouts[0].render,
        doc.model.unwrap().layouts[0].render
    );
}

#[test]
fn older_level_2_curve_segments_become_points() {
    // Polygons as written by the SBW layout viewer (2011).
    let xml = r##"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1">
  <model id="m">
    <annotation>
      <listOfLayouts xmlns="http://projects.eml.org/bcb/sbml/level2">
        <layout id="l">
          <annotation>
            <listOfRenderInformation xmlns="http://projects.eml.org/bcb/sbml/render/level2">
              <renderInformation id="r">
                <listOfStyles>
                  <style idList="sg">
                    <g>
                      <polygon fill="#FF0000">
                        <listOfCurveSegments>
                          <curveSegment xsi:type="LineSegment" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
                            <start x="0" y="0%" />
                            <end x="100%" y="0" />
                          </curveSegment>
                          <curveSegment xsi:type="CubicBezier" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
                            <start x="100%" y="0" />
                            <end x="0" y="100%" />
                            <basePoint1 x="100%" y="50%" />
                            <basePoint2 x="50%" y="100%" />
                          </curveSegment>
                        </listOfCurveSegments>
                      </polygon>
                    </g>
                  </style>
                </listOfStyles>
              </renderInformation>
            </listOfRenderInformation>
          </annotation>
          <dimensions width="10" height="10" />
        </layout>
      </listOfLayouts>
    </annotation>
  </model>
</sbml>"##;
    let doc = read(xml).unwrap().document;
    let l = &doc.model.unwrap().layouts[0];
    let st = &l.render.as_ref().unwrap().items[0].styles[0];
    let Primitive::Polygon(p) = &st.group.as_ref().unwrap().children[0] else {
        panic!("no polygon")
    };
    let xy: Vec<_> = p
        .elements
        .iter()
        .map(|e| (e.x.clone().unwrap().0, e.bezier.is_some()))
        .collect();
    assert_eq!(
        xy,
        [
            ("0".to_string(), false),
            ("100%".to_string(), false),
            ("0".to_string(), true)
        ]
    );
    let b = p.elements[2].bezier.as_ref().unwrap();
    assert_eq!((b.x1.clone(), b.y2.clone()), (rel("100%"), rel("100%")));
}

#[test]
fn a_level_2_layout_list_keeps_its_annotation() {
    // Before Render was read, an annotation on listOfLayouts was dropped.
    let xml = r##"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model id="m">
    <annotation>
      <listOfLayouts xmlns="http://projects.eml.org/bcb/sbml/level2">
        <annotation><x:keep xmlns:x="urn:x" a="1"/></annotation>
        <layout id="l">
          <dimensions width="10" height="10"/>
        </layout>
      </listOfLayouts>
    </annotation>
  </model>
</sbml>"##;
    let r = read(xml).unwrap();
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    let m = r.document.model.as_ref().unwrap();
    assert!(
        m.layouts_meta
            .annotation
            .as_deref()
            .is_some_and(|a| a.contains("<x:keep")),
        "{:?}",
        m.layouts_meta
    );
    let text = write(&r.document);
    assert!(text.contains("<x:keep"), "{text}");
    assert_eq!(read(&text).unwrap().document, r.document);
}

#[test]
fn relabs_values_read() {
    let v = |s: &str| RelAbs(s.to_string()).parse();
    assert_eq!(v("10"), Some((10.0, 0.0)));
    assert_eq!(v("50%"), Some((0.0, 50.0)));
    assert_eq!(v("-50%"), Some((0.0, -50.0)));
    assert_eq!(v("5 + 50%"), Some((5.0, 50.0)));
    assert_eq!(v("-2+50%"), Some((-2.0, 50.0)));
    assert_eq!(v("10 - 50%"), Some((10.0, -50.0)));
    assert_eq!(v("5+-50%"), Some((5.0, -50.0)));
    assert_eq!(v("1e-3%"), Some((0.0, 0.001)));
    assert_eq!(v("2e+1 + 1e1%"), Some((20.0, 10.0)));
    assert_eq!(v("wide"), None);
    assert_eq!(v("5 + x%"), None);
    assert_eq!(RelAbs("5 + 50%".into()).resolve(40.0), Some(25.0));
    assert_eq!(RelAbs::rel(50.0).0, "50%");
    assert_eq!(RelAbs::abs(8.0).0, "8");
}

#[test]
fn drawing_written_without_a_prefix_is_read() {
    // As SBMLNetwork's examples write it, and libSBML reads it: the drawing
    // inside a style with no render: prefix, so in the core namespace.
    let xml = r##"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:layout="http://www.sbml.org/sbml/level3/version1/layout/version1" xmlns:render="http://www.sbml.org/sbml/level3/version1/render/version1" level="3" version="1" layout:required="false" render:required="false">
  <model id="m">
    <layout:listOfLayouts>
      <layout:layout layout:id="l">
        <layout:dimensions layout:width="10" layout:height="10"/>
        <render:listOfRenderInformation>
          <render:renderInformation render:id="r">
            <render:listOfStyles>
              <render:style render:id="s" render:idList="sg">
                <g stroke="black" fill="#FF0000">
                  <rectangle x="0" y="0" width="100%" height="100%" rx="6"/>
                </g>
              </render:style>
            </render:listOfStyles>
          </render:renderInformation>
        </render:listOfRenderInformation>
      </layout:layout>
    </layout:listOfLayouts>
  </model>
</sbml>"##;
    let r = read(xml).unwrap();
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    let l = &r.document.model.unwrap().layouts[0];
    let g = l.render.as_ref().unwrap().items[0].styles[0]
        .group
        .as_ref()
        .expect("the group");
    assert_eq!(g.presentation.fill.as_deref(), Some("#FF0000"));
    let Primitive::Rectangle(rect) = &g.children[0] else {
        panic!("no rectangle")
    };
    assert_eq!(rect.rx, rel("6"));
}

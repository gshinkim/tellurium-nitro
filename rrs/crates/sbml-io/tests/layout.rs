//! The Layout package, read from the files in tests/data. Those were written
//! by python-libsbml (tests/scripts/make_layout_fixtures.py): the same model
//! and layout in Level 3 Version 1, and in Level 2 Version 4 where the layout
//! lives in the model's annotation.

use sbml_io::*;

fn load(name: &str) -> Document {
    let path = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
    let read = read(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert!(read.warnings.is_empty(), "{:?}", read.warnings);
    read.document
}

fn layout(doc: &Document) -> &Layout {
    let m = doc.model.as_ref().unwrap();
    assert_eq!(m.layouts.len(), 1);
    &m.layouts[0]
}

#[test]
fn level_3_layout_is_read_in_full() {
    let doc = load("layout-l3v1.xml");
    let m = doc.model.as_ref().unwrap();
    // Nothing of the layout is left as raw XML.
    assert!(m.meta.extra.is_empty(), "{:?}", m.meta.extra);
    let l = layout(&doc);
    assert_eq!(l.id, "layout1");
    assert_eq!(l.dimensions, Dimensions::new(400.0, 300.0));

    assert_eq!(l.compartment_glyphs.len(), 1);
    let cg = &l.compartment_glyphs[0];
    assert_eq!(cg.compartment.as_deref(), Some("cell"));
    let bb = cg.bounding_box.as_ref().unwrap();
    assert_eq!(bb.id.as_deref(), Some("bb_cell"));
    assert_eq!(bb.position, Point::new(5.0, 5.0));
    assert_eq!(bb.dimensions, Dimensions::new(390.0, 290.0));

    let ids: Vec<_> = l.species_glyphs.iter().map(|g| g.id.as_str()).collect();
    assert_eq!(ids, ["sg_Glc", "sg_G6P", "sg_ATP", "sg_ADP", "sg_HK"]);
    assert_eq!(l.text_glyphs.len(), 5);
    let tg = &l.text_glyphs[0];
    assert_eq!(tg.origin_of_text.as_deref(), Some("Glc"));
    assert_eq!(tg.graphical_object.as_deref(), Some("sg_Glc"));

    let rg = &l.reaction_glyphs[0];
    assert_eq!(rg.reaction.as_deref(), Some("J0"));
    assert_eq!(
        rg.curve.as_ref().unwrap().segments,
        [Segment::Line {
            start: Point::new(180.0, 100.0),
            end: Point::new(200.0, 100.0)
        }]
    );
    let roles: Vec<_> = rg
        .species_reference_glyphs
        .iter()
        .map(|s| (s.species_glyph.as_str(), s.role.clone().unwrap()))
        .collect();
    assert_eq!(
        roles,
        [
            ("sg_Glc", Role::Substrate),
            ("sg_ATP", Role::SideSubstrate),
            ("sg_G6P", Role::Product),
            ("sg_ADP", Role::SideProduct),
            ("sg_HK", Role::Activator),
        ]
    );
    let glc = &rg.species_reference_glyphs[0];
    assert_eq!(glc.species_reference.as_deref(), Some("sr_glc"));
    assert_eq!(
        glc.curve.as_ref().unwrap().segments,
        [Segment::Bezier {
            start: Point::new(190.0, 100.0),
            end: Point::new(50.0, 35.0),
            base1: Point::new(190.0, 60.5),
            base2: Point::new(50.0, 60.25),
        }]
    );

    let [
        AdditionalGlyph::Graphical(legend),
        AdditionalGlyph::General(gg),
    ] = &l.additional[..]
    else {
        panic!("{:?}", l.additional);
    };
    assert_eq!(legend.id, "legend");
    assert_eq!(gg.reference.as_deref(), Some("k1"));
    assert_eq!(gg.reference_glyphs[0].glyph, "rg_J0");
    assert_eq!(gg.reference_glyphs[0].role.as_deref(), Some("parameter"));
    let [AdditionalGlyph::Text(sub)] = &gg.sub_glyphs[..] else {
        panic!("{:?}", gg.sub_glyphs);
    };
    assert_eq!(sub.text.as_deref(), Some("k1"));

    assert!(validate(&doc).is_empty(), "{:?}", validate(&doc));
}

#[test]
fn level_2_layout_comes_out_of_the_annotation() {
    let l2 = load("layout-l2v4.xml");
    let l3 = load("layout-l3v1.xml");
    // The same layout, whichever level it was read from.
    assert_eq!(layout(&l2), layout(&l3));
    // The annotation held nothing else.
    assert_eq!(l2.model.as_ref().unwrap().meta.annotation, None);
}

#[test]
fn layouts_are_written_back_at_each_level() {
    for name in ["layout-l3v1.xml", "layout-l2v4.xml"] {
        let doc = load(name);
        let xml = write(&doc);
        let again = read(&xml).unwrap().document;
        assert_eq!(layout(&again), layout(&doc), "{name}");
        assert_eq!(write(&again), xml, "{name}");
        if doc.level == 2 {
            assert!(xml.contains(&format!("<listOfLayouts xmlns=\"{}\"", ns::LAYOUT_L2)));
        } else {
            assert!(xml.contains("<layout:listOfLayouts"));
        }
    }
}

#[test]
fn a_level_2_annotation_keeps_what_isnt_layout() {
    let mut doc = load("layout-l2v4.xml");
    let other =
        "\n      <rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"/>\n    ";
    doc.model.as_mut().unwrap().meta.annotation = Some(other.to_string());
    let again = read(&write(&doc)).unwrap().document;
    let m = again.model.as_ref().unwrap();
    assert_eq!(m.meta.annotation.as_deref(), Some(other));
    assert_eq!(m.layouts.len(), 1);
}

#[test]
fn layout_is_added_to_a_document_built_from_scratch() {
    let mut doc = Document::new(3, 1);
    doc.model = Some(Model {
        compartments: vec![Compartment {
            id: "c".into(),
            ..Default::default()
        }],
        species: vec![Species {
            id: "S".into(),
            compartment: "c".into(),
            ..Default::default()
        }],
        layouts: vec![Layout {
            id: "L".into(),
            dimensions: Dimensions::new(100.0, 50.0),
            species_glyphs: vec![SpeciesGlyph {
                id: "g".into(),
                species: Some("S".into()),
                bounding_box: Some(BoundingBox::new(1.0, 2.0, 30.0, 0.125)),
                ..Default::default()
            }],
            ..Default::default()
        }],
        ..Default::default()
    });
    let xml = write(&doc);
    assert!(
        xml.contains(&format!("xmlns:layout=\"{}\"", ns::LAYOUT)),
        "{xml}"
    );
    assert!(xml.contains("layout:required=\"false\""), "{xml}");
    let again = read(&xml).unwrap().document;
    assert_eq!(again.model, doc.model);
}

#[test]
fn validate_checks_layout_references() {
    let mut doc = load("layout-l3v1.xml");
    let l = &mut doc.model.as_mut().unwrap().layouts[0];
    l.species_glyphs[0].species = Some("nothing".into());
    l.species_glyphs[1].id = "sg_ATP".into();
    l.reaction_glyphs[0].species_reference_glyphs[0].species_glyph = "rg_J0".into();
    let issues: Vec<String> = validate(&doc).into_iter().map(|i| i.message).collect();
    assert!(
        issues
            .iter()
            .any(|i| i.contains("\"nothing\", which isn't in the model")),
        "{issues:#?}"
    );
    assert!(
        issues
            .iter()
            .any(|i| i.contains("glyph id \"sg_ATP\" is used more than once")),
        "{issues:#?}"
    );
    assert!(
        issues
            .iter()
            .any(|i| i.contains("\"rg_J0\", which isn't a species glyph")),
        "{issues:#?}"
    );
}

/// A layout without its size, as MINERVA 22.0.2 and Reactome write, is read
/// with a warning, its size taken from what it draws; one without an id, as
/// Reactome writes, gets one. Written back, both are complete. A glyph's box
/// without its place is still refused, with a reason that names it.
#[test]
fn an_incomplete_layout_is_completed_or_refused_saying_what_is_missing() {
    let path = format!("{}/tests/data/layout-l3v1.xml", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let first = |xml: &str| match read(xml) {
        Ok(_) => panic!("read, though incomplete"),
        Err(errors) => errors[0].message.clone(),
    };
    let complete = read(&text).unwrap().document;
    let start = text.find("<layout:dimensions").unwrap();
    let end = start + text[start..].find("/>").unwrap() + 2;
    let no_size = format!("{}{}", &text[..start], &text[end..]).replacen(
        r#"<layout:layout layout:id="layout1">"#,
        "<layout:layout>",
        1,
    );
    assert!(no_size.contains("<layout:layout>"));
    let r = read(&no_size).unwrap();
    let said: Vec<&str> = r.warnings.iter().map(|w| w.message.as_str()).collect();
    assert!(
        said.iter()
            .any(|m| m.contains("has no <dimensions>") && m.contains("SBML Layout requires")),
        "{said:?}"
    );
    assert!(
        said.iter().any(|m| m.contains("a layout has no id")),
        "{said:?}"
    );
    let layouts = &r.document.model.as_ref().unwrap().layouts;
    let l = &layouts[0];
    assert_eq!(l.id, "layout");
    // The compartment's box, at (5, 5) and 390 by 290, reaches furthest.
    assert_eq!(l.dimensions, Dimensions::new(395.0, 295.0));
    assert_eq!(
        complete.model.as_ref().unwrap().layouts[0].extent(),
        l.dimensions
    );
    // Written back and read again, nothing is missing.
    let again = read(&write(&r.document)).unwrap();
    assert!(again.warnings.is_empty(), "{:?}", again.warnings);
    // The id made up clashes with nothing.
    let taken = no_size.replacen(
        r#"<model id="glycolysis_fragment">"#,
        r#"<model id="layout">"#,
        1,
    );
    let r = read(&taken).unwrap();
    assert_eq!(r.document.model.as_ref().unwrap().layouts[0].id, "layout_2");

    let at = text.find("<layout:position").unwrap();
    let end = at + text[at..].find("/>").unwrap() + 2;
    let no_place = format!("{}{}", &text[..at], &text[end..]);
    let message = first(&no_place);
    assert!(message.contains("the bounding box of '"), "{message}");
}

/// A layout's extent reaches the furthest box or curve, wherever it is.
#[test]
fn the_extent_holds_every_box_and_curve() {
    let line = |x0, y0, x1, y1| Curve {
        segments: vec![Segment::Line {
            start: Point::new(x0, y0),
            end: Point::new(x1, y1),
        }],
    };
    let mut l = Layout::default();
    l.species_glyphs.push(SpeciesGlyph {
        bounding_box: Some(BoundingBox::new(10.0, 10.0, 50.0, 20.0)),
        ..SpeciesGlyph::default()
    });
    assert_eq!(l.extent(), Dimensions::new(60.0, 30.0));
    l.reaction_glyphs.push(ReactionGlyph {
        curve: Some(line(0.0, 0.0, 100.0, 5.0)),
        ..ReactionGlyph::default()
    });
    assert_eq!(l.extent(), Dimensions::new(100.0, 30.0));
    l.reaction_glyphs[0]
        .species_reference_glyphs
        .push(SpeciesReferenceGlyph {
            curve: Some(line(0.0, 0.0, 20.0, 80.0)),
            ..SpeciesReferenceGlyph::default()
        });
    assert_eq!(l.extent(), Dimensions::new(100.0, 80.0));
    l.text_glyphs.push(TextGlyph {
        bounding_box: Some(BoundingBox::new(150.0, 0.0, 10.0, 10.0)),
        ..TextGlyph::default()
    });
    l.additional
        .push(AdditionalGlyph::Graphical(GraphicalObject {
            bounding_box: Some(BoundingBox::new(0.0, 90.0, 10.0, 10.0)),
            ..GraphicalObject::default()
        }));
    assert_eq!(l.extent(), Dimensions::new(160.0, 100.0));
}

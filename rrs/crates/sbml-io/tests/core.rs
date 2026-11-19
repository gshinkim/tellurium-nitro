//! Core SBML: what is kept, MathML, errors, and validation.

use sbml_io::*;

/// A Level 3 Version 2 document with things `sbml-io` doesn't model: an fbc
/// attribute, an fbc element, a namespace declared below `<sbml>`, notes,
/// annotations and MathML outside the subset.
const KEPT: &str = r##"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" level="3" version="2" fbc:required="false">
  <model id="m" fbc:strict="true">
    <notes>
      <body xmlns="http://www.w3.org/1999/xhtml"><p>Some &amp; notes</p></body>
    </notes>
    <listOfCompartments>
      <compartment id="c" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species metaid="m_S" sboTerm="SBO:0000247" id="S" compartment="c" initialAmount="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" fbc:charge="-2" fbc:chemicalFormula="C6H12O6">
        <annotation>
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"><rdf:Description rdf:about="#m_S"/></rdf:RDF>
        </annotation>
      </species>
      <species id="T" compartment="c" xmlns:my="http://example.org/my" my:colour="red" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false">
        <annotation><my:info>kept</my:info></annotation>
      </species>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="1e-20" constant="true"/>
      <parameter id="x" constant="false"/>
    </listOfParameters>
    <listOfRules>
      <assignmentRule variable="x">
        <math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply><diff/><bvar><ci> t </ci></bvar><ci> S </ci></apply>
        </math>
      </assignmentRule>
    </listOfRules>
    <fbc:listOfObjectives fbc:activeObjective="obj">
      <fbc:objective fbc:id="obj" fbc:type="maximize"/>
    </fbc:listOfObjectives>
  </model>
</sbml>
"##;

#[test]
fn what_isnt_understood_is_kept() {
    let doc = read(KEPT).unwrap().document;
    let m = doc.model.as_ref().unwrap();

    let attr = |e: &Extra, name: &str| {
        e.attributes
            .iter()
            .find(|a| a.name == name)
            .map(|a| a.value.clone())
    };
    assert_eq!(
        attr(&doc.meta.extra, "fbc:required").as_deref(),
        Some("false")
    );
    assert_eq!(attr(&m.meta.extra, "fbc:strict").as_deref(), Some("true"));
    assert_eq!(m.meta.extra.elements.len(), 1);
    assert!(m.meta.extra.elements[0].starts_with("<fbc:listOfObjectives"));
    assert!(
        m.meta
            .notes
            .as_ref()
            .unwrap()
            .contains("<p>Some &amp; notes</p>")
    );

    let s = &m.species[0];
    assert_eq!(s.meta.metaid.as_deref(), Some("m_S"));
    assert_eq!(s.meta.sbo_term.as_deref(), Some("SBO:0000247"));
    assert_eq!(attr(&s.meta.extra, "fbc:charge").as_deref(), Some("-2"));
    assert!(
        s.meta
            .annotation
            .as_ref()
            .unwrap()
            .contains("<rdf:Description")
    );

    // A namespace declared on <species> is carried into the raw annotation.
    let t = &m.species[1];
    let my = t.meta.annotation.as_ref().unwrap();
    assert!(
        my.contains(r#"<my:info xmlns:my="http://example.org/my">kept</my:info>"#),
        "{my}"
    );

    assert_eq!(m.parameters[0].value, Some(1e-20));
    let Some(Math::Raw(diff)) = m.rules[0].math() else {
        panic!("{:?}", m.rules[0]);
    };
    assert!(diff.contains("<diff/>"));

    // Written and read again, nothing changes.
    let xml = write(&doc);
    let again = read(&xml).unwrap().document;
    assert_eq!(again, doc, "{xml}");
    // The prefix declared on <species> is declared again for its attribute.
    assert!(
        xml.contains(r#"xmlns:my="http://example.org/my" my:colour="red""#),
        "{xml}"
    );
    assert!(xml.contains(r#"value="1e-20""#));
}

#[test]
fn absent_attributes_stay_absent() {
    let doc = read(KEPT).unwrap().document;
    let x = &doc.model.as_ref().unwrap().parameters[1];
    assert_eq!(x.value, None);
    assert_eq!(x.units, None);
    assert!(write(&doc).contains(r#"<parameter id="x" constant="false"/>"#));
}

fn math_of(body: &str) -> Math {
    let xml = format!(
        r#"<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model><listOfParameters><parameter id="p" constant="false"/></listOfParameters>
    <listOfRules><assignmentRule variable="p">
      <math xmlns="http://www.w3.org/1998/Math/MathML" xmlns:sbml="http://www.sbml.org/sbml/level3/version2/core">{body}</math>
    </assignmentRule></listOfRules></model></sbml>"#
    );
    let doc = read(&xml).unwrap_or_else(|e| panic!("{e:?}")).document;
    let math = doc.model.as_ref().unwrap().rules[0].math().unwrap().clone();
    // It also survives a round trip.
    let again = read(&write(&doc)).unwrap().document;
    assert_eq!(again.model.unwrap().rules[0].math(), Some(&math));
    math
}

#[test]
fn mathml_is_read_as_a_tree() {
    use Math::*;
    let b = Box::new;
    assert_eq!(
        math_of("<apply><plus/><ci> a </ci><cn type='integer'> 2 </cn><cn> 0.5 </cn></apply>"),
        Apply(Op::Plus, vec![ident("a"), int(2), real(0.5)])
    );
    assert_eq!(
        math_of("<cn sbml:units='mole' type='e-notation'> 1.5 <sep/> -3 </cn>"),
        Number {
            value: sbml_io::Number::ENotation(1.5, -3),
            units: Some("mole".into())
        }
    );
    assert_eq!(
        math_of("<cn type='rational'> 1 <sep/> 3 </cn>"),
        Number {
            value: sbml_io::Number::Rational(1, 3),
            units: None
        }
    );
    assert_eq!(
        math_of("<apply><root/><degree><cn type='integer'>3</cn></degree><ci>x</ci></apply>"),
        Root {
            degree: Some(b(int(3))),
            arg: b(ident("x"))
        }
    );
    assert_eq!(
        math_of("<apply><log/><ci>x</ci></apply>"),
        Log {
            base: None,
            arg: b(ident("x"))
        }
    );
    assert_eq!(
        math_of("<apply><ci> f </ci><ci>x</ci><pi/></apply>"),
        Call(
            "f".into(),
            vec![ident("x"), Constant(sbml_io::Constant::Pi)]
        )
    );
    assert_eq!(
        math_of(
            "<apply><csymbol encoding='text' definitionURL='http://www.sbml.org/sbml/symbols/delay'>delay</csymbol>\
             <ci>x</ci><csymbol definitionURL='http://www.sbml.org/sbml/symbols/time'>t</csymbol></apply>"
        ),
        CallSymbol(
            sbml_io::Symbol::new(sbml_io::Symbol::DELAY, "delay"),
            vec![
                ident("x"),
                Symbol(sbml_io::Symbol::new(sbml_io::Symbol::TIME, "t"))
            ]
        )
    );
    assert_eq!(
        math_of(
            "<piecewise><piece><cn>1</cn><apply><lt/><ci>x</ci><cn>0</cn></apply></piece>\
             <otherwise><infinity/></otherwise></piecewise>"
        ),
        Piecewise {
            pieces: vec![(real(1.0), Apply(Op::Lt, vec![ident("x"), real(0.0)]))],
            otherwise: Some(b(Constant(sbml_io::Constant::Infinity)))
        }
    );
    let Semantics { math, annotations } =
        math_of("<semantics><ci>x</ci><annotation encoding='text'>note</annotation></semantics>")
    else {
        panic!()
    };
    assert_eq!(*math, ident("x"));
    // The annotation carries the namespaces it inherited from <math>.
    assert_eq!(annotations.len(), 1);
    assert!(
        annotations[0].starts_with(r#"<annotation xmlns="http://www.w3.org/1998/Math/MathML""#)
    );
    assert!(annotations[0].ends_with("encoding='text'>note</annotation>"));
    // Outside the subset: kept as raw XML.
    assert!(matches!(math_of("<apply><sum/><ci>x</ci></apply>"), Raw(_)));
    assert!(matches!(math_of("<ci type='real'>x</ci>"), Raw(_)));
}

fn ident(s: &str) -> Math {
    Math::ident(s)
}

fn int(i: i64) -> Math {
    Math::int(i)
}

fn real(r: f64) -> Math {
    Math::real(r)
}

#[test]
fn lambdas_and_function_definitions() {
    let xml = r#"<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">
  <model><listOfFunctionDefinitions><functionDefinition id="sq">
    <math xmlns="http://www.w3.org/1998/Math/MathML"><lambda><bvar><ci>x</ci></bvar>
      <apply><times/><ci>x</ci><ci>x</ci></apply></lambda></math>
  </functionDefinition></listOfFunctionDefinitions></model></sbml>"#;
    let doc = read(xml).unwrap().document;
    let f = &doc.model.as_ref().unwrap().function_definitions[0];
    assert_eq!(
        f.math,
        Some(Math::Lambda {
            args: vec!["x".into()],
            body: Box::new(Math::Apply(Op::Times, vec![ident("x"), ident("x")]))
        })
    );
    assert!(validate(&doc).is_empty());
    assert_eq!(read(&write(&doc)).unwrap().document, doc);
}

#[test]
fn errors_say_where() {
    // Not XML.
    let e = read("<sbml level='3'").unwrap_err();
    assert_eq!(e.len(), 1);
    assert_eq!(e[0].line, 1);

    // Level 1 isn't supported.
    let e = read(r#"<sbml xmlns="http://www.sbml.org/sbml/level1" level="1" version="2"/>"#)
        .unwrap_err();
    assert!(e[0].message.contains("Level 1 Version 2"), "{e:?}");

    // A missing id and a bad number, each with its place.
    let xml = "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\">\n\
               <model>\n<listOfSpecies>\n\
               <species compartment=\"c\"/>\n\
               <species id=\"S\" compartment=\"c\" initialAmount=\"lots\"/>\n\
               </listOfSpecies></model></sbml>";
    let e = read(xml).unwrap_err();
    assert_eq!(e.len(), 2, "{e:?}");
    assert_eq!((e[0].line, e[0].column), (4, 1));
    assert!(e[0].message.contains("needs a id"), "{e:?}");
    assert_eq!(e[1].line, 5);
    assert!(e[1].message.contains("initialAmount=\"lots\""), "{e:?}");

    // The namespace must match the level and version.
    let e = read(
        r#"<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="2"/>"#,
    )
    .unwrap_err();
    assert!(e[0].message.contains("namespace"), "{e:?}");
}

#[test]
fn validate_finds_duplicates_and_dangling_references() {
    let xml = r#"<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  <model>
    <listOfCompartments><compartment id="c" constant="true"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="S" compartment="nowhere" substanceUnits="furlongs" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="J" reversible="false">
        <listOfReactants><speciesReference species="X" constant="true"/></listOfReactants>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>S</ci><ci>q</ci><apply><ci>f</ci><ci>S</ci></apply></apply>
          </math>
          <listOfLocalParameters><localParameter id="k" value="1"/></listOfLocalParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"#;
    let doc = read(xml).unwrap().document;
    let issues: Vec<String> = validate(&doc).into_iter().map(|i| i.message).collect();
    let expect = [
        "the id \"S\" is used more than once",
        "species \"S\" refers to \"nowhere\"",
        "the units \"furlongs\" aren't defined",
        "reaction \"J\" refers to \"X\"",
        "the math refers to \"q\"",
        "the math refers to \"f\"",
    ];
    for e in expect {
        assert!(
            issues.iter().any(|i| i.contains(e)),
            "missing {e:?} in {issues:#?}"
        );
    }
    // k is a local parameter and S a species: neither is reported.
    assert!(!issues.iter().any(|i| i.contains("\"k\"")), "{issues:#?}");
    assert_eq!(issues.len(), expect.len(), "{issues:#?}");
}

#[test]
fn numbers_are_written_as_sbml_writes_them() {
    assert_eq!(fmt_f64(1.0), "1");
    assert_eq!(fmt_f64(-3.0), "-3");
    assert_eq!(fmt_f64(0.1), "0.1");
    assert_eq!(fmt_f64(1e-20), "1e-20");
    assert_eq!(fmt_f64(6.02214076e23), "6.02214076e23");
    assert_eq!(fmt_f64(f64::INFINITY), "INF");
    assert_eq!(fmt_f64(f64::NEG_INFINITY), "-INF");
    assert_eq!(fmt_f64(f64::NAN), "NaN");
}

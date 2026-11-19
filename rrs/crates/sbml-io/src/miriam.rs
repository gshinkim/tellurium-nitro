//! MIRIAM annotations: the RDF in an element's `<annotation>` that relates it
//! to database entries (controlled-vocabulary terms, `bqbiol:is` ChEBI 15422)
//! and, on a model, says who made it and when. [`read_annotation`] takes an
//! annotation, as `sbml-io` keeps it (raw XML), apart; [`write_annotation`]
//! writes terms and a history back as RDF, in the form libSBML writes.
//!
//! Anything else in an annotation, such as another tool's own XML, is
//! returned as it is, for the caller to keep or drop.

use crate::read::escape_attr;

/// The RDF namespace.
pub const RDF: &str = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";
/// Biology qualifiers (`bqbiol:is`, `bqbiol:isVersionOf`, …).
pub const BQBIOL: &str = "http://biomodels.net/biology-qualifiers/";
/// Model qualifiers (`bqmodel:is`, `bqmodel:isDescribedBy`, …).
pub const BQMODEL: &str = "http://biomodels.net/model-qualifiers/";
pub const DCTERMS: &str = "http://purl.org/dc/terms/";
pub const DC: &str = "http://purl.org/dc/elements/1.1/";
pub const VCARD: &str = "http://www.w3.org/2001/vcard-rdf/3.0#";
pub const VCARD4: &str = "http://www.w3.org/2006/vcard/ns#";

/// A qualifier and the resources (URIs, as written) it relates the element
/// to: `bqbiol:is` with `http://identifiers.org/chebi/CHEBI:15422`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CvTerm {
    /// A model qualifier (`bqmodel:`), else a biology one (`bqbiol:`).
    pub model: bool,
    /// Its local name: `is`, `isVersionOf`, `hasPart`, `isDescribedBy`, …
    pub qualifier: String,
    pub resources: Vec<String>,
}

/// A model's creator, from a vCard.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Creator {
    pub family: Option<String>,
    pub given: Option<String>,
    pub email: Option<String>,
    pub organization: Option<String>,
}

/// Who made a model and when: its creators, and its created and modified
/// dates (W3CDTF, as written).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct History {
    pub creators: Vec<Creator>,
    pub created: Option<String>,
    pub modified: Vec<String>,
}

impl History {
    pub fn is_empty(&self) -> bool {
        self.creators.is_empty() && self.created.is_none() && self.modified.is_empty()
    }
}

/// What [`read_annotation`] found.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Annotation {
    pub terms: Vec<CvTerm>,
    pub history: History,
    /// Every other element in the annotation: another tool's own, or RDF
    /// this reader doesn't take apart.
    pub other: Vec<Other>,
    /// The `rdf:Description`s in the RDF this reader doesn't take apart
    /// (statements about something else, or with what it doesn't know), each
    /// as raw XML standing on its own, to be written back inside the RDF.
    pub descriptions: Vec<String>,
}

/// One element of an annotation this reader doesn't take apart.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Other {
    /// The element's namespace; `None` when the annotation didn't read as
    /// XML, and `xml` is the whole of it.
    pub namespace: Option<String>,
    /// The element as raw XML, standing on its own: the namespaces it uses
    /// that were declared on the document's `<sbml>` element are declared
    /// on it, so it can be written into another document as it is.
    pub xml: String,
}

impl Other {
    /// RDF, which only [`read_annotation`]'s own reading makes sense of.
    pub fn is_rdf(&self) -> bool {
        self.namespace.as_deref() == Some(RDF)
    }
}

/// Takes an annotation's content apart. `namespaces` are those declared on
/// the document's `<sbml>` element (`Document::namespaces`), which the raw
/// XML may use without declaring. XML that doesn't read is returned whole
/// in `other`.
pub fn read_annotation(xml: &str, namespaces: &[(Option<String>, String)]) -> Annotation {
    let mut decl = String::new();
    for (prefix, uri) in namespaces {
        match prefix {
            Some(p) => decl.push_str(&format!(" xmlns:{p}=\"{}\"", escape_attr(uri))),
            None => decl.push_str(&format!(" xmlns=\"{}\"", escape_attr(uri))),
        }
    }
    let open = format!("<annotation{decl}>");
    let wrapped = format!("{open}{xml}</annotation>");
    let mut out = Annotation::default();
    let Ok(doc) = roxmltree::Document::parse(&wrapped) else {
        if !xml.trim().is_empty() {
            out.other.push(Other {
                namespace: None,
                xml: xml.to_owned(),
            });
        }
        return out;
    };
    for top in doc.root_element().children().filter(|n| n.is_element()) {
        if is(top, RDF, "RDF") && read_rdf(top, &wrapped, &mut out) {
            continue;
        }
        out.other.push(Other {
            namespace: top.tag_name().namespace().map(str::to_owned),
            xml: standing_alone(&wrapped[top.range()], namespaces),
        });
    }
    out
}

/// `element` (raw XML) with declarations added to its start tag for the
/// document's namespaces it uses but doesn't declare itself: a prefix
/// written in it, or the default namespace for an unprefixed tag.
fn standing_alone(element: &str, namespaces: &[(Option<String>, String)]) -> String {
    let name_end = element[1..]
        .find(|c: char| c.is_whitespace() || c == '>' || c == '/')
        .map_or(element.len(), |i| i + 1);
    let tag = &element[1..name_end];
    let mut decl = String::new();
    for (prefix, uri) in namespaces {
        let needed = match prefix {
            Some(p) => {
                let used = [
                    format!("<{p}:"),
                    format!(" {p}:"),
                    format!("\t{p}:"),
                    format!("\n{p}:"),
                ]
                .iter()
                .any(|u| element.contains(u.as_str()));
                used && !element.contains(&format!("xmlns:{p}="))
            }
            None => !tag.contains(':') && !element.contains("xmlns="),
        };
        if needed {
            match prefix {
                Some(p) => decl.push_str(&format!(" xmlns:{p}=\"{}\"", escape_attr(uri))),
                None => decl.push_str(&format!(" xmlns=\"{}\"", escape_attr(uri))),
            }
        }
    }
    format!("{}{decl}{}", &element[..name_end], &element[name_end..])
}

/// Reads an `rdf:RDF` into `out`; false, with `out` unchanged, if it holds
/// anything but descriptions. A description holding anything this reader
/// doesn't know is kept whole, as raw XML (`text` is the document it is in).
fn read_rdf(rdf: roxmltree::Node, text: &str, out: &mut Annotation) -> bool {
    let mut found = Annotation::default();
    for d in rdf.children().filter(|n| n.is_element()) {
        if !is(d, RDF, "Description") {
            return false;
        }
        let known = |c: &roxmltree::Node| {
            let ns = c.tag_name().namespace();
            let name = c.tag_name().name();
            ns == Some(BQBIOL)
                || ns == Some(BQMODEL)
                || ((ns == Some(DCTERMS) || ns == Some(DC)) && name == "creator")
                || (ns == Some(DCTERMS) && (name == "created" || name == "modified"))
        };
        if !d.children().filter(|n| n.is_element()).all(|c| known(&c)) {
            found.descriptions.push(with_scope(d, &text[d.range()]));
            continue;
        }
        for c in d.children().filter(|n| n.is_element()) {
            let ns = c.tag_name().namespace();
            let name = c.tag_name().name();
            if ns == Some(BQBIOL) || ns == Some(BQMODEL) {
                found.terms.push(CvTerm {
                    model: ns == Some(BQMODEL),
                    qualifier: name.to_owned(),
                    resources: resources(c),
                });
            } else if name == "creator" {
                found.history.creators.extend(creators(c));
            } else if name == "created" {
                found.history.created = date(c);
            } else {
                found.history.modified.extend(date(c));
            }
        }
    }
    out.descriptions.extend(found.descriptions);
    out.terms.extend(found.terms);
    out.history.creators.extend(found.history.creators);
    if found.history.created.is_some() {
        out.history.created = found.history.created;
    }
    out.history.modified.extend(found.history.modified);
    true
}

/// An element's raw XML with every namespace in scope that it uses declared
/// on its start tag, so it stands on its own wherever it is written.
fn with_scope(node: roxmltree::Node, element: &str) -> String {
    let name_end = element[1..]
        .find(|c: char| c.is_whitespace() || c == '>' || c == '/')
        .map_or(element.len(), |i| i + 1);
    let mut decl = String::new();
    for ns in node.namespaces() {
        let (used, declared) = match ns.name() {
            Some(p) => (
                element.contains(&format!("{p}:")),
                element[..element.find('>').unwrap_or(element.len())]
                    .contains(&format!("xmlns:{p}=")),
            ),
            None => (true, element.contains("xmlns=")),
        };
        if used && !declared && ns.name() != Some("xml") {
            match ns.name() {
                Some(p) => decl.push_str(&format!(" xmlns:{p}=\"{}\"", escape_attr(ns.uri()))),
                None => decl.push_str(&format!(" xmlns=\"{}\"", escape_attr(ns.uri()))),
            }
        }
    }
    format!("{}{decl}{}", &element[..name_end], &element[name_end..])
}

fn is(n: roxmltree::Node, ns: &str, name: &str) -> bool {
    n.is_element() && n.tag_name().namespace() == Some(ns) && n.tag_name().name() == name
}

/// The `rdf:resource` of every `rdf:li` in a qualifier (in a Bag, Seq or
/// Alt), in order.
fn resources(q: roxmltree::Node) -> Vec<String> {
    q.descendants()
        .filter(|n| is(*n, RDF, "li"))
        .filter_map(|n| n.attribute((RDF, "resource")).map(str::to_owned))
        .collect()
}

/// A date, as W3CDTF: `<dcterms:W3CDTF>2005-02-02T14:56:11Z</dcterms:W3CDTF>`.
fn date(n: roxmltree::Node) -> Option<String> {
    n.descendants()
        .find(|d| is(*d, DCTERMS, "W3CDTF"))
        .or(Some(n))
        .and_then(|d| d.text())
        .map(|t| t.trim().to_owned())
        .filter(|t| !t.is_empty())
}

/// The creators in a `dcterms:creator`, each an `rdf:li` with a vCard (3.0
/// or 4).
fn creators(c: roxmltree::Node) -> Vec<Creator> {
    c.descendants()
        .filter(|n| is(*n, RDF, "li"))
        .map(|li| {
            let text = |names: &[&str]| {
                li.descendants()
                    .find(|d| {
                        d.is_element()
                            && matches!(d.tag_name().namespace(), Some(VCARD) | Some(VCARD4))
                            && names.contains(&d.tag_name().name())
                    })
                    .and_then(|d| d.text())
                    .map(|t| t.trim().to_owned())
                    .filter(|t| !t.is_empty())
            };
            Creator {
                family: text(&["Family", "family-name"]),
                given: text(&["Given", "given-name"]),
                email: text(&["EMAIL", "hasEmail"]),
                organization: text(&["Orgname", "organization-name"]),
            }
        })
        .filter(|c| *c != Creator::default())
        .collect()
}

/// The RDF for `terms` and `history` on the element whose metaid is
/// `metaid`: creators, dates, then terms, as libSBML writes them, then any
/// `descriptions` kept as read; `None` when there is nothing to write.
pub fn write_annotation(
    metaid: &str,
    terms: &[CvTerm],
    history: &History,
    descriptions: &[String],
) -> Option<String> {
    let own = !(terms.iter().all(|t| t.resources.is_empty()) && history.is_empty());
    if !own && descriptions.is_empty() {
        return None;
    }
    if !own {
        let mut out = format!("<rdf:RDF xmlns:rdf=\"{RDF}\">\n");
        for d in descriptions {
            out.push_str(&format!("  {d}\n"));
        }
        out.push_str("</rdf:RDF>");
        return Some(out);
    }
    let esc = escape_text;
    let mut out = format!(
        "<rdf:RDF xmlns:rdf=\"{RDF}\" xmlns:dcterms=\"{DCTERMS}\" xmlns:vCard=\"{VCARD}\" \
         xmlns:bqbiol=\"{BQBIOL}\" xmlns:bqmodel=\"{BQMODEL}\">\n  <rdf:Description \
         rdf:about=\"#{}\">\n",
        escape_attr(metaid)
    );
    if !history.creators.is_empty() {
        out.push_str("    <dcterms:creator>\n      <rdf:Bag>\n");
        for c in &history.creators {
            out.push_str("        <rdf:li rdf:parseType=\"Resource\">\n");
            if c.family.is_some() || c.given.is_some() {
                out.push_str("          <vCard:N rdf:parseType=\"Resource\">\n");
                if let Some(f) = &c.family {
                    out.push_str(&format!(
                        "            <vCard:Family>{}</vCard:Family>\n",
                        esc(f)
                    ));
                }
                if let Some(g) = &c.given {
                    out.push_str(&format!(
                        "            <vCard:Given>{}</vCard:Given>\n",
                        esc(g)
                    ));
                }
                out.push_str("          </vCard:N>\n");
            }
            if let Some(e) = &c.email {
                out.push_str(&format!(
                    "          <vCard:EMAIL>{}</vCard:EMAIL>\n",
                    esc(e)
                ));
            }
            if let Some(o) = &c.organization {
                out.push_str(&format!(
                    "          <vCard:ORG rdf:parseType=\"Resource\">\n            \
                     <vCard:Orgname>{}</vCard:Orgname>\n          </vCard:ORG>\n",
                    esc(o)
                ));
            }
            out.push_str("        </rdf:li>\n");
        }
        out.push_str("      </rdf:Bag>\n    </dcterms:creator>\n");
    }
    let mut date = |tag: &str, d: &str| {
        out.push_str(&format!(
            "    <dcterms:{tag} rdf:parseType=\"Resource\">\n      \
             <dcterms:W3CDTF>{}</dcterms:W3CDTF>\n    </dcterms:{tag}>\n",
            esc(d)
        ));
    };
    if let Some(d) = &history.created {
        date("created", d);
    }
    for d in &history.modified {
        date("modified", d);
    }
    for t in terms.iter().filter(|t| !t.resources.is_empty()) {
        let prefix = if t.model { "bqmodel" } else { "bqbiol" };
        out.push_str(&format!(
            "    <{prefix}:{}>\n      <rdf:Bag>\n",
            t.qualifier
        ));
        for r in &t.resources {
            out.push_str(&format!(
                "        <rdf:li rdf:resource=\"{}\"/>\n",
                escape_attr(r)
            ));
        }
        out.push_str(&format!(
            "      </rdf:Bag>\n    </{prefix}:{}>\n",
            t.qualifier
        ));
    }
    out.push_str("  </rdf:Description>\n");
    for d in descriptions {
        out.push_str(&format!("  {d}\n"));
    }
    out.push_str("</rdf:RDF>");
    Some(out)
}

fn escape_text(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPECIES: &str = r##"
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
        <rdf:Description rdf:about="#S">
          <bqbiol:is>
            <rdf:Bag>
              <rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:17234"/>
              <rdf:li rdf:resource="urn:miriam:kegg.compound:C00031"/>
            </rdf:Bag>
          </bqbiol:is>
          <bqbiol:isPartOf><rdf:Bag><rdf:li rdf:resource="http://identifiers.org/go/GO:0005829"/></rdf:Bag></bqbiol:isPartOf>
        </rdf:Description>
      </rdf:RDF>
      <COPASI xmlns="http://www.copasi.org/static/sbml"><x/></COPASI>"##;

    #[test]
    fn terms_are_read_and_another_tools_xml_returned() {
        let a = read_annotation(SPECIES, &[]);
        assert_eq!(a.terms.len(), 2);
        assert_eq!(a.terms[0].qualifier, "is");
        assert!(!a.terms[0].model);
        assert_eq!(
            a.terms[0].resources,
            [
                "http://identifiers.org/chebi/CHEBI:17234",
                "urn:miriam:kegg.compound:C00031"
            ]
        );
        assert_eq!(a.terms[1].qualifier, "isPartOf");
        assert_eq!(a.other.len(), 1);
        assert!(a.other[0].xml.starts_with("<COPASI"), "{:?}", a.other);
        assert_eq!(
            a.other[0].namespace.as_deref(),
            Some("http://www.copasi.org/static/sbml")
        );
        assert!(!a.other[0].is_rdf());
    }

    /// Another tool's element that uses a prefix declared only on the
    /// document's `<sbml>` gets the declaration, so it stands on its own.
    #[test]
    fn another_tools_xml_carries_the_namespaces_it_uses() {
        let ns = [
            (
                None,
                "http://www.sbml.org/sbml/level3/version1/core".to_owned(),
            ),
            (Some("jd".to_owned()), "urn:jd".to_owned()),
            (Some("x".to_owned()), "urn:x".to_owned()),
        ];
        let a = read_annotation(
            r#"<jd:display x:size="2"><jd:at>1</jd:at></jd:display>"#,
            &ns,
        );
        assert_eq!(a.other.len(), 1);
        let xml = &a.other[0].xml;
        assert!(
            xml.starts_with(r#"<jd:display xmlns:jd="urn:jd" xmlns:x="urn:x" x:size"#),
            "{xml}"
        );
        assert!(roxmltree::Document::parse(xml).is_ok(), "{xml}");
        assert_eq!(a.other[0].namespace.as_deref(), Some("urn:jd"));
        // One that declares its own namespace is left as it is.
        let own = r#"<COPASI xmlns="urn:c"><x/></COPASI>"#;
        assert_eq!(read_annotation(own, &ns).other[0].xml, own);
    }

    #[test]
    fn prefixes_declared_on_the_sbml_element_are_used() {
        let xml = r##"<rdf:RDF><rdf:Description rdf:about="#S"><bqbiol:is><rdf:Bag><rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:1"/></rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF>"##;
        let ns = [
            (Some("rdf".to_owned()), RDF.to_owned()),
            (Some("bqbiol".to_owned()), BQBIOL.to_owned()),
        ];
        assert_eq!(read_annotation(xml, &ns).terms.len(), 1);
        // Without them it doesn't read, and comes back whole.
        let a = read_annotation(xml, &[]);
        assert!(a.terms.is_empty());
        assert_eq!(a.other.len(), 1);
        assert_eq!(a.other[0].xml, xml);
        assert_eq!(a.other[0].namespace, None);
    }

    #[test]
    fn a_model_history_is_read_and_written_back_the_same() {
        let terms = vec![CvTerm {
            model: true,
            qualifier: "isDescribedBy".into(),
            resources: vec!["http://identifiers.org/pubmed/10659856".into()],
        }];
        let history = History {
            creators: vec![Creator {
                family: Some("Elowitz".into()),
                given: Some("Michael".into()),
                email: Some("m@example.org".into()),
                organization: Some("Caltech & co".into()),
            }],
            created: Some("2005-02-02T14:56:11Z".into()),
            modified: vec!["2012-07-05T16:44:43Z".into()],
        };
        let xml = write_annotation("model1", &terms, &history, &[]).unwrap();
        let back = read_annotation(&xml, &[]);
        assert_eq!(back.terms, terms);
        assert_eq!(back.history, history);
        assert!(back.other.is_empty());
        assert_eq!(write_annotation("x", &[], &History::default(), &[]), None);
    }

    /// A description with something this reader doesn't know is kept whole,
    /// with the namespaces it uses from the RDF around it; the others are
    /// read as usual. Written back inside the RDF, it reads again the same
    /// (found by the BioModels round trip: COPASI's statements about
    /// qualifiers lost a model its references and history).
    #[test]
    fn a_description_with_something_unknown_is_kept_whole() {
        let xml = r##"<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:x="urn:x"><rdf:Description rdf:about="#S"><bqbiol:is><rdf:Bag><rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:1"/></rdf:Bag></bqbiol:is></rdf:Description><rdf:Description rdf:about="urn:other"><x:y>1</x:y></rdf:Description></rdf:RDF>"##;
        let a = read_annotation(xml, &[]);
        assert_eq!(a.terms.len(), 1);
        assert!(a.other.is_empty());
        assert_eq!(a.descriptions.len(), 1);
        let d = &a.descriptions[0];
        assert!(
            d.contains(r#"xmlns:x="urn:x""#) && d.contains("xmlns:rdf="),
            "{d}"
        );
        let back = write_annotation("S", &a.terms, &a.history, &a.descriptions).unwrap();
        let again = read_annotation(&back, &[]);
        assert_eq!(again.terms, a.terms);
        assert_eq!(again.descriptions.len(), 1);
        // With nothing else to say, the description alone.
        let alone = write_annotation("S", &[], &History::default(), &a.descriptions).unwrap();
        assert_eq!(read_annotation(&alone, &[]).descriptions.len(), 1);
        // Something other than descriptions in the RDF: returned whole.
        let odd = r##"<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"><rdf:Seq/></rdf:RDF>"##;
        let o = read_annotation(odd, &[]);
        assert_eq!(o.other.len(), 1);
        assert!(o.other[0].is_rdf());
    }
}

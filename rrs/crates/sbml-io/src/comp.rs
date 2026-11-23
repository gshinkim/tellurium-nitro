//! Flatten a hierarchical model (the SBML `comp` package) into plain SBML core,
//! as libsbml's "flatten comp" converter does: each submodel's elements are
//! copied into the containing model with ids prefixed `submodelId__`,
//! replaced and deleted elements are removed, references follow their
//! replacements, and conversion factors rewrite the math.
//!
//! Works on the XML tree, so anything this module doesn't know (notes,
//! annotations, other packages) is copied as it is.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::rc::Rc;

const MATHML: &str = "http://www.w3.org/1998/Math/MathML";
const XML_NS: &str = "http://www.w3.org/XML/1998/namespace";
const SIDREFS: &[&str] = &["compartment", "species", "variable", "symbol", "outside", "conversionFactor", "compartmentType", "speciesType"];
const UNITREFS: &[&str] = &["units", "substanceUnits", "timeUnits", "extentUnits", "volumeUnits", "areaUnits", "lengthUnits", "spatialSizeUnits"];
const ORDER: &[&str] = &["notes", "annotation", "listOfFunctionDefinitions", "listOfUnitDefinitions", "listOfCompartmentTypes", "listOfSpeciesTypes",
    "listOfCompartments", "listOfSpecies", "listOfParameters", "listOfInitialAssignments", "listOfRules", "listOfConstraints", "listOfReactions", "listOfEvents"];

fn is_comp(ns: &Option<String>) -> bool { ns.as_deref().is_some_and(|u| u.contains("/comp/")) }
fn is_core(ns: &Option<String>) -> bool {
    ns.as_deref().is_some_and(|u| u.starts_with("http://www.sbml.org/sbml/level3") && u.ends_with("/core") || u.starts_with("http://www.sbml.org/sbml/level2") || u.starts_with("http://www.sbml.org/sbml/level1"))
}

#[derive(Clone, Debug)]
struct At { ns: Option<String>, pre: Option<String>, name: String, val: String }
#[derive(Clone, Debug)]
enum Node { E(El), T(String) }
#[derive(Clone, Debug)]
struct El { ns: Option<String>, name: String, attrs: Vec<At>, kids: Vec<Node>, dead: bool }

impl El {
    fn new(ns: &str, name: &str, kids: Vec<Node>) -> El { El { ns: Some(ns.into()), name: name.into(), attrs: vec![], kids, dead: false } }
    fn is(&self, name: &str) -> bool { self.name == name }
    /// Core (no-namespace) attribute.
    fn get(&self, name: &str) -> Option<&str> { self.attrs.iter().find(|a| a.ns.is_none() && a.name == name).map(|a| a.val.as_str()) }
    fn comp(&self, name: &str) -> Option<&str> { self.attrs.iter().find(|a| is_comp(&a.ns) && a.name == name).map(|a| a.val.as_str()) }
    fn els(&self) -> impl Iterator<Item = &El> { self.kids.iter().filter_map(|n| if let Node::E(e) = n { (!e.dead).then_some(e) } else { None }) }
    fn els_mut(&mut self) -> impl Iterator<Item = &mut El> { self.kids.iter_mut().filter_map(|n| if let Node::E(e) = n { (!e.dead).then_some(e) } else { None }) }
    fn child(&self, name: &str) -> Option<&El> { self.els().find(|e| e.is(name)) }
    fn list<'a>(&'a self, list: &str, item: &'a str) -> impl Iterator<Item = &'a El> { self.child(list).into_iter().flat_map(move |l| l.els().filter(move |e| e.is(item))) }
    fn text(&self) -> String { self.kids.iter().filter_map(|n| if let Node::T(t) = n { Some(t.as_str()) } else { None }).collect::<String>().trim().to_string() }
}

fn convert(n: roxmltree::Node) -> El {
    let attrs = n.attributes().map(|a| {
        let pre = a.namespace().and_then(|u| if u == XML_NS { Some("xml") } else { n.namespaces().find(|x| x.uri() == u && x.name().is_some()).and_then(|x| x.name()) });
        At { ns: a.namespace().map(Into::into), pre: pre.map(Into::into), name: a.name().into(), val: a.value().into() }
    }).collect();
    let kids = n.children().filter_map(|c| if c.is_element() { Some(Node::E(convert(c))) } else if c.is_text() && !c.text().unwrap_or("").trim().is_empty() { Some(Node::T(c.text().unwrap().into())) } else { None }).collect();
    El { ns: n.tag_name().namespace().map(Into::into), name: n.tag_name().name().into(), attrs, kids, dead: false }
}

fn esc(s: &str, out: &mut String) {
    for c in s.chars() {
        match c { '&' => out.push_str("&amp;"), '<' => out.push_str("&lt;"), '>' => out.push_str("&gt;"), '"' => out.push_str("&quot;"), _ => out.push(c) }
    }
}

/// Writes without comp elements/attributes, dead elements, or emptied `listOf`s. `scope`: (prefix, uri) in force.
fn write(e: &El, scope: &mut Vec<(Option<String>, String)>, out: &mut String) {
    let n = scope.len();
    let bound = |scope: &Vec<(Option<String>, String)>, p: &Option<String>| scope.iter().rev().find(|(q, _)| q == p).map(|(_, u)| u.clone());
    out.push('<');
    out.push_str(&e.name);
    let decl = |p: Option<String>, u: String, scope: &mut Vec<(Option<String>, String)>, out: &mut String| {
        if bound(scope, &p).as_ref() != Some(&u) {
            out.push_str(&match &p { Some(p) => format!(" xmlns:{p}=\""), None => " xmlns=\"".into() });
            esc(&u, out);
            out.push('"');
            scope.push((p, u));
        }
    };
    decl(None, e.ns.clone().unwrap_or_default(), scope, out);
    for a in e.attrs.iter().filter(|a| !is_comp(&a.ns)) {
        let name = match (&a.ns, &a.pre) {
            (Some(u), Some(p)) if u != XML_NS => { decl(Some(p.clone()), u.clone(), scope, out); format!("{p}:{}", a.name) }
            (Some(_), Some(p)) => format!("{p}:{}", a.name),
            _ => a.name.clone(),
        };
        out.push_str(&format!(" {name}=\""));
        esc(&a.val, out);
        out.push('"');
    }
    out.push('>');
    for k in &e.kids {
        match k {
            Node::T(t) => esc(t, out),
            Node::E(c) if c.dead || is_comp(&c.ns) || c.name.starts_with("listOf") && !c.els().any(|x| !x.is("notes") && !x.is("annotation")) => {}
            Node::E(c) => write(c, scope, out),
        }
    }
    out.push_str(&format!("</{}>", e.name));
    scope.truncate(n);
}

/// What an SBaseRef points at, in some model's namespace.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Target { SId(String), Unit(String), Meta(String) }

struct Sub { meta: Option<String>, flat: Flat, map: HashMap<Target, Target> }
/// A flattened model, with what is needed to resolve references into it.
struct Flat { model: El, ports: HashMap<String, El>, subs: HashMap<String, Sub> }

/// A document's model definitions; `base`: directory external sources are relative to.
struct Scope { main: Option<El>, defs: Vec<El>, exts: Vec<El>, base: Option<PathBuf> }

fn scope_of(root: &El, base: Option<PathBuf>) -> Rc<Scope> {
    let comp_list = |l: &str| root.els().filter(|e| is_comp(&e.ns) && e.is(l)).flat_map(|l| l.els().cloned().collect::<Vec<_>>()).collect();
    Rc::new(Scope { main: root.child("model").cloned(), defs: comp_list("listOfModelDefinitions"), exts: comp_list("listOfExternalModelDefinitions"), base })
}

fn parse(xml: &str) -> Result<El, String> {
    let doc = roxmltree::Document::parse(xml).map_err(|e| e.to_string())?;
    Ok(convert(doc.root_element()))
}

fn find(scope: &Rc<Scope>, id: &str, depth: usize) -> Result<(El, Rc<Scope>), String> {
    if depth > 64 { return Err(format!("comp: model definitions nest too deep at '{id}'")); }
    let mid = |e: &El| e.get("id").or(e.comp("id")).map(str::to_string);
    if let Some(d) = scope.defs.iter().find(|d| mid(d).as_deref() == Some(id)) { return Ok((d.clone(), scope.clone())); }
    if let Some(x) = scope.exts.iter().find(|d| mid(d).as_deref() == Some(id)) {
        let src = x.comp("source").ok_or("comp: externalModelDefinition without source")?;
        let src = src.strip_prefix("file://").or(src.strip_prefix("file:")).unwrap_or(src);
        let path = scope.base.as_ref().map_or(PathBuf::from(src), |b| b.join(src));
        let text = std::fs::read_to_string(&path).map_err(|e| format!("comp: {}: {e}", path.display()))?;
        let ext = scope_of(&parse(&text)?, path.parent().map(Path::to_path_buf));
        return match x.comp("modelRef") {
            Some(r) => find(&ext, r, depth + 1),
            None => Ok((ext.main.clone().ok_or("comp: external document without a model")?, ext.clone())),
        };
    }
    match &scope.main {
        Some(m) if mid(m).as_deref() == Some(id) => Ok((m.clone(), scope.clone())),
        _ => Err(format!("comp: no model definition '{id}'")),
    }
}

/// Resolves an SBaseRef (port, deletion, replacement, sBaseRef) in `f`'s namespace.
fn resolve(f: &Flat, r: &El) -> Result<Target, String> {
    let head = if let Some(p) = r.comp("portRef") {
        resolve(f, f.ports.get(p).ok_or_else(|| format!("comp: no port '{p}'"))?)?
    } else if let Some(x) = r.comp("idRef") { Target::SId(x.into()) }
    else if let Some(x) = r.comp("unitRef") { Target::Unit(x.into()) }
    else if let Some(x) = r.comp("metaIdRef") { Target::Meta(x.into()) }
    else { return Err("comp: reference without a target".into()) };
    let Some(c) = r.els().find(|e| is_comp(&e.ns) && e.is("sBaseRef")) else { return Ok(head) };
    let sub = f.subs.iter().find(|(id, s)| match &head { Target::SId(x) => *id == x, Target::Meta(m) => s.meta.as_ref() == Some(m), _ => false })
        .map(|(_, s)| s).ok_or_else(|| format!("comp: {head:?} is not a submodel"))?;
    let t = resolve(&sub.flat, c)?;
    sub.map.get(&t).cloned().ok_or_else(|| format!("comp: {t:?} was deleted"))
}

/// The element `t` names (a global before a local parameter), as a path of child indices.
fn locate(e: &El, t: &Target) -> Option<Vec<usize>> {
    fn go(e: &El, t: &Target, path: &mut Vec<usize>, local: bool, out: &mut Option<(Vec<usize>, bool)>) {
        for (i, n) in e.kids.iter().enumerate() {
            let Node::E(c) = n else { continue };
            if !is_core(&c.ns) { continue; }
            path.push(i);
            let hit = match t {
                Target::SId(x) => !c.is("unitDefinition") && c.get("id") == Some(x),
                Target::Unit(x) => c.is("unitDefinition") && c.get("id") == Some(x),
                Target::Meta(x) => c.get("metaid") == Some(x),
            };
            let local = local || c.is("kineticLaw");
            if hit && out.as_ref().is_none_or(|(_, l)| *l && !local) { *out = Some((path.clone(), local)); }
            go(c, t, path, local, out);
            path.pop();
        }
    }
    let mut out = None;
    go(e, t, &mut vec![], false, &mut out);
    out.map(|(p, _)| p)
}

fn at_mut<'a>(e: &'a mut El, path: &[usize]) -> &'a mut El {
    path.iter().fold(e, |e, &i| match &mut e.kids[i] { Node::E(c) => c, _ => unreachable!() })
}

/// Renames ids, metaids and references (attributes, MathML `ci`, `sbml:units`) by `map`.
/// `local`: ids of the enclosing kinetic law's local parameters, which are not renamed.
fn rename(e: &mut El, map: &HashMap<Target, Target>, local: &[String]) {
    let get = |t: Target| if matches!(&t, Target::SId(x) if local.contains(x)) { None } else { map.get(&t) }.map(|t| match t { Target::SId(x) | Target::Unit(x) | Target::Meta(x) => x.clone() });
    if is_core(&e.ns) {
        let unitdef = e.is("unitDefinition");
        for a in e.attrs.iter_mut().filter(|a| a.ns.is_none()) {
            let v = a.val.clone();
            let new = match a.name.as_str() {
                "id" if unitdef => get(Target::Unit(v)),
                "id" => get(Target::SId(v)),
                "metaid" => get(Target::Meta(v)),
                n if SIDREFS.contains(&n) => get(Target::SId(v)),
                n if UNITREFS.contains(&n) => get(Target::Unit(v)),
                _ => None,
            };
            if let Some(n) = new { a.val = n; }
        }
    } else if e.ns.as_deref() == Some(MATHML) {
        if e.is("ci") {
            if let Some(n) = get(Target::SId(e.text())) { e.kids = vec![Node::T(n)]; }
        }
        for a in e.attrs.iter_mut().filter(|a| is_core(&a.ns) && a.name == "units") {
            if let Some(n) = get(Target::Unit(a.val.clone())) { a.val = n; }
        }
    }
    let mut local = local.to_vec();
    if e.is("kineticLaw") {
        local = e.els().filter(|l| l.name.starts_with("listOf")).flat_map(|l| l.els().filter_map(|p| p.get("id").map(String::from)).collect::<Vec<_>>()).collect();
    }
    for c in e.kids.iter_mut() {
        if let Node::E(c) = c { rename(c, map, &local); }
    }
}

fn ci(x: &str) -> El { El::new(MATHML, "ci", vec![Node::T(x.into())]) }
fn apply(op: &str, a: El, b: El) -> El { El::new(MATHML, "apply", vec![Node::E(El::new(MATHML, op, vec![])), Node::E(a), Node::E(b)]) }

/// Rewrites every MathML expression below `e` (not in function definitions) bottom-up with `f`.
fn map_math(e: &mut El, f: &dyn Fn(El) -> El) {
    if e.is("functionDefinition") { return; }
    for n in e.kids.iter_mut() {
        if let Node::E(c) = n {
            map_math(c, f);
            if c.ns.as_deref() == Some(MATHML) && !c.is("math") {
                let old = std::mem::replace(c, El::new(MATHML, "math", vec![]));
                *c = f(old);
            }
        }
    }
}

/// Applies `f` to the expression of the `<math>` of every element `pick` accepts.
fn map_top(e: &mut El, pick: &dyn Fn(&El) -> bool, f: &dyn Fn(El) -> El) {
    for c in e.els_mut() {
        if pick(c) {
            if let Some(m) = c.els_mut().find(|m| m.is("math") && m.ns.as_deref() == Some(MATHML)) {
                if let Some(Node::E(x)) = m.kids.iter_mut().find(|n| matches!(n, Node::E(_))) {
                    let old = std::mem::replace(x, ci(""));
                    *x = f(old);
                }
            }
        }
        map_top(c, pick, f);
    }
}

fn is_csymbol(e: &El, what: &str) -> bool { e.is("csymbol") && e.get("definitionURL").is_some_and(|u| u.ends_with(what)) }

/// Rewrites a submodel instance for a replacement's conversion factor: `x` -> `x / cf`, values assigned to x times cf.
fn convert_factor(inst: &mut El, x: &str, cf: &str) {
    map_math(inst, &|e| if e.is("ci") && e.text() == x { apply("divide", e, ci(cf)) } else { e });
    map_top(inst, &|e| matches!(e.name.as_str(), "assignmentRule" | "rateRule" | "eventAssignment") && e.get("variable") == Some(x)
        || e.is("initialAssignment") && e.get("symbol") == Some(x), &|m| apply("times", m, ci(cf)));
}

/// Time and extent conversion of a submodel instance (time -> time / tcf, rates / tcf, delays * tcf, kinetic laws * xcf / tcf).
fn convert_time(inst: &mut El, tcf: Option<&str>, xcf: Option<&str>) {
    let klmod = match (xcf, tcf) {
        (Some(x), Some(t)) => apply("divide", ci(x), ci(t)),
        (None, Some(t)) => apply("divide", El::new(MATHML, "cn", vec![Node::T("1".into())]), ci(t)),
        (Some(x), None) => ci(x),
        (None, None) => return,
    };
    // reaction ids in math keep the submodel's units; kinetic laws are converted
    fn reactions(e: &El, out: &mut Vec<String>) {
        for c in e.els() {
            if c.is("reaction") && is_core(&c.ns) { out.extend(c.get("id").map(String::from)); }
            reactions(c, out);
        }
    }
    let mut rids = vec![];
    reactions(inst, &mut rids);
    map_math(inst, &|e| if e.is("ci") && rids.contains(&e.text()) { apply("divide", e, klmod.clone()) } else { e });
    map_top(inst, &|e| e.is("kineticLaw"), &|m| apply("times", klmod.clone(), m));
    if let Some(t) = tcf {
        map_math(inst, &|e| {
            if is_csymbol(&e, "/time") { return apply("divide", e, ci(t)); }
            let mut e = e;
            if e.is("apply") && e.els().next().is_some_and(|h| is_csymbol(h, "/delay")) {
                if let Some(Node::E(d)) = e.kids.iter_mut().filter(|n| matches!(n, Node::E(_))).nth(2) {
                    let old = std::mem::replace(d, ci(""));
                    *d = apply("times", ci(t), old);
                }
            }
            e
        });
        map_top(inst, &|e| e.is("rateRule"), &|m| apply("divide", m, ci(t)));
        map_top(inst, &|e| e.is("delay"), &|m| apply("times", m, ci(t)));
    }
}

fn flatten_model(def: &El, scope: &Rc<Scope>, depth: usize) -> Result<Flat, String> {
    if depth > 64 { return Err("comp: submodels nest too deep".into()); }
    let mut m = def.clone();
    let comp_items = |m: &El, list: &str, item: &str| -> Vec<El> { m.els().filter(|e| is_comp(&e.ns) && e.is(list)).flat_map(|l| l.els().filter(|e| e.is(item)).cloned().collect::<Vec<_>>()).collect() };
    let ports = comp_items(&m, "listOfPorts", "port").into_iter().filter_map(|p| Some((p.comp("id")?.to_string(), p))).collect();
    let mut subs = HashMap::new();
    let mut insts: Vec<(String, El, El)> = vec![]; // (id, submodel element, instance)
    for s in comp_items(&m, "listOfSubmodels", "submodel") {
        let id = s.comp("id").ok_or("comp: submodel without id")?.to_string();
        let (d, sc) = find(scope, s.comp("modelRef").ok_or("comp: submodel without modelRef")?, 0)?;
        let flat = flatten_model(&d, &sc, depth + 1)?;
        let mut inst = flat.model.clone();
        for del in s.list("listOfDeletions", "deletion") {
            let t = resolve(&flat, del)?;
            let path = locate(&inst, &t).ok_or_else(|| format!("comp: nothing to delete at {t:?}"))?;
            at_mut(&mut inst, &path).dead = true;
        }
        subs.insert(id.clone(), Sub { meta: s.get("metaid").map(Into::into), flat, map: HashMap::new() });
        insts.push((id, s, inst));
    }

    // replacements on this model's own elements
    let mut cfs = vec![]; // (submodel, replacement id, conversion factor)
    fn walk(e: &mut El, subs: &mut HashMap<String, Sub>, insts: &mut [(String, El, El)], cfs: &mut Vec<(String, String, String)>) -> Result<(), String> {
        for c in e.els_mut().filter(|c| is_core(&c.ns)) {
            let (pid, pmeta, unit) = (c.get("id").map(String::from), c.get("metaid").map(String::from), c.is("unitDefinition"));
            let reps: Vec<El> = c.els().filter(|x| is_comp(&x.ns) && x.is("listOfReplacedElements")).flat_map(|l| l.els().cloned().collect::<Vec<_>>()).collect();
            let by: Vec<El> = c.els().filter(|x| is_comp(&x.ns) && x.is("replacedBy")).cloned().collect();
            for (r, replaced_by) in reps.iter().map(|r| (r, false)).chain(by.iter().map(|r| (r, true))) {
                if r.comp("deletion").is_some() { continue; } // replaces something already deleted
                let s = r.comp("submodelRef").ok_or("comp: replacement without submodelRef")?;
                let sub = subs.get_mut(s).ok_or_else(|| format!("comp: no submodel '{s}'"))?;
                let t = resolve(&sub.flat, r)?;
                let inst = &mut insts.iter_mut().find(|(id, _, _)| id == s).unwrap().2;
                let path = locate(inst, &t).ok_or_else(|| format!("comp: nothing to replace at {t:?}"))?;
                let q = at_mut(inst, &path);
                if let (Some(a), Some(b)) = (q.get("id"), &pid) {
                    let k = |x: &str| if unit { Target::Unit(x.into()) } else { Target::SId(x.into()) };
                    sub.map.insert(k(a), k(b));
                }
                if let (Some(a), Some(b)) = (q.get("metaid"), &pmeta) { sub.map.insert(Target::Meta(a.into()), Target::Meta(b.clone())); }
                if replaced_by {
                    c.dead = true;
                } else {
                    q.dead = true;
                    if let (Some(f), Some(p)) = (r.comp("conversionFactor"), &pid) { cfs.push((s.to_string(), p.clone(), f.to_string())); }
                }
            }
            walk(c, subs, insts, cfs)?;
        }
        Ok(())
    }
    walk(&mut m, &mut subs, &mut insts, &mut cfs)?;

    for (id, s, mut inst) in insts {
        let sub = subs.get_mut(&id).unwrap();
        fn ids(e: &El, map: &mut HashMap<Target, Target>, pre: &str, local: bool) {
            for c in e.els().filter(|c| is_core(&c.ns)) {
                let local = local || c.is("kineticLaw"); // local parameters keep their ids, as in libsbml
                if let Some(x) = c.get("id").filter(|_| !local) {
                    let (t, n) = if c.is("unitDefinition") { (Target::Unit(x.into()), Target::Unit(format!("{pre}{x}"))) } else { (Target::SId(x.into()), Target::SId(format!("{pre}{x}"))) };
                    map.entry(t).or_insert(n);
                }
                if let Some(x) = c.get("metaid") { map.entry(Target::Meta(x.into())).or_insert(Target::Meta(format!("{pre}{x}"))); }
                ids(c, map, pre, local);
            }
        }
        ids(&inst, &mut sub.map, &format!("{id}__"), false);
        rename(&mut inst, &sub.map, &[]);
        for (_, x, f) in cfs.iter().filter(|(s, _, _)| *s == id) { convert_factor(&mut inst, x, f); }
        convert_time(&mut inst, s.comp("timeConversionFactor"), s.comp("extentConversionFactor"));
        // merge the instance's lists after this model's own
        for l in inst.els().filter(|l| is_core(&l.ns) && l.name.starts_with("listOf")) {
            let items: Vec<Node> = l.els().map(|e| Node::E(e.clone())).collect();
            if let Some(x) = m.els_mut().find(|x| x.is(&l.name) && is_core(&x.ns)) { x.kids.extend(items); continue; }
            m.kids.push(Node::E(El { ns: l.ns.clone(), name: l.name.clone(), attrs: vec![], kids: items, dead: false }));
        }
    }
    m.kids.retain(|n| !matches!(n, Node::E(e) if is_comp(&e.ns)));
    let rank = |n: &Node| match n { Node::E(e) if is_core(&e.ns) => ORDER.iter().position(|o| e.is(o)).unwrap_or(ORDER.len()), _ => ORDER.len() };
    m.kids.sort_by_key(rank);
    Ok(Flat { model: m, ports, subs })
}

/// Flattens an SBML document that uses the `comp` package into SBML core text
/// (the comp namespace removed). `base_dir`: the directory external model
/// definitions are relative to (the current directory if `None`).
pub fn flatten(xml: &str, base_dir: Option<&Path>) -> Result<String, String> {
    let mut root = parse(xml)?;
    let scope = scope_of(&root, base_dir.map(Path::to_path_buf));
    let main = scope.main.as_ref().ok_or("comp: document without a model")?;
    let flat = flatten_model(main, &scope, 0)?;
    for n in root.kids.iter_mut() {
        if let Node::E(e) = n { if e.is("model") && is_core(&e.ns) { *e = flat.model.clone(); } }
    }
    let doc = roxmltree::Document::parse(xml).map_err(|e| e.to_string())?;
    let mut scope: Vec<(Option<String>, String)> = doc.root_element().namespaces().filter(|n| n.name().is_some() && !n.uri().contains("/comp/"))
        .map(|n| (n.name().map(Into::into), n.uri().into())).collect();
    let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    // declare the root's prefixes up front so annotations that use them in values stay valid
    let decls: String = scope.iter().map(|(p, u)| format!(" xmlns:{}=\"{}\"", p.as_deref().unwrap_or(""), u)).collect();
    write(&root, &mut scope, &mut out);
    let at = out.find("<sbml").map(|i| i + 5).unwrap_or(0);
    out.insert_str(at, &decls);
    Ok(out)
}

#[cfg(test)]
mod tests {
    #[test]
    fn submodel_with_replacement() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:comp="http://www.sbml.org/sbml/level3/version1/comp/version1" level="3" version="1" comp:required="true">
  <model id="top">
    <listOfParameters>
      <parameter id="k" value="2" constant="true">
        <comp:listOfReplacedElements><comp:replacedElement comp:submodelRef="A" comp:portRef="k_port"/></comp:listOfReplacedElements>
      </parameter>
    </listOfParameters>
    <comp:listOfSubmodels><comp:submodel comp:id="A" comp:modelRef="inner"/></comp:listOfSubmodels>
  </model>
  <comp:listOfModelDefinitions>
    <comp:modelDefinition id="inner">
      <listOfParameters>
        <parameter id="k" value="1" constant="true"/>
        <parameter id="x" value="0" constant="false"/>
      </listOfParameters>
      <listOfRules><assignmentRule variable="x"><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci> k </ci><csymbol encoding="text" definitionURL="http://www.sbml.org/sbml/symbols/time">t</csymbol></apply></math></assignmentRule></listOfRules>
      <comp:listOfPorts><comp:port comp:id="k_port" comp:idRef="k"/></comp:listOfPorts>
    </comp:modelDefinition>
  </comp:listOfModelDefinitions>
</sbml>"#;
        let out = super::flatten(xml, None).unwrap();
        assert!(!out.contains("comp"), "{out}");
        assert!(out.contains(r#"<parameter id="k" value="2""#), "{out}");
        assert!(out.contains(r#"<parameter id="A__x""#), "{out}");
        assert!(!out.contains(r#"value="1""#), "{out}");
        assert!(out.contains(r#"<assignmentRule variable="A__x">"#) && out.contains("<ci>k</ci>"), "{out}");
        crate::read(&out).unwrap();
    }
}

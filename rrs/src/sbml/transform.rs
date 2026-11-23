//! Model rewrites before compiling: fast reactions -> algebraic rules, algebraic rules -> assignment
//! rules (or DAE states), parameters read inside delay() -> states.
use super::math::*;
use super::{nimp, Error, Ev, Par, Rule, RK, M};
use crate::model::Rng;
use std::collections::{HashMap, HashSet};

fn unif(r: &mut Rng) -> f64 { 0.5 + 1.5 * r.unif() }

/// Id of a new non-constant parameter, named stem (plus underscores if taken).
fn fresh_param(m: &mut M, stem: &str) -> String {
    let taken = m.sids();
    let mut i = stem.to_string();
    while taken.contains(&i) { i.push('_'); }
    m.params.push(Par { id: i.clone(), value: None, constant: false });
    i
}

fn rule(kind: RK, var: &str, math: E) -> Rule { Rule { kind, var: var.to_string(), math } }

/// Rank of the rows x cols matrix a.
// ponytail: f64 elimination with a relative tolerance (Python used exact fractions); fine for stoichiometries
fn rank(mut a: Vec<Vec<f64>>) -> usize {
    let (rows, cols) = (a.len(), a.first().map_or(0, |r| r.len()));
    let tol = 1e-12 * a.iter().flatten().fold(1.0f64, |m, x| m.max(x.abs()));
    let mut r = 0;
    for c in 0..cols {
        let Some(p) = (r..rows).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs())).filter(|&p| a[p][c].abs() > tol) else { continue };
        a.swap(r, p);
        for i in 0..rows {
            if i != r {
                let f = a[i][c] / a[r][c];
                for k in 0..cols { a[i][k] -= f * a[r][k]; }
            }
        }
        r += 1;
    }
    r
}

/// SBML fast reactions are at equilibrium at every instant. Each fast species' amount (over its
/// conversion factor) is base + sum_j n_j * extent_j: base a rate-rule parameter the slow reactions
/// drive, extent_j fast reaction j's extent, determined by an algebraic rule v_j = 0. The fast species
/// become boundary species set by assignment rules. Returns the extents; `pre` gets the fast species'
/// values before they equilibrate (for delay() before the run).
pub fn fast_to_algebraic(m: &mut M, pre: &mut Vec<(String, E)>) -> Result<HashSet<String>, Error> {
    if !m.reactions.iter().any(|r| r.fast) { return Ok(HashSet::new()); }
    // a fast reaction whose rate is 0 for any values (e.g. k * A * ln(1)) never moves anything: drop it
    let fd = m.fdefs();
    let mut rng = Rng::new(2);
    let mut rnd: HashMap<String, f64> = HashMap::new();
    let mut dead = HashSet::new();
    for rx in m.reactions.iter().filter(|r| r.fast) {
        let Some(Some(k)) = rx.kl.as_ref().map(|k| &k.0) else { continue };
        let mut vals = vec![];
        for _ in 0..2 {
            match eval(k, &mut |x| Ok(*rnd.entry(x.to_string()).or_insert_with(|| unif(&mut rng))), &fd) {
                Ok(v) => vals.push(v),
                Err(Error::NotImplemented(_)) => break,
                Err(e) => return Err(e),
            }
        }
        if vals.len() == 2 && vals.iter().all(|&v| v == 0.0) { dead.insert(rx.id.clone()); }
    }
    m.reactions.retain(|r| !dead.contains(&r.id));
    let fast: Vec<usize> = (0..m.reactions.len()).filter(|&j| m.reactions[j].fast).collect();
    if fast.is_empty() { return Ok(HashSet::new()); }
    let spi: HashMap<String, usize> = m.species.iter().enumerate().map(|(i, s)| (s.id.clone(), i)).collect();
    let sp = |x: &str| &m.species[spi[x]];
    let mut fs: Vec<String> = vec![];
    for &j in &fast {
        for (_, sr) in m.reactions[j].refs() {
            let s = sp(&sr.species);
            if !s.boundary && !s.constant && !fs.contains(&sr.species) { fs.push(sr.species.clone()); }
        }
    }
    let ias: HashMap<String, E> = m.ias.iter().cloned().collect();
    let asg: HashMap<String, E> = m.rules.iter().filter(|r| r.kind == RK::Assign).map(|r| (r.var.clone(), r.math.clone())).collect();
    let targets: HashSet<String> = m.events.iter().flat_map(|e| e.asg.iter().map(|a| a.0.clone())).collect();
    let moving: HashSet<String> = m.rules.iter().map(|r| r.var.clone()).chain(targets.iter().cloned()).collect();

    // initial value of a symbol a stoichiometry uses; it may jump (events) but not drift
    fn val0(m: &M, fd: &Fdefs, asg: &HashMap<String, E>, ias: &HashMap<String, E>, moving: &HashSet<String>,
            targets: &HashSet<String>, x: &str) -> Result<f64, Error> {
        if moving.contains(x) && !asg.contains_key(x) && !targets.contains(x) {
            return Err(nimp(format!("fast reaction with stoichiometry {x} that changes continuously")));
        }
        let mut g = |y: &str| val0(m, fd, asg, ias, moving, targets, y);
        if let Some(e) = asg.get(x).or_else(|| ias.get(x)) { return eval(e, &mut g, fd); }
        if let Some(sr) = m.reactions.iter().flat_map(|r| r.refs()).map(|p| p.1).find(|sr| sr.id.as_deref() == Some(x)) {
            return Ok(sr.stoich.unwrap_or(1.0));
        }
        if let Some(v) = m.params.iter().find(|p| p.id == x).and_then(|p| p.value) { return Ok(v); }
        if let Some(v) = m.comps.iter().find(|c| c.id == x).and_then(|c| c.size) { return Ok(v); }
        Err(nimp(format!("fast reaction: no constant value for {x}")))
    }
    let v0 = |x: &str| val0(m, &fd, &asg, &ias, &moving, &targets, x);
    // net stoichiometry of species x in reaction j, as a formula (None if none)
    let ns = |j: usize, x: &str| -> Option<E> {
        let t: Vec<E> = m.reactions[j].refs().filter(|(_, sr)| sr.species == x).map(|(sg, sr)| {
            let st = if let Some(e) = &sr.smath { e.clone() }
                else if sr.id.as_ref().is_some_and(|i| moving.contains(i) || ias.contains_key(i)) { name(sr.id.as_ref().unwrap()) }
                else { num(sr.stoich.unwrap_or(1.0)) };
            mul(num(sg as f64), st)
        }).collect();
        if t.is_empty() { None } else { Some(sum(t)) }
    };
    let n = |j: usize, x: &str| -> Result<f64, Error> {
        match ns(j, x) { Some(e) => eval(&e, &mut |y| v0(y), &fd), None => Ok(0.0) }
    };
    let cf = |x: &str| sp(x).cf.clone().or_else(|| m.cf.clone());
    for x in &fs {
        if let Some(c) = cf(x).filter(|c| moving.contains(c)) {
            return Err(nimp(format!("fast reaction with conversion factor {c} that changes")));
        }
    }
    let nmat = fs.iter().map(|x| fast.iter().map(|&j| n(j, x)).collect::<Result<Vec<_>, _>>()).collect::<Result<Vec<_>, _>>()?;
    if rank(nmat) < fast.len() { return Err(nimp("dependent fast reactions".into())); }
    let over_cf = |x: &str, v: E| match cf(x) { Some(c) => div(v, name(&c)), None => v };
    // amounts over the conversion factor: reactions change these by n * rate
    let amt: HashMap<String, E> = fs.iter().map(|x| {
        let s = sp(x);
        (x.clone(), over_cf(x, if s.hosu { name(x) } else { mul(name(x), name(&s.comp)) }))
    }).collect();
    let amt0 = |x: &str| -> E { // x's initial amount (over its conversion factor), from its own initial value
        let s = sp(x);
        let c = name(&s.comp);
        let v = if let Some(e) = ias.get(x) { if s.hosu { e.clone() } else { mul(e.clone(), c) } }
            else if let Some(a) = s.amt { num(a) } else { mul(num(s.conc.unwrap_or(0.0)), c) };
        over_cf(x, v)
    };
    let conc_from = |x: &str, u: E| -> E { // species value (SBML units) from an amount over the conversion factor
        let a = match cf(x) { Some(c) => mul(u, name(&c)), None => u };
        if sp(x).hosu { a } else { div(a, name(&sp(x).comp)) }
    };
    let kinetic = |j: usize| -> Result<E, Error> { // kinetic law with local parameters inlined
        let rx = &m.reactions[j];
        let Some((Some(k), loc)) = &rx.kl else { return Err(nimp(format!("fast reaction {} without kinetic law", rx.id))) };
        Ok(loc.iter().fold(k.clone(), |a, (p, v)| a.sub1(p, &num(*v))))
    };
    // everything that needs `m` read-only, before mutating it
    let fids: Vec<String> = fast.iter().map(|&j| m.reactions[j].id.clone()).collect();
    let slow: Vec<usize> = (0..m.reactions.len()).filter(|&j| !m.reactions[j].fast).collect();
    let nsf: Vec<Vec<Option<E>>> = fast.iter().map(|&j| fs.iter().map(|x| ns(j, x)).collect()).collect();
    let nss: Vec<Vec<Option<E>>> = slow.iter().map(|&j| fs.iter().map(|x| ns(j, x)).collect()).collect();
    let slow_ids: Vec<String> = slow.iter().map(|&j| m.reactions[j].id.clone()).collect();
    let amt0s: Vec<E> = fs.iter().map(|x| amt0(x)).collect();
    let conc0: Vec<E> = fs.iter().zip(&amt0s).map(|(x, a)| conc_from(x, a.clone())).collect();
    let kins = fast.iter().map(|&j| kinetic(j)).collect::<Result<Vec<_>, _>>()?;
    let nvals = fast.iter().map(|&j| fs.iter().map(|x| n(j, x)).collect::<Result<Vec<_>, _>>()).collect::<Result<Vec<_>, _>>()?;
    let mut syms: HashSet<String> = fs.iter().cloned().collect();
    for row in &nsf { for e in row.iter().flatten() { e.names_into(&mut syms); } }
    let named: Vec<(String, f64)> = fast.iter().flat_map(|&j| m.reactions[j].refs()).filter_map(|(_, sr)| sr.id.clone().map(|i| (i, sr.stoich.unwrap_or(1.0)))).collect();
    // conc_from needs only species attributes, which do not change below
    let cfx: Vec<(Option<String>, bool, String)> = fs.iter().map(|x| (cf(x), sp(x).hosu, sp(x).comp.clone())).collect();
    let conc_of = |i: usize, u: E| -> E {
        let (c, h, comp) = &cfx[i];
        let a = match c { Some(c) => mul(u, name(c)), None => u };
        if *h { a } else { div(a, name(comp)) }
    };

    let xi: Vec<String> = fids.iter().map(|r| fresh_param(m, &format!("__fast_extent_{r}"))).collect();
    let base: Vec<String> = fs.iter().map(|x| fresh_param(m, &format!("__fast_base_{x}"))).collect();
    let mut now = vec![]; // fast species -> its amount (over the conversion factor): base + extents
    for (i, _) in fs.iter().enumerate() {
        m.ias.push((base[i].clone(), amt0s[i].clone()));
        let r: Vec<E> = slow_ids.iter().enumerate().filter_map(|(q, rid)| nss[q][i].clone().map(|e| mul(e, name(rid)))).collect();
        m.rules.push(rule(RK::Rate, &base[i], if r.is_empty() { num(0.0) } else { sum(r) }));
        let mut t = vec![name(&base[i])];
        t.extend((0..fids.len()).filter_map(|q| nsf[q][i].clone().map(|e| mul(e, name(&xi[q])))));
        now.push(sum(t));
    }
    for k in &kins {
        let mut e = k.clone();
        for (i, x) in fs.iter().enumerate() { e = e.sub1(x, &conc_of(i, now[i].clone())); }
        m.rules.push(Rule { kind: RK::Alg, var: String::new(), math: e });
    }
    // starting guess where an equilibrium stays a DAE constraint (e.g. k*A*B = 0 has two roots): every
    // fast reaction run to completion in the direction of its initial rate
    pre.extend(fs.iter().cloned().zip(conc0.iter().cloned())); // delay() before the run sees these
    for (q, k) in kins.iter().enumerate() {
        let mut v0 = k.clone();
        for (i, x) in fs.iter().enumerate() { v0 = v0.sub1(x, &conc0[i]); }
        let lim = |sg: f64| -> Vec<E> {
            (0..fs.len()).filter(|&i| sg * nvals[q][i] < 0.0).map(|i| div(amt0s[i].clone(), num(nvals[q][i].abs()))).collect()
        };
        let (fwd, bwd) = (lim(1.0), lim(-1.0));
        let fw = if fwd.is_empty() { num(0.0) } else { op("min", fwd) };
        let bw = if bwd.is_empty() { num(0.0) } else { neg(op("min", bwd)) };
        m.ias.push((xi[q].clone(), op("piecewise", vec![fw, op("gt", vec![v0.clone(), num(0.0)]), bw, op("lt", vec![v0, num(0.0)]), num(0.0)])));
    }
    // an event that assigns a fast species or a stoichiometry: the amounts jump to their assigned values
    // (base = amount - N_new * extent), then the fast reactions re-equilibrate from there
    for ev in m.events.iter_mut() {
        let eas: HashMap<String, E> = ev.asg.iter().filter(|a| syms.contains(&a.0)).cloned().collect();
        if eas.is_empty() { continue; }
        for (i, x) in fs.iter().enumerate() {
            let mut e = amt[x].subst(&eas);
            for q in 0..fids.len() {
                if let Some(nq) = &nsf[q][i] { e = sub(e, mul(nq.subst(&eas), name(&xi[q]))); }
            }
            ev.asg.push((base[i].clone(), e));
        }
        ev.asg.retain(|a| !fs.contains(&a.0));
    }
    for (i, x) in fs.iter().enumerate() {
        m.species[spi[x]].boundary = true;
        m.ias.retain(|a| &a.0 != x); // now set by a rule
        m.rules.push(rule(RK::Assign, x, conc_of(i, now[i].clone())));
    }
    m.reactions.retain(|r| !r.fast);
    for (sid, v) in named { // a named stoichiometry other math may use: kept as a parameter
        let constant = !moving.contains(&sid);
        m.params.push(Par { id: sid, value: Some(v), constant });
    }
    Ok(xi.into_iter().collect())
}

/// Symbols an algebraic rule may determine: not constant, not set by another rule, not a floating
/// species that reactions change.
fn free_symbols(m: &M) -> HashSet<String> {
    let taken: HashSet<&str> = m.rules.iter().filter(|r| r.kind != RK::Alg).map(|r| r.var.as_str())
        .chain(m.events.iter().flat_map(|e| e.asg.iter().map(|a| a.0.as_str()))).collect();
    let moved: HashSet<&str> = m.reactions.iter().flat_map(|r| r.refs()).map(|(_, sr)| sr.species.as_str()).collect();
    let mut free: HashSet<String> = HashSet::new();
    free.extend(m.params.iter().filter(|p| !p.constant).map(|p| p.id.clone()));
    free.extend(m.comps.iter().filter(|c| !c.constant).map(|c| c.id.clone()));
    free.extend(m.species.iter().filter(|s| !s.constant && (s.boundary || !moved.contains(s.id.as_str()))).map(|s| s.id.clone()));
    free.extend(m.reactions.iter().flat_map(|r| r.refs()).filter(|(_, sr)| !sr.constant).filter_map(|(_, sr)| sr.id.clone()));
    free.retain(|x| !taken.contains(x.as_str()));
    free
}

/// Fast-reaction extent z = f solved from v = 0, where an event can switch the reaction off (dv/dz = d
/// becomes 0, e.g. v = k * A with k := 0): then v = 0 holds for any z and the extent stays where it was.
/// A state h tracks f while the reaction runs (rate rateOf(f), reset to f when it switches on) and holds
/// still while it is off; z = (d == 0 ? h : f).
fn held(m: &mut M, z: &str, f: E, d: E) -> E {
    let s = fresh_param(m, &format!("__fast_solved_{z}"));
    let h = fresh_param(m, &format!("__fast_held_{z}"));
    let off = || op("eq", vec![d.clone(), num(0.0)]);
    m.rules.push(rule(RK::Assign, &s, f));
    m.ias.push((h.clone(), op("piecewise", vec![num(0.0), off(), name(&s)]))); // off at the start: no extent yet
    m.rules.push(rule(RK::Rate, &h, op("piecewise", vec![num(0.0), off(), op("rateOf", vec![name(&s)])])));
    m.events.push(Ev {
        id: None, uvftt: true, trigger: op("neq", vec![d.clone(), num(0.0)]), init: true, persistent: true,
        delay: None, priority: None, asg: vec![(h.clone(), name(&s))],
    });
    op("piecewise", vec![name(&h), off(), name(&s)])
}

fn isclose(a: f64, b: f64) -> bool { (a - b).abs() <= 1e-12 + 1e-9 * b.abs() }

/// Replace algebraic rules by assignment rules: match each rule to a symbol it determines, then
/// eliminate in turn. A rule affine in its symbol (checked numerically at random values) is solved
/// for it; any other stays a constraint: its symbol becomes an algebraic state (a rate rule holding
/// the residual, returned) and the model a DAE.
pub fn algebraic_to_assignment(m: &mut M, prefer: &HashSet<String>) -> Result<HashSet<String>, Error> {
    let alg: Vec<E> = m.rules.iter().filter(|r| r.kind == RK::Alg).map(|r| r.math.clone()).collect();
    if alg.is_empty() { return Ok(HashSet::new()); }
    let fd = m.fdefs();
    let free = free_symbols(m);
    let cand: Vec<Vec<String>> = alg.iter().map(|g| {
        let mut c: Vec<String> = g.names().into_iter().filter(|x| free.contains(x)).collect();
        c.sort_by_key(|x| (!prefer.contains(x), x.clone()));
        c
    }).collect();
    fn assign(i: usize, seen: &mut HashSet<String>, cand: &[Vec<String>], mat: &mut HashMap<String, usize>) -> bool {
        for x in &cand[i] {
            if seen.insert(x.clone()) && (!mat.contains_key(x) || assign(mat[x], seen, cand, mat)) {
                mat.insert(x.clone(), i);
                return true;
            }
        }
        false
    }
    let mut mat = HashMap::new(); // symbol -> rule index (maximum bipartite matching, augmenting paths)
    for i in 0..alg.len() {
        if !assign(i, &mut HashSet::new(), &cand, &mut mat) {
            return Err(nimp(format!("algebraic rule {} determines no free symbol", s(&alg[i]))));
        }
    }
    let var: HashMap<usize, String> = mat.into_iter().map(|(x, i)| (i, x)).collect();
    let targets: HashSet<String> = m.events.iter().flat_map(|e| e.asg.iter().map(|a| a.0.clone())).collect();
    let mut rng = Rng::new(1);
    let mut vals: HashMap<String, f64> = HashMap::new();
    let mut g = alg;
    let (mut sol, mut dae) = (vec![], vec![]);
    for i in 0..g.len() {
        let z = var[&i].clone();
        let (g0, g1) = (g[i].sub1(&z, &num(0.0)), g[i].sub1(&z, &num(1.0)));
        let d = sub(g1.clone(), g0.clone());
        let mut f = neg(div(g0.clone(), d.clone()));
        let mut linear = true;
        for zt in [2.0, -1.3, 3.7] { // affine in z?
            let mut val = |n: &str| Ok(*vals.entry(n.to_string()).or_insert_with(|| unif(&mut rng)));
            let ex = eval(&g[i], &mut |n| if n == z { Ok(zt) } else { val(n) }, &fd)?;
            let lin = eval(&add(g0.clone(), mul(num(zt), d.clone())), &mut val, &fd)?;
            linear = linear && isclose(ex, lin);
        }
        if !linear { // a true constraint: z stays a state of the DAE
            dae.push((z, g[i].clone()));
            continue;
        }
        if prefer.contains(&z) && d.names().iter().any(|x| targets.contains(x)) { f = held(m, &z, f, d); }
        for k in i + 1..g.len() { g[k] = g[k].sub1(&z, &f); } // eliminate z from the rules still to solve
        sol.push((z, f));
    }
    m.rules.retain(|r| r.kind != RK::Alg);
    for (z, f) in sol { m.rules.push(rule(RK::Assign, &z, f)); }
    for (z, f) in &dae {
        m.rules.push(rule(RK::Rate, z, f.clone())); // its "rate" is the residual; the engine holds it at 0
        if let Some(sp) = m.species.iter_mut().find(|s| &s.id == z) { sp.boundary = true; }
    }
    Ok(dae.into_iter().map(|d| d.0).collect())
}

/// delay() reads the run's state history: a parameter (or compartment, boundary species) that only
/// events change, used inside delay(), becomes a state with rate 0 so its history is kept.
pub fn delayed_params_to_states(m: &mut M) {
    fn walk(a: &E, inside: bool, found: &mut HashSet<String>) {
        if inside { if let E::Name(n) = a { found.insert(n.clone()); } }
        let d = matches!(a, E::Op(f, _) if f == "delay");
        for (i, k) in a.kids().iter().enumerate() { walk(k, inside || (d && i == 0), found); }
    }
    let mut found = HashSet::new();
    let mut all: Vec<&E> = m.rules.iter().map(|r| &r.math).chain(m.ias.iter().map(|a| &a.1)).chain(&m.constraints).collect();
    all.extend(m.reactions.iter().filter_map(|r| r.kl.as_ref().and_then(|k| k.0.as_ref())));
    for ev in &m.events {
        all.push(&ev.trigger);
        all.extend(ev.delay.iter().chain(&ev.priority).chain(ev.asg.iter().map(|a| &a.1)));
    }
    for e in all { walk(e, false, &mut found); }
    let ruled: HashSet<&str> = m.rules.iter().map(|r| r.var.as_str()).collect();
    let targets: HashSet<&str> = m.events.iter().flat_map(|e| e.asg.iter().map(|a| a.0.as_str())).collect();
    let mut add: Vec<String> = found.into_iter().filter(|x| targets.contains(x.as_str()) && !ruled.contains(x.as_str())
        && (m.params.iter().any(|p| &p.id == x) || m.comps.iter().any(|c| &c.id == x)
            || m.species.iter().any(|s| &s.id == x && s.boundary))).collect();
    add.sort();
    for x in add { m.rules.push(rule(RK::Rate, &x, num(0.0))); }
}

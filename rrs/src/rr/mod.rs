//! `RoadRunner`: libroadrunner's instance semantics on the Rust engine, for the C API and the Python
//! package. A port of tellurium_nitro's lazy RoadRunner (extended_roadrunner.py: `_lz_*`, `_RsModel`,
//! `_Cols`, `_rs_simulate`, `_rs_ssa`, `steadyState`, `_rs_mca`, `_rs_jac`, `getuCC`, ...).
//!
//! State: the values of `keys` = [global parameters | compartments | boundary species | extra states |
//! floating species] the engine holds (species as concentrations; boundary species whose engine value
//! is an amount stay amounts), plus the time. Everything else (rules, rates) is computed from them.
mod cons;
mod eig;
mod eig_aed;
mod ss;
#[cfg(test)]
mod regen_test;
pub use eig::eigenvalues;
mod sel;
mod settings;

pub use cons::Ls;
pub use sel::{parse as parse_selection, Sel};
pub use settings::{Solver, Value, INTEGRATORS, STEADY_STATE_SOLVERS};

use crate::model::{Ctl, Key, Model, Rows, Work};
use crate::sbml::Info;
use crate::sens::RateDerivs;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Row-major matrix with optional row and column labels.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Matrix { pub nr: usize, pub nc: usize, pub data: Vec<f64>, pub rows: Vec<String>, pub cols: Vec<String> }

impl Matrix {
    pub fn get(&self, i: usize, j: usize) -> f64 { self.data[i * self.nc + j] }
    fn new(rows: Vec<String>, cols: Vec<String>, data: Vec<f64>) -> Matrix {
        Matrix { nr: rows.len(), nc: cols.len(), data, rows, cols }
    }
}

/// One output column of a time course, from the engine's output registers `outs`.
#[derive(Debug, Clone)]
enum Col { Out(usize), Mul(usize, usize), Div(usize, usize), Lin(Vec<(usize, f64)>), Pre(String) }

struct Plan { outs: Vec<u32>, cols: Vec<Col> }

impl Plan {
    /// Engine rows y (len(outs) each) -> selection rows; Pre columns take `pre` (value before the run).
    fn apply(&self, y: &[f64], pre: &[f64]) -> Vec<f64> {
        let no = self.outs.len().max(1);
        let nt = if self.outs.is_empty() { 0 } else { y.len() / no };
        let mut out = Vec::with_capacity(nt * self.cols.len());
        for r in 0..nt {
            let y = &y[r * no..r * no + self.outs.len()];
            let mut k = 0;
            for c in &self.cols {
                out.push(match c {
                    Col::Out(i) => y[*i],
                    Col::Mul(i, j) => y[*i] * y[*j],
                    Col::Div(i, j) => y[*i] / y[*j],
                    Col::Lin(t) => t.iter().map(|&(i, c)| y[i] * c).sum(),
                    Col::Pre(_) => { k += 1; pre[k - 1] }
                });
            }
        }
        out
    }
}

/// The compiled model and its state.
struct Core {
    m: Model,
    info: Info,
    keys: Vec<String>,
    kres: Vec<Key>,
    kpos: HashMap<String, usize>,
    xo: Vec<u32>,
    /// values at load (origin), initial values (init(...), resetAll), current values
    x0: Vec<f64>, xi: Vec<f64>, x: Vec<f64>,
    t: f64,
    /// keys reset() restores: floating species and rate-rule targets
    tgt: Vec<usize>,
    /// initial assignments of reset targets, re-evaluated by reset() (None: reset to stored values)
    ia: Option<Vec<(String, String)>>,
    /// compartment volumes baked into the engine model
    cconst: HashSet<String>,
    sp: Vec<Option<usize>>, // per floating id (roadrunner order): engine state index
    rx: Vec<usize>,         // per reaction id: engine reaction index
    plain: bool,            // floating species are all engine species without rate rules
    nss: bool,              // events or time-dependent rates: no steady state
    /// libstruct's view: SBML species without boundaryCondition (constant ones too), N from the SBML
    ls_rows: Vec<String>, ls_n: Vec<f64>, ls_nr: usize,
    /// libstruct's analysis and csum_dep, computed on first use (a dense QR: seconds on big models)
    ls_cell: std::cell::OnceCell<(Ls, Vec<usize>)>,
    /// dependent species (index into ls's dependents) with a moiety total: those not set by rules
    /// initial concentrations of rule-set floating species as the SBML declares them (0 if unset)
    rule_init: HashMap<String, f64>,
    w: Work, fresh: bool, live: bool, // the running simulation (event state, delay history)
    q: Work, qok: bool,               // scratch work loaded with (x, t) when qok
    rd: Option<RateDerivs>,
    csum: Vec<f64>,
}

fn names_of(src: &str) -> Option<HashSet<String>> {
    let mut seen = HashSet::new();
    crate::expr::parse(src, &mut Default::default(), |n| { seen.insert(n.to_string()); Some(crate::expr::Op::Const(0)) }).ok()?;
    Some(seen)
}

/// Value of an infix formula with every name given by `get` (constant folding does the arithmetic).
fn eval(src: &str, get: &mut dyn FnMut(&str) -> Option<f64>) -> Result<f64, String> {
    let mut dag = crate::expr::Dag::default();
    let mut miss = None;
    let n = crate::expr::parse(src, &mut dag, |n| match get(n) {
        Some(v) => Some(crate::expr::Op::Const(v.to_bits())),
        None => { miss = Some(n.to_string()); None }
    })?;
    if let Some(n) = miss { return Err(format!("unknown symbol '{n}' in '{src}'")); }
    dag.val(n).ok_or_else(|| format!("cannot evaluate '{src}'"))
}

impl Core {
    fn new(m: Model, info: Info, nss: bool, text: &str) -> Result<Core, String> {
        let [fids, bids, gids, cids, rids] = &info.ids;
        let known = |k: &String| m.resolve(std::slice::from_ref(k)).is_ok();
        let mut keys: Vec<String> = gids.iter().chain(cids).chain(bids).filter(|k| known(k)).cloned().collect();
        // states in none of roadrunner's groups (e.g. a species reference with a rate rule)
        let extra: Vec<String> = m.state_ids.iter().filter(|k| !keys.contains(k) && !fids.contains(k)).cloned().collect();
        keys.extend(extra);
        keys.extend(fids.iter().filter(|k| known(k)).cloned());
        let outs: HashSet<&String> = m.outputs.iter().map(|o| &o.0).collect();
        if !m.state_ids.iter().all(|s| keys.contains(s)) || !info.ids.iter().flatten().all(|k| outs.contains(k)) {
            return Err("not all of the model is in the engine model".into());
        }
        let rx: Vec<usize> = rids.iter().map(|r| m.reactions.iter().position(|x| x == r))
            .collect::<Option<_>>().filter(|v: &Vec<usize>| v.len() == m.reactions.len()).ok_or("reactions differ from the engine model")?;
        let kres = m.resolve(&keys)?;
        let xo = m.output_regs(&keys)?;
        let mut q = m.work();
        q.atol_scale.fill(1.0);
        m.load(&mut q, &[], &[]);
        let mut x0 = vec![0.0; keys.len()];
        m.record(&mut q, 0.0, &xo, &mut x0);
        let kpos: HashMap<String, usize> = keys.iter().enumerate().map(|(i, k)| (k.clone(), i)).collect();
        let tgt: Vec<usize> = keys.iter().enumerate().filter(|(_, k)| fids.contains(k) || info.rate_rule_ids.contains(k)).map(|(i, _)| i).collect();
        // reset() re-evaluates initial assignments of the targets over species, compartments and
        // rule-free parameters with constant compartments (as the lazy Python path); else stored values
        let ias: HashMap<&String, &String> = info.ia.iter().map(|(k, v)| (k, v)).collect();
        let plain_p: HashSet<&String> = gids.iter().filter(|g| !info.rule_vars.contains(g)).collect();
        let mut ok = !cids.iter().any(|c| info.changing.contains(c) || ias.contains_key(c));
        let mut ia = vec![];
        for &j in &tgt {
            if let Some(f) = ias.get(&keys[j]) {
                ok = ok && names_of(f).is_some_and(|n| n.iter().all(|x| fids.contains(x) || cids.contains(x) || plain_p.contains(x)));
                ia.push((keys[j].clone(), f.to_string()));
            }
        }
        let cconst = fids.iter().filter(|s| !info.general.contains(*s)).map(|s| info.comp_of[s].clone()).collect();
        let sp: Vec<Option<usize>> = fids.iter().map(|s| m.state_ids.iter().position(|x| x == s)).collect();
        // assignment-rule species are not states: zero rows/columns in Jacobians and MCA (as roadrunner)
        let plain = fids.iter().all(|s| !info.general.contains(s) && (m.species.contains(s) || !m.state_ids.contains(s)));
        let doc = read_doc(text)?;
        let (mut ls_rows, mut ls_n) = (vec![], vec![]);
        if let Some(sm) = &doc.model {
            ls_rows = sm.species.iter().filter(|s| s.boundary_condition != Some(true)).map(|s| s.id.clone()).collect();
            let nr = sm.reactions.len();
            ls_n = vec![0.0; ls_rows.len() * nr];
            for (j, r) in sm.reactions.iter().enumerate() {
                for (sg, l) in [(-1.0, &r.reactants), (1.0, &r.products)] {
                    for sr in l {
                        if let Some(i) = ls_rows.iter().position(|x| *x == sr.species) { ls_n[i * nr + j] += sg * sr.stoichiometry.unwrap_or(1.0); }
                    }
                }
            }
        }
        let ls_nr = doc.model.as_ref().map_or(0, |m| m.reactions.len());
        let rule_init = doc.model.iter().flat_map(|sm| sm.species.iter().filter(|s| fids.contains(&s.id) && info.rule_vars.contains(&s.id)).map(|s| {
            let v = sm.compartments.iter().find(|c| c.id == s.compartment).and_then(|c| c.size).unwrap_or(1.0);
            (s.id.clone(), s.initial_concentration.or(s.initial_amount.map(|a| a / v)).unwrap_or(0.0))
        })).collect();
        let w = m.work();
        Ok(Core { xi: x0.clone(), x: x0.clone(), x0, t: 0.0, tgt, ia: ok.then_some(ia), cconst, sp, rx, plain, nss, ls_nr, ls_cell: Default::default(), w, fresh: true,
                  live: false, q, qok: false, rd: None, csum: vec![], ls_rows, ls_n, rule_init, keys, kres, kpos, xo, m, info })
    }

    fn fids(&self) -> &[String] { &self.info.ids[0] }
    fn ls_all(&self) -> &(Ls, Vec<usize>) {
        self.ls_cell.get_or_init(|| {
            let ls = cons::analyze(&self.ls_n, self.ls_rows.len(), self.ls_nr);
            let dep = (0..ls.nd()).filter(|&d| !self.info.rule_vars.contains(&self.ls_rows[ls.order[ls.ni + d]])).collect();
            (ls, dep)
        })
    }
    fn ls(&self) -> &Ls { &self.ls_all().0 }
    /// dependent species (index into ls's dependents) with a moiety total: those not set by rules
    fn csum_dep(&self) -> &[usize] { &self.ls_all().1 }
    fn comp_of(&self, s: &str) -> &str { &self.info.comp_of[s] }
    fn touch(&mut self) { self.live = false; self.qok = false; }

    /// Engine outputs (registers) at the current state and time.
    fn out(&mut self, regs: &[u32]) -> Vec<f64> {
        let mut row = vec![0.0; regs.len()];
        if self.live {
            self.m.record(&mut self.w, self.t, regs, &mut row);
        } else {
            if !self.qok { self.m.load(&mut self.q, &self.kres, &self.x); self.qok = true; }
            self.m.record(&mut self.q, self.t, regs, &mut row);
        }
        row
    }

    fn reg(&self, id: &str) -> Result<u32, String> {
        self.m.output_regs(&[id.to_string()]).map(|r| r[0]).map_err(|_| format!("unknown symbol '{id}'"))
    }
    fn is_species(&self, s: &str) -> bool { self.info.comp_of.contains_key(s) }

    /// (value, volume) of species s from the engine; value is an amount if s is in amount_params.
    fn species_raw(&mut self, s: &str) -> Result<(f64, f64), String> {
        let r = [self.reg(s)?, self.reg(self.comp_of(s))?];
        let v = self.out(&r);
        Ok((v[0], v[1]))
    }
    fn species(&mut self, s: &str, conc: bool) -> Result<f64, String> {
        let (v, vol) = self.species_raw(s)?;
        Ok(match (conc, self.info.amount_params.contains(s)) { (true, true) => v / vol, (false, false) => v * vol, _ => v })
    }

    /// Initial concentrations reset() gives the floating species (fids order).
    fn init_concs(&mut self) -> Result<Vec<f64>, String> { self.init_concs_v(true) }

    /// Initial volume of compartment `c`: its initial value, else (set by an assignment rule) the
    /// model's initial value of it.
    fn init_volume(&self, c: &str) -> Option<f64> {
        self.kpos.get(c).map(|&q| self.xi[q]).or_else(|| self.info.values.get(c).copied())
    }

    /// `keep_amounts`: as reset() (species of a compartment changed since its initial value keep their
    /// initial amount in the current volume); false: the initial concentrations themselves (init([S])).
    fn init_concs_v(&mut self, keep_amounts: bool) -> Result<Vec<f64>, String> {
        let mut vals: HashMap<String, f64> = self.tgt.iter().map(|&j| (self.keys[j].clone(), self.xi[j])).collect();
        if let Some(ia) = self.ia.clone().filter(|v| !v.is_empty()) {
            let cur: HashMap<&String, f64> = self.keys.iter().zip(&self.x).map(|(k, &v)| (k, v)).collect();
            let fids: HashSet<&String> = self.fids().iter().collect();
            let mut done: HashMap<String, f64> = HashMap::new();
            fn get(n: &str, c: &Core, ia: &[(String, String)], vals: &HashMap<String, f64>, cur: &HashMap<&String, f64>,
                   fids: &HashSet<&String>, done: &mut HashMap<String, f64>, depth: usize) -> Option<f64> {
                if let Some((_, f)) = ia.iter().find(|(k, _)| k == n) {
                    if let Some(&v) = done.get(n) { return Some(v); }
                    if depth > 100 { return None; }
                    let v = eval(f, &mut |m| get(m, c, ia, vals, cur, fids, done, depth + 1)).ok()?;
                    done.insert(n.to_string(), v);
                    return Some(v);
                }
                if fids.contains(&n.to_string()) {
                    let v = vals.get(n).copied()?;
                    return Some(if c.info.hosu.contains(n) { v * cur.get(&c.info.comp_of[n].to_string()).copied()? } else { v });
                }
                cur.get(&n.to_string()).copied()
            }
            for (k, _) in &ia {
                let v = get(k, self, &ia, &vals, &cur, &fids, &mut done, 0).ok_or_else(|| format!("initial assignment of {k}"))?;
                let v = if self.info.hosu.contains(k) && fids.contains(k) { v / cur[&self.info.comp_of[k]] } else { v };
                vals.insert(k.clone(), v);
            }
        }
        let fids = self.fids().to_vec();
        let mut out = vec![];
        for s in &fids {
            let v = match vals.get(s) {
                // amounts are kept when a compartment changed since its initial value (as roadrunner)
                // (a compartment reset() restores itself, e.g. a rate-rule one, keeps none)
                Some(&v) => match self.kpos.get(self.comp_of(s)) {
                    Some(&c) if keep_amounts && self.x[c] != 0.0 && !self.tgt.contains(&c) => v * self.xi[c] / self.x[c], _ => v },
                None => match self.rule_init.get(s) { Some(&v) => v, None => self.species(s, true)? },
            };
            out.push(v);
        }
        Ok(out)
    }

    fn reset(&mut self) -> Result<(), String> {
        let ic = self.init_concs()?;
        let vals: HashMap<String, f64> = self.fids().iter().cloned().zip(ic).collect();
        for &j in &self.tgt.clone() {
            self.x[j] = vals.get(&self.keys[j]).copied().unwrap_or(self.xi[j]);
        }
        self.t = 0.0;
        self.fresh = true;
        self.touch();
        Ok(())
    }

    /// Per-state atol multipliers: |initial concentration| of each floating species (1 if ~0), and
    /// |initial value| of other states (rate-rule variables; 1 if ~0), as roadrunner's CVODE scales a
    /// scalar atol by initial amounts (CVODEIntegrator::getAbsoluteToleranceVector).
    fn atol_scale(&mut self) -> Result<Vec<f64>, String> {
        let ic = self.init_concs()?;
        let sc = |v: f64| if v.abs() < 1e-30 || !v.is_finite() { 1.0 } else { v.abs() };
        let mut a: Vec<f64> = self.m.state_ids.iter().map(|id| self.kpos.get(id).map_or(1.0, |&k| sc(self.xi[k]))).collect();
        for (i, s) in self.sp.iter().enumerate() {
            if let Some(j) = *s { a[j] = sc(ic[i]); }
        }
        Ok(a)
    }

    /// advance() with libroadrunner's euler / rk4 / rk45 (crate::fixed); the model is left at the end.
    fn advance_fixed(&mut self, times: &[f64], outs: &[u32], method: crate::fixed::Method, step: f64, atol: f64)
                     -> Result<Vec<f64>, String> {
        let mut out = vec![0.0; times.len() * outs.len()];
        let Core { m, w, kres, x, fresh, .. } = self;
        if *fresh { m.events_reset(w); }
        m.load(w, kres, x);
        crate::fixed::integrate_fixed(m, w, method, times, step, 0.0, atol, outs, &mut out)?;
        let (tend, xo) = (times[times.len() - 1], self.xo.clone());
        self.m.record(&mut self.w, tend, &xo, &mut self.x);
        (self.t, self.fresh, self.qok, self.live) = (tend, false, false, true);
        Ok(out)
    }

    /// rk45 with variable_step_size (crate::fixed::integrate_rk45_variable) from t0 to t1: rows of
    /// `outs`, one per accepted step; the model is left at the last row's time.
    #[allow(clippy::too_many_arguments)]
    fn advance_rk45_variable(&mut self, t0: f64, t1: f64, outs: &[u32], hmin: f64, hmax: f64, eps: f64, max_rows: usize)
                             -> Result<Vec<f64>, String> {
        let Core { m, w, kres, x, fresh, .. } = self;
        if *fresh { m.events_reset(w); }
        m.load(w, kres, x);
        let (ts, out) = crate::fixed::integrate_rk45_variable(m, w, t0, t1, hmin, hmax, eps, outs, max_rows)?;
        let (tend, xo) = (ts[ts.len() - 1], self.xo.clone());
        self.m.record(&mut self.w, tend, &xo, &mut self.x);
        (self.t, self.fresh, self.qok, self.live) = (tend, false, false, true);
        Ok(out)
    }

    /// Integrate from the current state, carrying event state and retrying with tighter tolerances.
    /// `tol`: absolute tolerance per state (RoadRunner::cvode_tolerances); `ctl`: the cvode settings
    /// (with ctl.var: rows at the solver's steps from times[0] to the last time, model left at the last row).
    fn advance(&mut self, times: &[f64], outs: &[u32], rtol: f64, tol: &[f64], ctl: Ctl) -> Result<Vec<f64>, String> {
        let asc = self.atol_scale()?;
        // the retry: every state at most 1e-12 times its scale
        let retry = tol.iter().zip(&asc).map(|(&t, &a)| (1e-12 * a / t).min(1.0)).fold(1.0, f64::min);
        let (no, nt) = (outs.len(), times.len());
        let mut out = vec![0.0; nt * no];
        let Core { m, w, kres, x, fresh, .. } = self;
        w.atol_scale.copy_from_slice(tol);
        w.ctl = ctl;
        if *fresh { m.events_reset(w); }
        let saved = m.event_state(w);
        let hist = w.hist.save();
        let mut tried = f64::NAN;
        let mut ok = false;
        let atol = tol.iter().zip(&asc).map(|(&t, &a)| t / a).fold(0.0, f64::max);
        for a in [1.0, retry] {
            if a == tried { break; }
            tried = a;
            m.load(w, kres, x);
            if !*fresh { w.hist.restore(hist.clone()); }
            if m.simulate_cont(w, times, rtol, a, outs, &mut out) { ok = true; break; }
            let _ = m.set_event_state(w, saved.0.clone(), saved.1.clone());
        }
        w.ctl = Ctl::default();
        if !ok {
            w.hist.restore(hist);
            return Err(format!("integration failed (atol {atol:e} and 1e-12)"));
        }
        let (tend, out) = if ctl.var { (*w.rows.0.last().unwrap(), std::mem::take(&mut w.rows.1)) } else { (times[nt - 1], out) };
        if out.iter().any(|v| !v.is_finite()) && (0..out.len()).any(|i| !out[i].is_finite() && out[i % no.max(1)].is_finite()) {
            return Err("integration went non-finite".into());
        }
        let xo = self.xo.clone();
        self.m.record(&mut self.w, tend, &xo, &mut self.x);
        self.t = tend;
        self.fresh = false;
        self.qok = false;
        self.live = true;
        Ok(out)
    }

    /// Engine state after a solve (w.s) back into x.
    fn take_states(&mut self, s: &[f64]) {
        for (j, id) in self.m.state_ids.iter().enumerate() {
            if let Some(&k) = self.kpos.get(id) { self.x[k] = s[j]; }
        }
        self.touch();
    }

    /// Steady state from the current state into x with the solver settings (rr::ss); the solution
    /// is left in q.
    fn solve(&mut self, cfg: &ss::SsCfg) -> Result<Option<f64>, String> { self.solve_with(cfg) }

    /// d state/dt (engine units: concentrations) at the current state and time.
    fn state_rates(&mut self) -> Vec<f64> {
        let mut r = vec![0.0; self.m.nx];
        if self.live { self.m.state_rates(&mut self.w, self.t, &mut r); return r; }
        if !self.qok { self.m.load(&mut self.q, &self.kres, &self.x); self.qok = true; }
        self.m.state_rates(&mut self.q, self.t, &mut r);
        r
    }

    fn need_plain(&self) -> Result<(), String> {
        if self.plain { Ok(()) } else { Err("species with rate rules or outside the reaction network: not supported for this analysis".into()) }
    }
    fn need_mca(&self) -> Result<(), String> {
        self.need_plain()?;
        if self.nss { Err("MCA: the model has events or time-dependent rates (no steady state)".into()) } else { Ok(()) }
    }

    /// (v, d v / d [states | params]) at the current state, exact.
    fn rate_derivs(&mut self) -> (Vec<f64>, Vec<f64>) {
        if self.rd.is_none() { self.rd = Some(RateDerivs::new(&self.m)); }
        self.rd.as_ref().unwrap().eval(&self.m, &self.kres, &self.x, self.t)
    }
}

/// Engine model + info from SBML text; `sizes` override compartment sizes.
fn compile(text: &str, base: Option<&Path>, sizes: &[(String, f64)]) -> Result<(Model, Info, String), String> {
    let mut t = if text.contains("/comp/") { sbml_io::comp::flatten(text, base)? } else { text.to_string() };
    if !sizes.is_empty() {
        let mut doc = read_doc(&t)?;
        for c in doc.model.as_mut().map_or(&mut vec![], |m| &mut m.compartments) {
            if let Some((_, v)) = sizes.iter().find(|s| s.0 == c.id) { c.size = Some(*v); }
        }
        t = sbml_io::write(&doc);
    }
    let (m, i) = crate::sbml::from_sbml(&t).map_err(|e| e.to_string())?;
    Ok((m, i, t))
}

fn read_doc(t: &str) -> Result<sbml_io::Document, String> {
    sbml_io::read(t).map(|r| r.document).map_err(|e| format!("SBML errors: {}", e.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("; ")))
}

/// A libroadrunner-like simulator instance. Errors are `Err(String)` (the C API's getLastError).
pub struct RoadRunner {
    sbml: String,
    base: Option<PathBuf>,
    c: Option<Core>,
    pub start: f64,
    pub end: f64,
    pub points: usize,
    /// explicit output times (setTimes); overrides start/end/points when set
    pub times: Option<Vec<f64>>,
    sel: Vec<String>,
    ssel: Vec<String>,
    cma: bool,
    pub diff_step: f64,
    pub flux_threshold: f64,
    integ: Solver,
    integs: Vec<Solver>, // settings of the integrators not current
    ss: Solver,
    seed_ctr: (u64, u64),
    result: Option<Matrix>,
}

impl Default for RoadRunner { fn default() -> Self { Self::new() } }

impl RoadRunner {
    pub fn new() -> RoadRunner {
        RoadRunner { sbml: String::new(), base: None, c: None, start: 0.0, end: 5.0, points: 51, times: None, sel: vec![], ssel: vec![],
                     cma: false, diff_step: 0.02, flux_threshold: 1e-9, integ: Solver::integrator("cvode").unwrap(),
                     integs: INTEGRATORS[1..].iter().map(|n| Solver::integrator(n).unwrap()).collect(),
                     ss: Solver::steady_state("nleq2").unwrap(), seed_ctr: (u64::MAX, 0), result: None }
    }

    fn core(&self) -> Result<&Core, String> { self.c.as_ref().ok_or_else(|| "no model loaded".to_string()) }
    fn core_mut(&mut self) -> Result<&mut Core, String> { self.c.as_mut().ok_or_else(|| "no model loaded".to_string()) }

    // ------------------------------------------------------------ load / save

    /// Load SBML text, or a file path (its directory then resolves comp external models).
    pub fn load_sbml(&mut self, text_or_path: &str, base: Option<&Path>) -> Result<(), String> {
        let (text, base) = if !text_or_path.trim_start().starts_with('<') {
            let p = Path::new(text_or_path);
            let t = std::fs::read_to_string(p).map_err(|e| format!("{text_or_path}: {e}"))?;
            (t, base.map(Path::to_path_buf).or_else(|| p.canonicalize().ok().and_then(|p| p.parent().map(Path::to_path_buf))))
        } else {
            (text_or_path.to_string(), base.map(Path::to_path_buf))
        };
        let (old, olds) = (std::mem::take(&mut self.sbml), self.base.take());
        let prev = self.c.take();
        self.sbml = text;
        self.base = base;
        if let Err(e) = self.rebuild(None) {
            (self.sbml, self.base, self.c) = (old, olds, prev);
            return Err(e);
        }
        self.seed_ctr = (u64::MAX, 0);
        self.result = None;
        self.default_selections();
        Ok(())
    }

    /// Compile self.sbml; carry the state of `old` (by id; ids whose initial value changed take the
    /// new one), and bake changed compartment volumes into the engine model.
    fn rebuild(&mut self, old: Option<Core>) -> Result<(), String> {
        let compile1 = |sizes: &[(String, f64)]| -> Result<Core, String> {
            let (m, info, text) = compile(&self.sbml, self.base.as_deref(), sizes)?;
            let nss = self.sbml.contains("symbols/time") || info.n_events > 0;
            Core::new(m, info, nss, &text)
        };
        let mut c = compile1(&[])?;
        if let Some(o) = &old {
            // as libroadrunner's LLVMModelGenerator::regenerateModel: every current value (species,
            // compartments, parameters) and the time carry over; initial values are the edited SBML's
            // where the edit changed them, else the user's (libroadrunner keeps init() edits in its document)
            for (j, k) in c.keys.iter().enumerate() {
                if let Some(&i) = o.kpos.get(k) {
                    c.x[j] = o.x[i];
                    if o.x0[i] == c.x0[j] || (o.x0[i].is_nan() && c.x0[j].is_nan()) { c.xi[j] = o.xi[i]; }
                }
            }
            c.t = o.t;
        }
        let sizes: Vec<(String, f64)> = c.cconst.iter().filter_map(|k| c.kpos.get(k).filter(|&&j| c.x[j] != c.x0[j]).map(|&j| (k.clone(), c.x[j]))).collect();
        if !sizes.is_empty() {
            let mut c2 = compile1(&sizes)?;
            if c2.keys != c.keys { return Err("recompiled model differs".into()); }
            (c2.x0, c2.xi, c2.x, c2.t, c2.ia) = (c.x0, c.xi, c.x, c.t, c.ia);
            c = c2;
        }
        c.touch();
        self.c = Some(c);
        if self.cma { self.c.as_mut().unwrap().csum = self.csum_now()?; }
        Ok(())
    }

    pub fn is_model_loaded(&self) -> bool { self.c.is_some() }
    pub fn clear_model(&mut self) { self.c = None; self.sbml.clear(); self.result = None; self.sel.clear(); self.ssel.clear(); }

    /// The SBML as loaded (with edits).
    /// The model's SBML as libroadrunner's getSBML gives it: libsbml writes the stoichiometry of every
    /// species reference (default 1); the text is otherwise the one loaded.
    pub fn get_sbml(&self) -> Result<String, String> {
        self.core()?;
        let Ok(mut doc) = read_doc(&self.sbml) else { return Ok(self.sbml.clone()) };
        let Some(m) = doc.model.as_mut() else { return Ok(self.sbml.clone()) };
        let mut changed = false;
        for sr in m.reactions.iter_mut().flat_map(|r| r.reactants.iter_mut().chain(r.products.iter_mut())) {
            if sr.stoichiometry.is_none() && sr.stoichiometry_math.is_none() { sr.stoichiometry = Some(1.0); changed = true; }
        }
        Ok(if changed && !self.sbml.contains("/comp/") { sbml_io::write(&doc) } else { self.sbml.clone() })
    }

    /// The SBML with the current values as initial values (initial assignments removed).
    pub fn get_current_sbml(&mut self) -> Result<String, String> {
        let text = if self.sbml.contains("/comp/") { sbml_io::comp::flatten(&self.sbml, self.base.as_deref())? } else { self.sbml.clone() };
        let mut doc = read_doc(&text)?;
        let c = self.core_mut()?;
        let Some(m) = doc.model.as_mut() else { return Ok(text) };
        let mut val = |id: &str| -> Option<f64> { c.reg(id).ok().map(|r| c.out(&[r])[0]) };
        for p in &mut m.parameters { if let Some(v) = val(&p.id) { p.value = Some(v); } }
        for x in &mut m.compartments { if let Some(v) = val(&x.id) { x.size = Some(v); } }
        for s in &mut m.species {
            if !c.is_species(&s.id) || c.reg(&s.id).is_err() { continue; }
            if s.initial_concentration.is_some() { s.initial_concentration = Some(c.species(&s.id, true)?); }
            else { s.initial_amount = Some(c.species(&s.id, false)?); }
        }
        m.initial_assignments.clear();
        Ok(sbml_io::write(&doc))
    }

    /// Edit the SBML text; with `regen`, recompile carrying over the current values (regenerate()).
    pub fn edit(&mut self, f: impl FnOnce(&mut String) -> Result<(), String>, regen: bool) -> Result<(), String> {
        let mut s = self.sbml.clone();
        f(&mut s)?;
        self.sbml = s;
        if regen { self.regenerate() } else { Ok(()) }
    }

    /// Recompile the (edited) SBML, keeping current values, time and settings. The event state and
    /// delay history restart. Selections that no longer resolve fall back to the defaults.
    pub fn regenerate(&mut self) -> Result<(), String> { self.regenerate_opt(false) }

    /// regenerate(); with `reset`, then reset() unless the model ran since the last reset (as
    /// libroadrunner's regenerateModel(force, reset) after rule / initial-assignment / event edits).
    pub fn regenerate_opt(&mut self, reset: bool) -> Result<(), String> {
        let fresh = self.c.as_ref().is_some_and(|c| c.fresh);
        self.regenerate_inner()?;
        if reset && fresh { self.reset()?; }
        Ok(())
    }

    fn regenerate_inner(&mut self) -> Result<(), String> {
        let old = self.c.take();
        if let Err(e) = self.rebuild(old) { return Err(e); }
        if self.plan(&self.sel.clone()).is_err() || self.check(&self.ssel.clone()).is_err() { self.default_selections(); }
        Ok(())
    }

    /// Save the model, values, time and event state to a file (text format).
    pub fn save_state(&self, path: &Path) -> Result<(), String> {
        let c = self.core()?;
        let v = |x: &[f64]| x.iter().map(|v| format!("{v:?}")).collect::<Vec<_>>().join(" ");
        let (prev, pend) = c.m.event_state(&c.w);
        let mut s = format!("roadrunner_nitro state 1\nt {:?}\ncma {}\nfresh {}\nx {}\nxi {}\nprev {}\n", c.t, self.cma as u8, c.fresh as u8,
                            v(&c.x), v(&c.xi), prev.iter().map(|&b| if b { "1" } else { "0" }).collect::<Vec<_>>().join(" "));
        for (t, e, vals) in &pend { s += &format!("pend {t:?} {e} {}\n", v(vals)); }
        // the delay() history a continued run reads
        for (k, a) in ["hist_t", "hist_y", "hist_f"].iter().zip(c.w.hist.series()) { s += &format!("{k} {}\n", v(a)); }
        s += "sbml\n";
        s += &self.sbml;
        std::fs::write(path, s).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn load_state(&mut self, path: &Path) -> Result<(), String> {
        let s = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let bad = || format!("{}: not a roadrunner_nitro state file", path.display());
        let (head, sbml) = s.split_once("\nsbml\n").ok_or_else(bad)?;
        let mut f: HashMap<&str, &str> = HashMap::new();
        let mut pend = vec![];
        let nums = |s: &str| s.split_whitespace().map(|x| x.parse::<f64>()).collect::<Result<Vec<f64>, _>>().map_err(|_| bad());
        for l in head.lines().skip(1) {
            let (k, v) = l.split_once(' ').unwrap_or((l, ""));
            if k == "pend" {
                let mut it = v.splitn(3, ' ');
                let t = it.next().and_then(|x| x.parse().ok()).ok_or_else(bad)?;
                let e = it.next().and_then(|x| x.parse().ok()).ok_or_else(bad)?;
                pend.push((t, e, nums(it.next().unwrap_or(""))?));
            } else { f.insert(k, v); }
        }
        let get = |k: &str| f.get(k).copied().ok_or_else(bad);
        self.load_sbml(sbml, None)?;
        self.cma = get("cma")? == "1";
        let c = self.c.as_mut().unwrap();
        let (x, xi) = (nums(get("x")?)?, nums(get("xi")?)?);
        if x.len() != c.keys.len() || xi.len() != c.keys.len() { return Err(bad()); }
        (c.x, c.xi, c.t, c.fresh) = (x, xi, get("t")?.parse().map_err(|_| bad())?, get("fresh")? == "1");
        let prev = get("prev")?.split_whitespace().map(|b| b == "1").collect();
        c.m.set_event_state(&mut c.w, prev, pend)?;
        if let (Some(t), Some(y), Some(fv)) = (f.get("hist_t"), f.get("hist_y"), f.get("hist_f")) {
            c.w.hist.restore((nums(t)?, nums(y)?, nums(fv)?));
        }
        c.touch();
        if self.cma { self.c.as_mut().unwrap().csum = self.csum_now()?; }
        self.default_selections();
        Ok(())
    }

    // ------------------------------------------------------------ ids

    /// Floating species ids (conserved moiety analysis on: libstruct's order, independent first).
    pub fn floating_ids(&self) -> Vec<String> {
        let Some(c) = &self.c else { return vec![] };
        if self.cma { c.ls().order.iter().map(|&i| &c.ls_rows[i]).filter(|s| c.info.rows.contains(s)).cloned().collect() } else { c.fids().to_vec() }
    }
    pub fn independent_ids(&self) -> Vec<String> {
        self.c.as_ref().map_or(vec![], |c| c.ls().order[..c.ls().ni].iter().map(|&i| c.ls_rows[i].clone()).collect())
    }
    pub fn dependent_ids(&self) -> Vec<String> {
        self.c.as_ref().map_or(vec![], |c| c.ls().order[c.ls().ni..].iter().map(|&i| c.ls_rows[i].clone()).collect())
    }
    pub fn boundary_ids(&self) -> Vec<String> { self.ids(1) }
    /// Global parameters (+ conserved moiety totals `_CSUMi` with moiety analysis on).
    pub fn global_parameter_ids(&self) -> Vec<String> {
        let mut v = self.ids(2);
        let Some(c) = &self.c else { return v };
        let at = v.iter().position(|p| c.info.rule_vars.contains(p)).unwrap_or(v.len()); // before rule-set parameters
        v.splice(at..at, self.conserved_moiety_ids());
        v
    }
    pub fn compartment_ids(&self) -> Vec<String> { self.ids(3) }
    pub fn reaction_ids(&self) -> Vec<String> { self.ids(4) }
    pub fn event_ids(&self) -> Vec<String> { self.c.as_ref().map_or(vec![], |c| c.info.event_ids.clone()) }
    pub fn rate_rule_ids(&self) -> Vec<String> { self.c.as_ref().map_or(vec![], |c| c.info.rate_rule_ids.clone()) }
    pub fn conserved_moiety_ids(&self) -> Vec<String> {
        match &self.c { Some(c) if self.cma => (0..c.csum_dep().len()).map(|k| format!("_CSUM{k}")).collect(), _ => vec![] }
    }
    pub fn eigenvalue_ids(&self) -> Vec<String> {
        self.floating_ids().iter().flat_map(|s| ["eigen", "eigenReal", "eigenImag"].map(|f| format!("{f}({s})"))).collect()
    }
    /// Rate-rule targets, then the floating species not set by assignment rules (independent ones
    /// with moiety analysis on).
    pub fn rates_of_change_ids(&self) -> Vec<String> {
        let (fs, mut skip) = (self.floating_ids(), if self.cma { self.dependent_ids() } else { vec![] });
        let rr = self.rate_rule_ids();
        if let Some(c) = &self.c { skip.extend(c.info.rule_vars.iter().filter(|v| !rr.contains(v)).cloned()); } // assignment rules
        rr.iter().filter(|r| !fs.contains(r)).cloned().chain(fs.iter().filter(|s| !skip.contains(s)).cloned())
            .map(|s| format!("{s}'")).collect()
    }
    pub fn floating_init_ids(&self) -> Vec<String> { self.floating_ids().iter().map(|s| format!("init([{s}])")).collect() }
    fn ids(&self, k: usize) -> Vec<String> { self.c.as_ref().map_or(vec![], |c| c.info.ids[k].clone()) }

    /// Time-course symbol groups as the C API lists them (elasticities: `ec_ids`).
    pub fn available_time_course_symbols(&self) -> Vec<(&'static str, Vec<String>)> {
        let fs = self.floating_ids();
        let bs = self.boundary_ids();
        vec![("Time", vec!["time".into()]), ("Floating Species", fs.iter().map(|s| format!("[{s}]")).collect()),
             ("Boundary Species", bs.iter().map(|s| format!("[{s}]")).collect()), ("Floating Species (amount)", fs),
             ("Boundary Species (amount)", bs), ("Global Parameters", self.global_parameter_ids()), ("Fluxes", self.reaction_ids()),
             ("Rates of Change", self.rates_of_change_ids()), ("Volumes", self.compartment_ids()), ("Eigenvalues", self.eigenvalue_ids())]
    }
    /// Steady-state symbol groups as the C API lists them (control coefficients: `cc_ids`, `ec_ids`).
    pub fn available_steady_state_symbols(&self) -> Vec<(&'static str, Vec<String>)> {
        let mut v = self.available_time_course_symbols();
        v.retain(|g| !matches!(g.0, "Time" | "Rates of Change"));
        v
    }
    /// Per floating species (`flux`: per reaction): its control coefficient selections over global
    /// parameters, boundary species and moiety totals.
    pub fn cc_ids(&self, flux: bool, scaled: bool) -> Vec<(String, Vec<String>)> {
        let ps: Vec<String> = self.ids(2).into_iter().chain(self.boundary_ids()).chain(self.conserved_moiety_ids()).collect();
        let f = if scaled { "cc" } else { "ucc" };
        (if flux { self.reaction_ids() } else { self.floating_ids() }).into_iter()
            .map(|s| { let l = ps.iter().map(|p| format!("{f}({s}, {p})")).collect(); (s, l) }).collect()
    }
    /// Per reaction: elasticity selections over floating and boundary species, parameters, totals.
    pub fn ec_ids(&self, scaled: bool) -> Vec<(String, Vec<String>)> {
        let xs: Vec<String> = self.floating_ids().into_iter().chain(self.boundary_ids()).chain(self.ids(2)).chain(self.conserved_moiety_ids()).collect();
        let f = if scaled { "ec" } else { "uec" };
        self.reaction_ids().into_iter().map(|r| { let l = xs.iter().map(|x| format!("{f}({r}, {x})")).collect(); (r, l) }).collect()
    }

    // ------------------------------------------------------------ conserved moieties

    pub fn conservation_laws(&self) -> bool { self.cma }

    /// setComputeAndAssignConservationLaws: species reorder (libstruct), `_CSUMi` totals; selections
    /// go back to the defaults.
    pub fn set_conservation_laws(&mut self, on: bool) -> Result<(), String> {
        if let (true, Some(c)) = (on, &self.c) {
            if let Some(s) = c.m.event_src.iter().flat_map(|e| &e.assign).map(|a| &a.0).find(|s| c.info.rows.contains(*s)) {
                return Err(format!("cannot perform moiety conversion when floating species have events ({s})"));
            }
        }
        self.cma = on;
        if self.c.is_some() {
            if on { self.c.as_mut().unwrap().csum = self.csum_now()?; }
            self.default_selections();
        }
        Ok(())
    }

    /// Amounts of the libstruct species.
    fn amounts(&mut self) -> Result<Vec<f64>, String> {
        let c = self.core_mut()?;
        let rows = c.ls_rows.clone();
        rows.iter().map(|s| c.species(s, false)).collect()
    }

    /// Moiety totals G * amounts.
    fn csum_now(&mut self) -> Result<Vec<f64>, String> {
        let a = self.amounts()?;
        let c = self.core()?;
        let (ls, ni) = (&c.ls(), c.ls().ni);
        Ok(c.csum_dep().iter().map(|&d| a[ls.order[ni + d]] - (0..ni).map(|b| ls.l0[d * ni + b] * a[ls.order[b]]).sum::<f64>()).collect())
    }

    /// Dependent species from the totals and the independent species (moiety analysis on).
    fn fix_dependents(&mut self) -> Result<(), String> {
        let a = self.amounts()?;
        let c = self.core_mut()?;
        let ni = c.ls().ni;
        for (q, d) in c.csum_dep().to_vec().into_iter().enumerate() {
            let s = c.ls_rows[c.ls().order[ni + d]].clone();
            let amt = c.csum[q] + (0..ni).map(|b| c.ls().l0[d * ni + b] * a[c.ls().order[b]]).sum::<f64>();
            let Some(&j) = c.kpos.get(&s) else { continue };
            let (_, vol) = c.species_raw(&s)?;
            c.x[j] = amt / vol;
        }
        c.touch();
        Ok(())
    }

    // ------------------------------------------------------------ get / set

    pub fn time(&self) -> f64 { self.c.as_ref().map_or(0.0, |c| c.t) }

    pub fn get_value(&mut self, s: &str) -> Result<f64, String> {
        let sel = sel::parse(s)?;
        if let Sel::Id(k) = &sel {
            if let Some(k) = self.csum_index(k) { return Ok(self.core()?.csum[k]); }
        }
        let deps = if self.cma { self.dependent_ids() } else { vec![] };
        let c = self.core_mut()?;
        match sel {
            Sel::Id(k) if c.is_species(&k) => c.species(&k, false),
            Sel::Id(k) => { let r = c.reg(&k)?; Ok(c.out(&[r])[0]) }
            Sel::Conc(k) if c.is_species(&k) => c.species(&k, true),
            Sel::Init(k, conc) => self.init_value(&k, conc),
            Sel::Rate(k, conc) => self.rate_of(&k, conc),
            Sel::Eigen(k, part) => {
                let i = self.floating_ids().iter().position(|x| *x == k).ok_or_else(|| format!("'{k}' is not a floating species"))?;
                let e = self.full_eigenvalues()?[i];
                Ok(if part == 2 { e.1 } else { e.0 })
            }
            Sel::Cc(a, b, scaled) => if scaled { self.get_cc(&a, &b) } else { self.get_ucc(&a, &b) },
            Sel::Ec(a, b, scaled) => if scaled { self.get_ee(&a, &b, false) } else { self.get_uee(&a, &b, false) },
            Sel::Stoich(sp, _) if deps.contains(&sp) => {
                Err(format!("'{sp}' is a dependent species under conserved moiety analysis (as roadrunner)"))
            }
            Sel::Stoich(sp, r) => {
                let i = c.info.rows.iter().position(|x| *x == sp).ok_or_else(|| format!("'{sp}' is not a floating species"))?;
                let j = c.info.ids[4].iter().position(|x| *x == r).ok_or_else(|| format!("'{r}' is not a reaction"))?;
                let n = c.info.n[i * c.info.ids[4].len() + j];
                if n == 0.0 { return Err(format!("'{sp}' does not take part in reaction '{r}'")); } // as roadrunner
                Ok(n)
            }
            _ => Err(format!("invalid selection '{s}'")),
        }
    }

    fn csum_index(&self, k: &str) -> Option<usize> {
        let c = self.c.as_ref()?;
        k.strip_prefix("_CSUM").and_then(|n| n.parse::<usize>().ok()).filter(|&n| self.cma && n < c.csum_dep().len() && !c.kpos.contains_key(k))
    }

    fn init_value(&mut self, k: &str, conc: bool) -> Result<f64, String> {
        let c = self.core_mut()?;
        if let Some(i) = c.fids().iter().position(|x| x == k) {
            let ic = c.init_concs_v(false)?[i];
            if conc { return Ok(ic); }
            let comp = c.comp_of(k).to_string();
            let vol = match c.init_volume(&comp) { Some(v) => v, None => c.species_raw(k)?.1 };
            return Ok(ic * vol); // the initial amount: in the initial volume (as roadrunner)
        }
        let Some(&j) = c.kpos.get(k) else { // rule-set symbols: their current value (as roadrunner)
            let r = c.reg(k)?;
            return Ok(c.out(&[r])[0]);
        };
        if c.is_species(k) { // boundary species
            let comp = c.comp_of(k).to_string();
            let vol = c.init_volume(&comp).unwrap_or(1.0);
            return Ok(match (conc, c.info.amount_params.contains(k)) { (true, true) => c.xi[j] / vol, (false, false) => c.xi[j] * vol, _ => c.xi[j] });
        }
        if conc { return Err(format!("'{k}' is not a species")); }
        Ok(c.xi[j])
    }

    fn rate_of(&mut self, k: &str, conc: bool) -> Result<f64, String> {
        if self.cma && self.dependent_ids().iter().any(|d| d == k) {
            return Err(format!("no rate for '{k}': a dependent species under conserved moiety analysis (as roadrunner)"));
        }
        let c = self.core_mut()?;
        let i = c.m.state_ids.iter().position(|x| x == k).ok_or_else(|| format!("'{k}' has no rate of change"))?;
        let is_sp = c.is_species(k);
        if conc && !is_sp { return Err(format!("'{k}' is not a species")); }
        let rates = c.state_rates();
        let r = rates[i];
        if !is_sp || conc || c.info.amount_params.contains(k) { return Ok(r); }
        // d(C V)/dt = V dC/dt + C dV/dt
        let (v, vol) = c.species_raw(k)?;
        let comp = c.comp_of(k).to_string();
        let dv = if let Some(&reg) = c.info.assignment_volume_rates.get(&comp) {
            c.out(&[reg])[0]
        } else {
            c.m.state_ids.iter().position(|x| *x == comp).map_or(0.0, |j| rates[j])
        };
        Ok(r * vol + v * dv)
    }

    pub fn set_value(&mut self, s: &str, v: f64) -> Result<(), String> {
        let sel = sel::parse(s)?;
        if let Sel::Id(k) = &sel {
            if let Some(n) = self.csum_index(k) {
                self.core_mut()?.csum[n] = v;
                return self.fix_dependents();
            }
        }
        let conc = matches!(sel, Sel::Conc(_));
        let c = self.core_mut()?;
        match sel {
            Sel::Id(k) if k == "time" => { c.t = v; c.touch(); Ok(()) }
            Sel::Id(k) | Sel::Conc(k) if !conc || c.is_species(&k) => {
                self.set_key(&k, v, conc)?;
                if self.cma && self.core()?.ls_rows.contains(&k) {
                    match self.dependent_ids().iter().position(|d| *d == k) {
                        Some(d) => if let Some(q) = self.core()?.csum_dep().iter().position(|&x| x == d) {
                            let t = self.csum_now()?[q];
                            self.core_mut()?.csum[q] = t;
                        },
                        None => self.fix_dependents()?,
                    }
                }
                Ok(())
            }
            Sel::Init(k, conc) => {
                self.set_init(&k, conc, v)?;
                self.reset_all() // roadrunner reloads its model at the initial values
            }
            _ => Err(format!("cannot set '{s}'")),
        }
    }

    /// Write the initial value of `k` (`conc`: a species' initial concentration) and nothing else, as
    /// ExecutableModel::setFloatingSpeciesInitConcentrations; its initial assignment is dropped.
    /// setValue("init(...)") follows it with reset_all(); the C API's initial-condition calls do not.
    pub fn set_init(&mut self, k: &str, conc: bool, v: f64) -> Result<(), String> {
        let c = self.core_mut()?;
        let j = *c.kpos.get(k).ok_or_else(|| format!("no initial value for '{k}'"))?;
        let x = if c.is_species(k) {
            let comp = c.comp_of(k).to_string();
            let vol = c.init_volume(&comp).unwrap_or(1.0); // an initial amount: in the initial volume
            match (conc, c.info.amount_params.contains(k)) { (true, true) => v * vol, (false, false) => v / vol, _ => v }
        } else if conc { return Err(format!("'{k}' is not a species")) } else { v };
        c.xi[j] = x;
        if let Some(ia) = c.ia.as_mut() { ia.retain(|(n, _)| n != k); }
        Ok(())
    }

    /// Set the engine value of key k (species: amount, or concentration with `conc`).
    fn set_key(&mut self, k: &str, v: f64, conc: bool) -> Result<(), String> {
        let cma = self.cma;
        let c = self.core_mut()?;
        if cma && c.ls_rows.iter().any(|s| s == k) && !c.info.rows.iter().any(|s| s == k) {
            return Err(format!("cannot set '{k}': a constant species in a conserved moiety (as roadrunner)"));
        }
        let j = *c.kpos.get(k).ok_or_else(|| format!("cannot set '{k}': not a parameter, compartment or species value"))?;
        if c.info.conversion_factors.iter().any(|f| f == k) {
            return Err(format!("cannot set '{k}': conversion factors are compiled into the model"));
        }
        if c.is_species(k) {
            let (_, vol) = c.species_raw(k)?;
            let cv = if conc { v } else { v / vol };
            c.x[j] = if c.info.amount_params.contains(k) { cv * vol } else { cv };
        } else if c.cconst.contains(k) && c.x[j] != v {
            if c.info.ia.iter().any(|(n, _)| n == k) { return Err(format!("cannot set '{k}': it has an initial assignment")); }
            // species amounts stay, concentrations change (as roadrunner); the volume is baked in: recompile
            let old = c.x[j];
            for s in c.fids().to_vec() {
                if c.comp_of(&s) == k { if let Some(&q) = c.kpos.get(&s) { c.x[q] *= old / v; } }
            }
            c.x[j] = v;
            let c = self.c.take();
            self.rebuild(c)?;
            return Ok(());
        } else {
            c.x[j] = v;
        }
        c.touch();
        Ok(())
    }

    /// Values of several selections.
    pub fn values(&mut self, sels: &[String]) -> Result<Vec<f64>, String> { sels.iter().map(|s| self.get_value(s)).collect() }
    pub fn set_values(&mut self, sels: &[String], v: &[f64]) -> Result<(), String> {
        if sels.len() != v.len() { return Err("selections and values differ in length".into()); }
        sels.iter().zip(v).try_for_each(|(s, &x)| self.set_value(s, x))
    }
    fn fmt(ids: &[String], f: &str) -> Vec<String> { ids.iter().map(|s| f.replace("{}", s)).collect() }
    pub fn floating_concentrations(&mut self) -> Result<Vec<f64>, String> { let s = Self::fmt(&self.floating_ids(), "[{}]"); self.values(&s) }
    pub fn floating_amounts(&mut self) -> Result<Vec<f64>, String> { let s = self.floating_ids(); self.values(&s) }
    pub fn boundary_concentrations(&mut self) -> Result<Vec<f64>, String> { let s = Self::fmt(&self.boundary_ids(), "[{}]"); self.values(&s) }
    pub fn boundary_amounts(&mut self) -> Result<Vec<f64>, String> { let s = self.boundary_ids(); self.values(&s) }
    pub fn global_parameter_values(&mut self) -> Result<Vec<f64>, String> { let s = self.global_parameter_ids(); self.values(&s) }
    pub fn compartment_volumes(&mut self) -> Result<Vec<f64>, String> { let s = self.compartment_ids(); self.values(&s) }
    pub fn reaction_rates(&mut self) -> Result<Vec<f64>, String> { let s = self.reaction_ids(); self.values(&s) }
    /// Amount rates of change of the floating species.
    pub fn rates_of_change(&mut self) -> Result<Vec<f64>, String> { let s = self.rates_of_change_ids(); self.values(&s) }
    pub fn floating_init_concentrations(&mut self) -> Result<Vec<f64>, String> { let s = self.floating_init_ids(); self.values(&s) }
    pub fn rate_rule_values(&mut self) -> Result<Vec<f64>, String> { let s = self.rate_rule_ids(); self.values(&s) }

    // ------------------------------------------------------------ reset

    /// Floating species and rate-rule targets to their initial values (initial assignments
    /// re-evaluated), time 0; parameters keep their values.
    pub fn reset(&mut self) -> Result<(), String> {
        self.core_mut()?.reset()?;
        if self.cma { let t = self.csum_now()?; self.core_mut()?.csum = t; }
        Ok(())
    }
    /// reset() plus parameters, compartments and boundary species to their initial values.
    pub fn reset_all(&mut self) -> Result<(), String> {
        let c = self.core_mut()?;
        let vols_changed = c.cconst.iter().any(|k| c.kpos.get(k).is_some_and(|&j| c.x[j] != c.xi[j]));
        c.x = c.xi.clone();
        if vols_changed { let c = self.c.take(); self.rebuild(c)?; }
        self.reset()
    }
    /// Everything as loaded (initial values set since are dropped).
    /// Initial values set through init(...) stay: roadrunner writes them into its model.
    pub fn reset_to_origin(&mut self) -> Result<(), String> { self.reset_all() }
    /// Global parameters to their initial values.
    pub fn reset_parameter(&mut self) -> Result<(), String> {
        let c = self.core_mut()?;
        for g in c.info.ids[2].clone() { if let Some(&j) = c.kpos.get(&g) { c.x[j] = c.xi[j]; } }
        c.touch();
        Ok(())
    }

    // ------------------------------------------------------------ selections

    fn default_selections(&mut self) {
        let Some(c) = &self.c else { return };
        let fs = self.floating_ids();
        let mut sel = vec!["time".to_string()];
        sel.extend(fs.iter().map(|s| format!("[{s}]")));
        sel.extend(c.info.rate_rule_ids.iter().filter(|r| !fs.contains(r)).cloned());
        self.sel = sel;
        self.ssel = fs.iter().map(|s| format!("[{s}]")).collect();
    }

    pub fn selections(&self) -> &[String] { &self.sel }
    pub fn steady_state_selections(&self) -> &[String] { &self.ssel }
    pub fn set_selections(&mut self, sel: &[String]) -> Result<(), String> { self.plan(sel)?; self.sel = sel.to_vec(); Ok(()) }
    pub fn set_steady_state_selections(&mut self, sel: &[String]) -> Result<(), String> { self.check(sel)?; self.ssel = sel.to_vec(); Ok(()) }

    /// Every symbol a selection names exists.
    fn check(&self, sel: &[String]) -> Result<(), String> {
        let c = self.core()?;
        let known = |k: &str| c.reg(k).is_ok() || self.csum_index(k).is_some();
        for s in sel {
            let ok = match sel::parse(s)? {
                Sel::Id(k) | Sel::Init(k, false) | Sel::Rate(k, false) => known(&k),
                Sel::Conc(k) | Sel::Init(k, true) | Sel::Rate(k, true) | Sel::Eigen(k, _) => c.is_species(&k),
                Sel::Cc(a, b, _) | Sel::Ec(a, b, _) | Sel::Stoich(a, b) => known(&a) && known(&b),
            };
            if !ok { return Err(format!("selection '{s}': unknown symbol")); }
        }
        Ok(())
    }

    /// Time-course plan: engine outputs + per-column operations (as _RsModel._cols).
    fn plan(&self, sel: &[String]) -> Result<Plan, String> {
        self.check(sel)?;
        let c = self.core()?;
        let mut outs: Vec<u32> = vec![];
        let mut o = |r: u32| match outs.iter().position(|&x| x == r) { Some(i) => i, None => { outs.push(r); outs.len() - 1 } };
        let mut cols = vec![];
        for s in sel {
            cols.push(match sel::parse(s)? {
                Sel::Id(k) | Sel::Conc(k) if c.is_species(&k) => {
                    let br = s.trim().starts_with('[');
                    let (i, v) = (o(c.reg(&k)?), o(c.reg(c.comp_of(&k))?));
                    match (br, c.info.amount_params.contains(&k)) { (true, true) => Col::Div(i, v), (false, false) => Col::Mul(i, v), _ => Col::Out(i) }
                }
                Sel::Id(k) if self.csum_index(&k).is_some() => Col::Pre(s.clone()),
                Sel::Id(k) => Col::Out(o(c.reg(&k)?)),
                Sel::Init(..) => Col::Pre(s.clone()),
                Sel::Rate(k, conc) if c.plain && c.fids().iter().position(|x| *x == k).is_some_and(|i| c.sp[i].is_some()) && !(self.cma && self.dependent_ids().contains(&k)) => {
                    let i = c.fids().iter().position(|x| *x == k).unwrap();
                    let st = c.sp[i].unwrap();
                    let vol = c.m.vol[st];
                    let nv = c.m.nv;
                    let mut t = vec![];
                    for (q, rid) in c.info.ids[4].iter().enumerate() {
                        let n = c.m.stoich[st * nv + c.rx[q]];
                        if n != 0.0 { t.push((o(c.reg(rid)?), n / if conc { vol } else { 1.0 })); }
                    }
                    Col::Lin(t)
                }
                _ => return Err(format!("selection '{s}' is not available in time courses")),
            });
        }
        Ok(Plan { outs, cols })
    }

    // ------------------------------------------------------------ integrators

    pub fn integrator(&self) -> &Solver { &self.integ }
    /// The current integrator's settings (Solver::set refuses values the engine can't honour).
    pub fn integrator_mut(&mut self) -> &mut Solver { &mut self.integ }
    pub fn set_integrator(&mut self, name: &str) -> Result<(), String> {
        if name == self.integ.name { return Ok(()); }
        let i = self.integs.iter().position(|s| s.name == name).ok_or_else(|| format!("no integrator '{name}' (one of {INTEGRATORS:?})"))?;
        let new = self.integs.remove(i);
        self.integs.push(std::mem::replace(&mut self.integ, new));
        Ok(())
    }
    /// Settings of any integrator (current or not).
    pub fn integrator_named(&mut self, name: &str) -> Option<&mut Solver> {
        if self.integ.name == name { Some(&mut self.integ) } else { self.integs.iter_mut().find(|s| s.name == name) }
    }
    pub fn steady_state_solver(&self) -> &Solver { &self.ss }
    pub fn steady_state_solver_mut(&mut self) -> &mut Solver { &mut self.ss }
    pub fn set_steady_state_solver(&mut self, name: &str) -> Result<(), String> {
        if name != self.ss.name {
            self.ss = Solver::steady_state(name).ok_or_else(|| format!("no steady-state solver '{name}' (one of {STEADY_STATE_SOLVERS:?})"))?;
        }
        Ok(())
    }

    /// libroadrunner's CVODE state order for tolerances: rate-rule variables, then the independent
    /// floating species (all floating species with moiety analysis off).
    fn tolerance_ids(&self) -> Vec<String> {
        let fs = self.floating_ids();
        let ni = if self.cma { self.independent_ids().len() } else { fs.len() };
        self.rate_rule_ids().into_iter().filter(|r| !fs.contains(r)).chain(fs[..ni].iter().cloned()).collect()
    }

    /// CVODEIntegrator::getAbsoluteToleranceVector: absolute tolerances on amounts, one per
    /// tolerance_ids(). A scalar absolute_tolerance is scaled by each initial amount (initial value of
    /// a rate-rule variable; the compartment volume if ~0); a vector set on the cvode integrator is
    /// used as is (one of the wrong length: 1e-12, scaled, as roadrunner's CVODE_MIN_ABSOLUTE).
    pub fn absolute_tolerance_vector(&mut self) -> Result<Vec<f64>, String> {
        let ids = self.tolerance_ids();
        let a = match self.integrator_named("cvode").unwrap().get("absolute_tolerance")? {
            Value::Doubles(v) if v.len() == ids.len() => return Ok(v.clone()),
            Value::Doubles(_) => 1e-12,
            v => v.f64(),
        };
        let c = self.core_mut()?;
        let asc = c.atol_scale()?;
        Ok(ids.iter().map(|id| c.m.state_ids.iter().position(|s| s == id).map_or(a, |j| a * asc[j] * c.m.vol[j])).collect())
    }

    /// CVODEIntegrator::setIndividualTolerance: absolute tolerance `value` (on the amount) for a
    /// floating species (independent, with moiety analysis on) or rate-rule variable; the cvode
    /// absolute_tolerance becomes the vector of absolute_tolerance_vector().
    pub fn set_individual_tolerance(&mut self, sid: &str, value: f64) -> Result<(), String> {
        if self.integ.name != "cvode" { return Err(format!("{}: individual tolerances are a cvode setting", self.integ.name)); }
        let i = self.tolerance_ids().iter().position(|s| s == sid).ok_or_else(|| format!(
            "CVODEIntegrator::setIndividualTolerance failed, given sid {sid} is neither a floating species nor has a rate rule."))?;
        let mut v = self.absolute_tolerance_vector()?;
        v[i] = value;
        self.integ.set("absolute_tolerance", Value::Doubles(v))
    }

    /// Absolute tolerance per engine state for cvode: absolute_tolerance (scalar) times atol_scale, or
    /// the vector's amount tolerances in state units (states it does not list: 1e-12 scaled).
    fn cvode_tolerances(&mut self) -> Result<Vec<f64>, String> {
        let at = self.integrator_named("cvode").unwrap().get("absolute_tolerance")?.clone();
        let (ids, v) = match at {
            Value::Doubles(_) => (self.tolerance_ids(), self.absolute_tolerance_vector()?),
            _ => (vec![], vec![]),
        };
        let c = self.core_mut()?;
        let asc = c.atol_scale()?;
        let a = if let Value::Doubles(_) = at { 1e-12 } else { at.f64() };
        let mut tol: Vec<f64> = asc.iter().map(|s| a * s).collect();
        for (id, x) in ids.iter().zip(v) {
            if let Some(j) = c.m.state_ids.iter().position(|s| s == id) { tol[j] = x / c.m.vol[j]; }
        }
        if tol.iter().any(|t| !(*t > 0.0 && t.is_finite())) { return Err("absolute tolerances must be positive".into()); }
        Ok(tol)
    }

    /// The cvode settings the engine honours; maximum_num_steps 0: CVODE's default 20000 (as roadrunner).
    fn cvode_ctl(&self) -> Ctl {
        let c = self.integ.ctl();
        Ctl { max_steps: if c.max_steps > 0 { c.max_steps } else { 20000 }, ..c }
    }

    pub fn seed(&mut self) -> u64 { self.integrator_named("gillespie").map_or(0, |s| s.num("seed") as u64) }
    pub fn set_seed(&mut self, seed: u64) -> Result<(), String> {
        self.integrator_named("gillespie").unwrap().set("seed", Value::UInt(seed))?;
        self.seed_ctr = (seed, 0); // a seed set again restarts its stream (roadrunner re-seeds)
        Ok(())
    }

    // ------------------------------------------------------------ time courses

    fn grid(&self) -> Result<Vec<f64>, String> {
        if let Some(t) = &self.times {
            if t.len() < 2 || t.windows(2).any(|w| !(w[1] >= w[0])) { return Err("times must be increasing, at least 2".into()); }
            return Ok(t.clone());
        }
        if self.points < 2 || !(self.end > self.start) { return Err("need end > start and at least 2 points".into()); }
        let n = self.points - 1;
        Ok((0..=n).map(|i| self.start + (self.end - self.start) * i as f64 / n as f64).collect())
    }

    /// Time course with the current integrator and selections; the model is left at the end.
    pub fn simulate(&mut self) -> Result<&Matrix, String> {
        let times = self.grid()?;
        let sel = self.sel.clone();
        let plan = self.plan(&sel)?;
        let pre: Vec<f64> = plan.cols.iter().filter_map(|c| if let Col::Pre(s) = c { Some(s.clone()) } else { None })
            .collect::<Vec<_>>().iter().map(|s| self.get_value(s)).collect::<Result<_, _>>()?;
        let y = match self.integ.name {
            "cvode" => {
                let (rtol, tol, ctl) = (self.integ.num("relative_tolerance"), self.cvode_tolerances()?, self.cvode_ctl());
                // An explicit sampling grid takes precedence over the output mode.
                let ctl = Ctl { var: ctl.var && self.times.is_none(), ..ctl };
                let ts = if ctl.var { vec![times[0], times[times.len() - 1]] } else { times };
                self.core_mut()?.advance(&ts, &plan.outs, rtol, &tol, ctl)?
            }
            "gillespie" => {
                let var = self.integ.get("variable_step_size")?.bool();
                let (y, _) = self.ssa(&plan.outs, if var { None } else { Some(&times) })?;
                y
            }
            "rk45" if self.integ.get("variable_step_size")?.bool() => {
                let i = &self.integ;
                let (hmin, hmax, eps, rows) = (i.num("minimum_time_step"), i.num("maximum_time_step"), i.num("epsilon"), i.num("max_output_rows").max(0.0) as usize);
                let (t0, t1) = (times[0], times[times.len() - 1]);
                self.core_mut()?.advance_rk45_variable(t0, t1, &plan.outs, hmin, hmax, eps, rows)?
            }
            name => {
                use crate::fixed::Method;
                let h = times[1] - times[0];
                let (method, step, atol) = match name {
                    "euler" => (Method::Euler, h / self.integ.num("subdivision_steps").max(1.0), 0.0),
                    "rk4" => (Method::Rk4, 0.0, 0.0),
                    "rk45" => (Method::Rk45, self.integ.num("maximum_time_step"), self.integ.num("epsilon")),
                    _ => return Err(format!("integrator {name} cannot simulate")),
                };
                self.core_mut()?.advance_fixed(&times, &plan.outs, method, step, atol)?
            }
        };
        self.result = Some(Matrix { nr: 0, nc: sel.len(), data: plan.apply(&y, &pre), rows: vec![], cols: sel });
        let r = self.result.as_mut().unwrap();
        r.nr = r.data.len() / r.nc.max(1);
        Ok(r)
    }

    /// The last simulate()/gillespie() result.
    pub fn result(&self) -> Option<&Matrix> { self.result.as_ref() }

    /// One integration step from t0 to t0 + h with the current integrator (as Integrator::integrate);
    /// returns the new time. Variable-step integrators may end before t0 + h.
    pub fn one_step(&mut self, t0: f64, h: f64) -> Result<f64, String> {
        use crate::fixed::Method;
        let var = self.integ.get("variable_step_size").is_ok_and(Value::bool);
        let i = &self.integ;
        let (hmin, hmax, eps) = (i.num("minimum_time_step"), i.num("maximum_time_step"), i.num("epsilon"));
        let (method, step) = match self.integ.name {
            "cvode" => {
                let (rtol, tol, ctl) = (self.integ.num("relative_tolerance"), self.cvode_tolerances()?, self.cvode_ctl());
                self.core_mut()?.advance(&[t0, t0 + h], &[], rtol, &tol, Ctl { max_rows: 2, ..ctl })?;
                return Ok(self.time());
            }
            "gillespie" => {
                self.core_mut()?.t = t0;
                if !var { return Ok(self.ssa(&[], Some(&[t0, t0 + h]))?.1); }
                let (s0, e0) = (self.start, self.end);
                (self.start, self.end) = (t0, t0 + h);
                let g = self.integrator_named("gillespie").unwrap();
                let rows = g.get("max_output_rows")?.clone();
                g.set("max_output_rows", Value::Int(2))?;
                let r = self.ssa(&[], None);
                self.integrator_named("gillespie").unwrap().set("max_output_rows", rows)?;
                (self.start, self.end) = (s0, e0);
                return Ok(r?.1);
            }
            "rk45" if var => {
                self.core_mut()?.advance_rk45_variable(t0, t0 + h, &[], hmin, hmax, eps, 2)?;
                return Ok(self.time());
            }
            "rk45" => (Method::Rk45, hmax),
            "rk4" => (Method::Rk4, 0.0),
            "euler" => (Method::Euler, h / self.integ.num("subdivision_steps").max(1.0)),
            name => return Err(format!("oneStep: integrator {name}")),
        };
        self.core_mut()?.advance_fixed(&[t0, t0 + h], &[], method, step, eps.max(0.0))?;
        Ok(t0 + h)
    }

    /// Next run's seed from the gillespie seed: successive runs draw from one seeded stream.
    fn run_seed(&mut self) -> u64 {
        let seed = self.seed();
        if self.seed_ctr.0 != seed { self.seed_ctr = (seed, 0); }
        self.seed_ctr.1 += 1;
        seed.wrapping_add(self.seed_ctr.1.wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }

    /// Gillespie from the current state: on `grid`, or a row per firing (variable steps, from
    /// start to end). Returns (rows of outs, end time); the model is left at the end.
    fn ssa(&mut self, outs: &[u32], grid: Option<&[f64]>) -> Result<(Vec<f64>, f64), String> {
        let seed = self.run_seed();
        let ctl = self.integrator_named("gillespie").unwrap().ctl();
        let max_rows = if ctl.max_rows > 0 { ctl.max_rows } else { usize::MAX };
        let (start, end) = (self.start, self.end);
        // variable steps with minimum_time_step: libroadrunner's rows every minimum_time_step
        let hgrid: Option<Vec<f64>> = (grid.is_none() && ctl.hmin > 0.0).then(|| {
            let mut g: Vec<f64> = (0..).map(|k| start + k as f64 * ctl.hmin).take_while(|&t| t < end).take(max_rows).collect();
            if g.len() < max_rows { g.push(end); }
            g
        });
        let grid = grid.or(hgrid.as_deref());
        let c = self.core_mut()?;
        if c.m.has_algebraic() { return Err("gillespie: algebraic rules are not supported".into()); }
        let all: Vec<u32> = outs.iter().chain(&c.xo).copied().collect();
        let no = all.len();
        c.m.load(&mut c.w, &c.kres, &c.x);
        c.w.ctl = ctl;
        let r = match grid {
            Some(t) => {
                let mut y = vec![0.0; t.len() * no];
                c.m.gillespie(&mut c.w, t, seed, &all, &mut y).then(|| (y, t[t.len() - 1], no))
            }
            None => c.m.gillespie_steps(&mut c.w, start, end, seed, &all, max_rows).map(|r| {
                let cols = no + 1;
                let last = r[r.len() - cols];
                (r, last, cols)
            }),
        };
        c.w.ctl = Ctl::default();
        let (y, tend, cols) = r.ok_or_else(|| if ctl.max_steps > 0 { format!("gillespie failed (maximum_num_steps {} reached?)", ctl.max_steps) } else { "gillespie failed".into() })?;
        let nt = y.len() / cols;
        let off = cols - no; // variable steps: time first
        let last = &y[(nt - 1) * cols + off..nt * cols];
        let nk = c.keys.len();
        c.x.copy_from_slice(&last[outs.len()..outs.len() + nk]);
        c.t = tend;
        c.fresh = true;
        c.touch();
        let rows: Vec<f64> = (0..nt).flat_map(|r| y[r * cols + off..r * cols + off + outs.len()].to_vec()).collect();
        if rows.iter().any(|v| !v.is_finite()) { return Err("gillespie went non-finite".into()); }
        Ok((rows, c.t))
    }

    /// Gillespie time course with the current selections: on `points` from start to end, or
    /// (None) a row per reaction firing. Seed: the gillespie integrator's.
    pub fn gillespie(&mut self, start: f64, end: f64, points: Option<usize>) -> Result<&Matrix, String> {
        let (s0, e0) = (self.start, self.end);
        (self.start, self.end) = (start, end);
        let r = (|| {
            let sel = self.sel.clone();
            let plan = self.plan(&sel)?;
            let pre: Vec<f64> = plan.cols.iter().filter_map(|c| if let Col::Pre(s) = c { Some(s.clone()) } else { None })
                .collect::<Vec<_>>().iter().map(|s| self.get_value(s)).collect::<Result<_, _>>()?;
            let grid = match points {
                Some(n) if n >= 2 => Some((0..n).map(|i| start + (end - start) * i as f64 / (n - 1) as f64).collect::<Vec<_>>()),
                Some(_) => return Err("need at least 2 points".to_string()),
                None => None,
            };
            let (y, _) = self.ssa(&plan.outs, grid.as_deref())?;
            let data = plan.apply(&y, &pre);
            Ok(Matrix { nr: data.len() / sel.len().max(1), nc: sel.len(), data, rows: vec![], cols: sel })
        })();
        (self.start, self.end) = (s0, e0);
        self.result = Some(r?);
        Ok(self.result.as_ref().unwrap())
    }

    /// n Gillespie runs from the current state on `points` from start to end (all CPU cores), with
    /// the current selections -> (n x points x selections, flat). Run i uses seed + i; the model is
    /// not changed. Models without reactions run their ODE (events, distrib draws) once per seed.
    pub fn gillespie_runs(&mut self, n: usize, start: f64, end: f64, points: usize, seed: Option<u64>) -> Result<Vec<f64>, String> {
        if points < 2 { return Err("need at least 2 points".into()); }
        let seed = match seed { Some(s) => s, None => self.run_seed() };
        let times: Vec<f64> = (0..points).map(|i| start + (end - start) * i as f64 / (points - 1) as f64).collect();
        let sel = self.sel.clone();
        let plan = self.plan(&sel)?;
        let pre: Vec<f64> = plan.cols.iter().filter_map(|c| if let Col::Pre(s) = c { Some(s.clone()) } else { None })
            .collect::<Vec<_>>().iter().map(|s| self.get_value(s)).collect::<Result<_, _>>()?;
        let rtol = self.integrator_named("cvode").unwrap().num("relative_tolerance");
        let (atol, asc) = (1.0, self.cvode_tolerances()?); // per-state absolute tolerances
        let g = self.integrator_named("gillespie").unwrap().ctl();
        let ctl = Ctl { var: false, ..g }; // ensembles always sample a fixed grid
        let c = self.core_mut()?;
        if c.m.has_algebraic() { return Err("gillespie: algebraic rules are not supported".into()); }
        let (m, kres, x, outs) = (&c.m, &c.kres, &c.x, &plan.outs);
        let no = outs.len();
        let mut y = vec![f64::NAN; n * points * no];
        let ode = m.nr == 0;
        let rows = Rows::new(&mut y, points * no);
        m.par(n, |w, i| {
            let r = unsafe { rows.row(i) }; // SAFETY: one row per run
            m.load(w, kres, x);
            let s = seed.wrapping_add(i as u64);
            let ok = if ode { w.atol_scale.copy_from_slice(&asc); m.reseed(w, s); m.simulate(w, &times, rtol, atol, outs, r) }
                     else { w.ctl = ctl; m.gillespie(w, &times, s, outs, r) };
            if !ok { r.fill(f64::NAN); }
        });
        let per = points * no;
        Ok((0..n).flat_map(|i| plan.apply(&y[i * per..(i + 1) * per], &pre)).collect())
    }

    /// Mean and standard deviation (n - 1) over n Gillespie runs on the grid.
    pub fn gillespie_mean_sd(&mut self, n: usize, start: f64, end: f64, points: usize) -> Result<(Matrix, Matrix), String> {
        if n < 2 { return Err("need at least 2 simulations".into()); }
        let y = self.gillespie_runs(n, start, end, points, None)?;
        let nc = self.sel.len();
        let per = points * nc;
        let mut mean = vec![0.0; per];
        for i in 0..n { for (a, b) in mean.iter_mut().zip(&y[i * per..(i + 1) * per]) { *a += b / n as f64; } }
        let mut sd = vec![0.0; per];
        for i in 0..n { for (k, b) in sd.iter_mut().enumerate() { *b += (y[i * per + k] - mean[k]).powi(2) / (n - 1) as f64; } }
        sd.iter_mut().for_each(|v| *v = v.sqrt());
        let mk = |d| Matrix { nr: points, nc, data: d, rows: vec![], cols: self.sel.clone() };
        Ok((mk(mean), mk(sd)))
    }

    // ------------------------------------------------------------ steady state

    /// Steady state from the current state (the model is left there); returns |d[S]/dt| of the
    /// floating species.
    pub fn steady_state(&mut self) -> Result<f64, String> {
        let cfg = ss::SsCfg::of(&self.ss)?;
        let c = self.core_mut()?;
        if let Some(tol) = c.solve(&cfg)? { return Ok(tol); }
        let r = c.state_rates();
        Ok(c.sp.iter().flatten().map(|&j| r[j] * r[j]).sum::<f64>().sqrt())
    }

    /// Steady state, then the steady-state selections.
    pub fn steady_state_values(&mut self) -> Result<Vec<f64>, String> {
        self.steady_state()?;
        let s = self.ssel.clone();
        self.values(&s)
    }

    // ------------------------------------------------------------ structure

    fn pick(&self, ids: &[String]) -> Result<Vec<usize>, String> {
        let c = self.core()?;
        Ok(ids.iter().map(|s| c.ls_rows.iter().position(|r| r == s).unwrap()).collect())
    }

    /// Full stoichiometry (floating species x reactions; moiety analysis on: libstruct's order).
    pub fn stoichiometry_matrix(&self) -> Result<Matrix, String> {
        let c = self.core()?;
        let nr = c.info.ids[4].len();
        if !self.cma { return Ok(Matrix::new(c.ls_rows.clone(), c.info.ids[4].clone(), c.ls_n.clone())); }
        let rows: Vec<String> = c.ls().order.iter().map(|&i| c.ls_rows[i].clone()).collect();
        let data = self.pick(&rows)?.iter().flat_map(|&i| c.ls_n[i * nr..(i + 1) * nr].to_vec()).collect();
        Ok(Matrix::new(rows, c.info.ids[4].clone(), data))
    }
    /// Nr: rows of the independent species.
    pub fn reduced_stoichiometry_matrix(&self) -> Result<Matrix, String> {
        let c = self.core()?;
        let rows = self.independent_ids();
        let nr = c.info.ids[4].len();
        let data = self.pick(&rows)?.iter().flat_map(|&i| c.ls_n[i * nr..(i + 1) * nr].to_vec()).collect();
        Ok(Matrix::new(rows, c.info.ids[4].clone(), data))
    }
    /// L0: dependent x independent.
    pub fn l0_matrix(&self) -> Result<Matrix, String> {
        let c = self.core()?;
        Ok(Matrix::new(self.dependent_ids(), self.independent_ids(), c.ls().l0.clone()))
    }
    /// L = [I; L0]: all floating species (libstruct's order) x independent species.
    pub fn link_matrix(&self) -> Result<Matrix, String> {
        let c = self.core()?;
        let (ni, nd) = (c.ls().ni, c.ls().nd());
        let mut d = vec![0.0; (ni + nd) * ni];
        for b in 0..ni { d[b * ni + b] = 1.0; }
        d[ni * ni..].copy_from_slice(&c.ls().l0);
        let rows = self.independent_ids().into_iter().chain(self.dependent_ids()).collect();
        Ok(Matrix::new(rows, self.independent_ids(), d))
    }
    /// Gamma = [-L0 | I]: moiety totals x floating species (moiety analysis on only).
    pub fn conservation_matrix(&self) -> Result<Matrix, String> {
        if !self.cma { return Err("conservation matrix: conserved moiety analysis is off".into()); }
        let c = self.core()?;
        let (ni, nd) = (c.ls().ni, c.ls().nd());
        // every dependent species has a row (as roadrunner); those set by rules have no total and
        // are labelled by their id
        let mut d = vec![0.0; nd * (ni + nd)];
        for k in 0..nd {
            for b in 0..ni { let v = c.ls().l0[k * ni + b]; d[k * (ni + nd) + b] = if v != 0.0 { -v } else { 0.0 }; }
            d[k * (ni + nd) + ni + k] = 1.0;
        }
        let deps = self.dependent_ids();
        let rows = (0..nd).map(|k| match c.csum_dep().iter().position(|&x| x == k) { Some(q) => format!("_CSUM{q}"), None => deps[k].clone() }).collect();
        let cols = c.ls().order.iter().map(|&i| c.ls_rows[i].clone()).collect();
        Ok(Matrix::new(rows, cols, d))
    }

    // ------------------------------------------------------------ Jacobian, eigenvalues

    /// d(d[S]/dt)/d[S] over the floating species (current order) at the current state and time.
    pub fn full_jacobian(&mut self) -> Result<Matrix, String> {
        let ids = self.floating_ids();
        let c = self.core_mut()?;
        c.need_plain()?;
        let nx = c.m.nx;
        let mut j = vec![0.0; nx * nx];
        if !c.qok { c.m.load(&mut c.q, &c.kres, &c.x); c.qok = true; }
        c.m.jacobian(&mut c.q, c.t, &mut j);
        c.qok = false; // jacobian evaluates rates into q's registers only, but keep it simple
        let st: Vec<Option<usize>> = ids.iter().map(|s| c.sp[c.fids().iter().position(|x| x == s).unwrap()]).collect();
        let j = &j;
        let data = st.iter().flat_map(|&a| st.iter().map(move |&b| match (a, b) { (Some(a), Some(b)) => j[a * nx + b], _ => 0.0 })).collect();
        Ok(Matrix::new(ids.clone(), ids, data))
    }

    /// Jacobian of the independent species (moiety analysis on); the full Jacobian otherwise.
    pub fn reduced_jacobian(&mut self) -> Result<Matrix, String> {
        let j = self.full_jacobian()?;
        if !self.cma { return Ok(j); }
        let l = self.link_matrix()?; // rows = floating_ids order (cma), amounts
        let ind = self.independent_ids();
        let c = self.core_mut()?;
        let vols: Vec<f64> = j.rows.iter().map(|s| c.species_raw(s).map(|v| v.1)).collect::<Result<_, _>>()?;
        let (n, ni) = (j.nr, ind.len());
        let mut out = vec![0.0; ni * ni];
        for (a, _) in ind.iter().enumerate() {
            for b in 0..ni {
                // (J Lc)[ind a, b], Lc = L / V[row] * V[ind b]
                out[a * ni + b] = (0..n).map(|k| j.get(a, k) * l.get(k, b) / vols[k] * vols[b]).sum();
            }
        }
        Ok(Matrix::new(ind.clone(), ind, out))
    }

    pub fn full_eigenvalues(&mut self) -> Result<Vec<(f64, f64)>, String> { let j = self.full_jacobian()?; eig::eigenvalues(&j.data, j.nr) }
    pub fn reduced_eigenvalues(&mut self) -> Result<Vec<(f64, f64)>, String> { let j = self.reduced_jacobian()?; eig::eigenvalues(&j.data, j.nr) }

    // ------------------------------------------------------------ MCA

    /// Elasticities d v / d[S] (reactions x floating species) at the current state; scaled: * [S] / v.
    pub fn elasticity_matrix(&mut self, scaled: bool) -> Result<Matrix, String> {
        let ids = self.floating_ids();
        let c = self.core_mut()?;
        c.need_plain()?;
        let (v, d) = c.rate_derivs();
        let nz = c.m.nx + c.m.pnames.len();
        let st: Vec<Option<usize>> = ids.iter().map(|s| c.sp[c.fids().iter().position(|x| x == s).unwrap()]).collect();
        let s: Vec<f64> = st.iter().map(|&j| j.map_or(0.0, |j| c.x[c.kpos[&c.m.state_ids[j]]])).collect();
        let mut data = vec![];
        for &r in &c.rx {
            for (k, &j) in st.iter().enumerate() {
                let e = j.map_or(0.0, |j| d[r * nz + j]);
                data.push(if scaled && e != 0.0 { e * s[k] / v[r] } else { e }); // no dependence: 0 even at zero flux
            }
        }
        Ok(Matrix::new(c.info.ids[4].clone(), ids, data))
    }

    /// (unscaled CCC species x reactions, unscaled FCC, v, [S]) at steady state (model left there).
    fn mca(&mut self) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<String>), String> {
        let ids = self.floating_ids();
        let cfg = ss::SsCfg::of(&self.ss)?;
        let c = self.core_mut()?;
        c.need_mca()?;
        c.solve(&cfg)?;
        let (ns, nr) = (c.m.ns, c.m.nr);
        let (mut e, mut uc, mut uf, mut cc, mut fc, mut el) = (vec![0.0; nr * ns], vec![0.0; ns * nr], vec![0.0; nr * nr], vec![0.0; ns * nr], vec![0.0; nr * nr], vec![0.0; nr * ns]);
        if !c.m.mca(&mut c.q, &mut e, &mut uc, &mut uf, &mut cc, &mut fc, &mut el) { return Err("MCA: singular Jacobian".into()); }
        let st: Vec<Option<usize>> = ids.iter().map(|s| c.sp[c.fids().iter().position(|x| x == s).unwrap()]).collect();
        let v: Vec<f64> = c.rx.iter().map(|&r| c.q.v[r]).collect();
        let s: Vec<f64> = st.iter().map(|&j| j.map_or(1.0, |j| c.q.s[j])).collect();
        let (uc, uf) = (&uc, &uf);
        let ucc = st.iter().flat_map(|&i| c.rx.iter().map(move |&r| i.map_or(0.0, |i| uc[i * nr + r]))).collect();
        let ufc = c.rx.iter().flat_map(|&i| c.rx.iter().map(move |&r| uf[i * nr + r])).collect();
        Ok((ucc, ufc, v, s, ids))
    }

    /// Concentration control coefficients (floating species x reactions) at steady state.
    pub fn cc_matrix(&mut self, scaled: bool) -> Result<Matrix, String> {
        let (mut uc, _, v, s, ids) = self.mca()?;
        let nr = v.len();
        if scaled { for i in 0..s.len() { for j in 0..nr { if uc[i * nr + j] != 0.0 { uc[i * nr + j] *= v[j] / s[i]; } } } }
        Ok(Matrix::new(ids, self.reaction_ids(), uc))
    }

    /// Flux control coefficients (reactions x reactions) at steady state; scaled rows of fluxes
    /// below the flux threshold are 0.
    pub fn fcc_matrix(&mut self, scaled: bool) -> Result<Matrix, String> {
        let (_, mut uf, v, _, _) = self.mca()?;
        let nr = v.len();
        if scaled {
            for i in 0..nr { for j in 0..nr {
                uf[i * nr + j] = if v[i].abs() < self.flux_threshold { 0.0 } else { uf[i * nr + j] * v[j] / v[i] };
            } }
        }
        let r = self.reaction_ids();
        Ok(Matrix::new(r.clone(), r, uf))
    }

    /// Unscaled control coefficient d var / d par at steady state (5-point stencil of steady states,
    /// step diff_step * par), var a floating species or reaction, par a parameter or boundary species.
    pub fn get_ucc(&mut self, var: &str, par: &str) -> Result<f64, String> { Ok(self.ucc(var, par)?.0) }

    pub fn get_cc(&mut self, var: &str, par: &str) -> Result<f64, String> {
        let (u, p, now, flux) = self.ucc(var, par)?;
        if flux && now.abs() < self.flux_threshold { return Ok(0.0); }
        Ok(u * p / now)
    }

    fn ucc(&mut self, var: &str, par: &str) -> Result<(f64, f64, f64, bool), String> {
        let ds = self.diff_step;
        let cfg = ss::SsCfg::of(&self.ss)?;
        let c = self.core_mut()?;
        c.need_mca()?;
        let var = var.trim_matches(['[', ']']).to_string();
        let off = c.keys.len() - c.fids().iter().filter(|s| c.kpos.contains_key(*s)).count();
        let j = *c.kpos.get(par).filter(|&&j| j < off).ok_or_else(|| format!("cc: '{par}' is not a parameter or boundary species"))?;
        if c.cconst.contains(par) { return Err(format!("cc: compartment '{par}' is compiled into the model")); }
        let (flux, idx) = match c.fids().iter().position(|x| *x == var) {
            Some(i) => (false, c.sp[i].unwrap()),
            None => (true, c.rx[c.info.ids[4].iter().position(|x| *x == var).ok_or_else(|| format!("cc: '{var}' is not a floating species or reaction"))?]),
        };
        let p0 = c.x[j];
        let mut h = ds * p0;
        if h.abs() < 1e-12 { h = ds; }
        let mut f = [0.0; 4];
        let (x0, t0, fresh) = (c.x.clone(), c.t, c.fresh);
        for (k, dx) in [h, 2.0 * h, -h, -2.0 * h].into_iter().enumerate() {
            c.x = x0.clone();
            c.x[j] += dx;
            c.touch();
            let r = c.solve(&cfg).map_err(|e| format!("cc: no steady state at a stencil point ({e})"));
            if r.is_ok() { f[k] = if flux { c.q.v[idx] } else { c.q.s[idx] }; }
            if let Err(e) = r { (c.x, c.t, c.fresh) = (x0, t0, fresh); c.touch(); return Err(e); }
        }
        (c.x, c.t, c.fresh) = (x0, t0, fresh);
        c.touch();
        let now = if flux { let r = c.reg(&var)?; c.out(&[r])[0] } else { c.species(&var, true)? };
        Ok(((f[3] + 8.0 * f[0] - 8.0 * f[2] - f[1]) / (12.0 * h), p0, now, flux))
    }

    /// Unscaled elasticity d v_reaction / d x (x: floating or boundary species concentration, or
    /// parameter), exact, at the current state (after a steady state with `steady`).
    pub fn get_uee(&mut self, reaction: &str, x: &str, steady: bool) -> Result<f64, String> { Ok(self.uee(reaction, x, steady)?.0) }

    pub fn get_ee(&mut self, reaction: &str, x: &str, steady: bool) -> Result<f64, String> {
        let (u, v) = self.uee(reaction, x, steady)?;
        let x = x.trim_matches(['[', ']']);
        let xv = if self.core()?.is_species(x) { self.get_value(&format!("[{x}]"))? } else { self.get_value(x)? };
        Ok(u * xv / v)
    }

    fn uee(&mut self, reaction: &str, x: &str, steady: bool) -> Result<(f64, f64), String> {
        if steady { self.steady_state()?; }
        let c = self.core_mut()?;
        let x = x.trim_matches(['[', ']']);
        let r = c.rx[c.info.ids[4].iter().position(|q| q == reaction).ok_or_else(|| format!("'{reaction}' is not a reaction"))?];
        let (v, d) = c.rate_derivs();
        let nz = c.m.nx + c.m.pnames.len();
        let col = match (c.m.state_ids.iter().position(|q| q == x), c.m.pnames.iter().position(|q| q == x)) {
            (Some(j), _) => j,
            (None, Some(k)) => c.m.nx + k,
            _ => return Err(format!("ee: '{x}' is not a species or parameter of the engine model")),
        };
        let mut e = d[r * nz + col];
        if c.info.amount_params.contains(x) { e *= c.species_raw(x)?.1; } // engine value is an amount: d/d[S]
        Ok((e, v[r]))
    }
}


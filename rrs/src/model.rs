//! Reaction-network model: conservation analysis, damped Newton steady state, BDF and Rodas4 integrators
//! with root-located events, Gillespie SSA and per-sample MCA. All per-sample work happens in a
//! preallocated `Work`; nothing allocates in the inner loops.
//!
//! State x = [floating species (concentrations) | parameters promoted to state by rate rules].
//! Rate rules are pseudo-reactions: v = [reaction rates | rate-rule rhs], and the generalized
//! stoichiometry (nx x nv) carries a unit entry per rate rule, so steady state, conservation and
//! the Jacobian treat them uniformly. Every tape uses the register layout [x | params | time | ...].
use crate::expr::{children, parse, Dag, Op, Tape};
use rayon::prelude::*;
use std::cell::RefCell;
use std::sync::Mutex;
use std::collections::HashMap;

/// Rules and events on top of the reaction network. Formulas use species/parameter ids and `time`.
#[derive(Default)]
pub struct Rules {
    /// (target, d/dt formula). Target is a reaction-less species (formula in concentration/time)
    /// or a parameter, which then becomes ODE state.
    pub rate: Vec<(String, String)>,
    /// (id, formula): output-only variables.
    pub assign: Vec<(String, String)>,
    pub events: Vec<EventDef>,
    /// delay(expr, tau) calls: (parameter holding the delayed value, expr, tau, expr before the run).
    /// The parameter is refilled with expr at time t - tau before every evaluation, from the run's
    /// history; before the run starts, with the last formula at the initial state and time t - tau.
    /// Inner delays come first (an expr or tau may use an earlier delay's parameter).
    pub delays: Vec<(String, String, String, String)>,
    /// rate-rule targets that are algebraic: their formula is a residual held at 0 (mass 0), so the
    /// model is a DAE and time courses go through diffsol's BDF
    pub algebraic: Vec<String>,
}

/// State history of the current run for delay(): (t, state, d state/dt) at accepted steps, and
/// scratch registers of the delay tape.
#[derive(Default)]
pub struct Hist { t: Vec<f64>, y: Vec<f64>, f: Vec<f64>, rd: Mutex<Vec<f64>>, yt: Mutex<Vec<f64>> }

impl Hist {
    pub fn clear(&mut self) { self.t.clear(); self.y.clear(); self.f.clear(); }

    /// The recorded series (t, y, f), taken out (a reload clears them; a continued run puts them back).
    pub fn save(&mut self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (std::mem::take(&mut self.t), std::mem::take(&mut self.y), std::mem::take(&mut self.f))
    }
    pub fn restore(&mut self, s: (Vec<f64>, Vec<f64>, Vec<f64>)) { (self.t, self.y, self.f) = s; }
    /// The recorded series (t, y, f), borrowed.
    pub fn series(&self) -> [&[f64]; 3] { [&self.t, &self.y, &self.f] }

    pub(crate) fn snapshot(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (self.t.clone(), self.y.clone(), self.f.clone())
    }

    /// State at tp into out (Hermite between recorded points; before the first: the first; after the
    /// last: linear to (t_now, s_now)).
    fn state_at(&self, tp: f64, s_now: &[f64], t_now: f64, out: &mut [f64]) {
        let (n, nx) = (self.t.len(), s_now.len());
        if n == 0 || tp >= t_now { out.copy_from_slice(s_now); return; }
        if tp <= self.t[0] { out.copy_from_slice(&self.y[..nx]); return; }
        let last = n - 1;
        if tp >= self.t[last] {
            let (ta, ya) = (self.t[last], &self.y[last * nx..]);
            let th = if t_now > ta { (tp - ta) / (t_now - ta) } else { 1.0 };
            for i in 0..nx { out[i] = ya[i] + th * (s_now[i] - ya[i]); }
            return;
        }
        let i = self.t.partition_point(|&x| x <= tp) - 1; // t[i] <= tp < t[i+1]
        let (ta, tb) = (self.t[i], self.t[i + 1]);
        let (ya, yb) = (&self.y[i * nx..(i + 1) * nx], &self.y[(i + 1) * nx..(i + 2) * nx]);
        let (fa, fb) = (&self.f[i * nx..(i + 1) * nx], &self.f[(i + 1) * nx..(i + 2) * nx]);
        let h = tb - ta;
        let u = (tp - ta) / h;
        let (h00, h10, h01, h11) = ((1.0 + 2.0 * u) * (1.0 - u).powi(2), u * (1.0 - u).powi(2), u * u * (3.0 - 2.0 * u), u * u * (u - 1.0));
        for k in 0..nx { out[k] = h00 * ya[k] + h10 * h * fa[k] + h01 * yb[k] + h11 * h * fb[k]; }
    }
}

#[derive(Clone)]
pub struct EventDef {
    pub trigger: String,
    pub assign: Vec<(String, String)>, // (species/parameter id, formula); species in concentration
    pub delay: Option<String>,
    pub priority: Option<String>,
    pub initial_value: bool,
    pub persistent: bool,
    pub use_trigger_values: bool,
    /// per assignment (or empty): a formula the value is divided by at execution time (an amount set
    /// with trigger-time values becomes a concentration with the volume when the event executes)
    pub div: Vec<Option<String>>,
    /// (state, its compartment's volume formula): concentrations scaled by V_before / V_after when the
    /// event executes, so amounts survive a compartment the event resizes
    pub rescale: Vec<(String, String)>,
}

/// Event compiled to registers of the aux tape; assignment targets are registers (< nx: state).
struct Event { trig: u32, delay: Option<u32>, prio: Option<u32>, asg: Vec<(u32, u32)>, div: Vec<Option<u32>>, rescale: Vec<(u32, u32)>, init: bool, persistent: bool, uvftt: bool }

pub struct Model {
    pub event_src: Vec<EventDef>, // events as given (sensitivity models re-use them)
    pub species: Vec<String>,
    pub reactions: Vec<String>,
    pub pnames: Vec<String>,
    pub state_ids: Vec<String>,
    pub outputs: Vec<(String, u32)>, // output id -> aux register
    keymap: HashMap<String, Key>,
    outmap: HashMap<String, u32>,
    pub init: Vec<f64>,   // nx
    pub vol: Vec<f64>,    // nx (1 for promoted parameters)
    pub pvals: Vec<f64>,
    pub stoich: Vec<f64>, // generalized, nx x nv row-major
    pub ns: usize,
    pub nr: usize,
    pub nx: usize,
    pub nv: usize,
    rates: Tape,          // outs: v
    full: Tape,           // outs: v, jacobian entries, dv/dt entries
    aux: Tape,            // outs: assignment rules, rates, event expressions
    dtape: Tape,          // outs: per delay, tau then expr
    alg: Vec<bool>,       // per state: algebraic (mass 0)
    dly: Vec<(u32, u32, u32, u32)>, // per delay: (parameter index, tau, expr, expr-before-start registers) in dtape
    events: Vec<Event>,
    comps: Vec<u32>, // aux registers of the comparisons inside triggers: a step is cut where any flips
    col: Vec<Vec<(u32, f64)>>, // per reaction: (species, N_ij), for SSA
    jac: Vec<(u32, u32)>, // (rate j, state k) for each jacobian entry
    dt: Vec<u32>,         // reactions with explicit time dependence
    nz: Vec<(u32, u32, f64)>, // (state i, rate j, N_ij / V_i)
    jfull: Vec<(u32, u32, f64)>, // J = diag(1/V) N E as (dst, jac entry, coef)
    sp: SparseLu,
    // conservation: S_dep = sum_b lc[d] * x_b + tc_d (concentrations)
    pub indep: Vec<usize>,
    pub dep: Vec<usize>,
    pub l0: Vec<f64>, // ndep x ni amount link coefficients
    lc: Vec<Vec<(u32, f64)>>,
    nzr: Vec<(u32, u32, f64)>,   // reduced residual f_a = sum N_ij v_j / V_i, i = indep[a]
    jred: Vec<(u32, u32, f64)>,  // reduced jacobian (dst, jac entry, coef)
    mca_c: Option<Cons>, // species-only conservation for MCA when there are rate rules
}

struct Cons { indep: Vec<usize>, dep: Vec<usize>, l0: Vec<f64>, lc: Vec<Vec<(u32, f64)>>,
              nzr: Vec<(u32, u32, f64)>, jred: Vec<(u32, u32, f64)> }

/// Conservation over the first `rows` states (RREF of N^T, pivot columns independent) and the
/// reduced residual/Jacobian maps for the entries `nz` (derivatives w.r.t. states >= rows dropped).
fn conservation(stoich: &[f64], nx: usize, nv: usize, rows: usize, vol: &[f64], nz: &[(u32, u32, f64)],
                jac: &[(u32, u32)]) -> Cons {
    let tol = 1e-9 * stoich.iter().fold(1.0f64, |m, x| m.max(x.abs()));
    let used: Vec<bool> = (0..nv).map(|j| nz.iter().any(|e| e.1 as usize == j)).collect();
    let mut a: Vec<f64> = (0..nv * rows).map(|t| if used[t / rows] { stoich[(t % rows) * nv + t / rows] } else { 0.0 }).collect();
    let (mut indep, mut dep, mut row) = (vec![], vec![], 0);
    for c in 0..rows {
        let p = (row..nv).max_by(|&x, &y| a[x * rows + c].abs().total_cmp(&a[y * rows + c].abs()));
        match p {
            Some(p) if a[p * rows + c].abs() > tol => {
                for j in 0..rows { a.swap(row * rows + j, p * rows + j); }
                let inv = 1.0 / a[row * rows + c];
                for j in 0..rows { a[row * rows + j] *= inv; }
                for r in 0..nv {
                    let m = a[r * rows + c];
                    if r != row && m != 0.0 { for j in 0..rows { a[r * rows + j] -= m * a[row * rows + j]; } }
                }
                indep.push(c);
                row += 1;
            }
            _ => dep.push(c),
        }
    }
    let ni = indep.len();
    let mut l0 = vec![0.0; dep.len() * ni];
    for (di, &d) in dep.iter().enumerate() {
        for b in 0..ni {
            let x = a[b * rows + d];
            l0[di * ni + b] = if x.abs() < 1e-12 { 0.0 } else { x };
        }
    }
    let lc: Vec<Vec<(u32, f64)>> = dep.iter().enumerate().map(|(di, &d)| (0..ni)
        .filter(|&b| l0[di * ni + b] != 0.0)
        .map(|b| (b as u32, l0[di * ni + b] * vol[indep[b]] / vol[d])).collect()).collect();
    let mut pos = vec![None; nx];
    for (b, &i) in indep.iter().enumerate() { pos[i] = Some((true, b)); }
    for (d, &i) in dep.iter().enumerate() { pos[i] = Some((false, d)); }
    let (mut nzr, mut jred) = (vec![], vec![]);
    for &(i, j, c) in nz {
        let Some((true, a_)) = pos[i as usize] else { continue };
        nzr.push((a_ as u32, j, c));
        for (e, &(_, k)) in jac.iter().enumerate().filter(|(_, &(jj, _))| jj == j) {
            match pos[k as usize] {
                Some((true, b)) => jred.push((a_ as u32 * ni as u32 + b as u32, e as u32, c)),
                Some((false, d)) => for &(b, l) in &lc[d] { jred.push((a_ as u32 * ni as u32 + b, e as u32, c * l)) },
                None => {}
            }
        }
    }
    Cons { indep, dep, l0, lc, nzr, jred }
}

/// Newton settings (crate::rr maps the steady-state solver registry onto these).
#[derive(Debug, Clone)]
pub struct NewtonOpts {
    pub max_iter: usize,
    /// converged when the full step is <= rtol * |x| (and the residual small)
    pub rtol: f64,
    /// smallest line-search damping factor
    pub min_damping: f64,
    /// backtracking line search (false: full Newton steps)
    pub linesearch: bool,
    /// singular Jacobian: least-squares step instead of failing
    pub rank_deficient: bool,
    pub allow_negative: bool,
    /// > 0: also converged when |f|_inf <= fnorm_tol / when the step <= step_tol
    pub fnorm_tol: f64,
    pub step_tol: f64,
}
impl Default for NewtonOpts {
    fn default() -> Self {
        NewtonOpts { max_iter: 100, rtol: 1e-12, min_damping: 1e-10, linesearch: true, rank_deficient: false, allow_negative: false, fnorm_tol: 0.0, step_tol: 0.0 }
    }
}

#[derive(Clone, Copy)]
pub enum Key { Par(usize), Spec(usize) }

pub struct Work {
    rr: Vec<f64>, // registers of the rates-only tape
    rf: Vec<f64>, // registers of the full tape
    pub v: Vec<f64>,
    ev: Vec<f64>,
    pub s: Vec<f64>,
    tc: Vec<f64>,
    x: Vec<f64>, xt: Vec<f64>, dx: Vec<f64>, f: Vec<f64>,
    jr: Vec<f64>, piv: Vec<usize>,
    // rosenbrock
    yn: Vec<f64>, dy: Vec<f64>, fx: Vec<f64>, f0: Vec<f64>, k: [Vec<f64>; 6], jf: Vec<f64>, pivf: Vec<usize>,
    lo: Vec<(u32, u32, f64)>, up: Vec<(u32, u32, f64)>,
    // mca
    xm: Vec<f64>,
    // events / outputs / ssa
    ra: Vec<f64>, // registers of the aux tape
    y0: Vec<f64>, d1: Vec<f64>, d2: Vec<f64>, // last step start + dense output coefficients
    prev: Vec<bool>, // trigger value at the current time
    prevc: Vec<bool>, // value of each comparison inside the triggers at the current time
    pend: Vec<(f64, usize, Vec<f64>)>, // (execution time, event, values if fixed at trigger time)
    tmp: Vec<f64>,
    amt: Vec<f64>, // SSA amounts
    bdf: Bdf,
    pub hist: Hist,
    rng: Rng, // ties between simultaneous events of equal priority are broken at random (SBML L3)
    /// per-state multiplier of the absolute tolerance (roadrunner scales atol by initial values)
    pub atol_scale: Vec<f64>,
    /// solver counters since the last reset: [BDF steps, Rodas4 steps, rejected steps, BDF->Rodas4 switches]
    pub stats: [usize; 4],
    pub ctl: Ctl,
    /// variable-step output (ctl.var): (times, rows of outs)
    pub rows: (Vec<f64>, Vec<f64>),
}

/// libroadrunner integrator settings the solvers honour (0 / false: unset, as upstream). CVODE:
/// initial / minimum / maximum step size, maximum steps between outputs (mxstep), maximum BDF
/// order, a row per step (variable_step_size, at most max_rows rows), multi (multiple_steps: an
/// output time between steps reports the state at the end of the step past it). Gillespie: max_steps
/// (firings between outputs), nonneg (skip firings that would make an amount negative).
#[derive(Default, Clone, Copy, Debug)]
pub struct Ctl {
    pub h0: f64,
    pub hmin: f64,
    pub hmax: f64,
    pub max_steps: usize,
    pub max_order: usize,
    pub nonstiff: bool,
    pub var: bool,
    pub max_rows: usize,
    pub multi: bool,
    pub nonneg: bool,
}

impl Ctl {
    /// Step size h within [hmin, hmax].
    fn lim(&self, h: f64) -> f64 { let h = h.max(self.hmin); if self.hmax > 0.0 { h.min(self.hmax) } else { h } }
    // CV_ONE_STEP calls have their own step budget: each call takes at most one step.
    fn step_limit(&self) -> usize {
        if self.var || self.multi || self.max_steps == 0 { usize::MAX } else { self.max_steps }
    }
    /// Factor of a step size reduction from h, stopping at hmin; None (failure) if h is already there.
    fn shrink(&self, h: f64, fac: f64) -> Option<f64> { if h <= self.hmin { None } else { Some(fac.max(self.hmin / h)) } }
}

/// Variable-order (1-5) BDF state: scaled backward differences D (rows of nx), as scipy's BDF /
/// ode15s. order 0 = restart at the next step.
#[derive(Default)]
struct Bdf { d: Vec<f64>, order: usize, h: f64, t: f64, neq: usize, lu: bool, lufac: f64, crate_: f64, jev: Vec<f64>, sc: [Vec<f64>; 7],
             nlu: usize, jage: usize, defer: bool } // steps since LU / Jacobian; error failure in this step

const BDF_MAX: usize = 5;
// NDF coefficients (kappa = [0, -0.185, -1/9, -0.0823, -0.0415, 0]): gamma_k = sum_{j<=k} 1/j,
// alpha_k = (1 - kappa_k) gamma_k, error constant kappa_k gamma_k + 1/(k+1)
const BDF_GAMMA: [f64; 6] = [0.0, 1.0, 1.5, 1.8333333333333333, 2.0833333333333335, 2.283333333333333];
const BDF_ALPHA: [f64; 6] = [0.0, 1.185, 1.6666666666666667, 1.9842166666666667, 2.1697916666666672, 2.283333333333333];
const BDF_ERR: [f64; 6] = [1.0, 0.315, 0.16666666666666666, 0.09911666666666669, 0.11354166666666667, 0.16666666666666666];

/// R(order, factor) of scipy's BDF: re-interpolates differences to step size h * factor.
fn bdf_r(order: usize, factor: f64) -> [[f64; 6]; 6] {
    let mut r = [[0.0; 6]; 6];
    for j in 0..=order { r[0][j] = 1.0; }
    for i in 1..=order {
        for j in 0..=order {
            let m = if j == 0 { 0.0 } else { (i as f64 - 1.0 - factor * j as f64) / i as f64 };
            r[i][j] = r[i - 1][j] * m;
        }
    }
    r
}

/// D[:order+1] <- (R U)^T D[:order+1] for a step size change by `factor`.
fn bdf_change(d: &mut [f64], tmp: &mut [f64], nx: usize, order: usize, factor: f64) {
    let (r, u) = (bdf_r(order, factor), bdf_r(order, 1.0));
    let n = order + 1;
    let mut ru = [[0.0; 6]; 6];
    for i in 0..n { for j in 0..n { ru[i][j] = (0..n).map(|k| r[i][k] * u[k][j]).sum(); } }
    for c in 0..nx {
        for j in 0..n { tmp[j] = (0..n).map(|i| ru[i][j] * d[i * nx + c]).sum(); }
        for j in 0..n { d[j * nx + c] = tmp[j]; }
    }
}

#[inline]
fn rms(n: usize, f: impl Fn(usize) -> f64) -> f64 { ((0..n).map(|i| f(i).powi(2)).sum::<f64>() / n.max(1) as f64).sqrt() }

// ---------------------------------------------------------------- dense LU (row-major, in place)

fn lu(a: &mut [f64], n: usize, piv: &mut [usize], rel_tol: f64) -> bool {
    let tol = if rel_tol > 0.0 { rel_tol * a[..n * n].iter().fold(0.0f64, |m, x| m.max(x.abs())) } else { 0.0 };
    for c in 0..n {
        let mut p = c;
        let mut best = a[c * n + c].abs();
        for r in c + 1..n {
            let x = a[r * n + c].abs();
            if x > best { best = x; p = r; }
        }
        if !(best > tol) { return false; }
        piv[c] = p;
        if p != c { for j in 0..n { a.swap(c * n + j, p * n + j); } }
        let inv = 1.0 / a[c * n + c];
        let (top, bot) = a.split_at_mut((c + 1) * n);
        let prow = &top[c * n + c + 1..c * n + n];
        for r in 0..n - c - 1 {
            let row = &mut bot[r * n..r * n + n];
            let m = row[c] * inv;
            if m != 0.0 {
                row[c] = m;
                for (x, &y) in row[c + 1..].iter_mut().zip(prow) { *x -= m * y; }
            } else { row[c] = 0.0; }
        }
    }
    true
}

/// Solve A X = B for B row-major n x m (m right-hand sides), in place.
fn lu_solve(a: &[f64], n: usize, piv: &[usize], b: &mut [f64], m: usize) {
    for c in 0..n {
        if piv[c] != c { for j in 0..m { b.swap(c * m + j, piv[c] * m + j); } }
    }
    for i in 1..n {
        let (top, rest) = b.split_at_mut(i * m);
        let bi = &mut rest[..m];
        for k in 0..i {
            let l = a[i * n + k];
            if l != 0.0 { for (x, &y) in bi.iter_mut().zip(&top[k * m..k * m + m]) { *x -= l * y; } }
        }
    }
    for i in (0..n).rev() {
        let (top, rest) = b.split_at_mut((i + 1) * m);
        let bi = &mut top[i * m..];
        for k in i + 1..n {
            let u = a[i * n + k];
            if u != 0.0 { for (x, &y) in bi.iter_mut().zip(&rest[(k - i - 1) * m..(k - i) * m]) { *x -= u * y; } }
        }
        let d = 1.0 / a[i * n + i];
        for x in bi.iter_mut() { *x *= d; }
    }
}

/// Nonzeros of packed LU factors in solve order: `lo` = L (row ascending); `up` = U rows descending,
/// each row's off-diagonals then (i, i, 1/u_ii). Sparse Jacobians then solve in O(nnz).
fn lu_pattern(a: &[f64], n: usize, lo: &mut Vec<(u32, u32, f64)>, up: &mut Vec<(u32, u32, f64)>) {
    lo.clear();
    up.clear();
    for i in 0..n {
        for k in 0..i { let x = a[i * n + k]; if x != 0.0 { lo.push((i as u32, k as u32, x)); } }
    }
    for i in (0..n).rev() {
        for k in i + 1..n { let x = a[i * n + k]; if x != 0.0 { up.push((i as u32, k as u32, x)); } }
        up.push((i as u32, i as u32, 1.0 / a[i * n + i]));
    }
}

/// Static-pattern LU of I/(h gamma) - J (no pivoting): symbolic fill computed once per model.
#[derive(Default)]
struct SparseLu { lcol: Vec<Vec<u32>>, urow: Vec<Vec<u32>>, lpat: Vec<(u32, u32)>, upat: Vec<(u32, u32)>, pos: Vec<u32>, perm: Vec<u32> }

/// Minimum-degree elimination order of the symmetrized pattern (fill-reducing).
fn min_degree(n: usize, nz: &[bool]) -> Vec<u32> {
    let mut adj: Vec<std::collections::BTreeSet<usize>> = (0..n).map(|i| (0..n).filter(|&j| j != i && (nz[i * n + j] || nz[j * n + i])).collect()).collect();
    let mut done = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for _ in 0..n {
        let v = (0..n).filter(|&i| !done[i]).min_by_key(|&i| adj[i].len()).unwrap();
        done[v] = true;
        order.push(v as u32);
        let nb: Vec<usize> = std::mem::take(&mut adj[v]).into_iter().collect();
        for &a in &nb {
            adj[a].remove(&v);
            for &b in &nb { if a != b { adj[a].insert(b); } }
        }
    }
    order
}

impl SparseLu {
    /// Symbolic LU of P A P^T, P the minimum-degree order; indices below are permuted ones.
    fn new(n: usize, entries: impl Iterator<Item = usize>) -> SparseLu {
        let mut orig = vec![false; n * n];
        for d in entries { orig[d] = true; }
        let perm = min_degree(n, &orig);
        let mut nz = vec![false; n * n];
        for i in 0..n { for j in 0..n { nz[i * n + j] = i == j || orig[perm[i] as usize * n + perm[j] as usize]; } }
        let mut sp = SparseLu { perm, ..Default::default() };
        for k in 0..n {
            let rows: Vec<u32> = (k + 1..n).filter(|&i| nz[i * n + k]).map(|i| i as u32).collect();
            let cols: Vec<u32> = (k + 1..n).filter(|&j| nz[k * n + j]).map(|j| j as u32).collect();
            for &i in &rows { for &j in &cols { nz[i as usize * n + j as usize] = true; } }
            sp.lcol.push(rows);
            sp.urow.push(cols);
        }
        for i in 0..n { for k in 0..i { if nz[i * n + k] { sp.lpat.push((i as u32, k as u32)); } } }
        for i in (0..n).rev() {
            sp.upat.extend(sp.urow[i].iter().map(|&j| (i as u32, j)));
            sp.upat.push((i as u32, i as u32));
        }
        let pm = &sp.perm;
        sp.pos = (0..n * n).filter(|&q| nz[q]).map(|q| pm[q / n] * n as u32 + pm[q % n]).collect();
        sp
    }

    /// Factor in place and emit solve lists as `lu_pattern` does. False if a pivot fails the
    /// threshold test |a_kk| > 0.1 max|a_ik| (caller falls back to partial pivoting).
    /// Solve lists use original indices, so `sparse_solve` needs no permutation.
    fn factor(&self, a: &mut [f64], n: usize, lo: &mut Vec<(u32, u32, f64)>, up: &mut Vec<(u32, u32, f64)>) -> bool {
        let pm = &self.perm;
        let at = |i: usize, j: usize| pm[i] as usize * n + pm[j] as usize;
        for k in 0..n {
            let p = a[at(k, k)];
            let cmax = self.lcol[k].iter().fold(0.0f64, |m, &i| m.max(a[at(i as usize, k)].abs()));
            if !(p.abs() > 0.1 * cmax) { return false; }
            let inv = 1.0 / p;
            for &i in &self.lcol[k] {
                let i = i as usize;
                let l = a[at(i, k)] * inv;
                a[at(i, k)] = l;
                for &j in &self.urow[k] { a[at(i, j as usize)] -= l * a[at(k, j as usize)]; }
            }
        }
        lo.clear();
        up.clear();
        lo.extend(self.lpat.iter().map(|&(i, k)| (pm[i as usize], pm[k as usize], a[at(i as usize, k as usize)])));
        up.extend(self.upat.iter().map(|&(i, j)| {
            let x = a[at(i as usize, j as usize)];
            (pm[i as usize], pm[j as usize], if i == j { 1.0 / x } else { x })
        }));
        true
    }
}

fn sparse_solve(piv: &[usize], lo: &[(u32, u32, f64)], up: &[(u32, u32, f64)], b: &mut [f64]) {
    for (c, &p) in piv.iter().enumerate() { if p != c { b.swap(c, p); } }
    // each row accumulates in a register (entries come row by row, rows in solve order)
    let (mut r, mut acc) = (u32::MAX, 0.0);
    for &(i, k, l) in lo {
        if i != r { if r != u32::MAX { b[r as usize] = acc; } r = i; acc = b[i as usize]; }
        acc -= l * b[k as usize];
    }
    if r != u32::MAX { b[r as usize] = acc; }
    r = u32::MAX;
    for &(i, k, u) in up { // row: off-diagonals, then (i, i, 1/u_ii)
        if i != r { r = i; acc = b[i as usize]; }
        if i == k { b[i as usize] = acc * u } else { acc -= u * b[k as usize] }
    }
}

fn inf_norm(x: &[f64]) -> f64 { x.iter().fold(0.0f64, |m, v| m.max(v.abs())) }

/// simulate uses BDF below this relative tolerance, Rodas4 otherwise. BDF is 1.0-4.8x faster from
/// rtol 1e-5 down (roadrunner's default is 1e-6) with roadrunner-level error; Rodas4 wins on small
/// or event-heavy models at looser tolerances.
const BDF_BELOW_RTOL: f64 = 2e-5;
// Rodas4 (Hairer & Wanner), transformed coefficients.
const GAM: f64 = 0.25;
const A: [[f64; 5]; 5] = [
    [1.544, 0.0, 0.0, 0.0, 0.0],
    [0.9466785280815826, 0.2557011698983284, 0.0, 0.0, 0.0],
    [3.314825187068521, 2.896124015972201, 0.9986419139977817, 0.0, 0.0],
    [1.221224509226641, 6.019134481288629, 12.53708332932087, -0.6878860361058950, 0.0],
    [1.221224509226641, 6.019134481288629, 12.53708332932087, -0.6878860361058950, 1.0],
];
const C: [[f64; 5]; 5] = [
    [-5.6688, 0.0, 0.0, 0.0, 0.0],
    [-2.430093356833875, -0.2063599157091915, 0.0, 0.0, 0.0],
    [-0.1073529058151375, -9.594562251023355, -20.47028614809616, 0.0, 0.0],
    [7.496443313967647, -10.24680431464352, -33.99990352819905, 11.70890893206160, 0.0],
    [8.083246795921522, -7.981132988064893, -31.52159432874371, 16.31930543123136, -6.058818238834054],
];
const CT: [f64; 6] = [0.0, 0.386, 0.21, 0.63, 1.0, 1.0];
const D: [f64; 6] = [0.25, -0.1043, 0.1035, -0.0362, 0.0, 0.0];
// Rodas4 dense output: y(t + th h) = (1-th) y0 + th (y1 + (1-th) (d1 + th d2)), d = H k[0..5].
const H: [[f64; 5]; 2] = [
    [10.12623508344586, -7.487995877610167, -34.80091861555747, -7.992771707568823, 1.025137723295662],
    [-0.6762803392801253, 6.087714651680015, 16.43084320892478, 24.76722511418386, -6.594389125716872],
];
const MAX_STEPS: usize = 200_000;

/// x rounded to 12 significant digits (event times located by bisection snap to round times).
fn round12(x: f64) -> f64 { let sc = 10f64.powi(12 - x.abs().max(1e-300).log10().ceil() as i32); (x * sc).round() / sc }

/// xoshiro256++ seeded through splitmix64.
pub struct Rng([u64; 4]);
impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut z = seed;
        Rng(std::array::from_fn(|_| {
            z = z.wrapping_add(0x9E3779B97F4A7C15);
            let x = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            let x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
            x ^ (x >> 31)
        }))
    }
    fn next(&mut self) -> u64 {
        let s = &mut self.0;
        let r = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = s[3].rotate_left(45);
        r
    }
    /// Uniform in (0, 1].
    pub fn unif(&mut self) -> f64 { ((self.next() >> 11) + 1) as f64 * (1.0 / (1u64 << 53) as f64) }
    pub fn next_u64(&mut self) -> u64 { self.next() }
}

impl Model {
    #[allow(clippy::too_many_arguments)]
    pub fn new(species: Vec<String>, mut init: Vec<f64>, mut vol: Vec<f64>, reactions: Vec<String>,
               rates: Vec<String>, stoich: Vec<f64>, pnames: Vec<String>, pvals: Vec<f64>, x: Rules) -> Result<Model, String> {
        let (ns, nr, np) = (species.len(), reactions.len(), pnames.len());
        if init.len() != ns || vol.len() != ns { return Err("init/volumes length must equal number of species".into()); }
        if rates.len() != nr { return Err("rates length must equal number of reactions".into()); }
        if stoich.len() != ns * nr { return Err(format!("stoich must be {ns} x {nr}")); }
        if pvals.len() != np { return Err("param_names/param_values length mismatch".into()); }
        // rate-rule targets: reaction-less species rows, or parameters promoted to extra state
        let mut state_ids = species.clone();
        let mut rr_row = vec![];
        for (name, _) in &x.rate {
            let i = if let Some(i) = species.iter().position(|s| s == name) {
                if (0..nr).any(|j| stoich[i * nr + j] != 0.0) { return Err(format!("species '{name}' has a rate rule and reactions")); }
                i
            } else {
                let k = pnames.iter().position(|p| p == name).ok_or(format!("rate rule target '{name}' is not a species or parameter"))?;
                state_ids.push(name.clone()); init.push(pvals[k]); vol.push(1.0);
                state_ids.len() - 1
            };
            rr_row.push(i);
        }
        let (nx, nv) = (state_ids.len(), nr + x.rate.len());
        for a in &x.algebraic {
            if !x.rate.iter().any(|(n, _)| n == a) { return Err(format!("algebraic state '{a}' has no residual")); }
        }
        let alg: Vec<bool> = state_ids.iter().map(|s| x.algebraic.contains(s)).collect();
        let mut st = vec![0.0; nx * nv];
        for i in 0..ns { st[i * nv..i * nv + nr].copy_from_slice(&stoich[i * nr..(i + 1) * nr]); }
        for (q, &i) in rr_row.iter().enumerate() { st[i * nv + nr + q] = 1.0; }
        let stoich = st;

        let mut sym: HashMap<&str, Op> = HashMap::new();
        for (i, p) in pnames.iter().enumerate() { sym.insert(p, Op::Par(i as u32)); }
        for (i, s) in state_ids.iter().enumerate() {
            if sym.insert(s, Op::Var(i as u32)).is_some() && i < ns { return Err(format!("'{s}' is both a species and a parameter")); }
        }
        sym.entry("time").or_insert(Op::Time); // a model id named 'time' wins; the clock stays as __time__
        sym.insert("__time__", Op::Time);
        let mut dag = Dag::default();
        let p = |src: &str, what: &dyn std::fmt::Display, dag: &mut Dag| parse(src, dag, |n| sym.get(n).copied())
            .map_err(|e| format!("{what}: {e}"));
        let mut roots = Vec::with_capacity(nv);
        for (j, r) in rates.iter().enumerate() { roots.push(p(r, &format!("reaction {}", reactions[j]), &mut dag)?); }
        for (n, f) in &x.rate { roots.push(p(f, &format!("rate rule {n}"), &mut dag)?); }
        let (mut jac, mut dt, mut jroots, mut troots) = (vec![], vec![], vec![], vec![]);
        for (j, &r) in roots.iter().enumerate() {
            let mut has_t = false;
            for k in dag.deps(r, &mut has_t) {
                let d = dag.diff(r, Op::Var(k));
                if dag.val(d) != Some(0.0) { jac.push((j as u32, k)); jroots.push(d); }
            }
            if has_t {
                let d = dag.diff(r, Op::Time);
                if dag.val(d) != Some(0.0) { dt.push(j as u32); troots.push(d); }
            }
        }
        // delayed values: parameters refilled before every evaluation (never hoisted)
        let mut dly = vec![];
        let mut droots = vec![];
        for (pn, e, tau, pre) in &x.delays {
            let k = pnames.iter().position(|p| p == pn).ok_or(format!("delay parameter '{pn}' missing"))? as u32;
            droots.push(p(tau, &format!("delay time of {pn}"), &mut dag)?);
            droots.push(p(e, &format!("delayed expression of {pn}"), &mut dag)?);
            droots.push(p(pre, &format!("delayed expression of {pn}"), &mut dag)?);
            dly.push(k);
        }
        let dynp = dly.clone();
        let dtape = Tape::new_dyn(&dag, &droots, nx, np, &dynp);
        let o3 = |q: usize, i: usize| dtape.outs[3 * q + i];
        let dly: Vec<(u32, u32, u32, u32)> = dly.iter().enumerate().map(|(q, &k)| (k, o3(q, 0), o3(q, 1), o3(q, 2))).collect();
        let rates_tape = Tape::new_dyn(&dag, &roots, nx, np, &dynp);
        let all: Vec<_> = roots.iter().chain(&jroots).chain(&troots).copied().collect();
        let full = Tape::new_dyn(&dag, &all, nx, np, &dynp);

        // aux tape: [assignment rules | reaction rates | per event: trigger, delay?, priority?, values]
        let mut aroots = vec![];
        for (n, f) in &x.assign { aroots.push(p(f, &format!("assignment rule {n}"), &mut dag)?); }
        aroots.extend_from_slice(&roots[..nr]);
        let mut evs = vec![];
        for (e, d) in x.events.iter().enumerate() {
            let what = format!("event {e}");
            let at = |f: &str, dag: &mut Dag, aroots: &mut Vec<u32>| -> Result<u32, String> {
                aroots.push(p(f, &what, dag)?); Ok(aroots.len() as u32 - 1)
            };
            let trig = at(&d.trigger, &mut dag, &mut aroots)?;
            let delay = d.delay.as_deref().map(|f| at(f, &mut dag, &mut aroots)).transpose()?;
            let prio = d.priority.as_deref().map(|f| at(f, &mut dag, &mut aroots)).transpose()?;
            let mut asg = vec![];
            for (tg, f) in &d.assign {
                let r = match sym.get(tg.as_str()) {
                    Some(Op::Var(i)) => *i,
                    Some(Op::Par(k)) => (nx + *k as usize) as u32,
                    _ => return Err(format!("event {e}: unknown assignment target '{tg}'")),
                };
                asg.push((r, at(f, &mut dag, &mut aroots)?));
            }
            let div = (0..asg.len()).map(|q| d.div.get(q).cloned().flatten().map(|f| at(&f, &mut dag, &mut aroots)).transpose())
                .collect::<Result<Vec<_>, _>>()?;
            let mut rescale = vec![];
            for (sn, vf) in &d.rescale {
                let Some(Op::Var(i)) = sym.get(sn.as_str()) else { return Err(format!("event {e}: '{sn}' is not a state")) };
                rescale.push((*i, at(vf, &mut dag, &mut aroots)?));
            }
            evs.push(Event { trig, delay, prio, asg, div, rescale, init: d.initial_value, persistent: d.persistent, uvftt: d.use_trigger_values });
        }
        // every comparison inside a trigger is watched on its own (as CVODE root functions): a step
        // can cross a whole window (t >= a && t <= b) with the trigger false at both ends
        let mut cmp = vec![];
        let mut stack: Vec<u32> = evs.iter().map(|e: &Event| aroots[e.trig as usize]).collect();
        let mut seen = std::collections::HashSet::new();
        while let Some(nd) = stack.pop() {
            if !seen.insert(nd) { continue; }
            let op = dag.nodes[nd as usize];
            if matches!(op, Op::Lt(..) | Op::Le(..) | Op::Eq(..) | Op::Ne(..)) { cmp.push(aroots.len() as u32); aroots.push(nd); }
            stack.extend(children(&op).iter().flatten());
        }
        let aux = Tape::new_dyn(&dag, &aroots, nx, np, &dynp);
        let comps: Vec<u32> = cmp.iter().map(|&i| aux.outs[i as usize]).collect();
        let o = |i: u32| aux.outs[i as usize];
        for ev in &mut evs {
            ev.trig = o(ev.trig);
            ev.delay = ev.delay.map(o);
            ev.prio = ev.prio.map(o);
            for a in &mut ev.asg { a.1 = o(a.1); }
            for r in ev.div.iter_mut().flatten() { *r = o(*r); }
            for r in &mut ev.rescale { r.1 = o(r.1); }
        }
        let na = x.assign.len();
        let mut outputs: Vec<(String, u32)> = state_ids.iter().enumerate().map(|(i, s)| (s.clone(), i as u32)).collect();
        outputs.extend(x.assign.iter().enumerate().map(|(q, (n, _))| (n.clone(), o(q as u32))));
        outputs.extend(reactions.iter().enumerate().map(|(j, n)| (n.clone(), o((na + j) as u32))));
        outputs.extend(pnames.iter().enumerate().filter(|(_, n)| !state_ids.contains(n))
            .map(|(k, n)| (n.clone(), (nx + k) as u32)));
        if !outputs.iter().any(|o| o.0 == "time") { outputs.push(("time".into(), (nx + np) as u32)); }
        outputs.push(("__time__".into(), (nx + np) as u32));
        let col = (0..nr).map(|j| (0..ns).filter(|&i| stoich[i * nv + j] != 0.0)
            .map(|i| (i as u32, stoich[i * nv + j])).collect()).collect();

        let mut nz = vec![];
        for i in 0..nx { for j in 0..nv {
            let c = stoich[i * nv + j];
            if c != 0.0 { nz.push((i as u32, j as u32, if j < nr { c / vol[i] } else { c })); } // rate rules: already d(conc)/dt
        }}
        let by_rxn = |j: u32| jac.iter().enumerate().filter(move |(_, &(jj, _))| jj == j);
        let mut jfull = vec![];
        for &(i, j, c) in &nz {
            for (e, &(_, k)) in by_rxn(j) { jfull.push((i * nx as u32 + k, e as u32, c)); }
        }

        let Cons { indep, dep, l0, lc, nzr, jred } = conservation(&stoich, nx, nv, nx, &vol, &nz, &jac);
        // MCA as roadrunner: species x reactions only, rate-rule states held at their steady value
        let mca_c = (nv != nr).then(|| conservation(&stoich, nx, nv, ns, &vol,
            &nz.iter().filter(|e| (e.0 as usize) < ns && (e.1 as usize) < nr).copied().collect::<Vec<_>>(), &jac));
        let mut keymap: HashMap<String, Key> = pnames.iter().enumerate().map(|(i, p)| (p.clone(), Key::Par(i))).collect();
        keymap.extend(state_ids.iter().enumerate().map(|(i, p)| (p.clone(), Key::Spec(i))));
        let outmap = outputs.iter().cloned().collect();
        Ok(Model { event_src: x.events, species, reactions, pnames, state_ids, outputs, keymap, outmap, init, vol, pvals, stoich, ns, nr, nx, nv,
                   sp: SparseLu::new(nx, jfull.iter().map(|e| e.0 as usize)),
                   rates: rates_tape, full, aux, dtape, dly, alg, events: evs, comps, col, jac, dt, nz, jfull, indep, dep, l0, lc, nzr, jred, mca_c })
    }

    pub fn ni(&self) -> usize { self.indep.len() }

    /// Keys are state ids (initial value) or parameters.
    pub fn resolve(&self, keys: &[String]) -> Result<Vec<Key>, String> {
        keys.iter().map(|k| self.keymap.get(k).copied().ok_or(format!("unknown parameter or species '{k}'"))).collect()
    }

    /// Output ids -> aux registers.
    pub fn output_regs(&self, ids: &[String]) -> Result<Vec<u32>, String> {
        ids.iter().map(|k| self.outmap.get(k).copied().ok_or(format!("unknown output '{k}'"))).collect()
    }

    pub fn has_events(&self) -> bool { !self.events.is_empty() }

    /// Any algebraic state (DAE): the SSA does not support these.
    pub fn has_algebraic(&self) -> bool { self.alg.contains(&true) }

    pub fn work(&self) -> Work {
        let (nx, nv, nr, ni) = (self.nx, self.nv, self.nr, self.ni());
        let z = |n| vec![0.0; n];
        Work {
            rr: self.rates.regs(), rf: self.full.regs(),
            v: z(nv), ev: z(self.jac.len() + self.dt.len()), s: z(nx), tc: z(self.dep.len()),
            x: z(ni), xt: z(ni), dx: z(ni), f: z(ni), jr: z(ni * ni), piv: vec![0; ni.max(nx)],
            yn: z(nx), dy: z(nx), fx: z(nx), f0: z(nx), k: std::array::from_fn(|_| z(nx)), jf: z(nx * nx), pivf: vec![0; nx],
            lo: vec![], up: vec![], xm: z(ni * nr),
            ra: self.aux.regs(), y0: z(nx), d1: z(nx), d2: z(nx), prev: vec![false; self.events.len()], prevc: vec![false; self.comps.len()],
            pend: vec![], tmp: vec![], amt: z(nx), atol_scale: vec![1.0; nx], stats: [0; 4], ctl: Ctl::default(), rows: (vec![], vec![]),
            bdf: Bdf { d: z((BDF_MAX + 3) * nx), jev: z(self.jac.len() + self.dt.len()), sc: std::array::from_fn(|_| z(nx)), ..Default::default() },
            hist: Hist { rd: Mutex::new(self.dtape.regs()), yt: Mutex::new(z(nx)), ..Default::default() },
            rng: Rng::new(0x5eed_e7e7),
        }
    }

    /// Load base values + one sample row into the work registers and initial state.
    pub fn load(&self, w: &mut Work, keys: &[Key], row: &[f64]) {
        let (nx, np) = (self.nx, self.pnames.len());
        w.rr[nx..nx + np].copy_from_slice(&self.pvals);
        w.s.copy_from_slice(&self.init);
        for (k, &x) in keys.iter().zip(row) {
            match *k { Key::Par(i) => w.rr[nx + i] = x, Key::Spec(i) => w.s[i] = x }
        }
        w.hist.clear();
        let a = &w.rr[nx..nx + np];
        w.rf[nx..nx + np].copy_from_slice(a);
        w.ra[nx..nx + np].copy_from_slice(a);
        self.prep(w);
    }

    /// Re-evaluate the parameter-only parts of the tapes (after parameters change).
    fn prep(&self, w: &mut Work) {
        self.rates.prep(&mut w.rr);
        self.full.prep(&mut w.rf);
        self.aux.prep(&mut w.ra);
    }

    /// Seed the run's random stream (event ties, distrib draws): batches give each sample its own.
    pub fn reseed(&self, w: &mut Work, seed: u64) { w.rng = Rng::new(seed); }

    /// Current value of parameter k (events may have changed it).
    pub fn par_value(&self, w: &Work, k: usize) -> f64 { w.rr[self.nx + k] }

    // ------------------------------------------------------------ evaluation

    /// Rates only at (s, t) into v.
    #[inline]
    fn rates_at(&self, hs: &Hist, regs: &mut [f64], v: &mut [f64], s: &[f64], t: f64) {
        let (nx, np) = (self.nx, self.pnames.len());
        self.fill_delays(hs, regs, s, t);
        regs[..nx].copy_from_slice(s);
        regs[nx + np] = t;
        self.rates.run(regs);
        for (o, &r) in v.iter_mut().zip(&self.rates.outs) { *o = regs[r as usize]; }
    }

    /// Rates, jacobian entries and dv/dt at (s, t).
    #[inline]
    fn full_at(&self, hs: &Hist, regs: &mut [f64], v: &mut [f64], ev: &mut [f64], s: &[f64], t: f64) {
        let (nx, np, nv) = (self.nx, self.pnames.len(), self.nv);
        self.fill_delays(hs, regs, s, t);
        regs[..nx].copy_from_slice(s);
        regs[nx + np] = t;
        self.full.run(regs);
        for (o, &r) in v.iter_mut().zip(&self.full.outs[..nv]) { *o = regs[r as usize]; }
        for (o, &r) in ev.iter_mut().zip(&self.full.outs[nv..]) { *o = regs[r as usize]; }
    }

    /// Aux tape (rules, rates, event expressions) at (s, t).
    #[inline]
    fn aux_at(&self, hs: &Hist, regs: &mut [f64], s: &[f64], t: f64) {
        let (nx, np) = (self.nx, self.pnames.len());
        self.fill_delays(hs, regs, s, t);
        regs[..nx].copy_from_slice(s);
        regs[nx + np] = t;
        self.aux.run(regs);
    }

    /// delay(): write each delayed value (expr at t - tau, from the run's history) into its parameter.
    fn fill_delays(&self, hs: &Hist, regs: &mut [f64], s: &[f64], t: f64) {
        if self.dly.is_empty() { return; }
        let (nx, np) = (self.nx, self.pnames.len());
        let (mut rd, mut yt) = (hs.rd.lock().unwrap(), hs.yt.lock().unwrap());
        rd[nx..nx + np].copy_from_slice(&regs[nx..nx + np]);
        self.dtape.prep(&mut rd);
        for &(k, tr, er, pr) in &self.dly { // in order: later delays see the earlier ones' new values
            rd[..nx].copy_from_slice(s);
            rd[nx + np] = t;
            self.dtape.run(&mut rd);
            let tp = t - rd[tr as usize];
            let before = tp < hs.t.first().copied().unwrap_or(0.0);
            hs.state_at(tp, s, t, &mut yt);
            rd[..nx].copy_from_slice(&yt);
            rd[nx + np] = tp.min(t);
            self.dtape.run(&mut rd);
            let x = rd[if before { pr } else { er } as usize];
            regs[nx + k as usize] = x;
            rd[nx + k as usize] = x;
        }
    }

    /// Record (t, w.s, d w.s/dt) in the delay history (delay models only).
    fn hist_push(&self, w: &mut Work, t: f64) {
        if self.dly.is_empty() { return; }
        let Work { rr, v, s, f0, hist, .. } = w;
        self.rates_at(hist, rr, v, s, t);
        self.rhs(v, f0);
        // same time twice (an event): both kept, the value just before the jump and the one after;
        // state_at brackets strictly-before times with the first and later times with the second
        hist.t.push(t);
        hist.y.extend_from_slice(s);
        hist.f.extend_from_slice(f0);
    }

    /// Smallest positive delay time at the current state (step-size cap for delay models).
    fn delay_cap(&self, w: &mut Work, t: f64) -> f64 {
        if self.dly.is_empty() { return f64::INFINITY; }
        let (nx, np) = (self.nx, self.pnames.len());
        let Work { rr, s, hist, .. } = w;
        let mut rd = hist.rd.lock().unwrap();
        rd[nx..nx + np].copy_from_slice(&rr[nx..nx + np]);
        rd[..nx].copy_from_slice(s);
        rd[nx + np] = t;
        self.dtape.run_all(&mut rd);
        self.dly.iter().map(|&(_, tr, _, _)| rd[tr as usize]).filter(|&x| x > 0.0).fold(f64::INFINITY, f64::min)
    }

    /// Full Jacobian d(dx/dt)/dx (nx x nx, row-major) at (w.s, t).
    pub fn jacobian(&self, w: &mut Work, t: f64, out: &mut [f64]) {
        let Work { rf, v, ev, s, hist, .. } = w;
        self.full_at(hist, rf, v, ev, s, t);
        out.fill(0.0);
        for &(d, e, c) in &self.jfull { out[d as usize] += c * ev[e as usize]; }
    }

    /// d state/dt (concentrations, state_ids order) at (w.s, t) from the loaded parameters.
    pub fn state_rates(&self, w: &mut Work, t: f64, out: &mut [f64]) {
        let Work { rr, v, s, hist, .. } = w;
        self.rates_at(hist, rr, v, s, t);
        self.rhs(v, out);
    }

    #[inline]
    fn rhs(&self, v: &[f64], out: &mut [f64]) {
        out.fill(0.0);
        let (mut r, mut acc) = (u32::MAX, 0.0); // nz is row-major: accumulate each row in a register
        for &(i, j, c) in &self.nz {
            if i != r { if r != u32::MAX { out[r as usize] = acc; } r = i; acc = 0.0; }
            acc += c * v[j as usize];
        }
        if r != u32::MAX { out[r as usize] = acc; }
    }

    /// Reduced residual dS_indep/dt into f, returns (inf norm, term scale).
    #[inline]
    fn resid(&self, v: &[f64], f: &mut [f64]) -> (f64, f64) {
        f.fill(0.0);
        let mut sc = 0.0f64;
        for &(a, j, c) in &self.nzr {
            let t = c * v[j as usize];
            f[a as usize] += t;
            sc = sc.max(t.abs());
        }
        (inf_norm(f), sc)
    }

    #[inline]
    fn expand(&self, x: &[f64], tc: &[f64], s: &mut [f64]) {
        for (b, &i) in self.indep.iter().enumerate() { s[i] = x[b]; }
        for (d, &i) in self.dep.iter().enumerate() {
            s[i] = tc[d] + self.lc[d].iter().map(|&(b, l)| l * x[b as usize]).sum::<f64>();
        }
    }

    // ------------------------------------------------------------ steady state

    /// Newton on the reduced system from w.x (independent states) as `o` configures it.
    fn newton(&self, w: &mut Work, o: &NewtonOpts) -> bool {
        let ni = self.ni();
        let Work { rr, rf, v, ev, s, tc, x, xt, dx, f, jr, piv, hist, .. } = w;
        self.expand(x, tc, s);
        self.rates_at(hist, rr, v, s, 0.0);
        let (mut fn0, mut sc) = self.resid(v, f);
        let done = |fnv: f64, sc: f64| fnv <= 1e-14 * sc || (o.fnorm_tol > 0.0 && fnv <= o.fnorm_tol);
        for _ in 0..o.max_iter {
            if !fn0.is_finite() { return false; }
            if done(fn0, sc) { return true; }
            self.full_at(hist, rf, v, ev, s, 0.0);
            jr.fill(0.0);
            for &(d, e, c) in &self.jred { jr[d as usize] += c * ev[e as usize]; }
            for (d, &y) in dx.iter_mut().zip(f.iter()) { *d = -y; }
            if lu(jr, ni, piv, 1e-15) {
                lu_solve(jr, ni, piv, dx, 1);
            } else if o.rank_deficient {
                // NLEQ2's rank strategy: the minimum-norm least-squares step
                jr.fill(0.0);
                for &(d, e, c) in &self.jred { jr[d as usize] += c * ev[e as usize]; }
                let j = nalgebra::DMatrix::from_row_slice(ni, ni, jr);
                let eps = 1e-12 * inf_norm(jr).max(f64::MIN_POSITIVE);
                let Ok(st) = j.svd(true, true).solve(&nalgebra::DVector::from_column_slice(dx), eps) else { return false };
                dx.copy_from_slice(st.as_slice());
            } else {
                return false;
            }
            let mut lam = 1.0;
            loop {
                for b in 0..ni { xt[b] = x[b] + lam * dx[b]; }
                self.expand(xt, tc, s);
                self.rates_at(hist, rr, v, s, 0.0);
                let (fnt, sct) = self.resid(v, f);
                let decrease = fnt <= (1.0 - 1e-4 * lam) * fn0 || fnt <= 1e-14 * sct;
                if fnt.is_finite() && (decrease || !o.linesearch) {
                    x.copy_from_slice(xt);
                    let small = lam == 1.0 && inf_norm(dx) <= o.rtol * inf_norm(x) && fnt <= 1e-8 * sct;
                    let step = o.step_tol > 0.0 && lam * inf_norm(dx) <= o.step_tol;
                    fn0 = fnt; sc = sct;
                    if small || step || done(fnt, sct) { return true; }
                    break;
                }
                if !o.linesearch { return false; }
                lam *= 0.5;
                if lam < o.min_damping { return false; }
            }
        }
        false
    }

    fn valid(&self, w: &Work, allow_negative: bool) -> bool {
        w.s.iter().enumerate().all(|(i, &c)| c.is_finite() && (allow_negative || i >= self.ns || c >= -1e-10)) && w.v.iter().all(|r| r.is_finite())
    }

    /// Moiety totals of the dependent states from w.s (into w.tc).
    fn totals(&self, w: &mut Work) {
        for (d, &i) in self.dep.iter().enumerate() {
            w.tc[d] = w.s[i] - self.lc[d].iter().map(|&(b, l)| l * w.s[self.indep[b as usize]]).sum::<f64>();
        }
    }

    /// Newton from w.s (totals in w.tc). On success w.s, w.v, w.ev hold the solution; on failure
    /// w.s holds Newton's last iterate.
    fn newton_from_s(&self, w: &mut Work, o: &NewtonOpts) -> bool {
        for (b, &i) in self.indep.iter().enumerate() { w.x[b] = w.s[i]; }
        if self.newton(w, o) && self.valid(w, o.allow_negative) {
            let Work { rf, v, ev, s, hist, .. } = w;
            self.full_at(hist, rf, v, ev, s, 0.0);
            return self.valid(w, o.allow_negative);
        }
        self.expand(&w.x.clone(), &w.tc.clone(), &mut w.s); // ponytail: alloc only on the failure path
        false
    }

    /// Newton alone from w.s (params loaded), as `o` configures it.
    pub fn newton_solve(&self, w: &mut Work, o: &NewtonOpts) -> bool {
        self.totals(w);
        self.newton_from_s(w, o)
    }

    /// Steady state from w.s (initial state, params loaded): Newton, else integrate towards it over
    /// growing horizons. On success w.s, w.v, w.ev hold the solution.
    pub fn steady_state(&self, w: &mut Work) -> bool {
        self.totals(w);
        let o = NewtonOpts::default();
        let (mut t, mut h) = (0.0, 0.0);
        for horizon in [0.0, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6] {
            if horizon > 0.0 && !self.integrate(w, t, horizon, &mut h, 1e-6, 1e-12) { break; }
            t = horizon;
            if self.newton_from_s(w, &o) { return true; }
        }
        false
    }

    // ------------------------------------------------------------ Rodas4 integrator

    /// Integrate w.s from t0 to t1 (full system). `h` carries the step size between calls (0 = pick).
    pub fn integrate(&self, w: &mut Work, t0: f64, t1: f64, h: &mut f64, rtol: f64, atol: f64) -> bool {
        let mut t = t0;
        for _ in 0..MAX_STEPS {
            if t >= t1 { return true; }
            match self.step(w, t, t1, h, rtol, atol) { Some(tn) => t = tn, None => return false }
        }
        false
    }

    /// One accepted step from (t, w.s) towards t1 (never past it). Keeps the start state in w.y0 and
    /// the stages in w.k for dense output. Returns the new time.
    fn step(&self, w: &mut Work, t: f64, t1: f64, h: &mut f64, rtol: f64, atol: f64) -> Option<f64> {
        let nx = self.nx;
        let nj = self.jac.len();
        let Work { rr, rf, v, ev, s: y, yn, dy, fx, f0, k, jf, pivf, y0, lo, up, atol_scale: asc, stats: w_stats, hist, ctl, .. } = w;
        if !(*h > 0.0) { *h = if ctl.h0 > 0.0 { ctl.h0 } else { (1e-6 * (t1 - t).abs()).max(1e-12) }; }
        *h = ctl.lim(*h);
        self.full_at(hist, rf, v, ev, y, t);
        self.rhs(v, f0);
        if !self.dt.is_empty() {
            fx.fill(0.0);
            for &(i, j, c) in &self.nz {
                if let Some(p) = self.dt.iter().position(|&r| r == j) { fx[i as usize] += c * ev[nj + p]; }
            }
        }
        let mut rejected = false;
        for _ in 0..50 {
            let clipped = t + *h >= t1;
            let hh = if clipped { t1 - t } else { *h };
            if hh <= 1e-14 * t.abs().max(1.0) && !clipped { return None; }
            let fac = 1.0 / (hh * GAM);
            let fill = |jf: &mut [f64], all: bool| {
                if all { jf.fill(0.0) } else { for &q in &self.sp.pos { jf[q as usize] = 0.0; } }
                for &(d, e, c) in &self.jfull { jf[d as usize] -= c * ev[e as usize]; }
                for i in 0..nx { jf[i * nx + i] += fac; }
            };
            // sparse LU on the precomputed fill pattern; dense partial pivoting if a pivot is weak
            fill(jf, false);
            if self.sp.factor(jf, nx, lo, up) {
                for (i, p) in pivf.iter_mut().enumerate() { *p = i; }
            } else {
                fill(jf, true);
                if !lu(jf, nx, pivf, 0.0) { *h = hh * ctl.shrink(hh, 0.25)?; rejected = true; w_stats[2] += 1; continue; }
                lu_pattern(jf, nx, lo, up);
            }
            // stage 1
            for i in 0..nx { k[0][i] = f0[i] + hh * D[0] * fx[i]; }
            sparse_solve(pivf, lo, up, &mut k[0]);
            // stages 2..6
            let mut ok = true;
            for st in 1..6 {
                yn.copy_from_slice(y);
                for j in 0..st {
                    let a = A[st - 1][j];
                    for (o, &x) in yn.iter_mut().zip(&k[j]) { *o += a * x; }
                }
                self.rates_at(hist, rr, v, yn, t + CT[st] * hh);
                self.rhs(v, dy);
                let (head, tail) = k.split_at_mut(st);
                let ks = &mut tail[0];
                let hd = hh * D[st];
                for ((o, &a), &b) in ks.iter_mut().zip(dy.iter()).zip(fx.iter()) { *o = a + hd * b; }
                for j in 0..st {
                    let c = C[st - 1][j] / hh;
                    for (o, &x) in ks.iter_mut().zip(&head[j]) { *o += c * x; }
                }
                sparse_solve(pivf, lo, up, ks);
                if !ks.iter().all(|x| x.is_finite()) { ok = false; break; }
            }
            let mut err = 0.0;
            if ok {
                for i in 0..nx {
                    let sol = yn[i] + k[5][i]; // yn is the stage-6 point (embedded solution)
                    let sc = atol * asc[i] + rtol * y[i].abs().max(sol.abs());
                    err += (k[5][i] / sc).powi(2);
                    dy[i] = sol;
                }
                err = (err / nx.max(1) as f64).sqrt();
            }
            if ok && err <= 1.0 {
                y0.copy_from_slice(y);
                y.copy_from_slice(dy);
                w_stats[1] += 1;
                let mut hn = hh * (0.9 * err.powf(-0.25)).clamp(0.2, 6.0);
                if rejected { hn = hn.min(hh); }
                if !(clipped && hn > hh) { *h = ctl.lim(hn); }
                return Some(if clipped { t1 } else { t + hh });
            }
            *h = hh * ctl.shrink(hh, if ok { (0.9 * err.powf(-0.25)).clamp(0.2, 1.0) } else { 0.25 })?;
            rejected = true;
            w_stats[2] += 1;
        }
        None
    }

    // ------------------------------------------------------------ BDF integrator

    /// LU of fac*I - J (J given by its jacobian entries `jev`) into jf/pivf/lo/up.
    #[allow(clippy::too_many_arguments)]
    fn factor(&self, jev: &[f64], fac: f64, jf: &mut [f64], pivf: &mut [usize],
              lo: &mut Vec<(u32, u32, f64)>, up: &mut Vec<(u32, u32, f64)>) -> bool {
        let nx = self.nx;
        let fill = |jf: &mut [f64], all: bool| {
            if all { jf.fill(0.0) } else { for &q in &self.sp.pos { jf[q as usize] = 0.0; } }
            for &(d, e, c) in &self.jfull { jf[d as usize] -= c * jev[e as usize]; }
            for i in 0..nx { jf[i * nx + i] += fac; }
        };
        fill(jf, false);
        if self.sp.factor(jf, nx, lo, up) {
            for (i, p) in pivf.iter_mut().enumerate() { *p = i; }
            return true;
        }
        fill(jf, true);
        if !lu(jf, nx, pivf, 0.0) { return false; }
        lu_pattern(jf, nx, lo, up);
        true
    }

    /// One accepted BDF step from (t, w.s) towards t1 (never past it): variable-order (1-5) NDF in
    /// scipy's difference form (after Shampine & Reichelt's ode15s) with CVODE's Newton and step
    /// control: convergence-rate test at 0.1 of the error test (max 3 iterations), Jacobian and LU
    /// reused (refactor when h/alpha moved > 30% or after 20 steps, fresh Jacobian after 50 steps
    /// or a convergence failure), step size changed only for a gain >= 1.5 with CVODE's error bias.
    /// Returns the new time.
    fn bdf_step(&self, w: &mut Work, t: f64, t1: f64, rtol: f64, atol: f64) -> Option<f64> {
        let nx = self.nx;
        let Work { rr, rf, v, ev, s: y, jf, pivf, lo, up, tmp, bdf: b, atol_scale: asc, stats: w_stats, hist, ctl, .. } = w;
        let max_order = if ctl.max_order > 0 { ctl.max_order.min(BDF_MAX) } else { BDF_MAX };
        let Bdf { d, order, h, t: bt, neq, lu: has_lu, lufac, crate_, jev, sc, nlu, jage, defer } = b;
        let [yp, psi, dd, yn, f, dy, wt] = sc; // wt: inverse error weights 1 / (atol * scale + rtol |y|)
        tmp.resize(BDF_MAX + 3, 0.0);
        let newton_tol = 10.0 * f64::EPSILON / rtol;
        if *order == 0 { // (re)start: order 1, Hairer's initial step, fresh Jacobian
            self.rates_at(hist, rr, v, y, t);
            self.rhs(v, f);
            for i in 0..nx { wt[i] = 1.0 / (atol * asc[i] + rtol * y[i].abs()); }
            let (d0, d1) = (rms(nx, |i| y[i] * wt[i]), rms(nx, |i| f[i] * wt[i]));
            let h0 = if d0 < 1e-5 || d1 < 1e-5 { 1e-6 } else { 0.01 * d0 / d1 }.min(t1 - t);
            for i in 0..nx { yn[i] = y[i] + h0 * f[i]; }
            self.rates_at(hist, rr, v, yn, t + h0);
            self.rhs(v, dy);
            let d2 = rms(nx, |i| (dy[i] - f[i]) * wt[i]) / h0;
            let h1 = if d1 <= 1e-15 && d2 <= 1e-15 { (h0 * 1e-3).max(1e-6) } else { (0.01 / d1.max(d2)).sqrt() };
            *h = ctl.lim(if ctl.h0 > 0.0 { ctl.h0 } else { (100.0 * h0).min(h1) }).min(t1 - t);
            d.fill(0.0);
            d[..nx].copy_from_slice(y);
            for i in 0..nx { d[nx + i] = f[i] * *h; }
            *order = 1;
            *neq = 0;
            *has_lu = false;
            *crate_ = 1.0;
            *jage = 0;
            *defer = false;
            self.full_at(hist, rf, v, ev, y, t);
            jev.copy_from_slice(ev);
        }
        // step size change by `fac`: re-interpolate the differences; the last correction (used only
        // for the order k+1 error estimate) scales as h^(k+1)
        let rescale = |d: &mut [f64], tmp: &mut [f64], k: usize, fac: f64| {
            bdf_change(d, tmp, nx, k, fac);
            let s = fac.powi(k as i32 + 1);
            for x in &mut d[(k + 1) * nx..(k + 2) * nx] { *x *= s; }
        };
        let min_step = 10.0 * f64::EPSILON * t.abs().max(f64::MIN_POSITIVE);
        for ((w, &x), &a) in wt.iter_mut().zip(y.iter()).zip(asc.iter()) { *w = 1.0 / (atol * a + rtol * x.abs()); } // per step, as CVODE
        let (mut fresh_jac, mut nef) = (false, 0);
        let (ynew_t, err_norm) = loop {
            if *h < min_step { return None; }
            let mut tn = t + *h;
            if tn >= t1 {
                tn = t1;
                rescale(d, tmp, *order, (t1 - t) / *h);
                *h = t1 - t;
            }
            let k = *order;
            yp.copy_from_slice(&d[..nx]);
            psi.fill(0.0);
            for j in 1..=k {
                let g = BDF_GAMMA[j] / BDF_ALPHA[k];
                for ((y, p), &x) in yp.iter_mut().zip(psi.iter_mut()).zip(&d[j * nx..(j + 1) * nx]) { *y += x; *p += x * g; }
            }
            let c = *h / BDF_ALPHA[k];
            let converged = loop {
                if !*has_lu || (c * *lufac - 1.0).abs() > 0.3 || *nlu >= 20 {
                    if *jage >= 50 && !fresh_jac {
                        self.full_at(hist, rf, v, ev, yp, tn);
                        jev.copy_from_slice(ev);
                        *jage = 0;
                        fresh_jac = true;
                    }
                    *lufac = 1.0 / c;
                    *crate_ = 1.0;
                    *nlu = 0;
                    if !self.factor(jev, *lufac, jf, pivf, lo, up) { break false; }
                    *has_lu = true;
                }
                // Newton on (I - c J) dy = c f - psi - d with the LU of I - c' J: CVODE scales the
                // solution by 2 / (1 + c / c') and tests convergence with the rate estimate kept
                // across steps against 0.1 of the error test
                let gs = 2.0 / (1.0 + c * *lufac);
                yn.copy_from_slice(yp);
                dd.fill(0.0);
                let ntol = newton_tol.max(0.1 / BDF_ERR[k]);
                let (mut old, mut conv) = (0.0, false);
                for kk in 0..3 {
                    self.rates_at(hist, rr, v, yn, tn);
                    self.rhs(v, f);
                    if !f.iter().all(|x| x.is_finite()) { break; }
                    for (((o, &x), &p), &e) in dy.iter_mut().zip(f.iter()).zip(psi.iter()).zip(dd.iter()) { *o = (c * x - p - e) * *lufac; }
                    sparse_solve(pivf, lo, up, dy);
                    if gs != 1.0 { for x in dy.iter_mut() { *x *= gs; } }
                    let dn = rms(nx, |i| dy[i] * wt[i]);
                    if kk > 0 {
                        *crate_ = (0.3 * *crate_).max(dn / old);
                        if dn > 2.0 * old { break; }
                    }
                    for ((y, e), &x) in yn.iter_mut().zip(dd.iter_mut()).zip(dy.iter()) { *y += x; *e += x; }
                    if dn * crate_.min(1.0) <= ntol { conv = true; break; }
                    old = dn;
                }
                if conv || fresh_jac { break conv; }
                self.full_at(hist, rf, v, ev, yp, tn); // stale Jacobian: refresh once, then retry
                jev.copy_from_slice(ev);
                *jage = 0;
                *has_lu = false;
                fresh_jac = true;
            };
            if !converged {
                w_stats[2] += 1;
                let fac = ctl.shrink(*h, 0.25)?;
                *h *= fac;
                rescale(d, tmp, *order, fac);
                *neq = 0;
                continue;
            }
            let en = rms(nx, |i| BDF_ERR[k] * dd[i] * wt[i]);
            if en > 1.0 { // CVODE: h ratio from the error with bias 6, in [0.1, 0.2 after 2 failures]
                w_stats[2] += 1;
                nef += 1;
                *defer = true;
                let mut fac = if nef > 3 { 0.1 } else { (1.0 / ((6.0 * en).powf(1.0 / (k as f64 + 1.0)) + 1e-6)).max(0.1) };
                if nef >= 2 { fac = fac.min(0.2); }
                if nef > 3 && k > 1 { *order = k - 1; }
                let fac = ctl.shrink(*h, fac)?;
                *h *= fac;
                rescale(d, tmp, *order, fac);
                *neq = 0;
                continue;
            }
            break (tn, en);
        };
        let k = *order;
        w_stats[0] += 1;
        y.copy_from_slice(yn);
        *bt = ynew_t;
        *neq += 1;
        *nlu += 1;
        *jage += 1;
        let (a, b) = d.split_at_mut((k + 2) * nx);
        for ((p2, p1), &c) in b[..nx].iter_mut().zip(&mut a[(k + 1) * nx..]).zip(dd.iter()) { *p2 = c - *p1; *p1 = c; }
        for j in (0..=k).rev() {
            let (a, b) = d.split_at_mut((j + 1) * nx);
            for (x, &y) in a[j * nx..].iter_mut().zip(&b[..nx]) { *x += y; }
        }
        if std::mem::replace(defer, false) { return Some(ynew_t); } // no growth right after an error failure
        // CVODE's step/order choice: h ratios at orders k-1, k, k+1 with error biases 6, 6, 10
        let eta = |e: f64, q: usize, bias: f64| 1.0 / ((bias * e).powf(1.0 / q as f64) + 1e-6);
        let (mut fac, mut nk) = (eta(err_norm, k + 1, 6.0), k);
        if *neq >= k + 1 {
            let em = if k > 1 { let e = BDF_ERR[k - 1]; eta(rms(nx, |i| e * d[k * nx + i] * wt[i]), k, 6.0) } else { 0.0 };
            let ep = if k < max_order { let e = BDF_ERR[k + 1]; eta(rms(nx, |i| e * d[(k + 2) * nx + i] * wt[i]), k + 2, 10.0) } else { 0.0 };
            if em > fac && em >= ep { (fac, nk) = (em, k - 1) } else if ep > fac { (fac, nk) = (ep, k + 1) }
        }
        if fac < 1.5 { return Some(ynew_t); } // small gain: keep h and order (and the LU)
        fac = ctl.lim(*h * fac.min(10.0)) / *h;
        if nk != k { *order = nk; *neq = 0; }
        *h *= fac;
        rescale(d, tmp, *order, fac);
        Some(ynew_t)
    }

    /// BDF interpolant of the last step at time te into w.yn.
    fn bdf_interp(&self, w: &mut Work, te: f64) {
        let nx = self.nx;
        let Work { yn, bdf: b, .. } = w;
        yn.copy_from_slice(&b.d[..nx]);
        let mut p = 1.0;
        for j in 0..b.order {
            p *= (te - (b.t - b.h * j as f64)) / (b.h * (j + 1) as f64);
            for i in 0..nx { yn[i] += b.d[(j + 1) * nx + i] * p; }
        }
    }

    // ------------------------------------------------------------ events

    /// Any trigger value at (s, t) differs from w.prev.
    fn changed(&self, hs: &Hist, ra: &mut [f64], prev: &[bool], prevc: &[bool], s: &[f64], t: f64) -> bool {
        self.aux_at(hs, ra, s, t);
        self.events.iter().zip(prev).any(|(e, &p)| (ra[e.trig as usize] != 0.0) != p)
            || self.comps.iter().zip(prevc).any(|(&r, &p)| (ra[r as usize] != 0.0) != p)
    }

    /// Earliest trigger change inside the last step (ta, tb]: bisection on the active solver's dense output
    /// to ~1e-12 relative in time. Returns the first time at which a trigger has changed.
    fn locate(&self, w: &mut Work, ta: f64, tb: f64, bdf: bool) -> f64 {
        if !bdf { self.dense_coef(w); }
        let tol = 1e-12 * tb.abs().max(tb - ta) / (tb - ta);
        let (mut lo, mut hi) = (0.0, 1.0);
        while hi - lo > tol {
            let th = 0.5 * (lo + hi);
            if bdf { self.bdf_interp(w, ta + th * (tb - ta)) } else { self.interp(w, th) }
            let Work { yn, ra, prev, prevc, hist, .. } = w;
            if self.changed(hist, ra, prev, prevc, yn, ta + th * (tb - ta)) { hi = th } else { lo = th }
        }
        if hi == 1.0 { tb } else { ta + hi * (tb - ta) }
    }

    /// Dense-output coefficients of the last step.
    fn dense_coef(&self, w: &mut Work) {
        let Work { k, d1, d2, .. } = w;
        for i in 0..self.nx {
            d1[i] = (0..5).map(|j| H[0][j] * k[j][i]).sum();
            d2[i] = (0..5).map(|j| H[1][j] * k[j][i]).sum();
        }
    }

    /// Last step's dense output at fraction th into w.yn (needs dense_coef).
    #[inline]
    fn interp(&self, w: &mut Work, th: f64) {
        let Work { s: y1, y0, d1, d2, yn, .. } = w;
        for i in 0..self.nx { yn[i] = (1.0 - th) * y0[i] + th * (y1[i] + (1.0 - th) * (d1[i] + th * d2[i])); }
    }

    /// Event state: (trigger value per event, pending (time, event, values fixed at trigger time)).
    pub fn event_state(&self, w: &Work) -> (Vec<bool>, Vec<(f64, usize, Vec<f64>)>) { (w.prev.clone(), w.pend.clone()) }

    pub fn set_event_state(&self, w: &mut Work, prev: Vec<bool>, pend: Vec<(f64, usize, Vec<f64>)>) -> Result<(), String> {
        let ne = self.events.len();
        if prev.len() != ne || pend.iter().any(|p| p.1 >= ne || !(p.2.is_empty() || p.2.len() == self.events[p.1].asg.len())) {
            return Err("event state does not match this model".into());
        }
        w.prev = prev;
        w.pend = pend;
        Ok(())
    }

    /// Reset event bookkeeping for a run starting at t0 (trigger state before t0 = initialValue).
    pub fn events_reset(&self, w: &mut Work) {
        w.hist.clear(); // a fresh run: delay() sees the initial state before it
        w.pend.clear();
        for (p, e) in w.prev.iter_mut().zip(&self.events) { *p = e.init; }
    }

    /// Process trigger transitions and execute due events at time t, one at a time by priority
    /// (ties: document order), re-evaluating triggers after each (cascades, persistence).
    /// Returns None if a cascade does not terminate, else whether any event executed.
    fn events_at(&self, w: &mut Work, t: f64) -> Option<bool> {
        let nx = self.nx;
        let Work { rr, rf, ra, s, prev, prevc, pend, tmp, hist, rng, .. } = w;
        let mut fired = false;
        for _ in 0..10_000 {
            self.aux_at(hist, ra, s, t);
            for (p, &r) in prevc.iter_mut().zip(&self.comps) { *p = ra[r as usize] != 0.0; }
            for (e, ev) in self.events.iter().enumerate() {
                let cur = ra[ev.trig as usize] != 0.0;
                if cur && !prev[e] {
                    let d = ev.delay.map_or(0.0, |r| ra[r as usize]);
                    let vals = if ev.uvftt { ev.asg.iter().map(|&(_, r)| ra[r as usize]).collect() } else { vec![] };
                    pend.push((t + d, e, vals));
                } else if !cur && prev[e] && !ev.persistent {
                    pend.retain(|p| p.1 != e);
                }
                prev[e] = cur;
            }
            let pr = |e: usize| self.events[e].prio.map_or(0.0, |r| ra[r as usize]);
            let due: Vec<usize> = (0..pend.len()).filter(|&i| pend[i].0 <= t).collect();
            if due.is_empty() { return Some(fired) }
            let top = due.iter().map(|&i| pr(pend[i].1)).fold(f64::NEG_INFINITY, f64::max);
            let tied: Vec<usize> = due.into_iter().filter(|&i| pr(pend[i].1) == top || top.is_nan()).collect();
            // equal explicit priorities: uniformly random (SBML L3); otherwise document order
            let i = if tied.len() > 1 && tied.iter().all(|&i| self.events[pend[i].1].prio.is_some()) {
                tied[((rng.unif() * tied.len() as f64) as usize).min(tied.len() - 1)]
            } else {
                *tied.iter().min_by(|&&a, &&b| pend[a].1.cmp(&pend[b].1).then(a.cmp(&b))).unwrap()
            };
            let (_, e, vals) = pend.remove(i);
            let ev = &self.events[e];
            tmp.clear();
            if ev.uvftt { tmp.extend(vals) } else { tmp.extend(ev.asg.iter().map(|&(_, r)| ra[r as usize])) }
            let vold: Vec<f64> = ev.rescale.iter().map(|&(_, r)| ra[r as usize]).collect(); // ponytail: allocates per event
            let mut par = false;
            for ((&(tg, _), &x), dv) in ev.asg.iter().zip(tmp.iter()).zip(&ev.div) {
                let x = dv.map_or(x, |r| x / ra[r as usize]);
                let tg = tg as usize;
                if tg < nx { s[tg] = x } else { rr[tg] = x; rf[tg] = x; ra[tg] = x; par = true }
            }
            if par { self.rates.prep(rr); self.full.prep(rf); self.aux.prep(ra); }
            if !ev.rescale.is_empty() {
                self.aux_at(hist, ra, s, t);
                for (&(i, r), &v0) in ev.rescale.iter().zip(&vold) { s[i as usize] *= v0 / ra[r as usize]; }
            }
            fired = true;
        }
        None
    }

    /// dy/dt at (t, y) into out (full state, same units as `simulate`), via the rates-only tape.
    pub fn deriv(&self, w: &mut Work, t: f64, y: &[f64], out: &mut [f64]) {
        let Work { rr, v, hist, .. } = w;
        self.rates_at(hist, rr, v, y, t);
        self.rhs(v, out);
    }

    /// Any delay() in the model.
    pub fn has_delays(&self) -> bool { !self.dly.is_empty() }

    /// Does any event trigger or tracked comparison differ from its last accepted value?
    pub(crate) fn trigger_changed(&self, w: &mut Work, y: &[f64], t: f64) -> bool {
        self.changed(&w.hist, &mut w.ra, &w.prev, &w.prevc, y, t)
    }

    /// Limit a step to the next scheduled event execution.
    pub(crate) fn next_event_time(&self, w: &Work, end: f64) -> f64 {
        w.pend.iter().fold(end, |t, p| t.min(p.0))
    }

    pub(crate) fn event_roots(&self, w: &mut Work, t: f64) -> Vec<f64> {
        self.aux_at(&w.hist, &mut w.ra, &w.s, t);
        self.events.iter().map(|e| if w.ra[e.trig as usize] != 0.0 { 1.0 } else { -1.0 }).collect()
    }

    /// Is any event trigger true at (w.s, t)?
    pub fn any_trigger(&self, w: &mut Work, t: f64) -> bool {
        let Work { ra, s, hist, .. } = w;
        self.aux_at(hist, ra, s, t);
        self.events.iter().any(|e| ra[e.trig as usize] != 0.0)
    }

    /// Process trigger transitions and execute due events at t (as `simulate` does between steps).
    /// None if a cascade does not terminate, else whether any event executed.
    pub fn apply_events(&self, w: &mut Work, t: f64) -> Option<bool> { self.events_at(w, t) }

    /// Write the requested output registers at (w.s, t) into row.
    pub fn record(&self, w: &mut Work, t: f64, outs: &[u32], row: &mut [f64]) {
        self.record_y(&w.hist, &mut w.ra, &w.s, t, outs, row);
    }

    fn record_y(&self, hs: &Hist, ra: &mut [f64], y: &[f64], t: f64, outs: &[u32], row: &mut [f64]) {
        let (nx, tr) = (self.nx, self.nx + self.pnames.len()); // time register; params sit in between
        if self.dly.is_empty() && outs.iter().all(|&r| (r as usize) <= tr) { // states, parameters, time: no tape run
            for (o, &r) in row.iter_mut().zip(outs) {
                let r = r as usize;
                *o = if r < nx { y[r] } else if r == tr { t } else { ra[r] };
            }
        } else {
            self.aux_at(hs, ra, y, t);
            for (o, &r) in row.iter_mut().zip(outs) { *o = ra[r as usize]; }
        }
    }

    // ------------------------------------------------------------ MCA

    /// Unscaled + scaled control coefficients at the steady state held in w (after steady_state).
    /// Outputs: e (nr x ns), uc (ns x nr), uf (nr x nr), cc, fc, el. Returns false if singular.
    pub fn mca(&self, w: &mut Work, e: &mut [f64], uc: &mut [f64], uf: &mut [f64],
               cc: &mut [f64], fc: &mut [f64], el: &mut [f64]) -> bool {
        let (ns, nr) = (self.ns, self.nr);
        let (indep, dep, lc, nzr, jred) = match &self.mca_c {
            Some(c) => (&c.indep, &c.dep, &c.lc, &c.nzr, &c.jred),
            None => (&self.indep, &self.dep, &self.lc, &self.nzr, &self.jred),
        };
        let ni = indep.len();
        let Work { v, ev, s, jr, piv, xm, .. } = w;
        let (jr, xm) = (&mut jr[..ni * ni], &mut xm[..ni * nr]);
        jr.fill(0.0);
        for &(d, e_, c) in jred { jr[d as usize] += c * ev[e_ as usize]; }
        if !lu(jr, ni, piv, 1e-13) { return false; }
        xm.fill(0.0);
        for &(a, j, c) in nzr { xm[a as usize * nr + j as usize] += c; }
        lu_solve(jr, ni, piv, xm, nr);
        for (b, &i) in indep.iter().enumerate() {
            for j in 0..nr { uc[i * nr + j] = -xm[b * nr + j]; }
        }
        for (d, &i) in dep.iter().enumerate() {
            let row = &mut uc[i * nr..(i + 1) * nr];
            row.fill(0.0);
            for &(b, l) in &lc[d] {
                for (o, &x) in row.iter_mut().zip(&xm[b as usize * nr..(b as usize + 1) * nr]) { *o -= l * x; }
            }
        }
        uf.fill(0.0);
        for j in 0..nr { uf[j * nr + j] = 1.0; }
        e.fill(0.0);
        for (q, &(j, kk)) in self.jac.iter().enumerate() {
            let (j, kk) = (j as usize, kk as usize);
            if j >= nr || kk >= ns { continue; }
            e[j * ns + kk] = ev[q];
            for c in 0..nr { uf[j * nr + c] += ev[q] * uc[kk * nr + c]; }
        }
        for i in 0..ns { for j in 0..nr { cc[i * nr + j] = uc[i * nr + j] * v[j] / s[i]; } }
        for i in 0..nr { for j in 0..nr { fc[i * nr + j] = uf[i * nr + j] * v[j] / v[i]; } }
        for i in 0..nr { for j in 0..ns { el[i * ns + j] = e[i * ns + j] * s[j] / v[i]; } }
        true
    }

    /// d v_reactions / d state (nr x nx).
    pub fn elasticity_into(&self, w: &Work, e: &mut [f64]) {
        e.fill(0.0);
        for (q, &(j, kk)) in self.jac.iter().enumerate() {
            if (j as usize) < self.nr { e[j as usize * self.nx + kk as usize] = w.ev[q]; }
        }
    }

    // ------------------------------------------------------------ batches

    /// Run `f(work, sample index)` over n samples in parallel with per-thread workspaces.
    pub fn par<F: Fn(&mut Work, usize) + Sync + Send>(&self, n: usize, f: F) {
        (0..n).into_par_iter().with_min_len(4).for_each_init(|| self.work(), |w, i| f(w, i));
    }

    /// Time course over `times` from the loaded state, writing outputs `outs` (registers) per time
    /// into out (len(times) x len(outs)). BDF below BDF_BELOW_RTOL, Rodas4 otherwise (w.stats counts
    /// steps). Adaptive steps are clipped only at the end time and at event times; output points in
    /// between come from the solver's dense output. Events fire with
    /// roadrunner semantics; output points at an event time report the post-event state.
    /// On success w.s / params hold the final state.
    pub fn simulate(&self, w: &mut Work, times: &[f64], rtol: f64, atol: f64, outs: &[u32], out: &mut [f64]) -> bool {
        self.events_reset(w);
        self.simulate_cont(w, times, rtol, atol, outs, out)
    }

    /// `simulate` keeping the event state (trigger values, pending delayed events) in `w`, as a
    /// roadrunner simulate call that continues from the previous one.
    pub fn simulate_cont(&self, w: &mut Work, times: &[f64], rtol: f64, atol: f64, outs: &[u32], out: &mut [f64]) -> bool {
        self.simulate_until(w, times, rtol, atol, outs, out, None)
    }

    /// One accepted native integration step, stopping early at an event.
    pub(crate) fn integrate_one(&self, w: &mut Work, t: f64, end: f64, rtol: f64, h: &mut f64) -> Option<f64> {
        let mut reached = t;
        self.simulate_until(w, &[t, end], rtol, 1.0, &[], &mut [], Some((&mut reached, h))).then_some(reached)
    }

    #[allow(clippy::too_many_arguments)]
    fn simulate_until(&self, w: &mut Work, times: &[f64], rtol: f64, atol: f64, outs: &[u32], out: &mut [f64],
                      step_out: Option<(&mut f64, &mut f64)>) -> bool {
        crate::rand::seed(w.rng.next_u64()); // distrib draws
        w.rows.0.clear();
        w.rows.1.clear();
        if self.alg.contains(&true) { return self.simulate_dae(w, times, rtol, atol, outs, out, step_out.map(|(t, _)| t)); }
        let (no, nt) = (outs.len(), times.len());
        let ev = self.has_events();
        let tend = times[nt - 1];
        let mut t = times[0];
        // Prefer the one-step method in nonstiff mode. The stiff default uses BDF at tight
        // tolerances, or whenever a lower BDF order was explicitly requested.
        let mut bdf = !w.ctl.nonstiff && (rtol < BDF_BELOW_RTOL || (1..BDF_MAX).contains(&w.ctl.max_order));
        let dl = !self.dly.is_empty(); // delay(): one-step Rodas4, steps no longer than the shortest delay
        if dl { bdf = false; }
        let mut h = step_out.as_ref().map_or(0.0, |(_, h)| **h);
        if h == 0.0 { w.bdf.order = 0; }
        if ev {
            match self.events_at(w, t) {
                None => return false,
                Some(true) => w.bdf.order = 0,
                Some(false) => {}
            }
        }
        if dl && w.hist.t.last().is_none_or(|&x| x < t) { self.hist_push(w, t); }
        let (mut q, mut n) = (0, 0); // n: steps since the last output
        if !self.push_row(w, t, outs) { return true; }
        for _ in 0..MAX_STEPS * 10 {
            while q < nt && times[q] <= t { self.record(w, t, outs, &mut out[q * no..(q + 1) * no]); q += 1; n = 0; }
            if q == nt { return true; }
            n += 1;
            if n > w.ctl.step_limit() { return false; }
            let mut t1 = self.next_event_time(w, tend);
            // Clip at the delay boundary even when the solver's requested minimum is larger.
            if dl {
                let cap = self.delay_cap(w, t);
                if h > cap { h = cap; }
                if w.ctl.hmin > cap || (h == 0.0 && w.ctl.h0 > cap) { t1 = t1.min(t + cap); }
            }
            let mut step = if bdf { self.bdf_step(w, t, t1, rtol, atol) } else { self.step(w, t, t1, &mut h, rtol, atol) };
            if step.is_none() && bdf { // BDF gave up (e.g. loose atol on a near-singular system): Rodas4 from here
                bdf = false;
                w.stats[3] += 1;
                step = self.step(w, t, t1, &mut h, rtol, atol);
            }
            let Some(tn) = step else { return false };
            let mut te = tn;
            if ev {
                let Work { ra, s, prev, prevc, hist, .. } = w;
                if self.changed(hist, ra, prev, prevc, s, tn) { te = self.locate(w, t, tn, bdf); }
            }
            if q < nt && times[q] < te && !w.ctl.multi { // interior output points from the dense output
                if !bdf { self.dense_coef(w); }
                while q < nt && times[q] < te {
                    if bdf { self.bdf_interp(w, times[q]) } else { self.interp(w, (times[q] - t) / (tn - t)) }
                    let Work { ra, yn, hist, .. } = w;
                    self.record_y(hist, ra, yn, times[q], outs, &mut out[q * no..(q + 1) * no]);
                    q += 1;
                    n = 0;
                }
            }
            if te < tn && bdf { // state at the event from the interpolant; restart after it
                self.bdf_interp(w, te);
                w.s.copy_from_slice(&w.yn);
            } else if te < tn { // redo the step exactly up to the event
                w.s.copy_from_slice(&w.y0);
                let mut h2 = te - t;
                if !self.integrate(w, t, te, &mut h2, rtol, atol) { return false; }
            }
            t = te;
            self.hist_push(w, t);
            if !self.push_row(w, t, outs) { return true; } // variable steps: the row before the events
            if ev {
                let fired = match self.events_at(w, t) { None => return false, Some(f) => f };
                if fired || te < tn { w.bdf.order = 0; n = 0; }
                if dl { self.hist_push(w, t); } // the state after the events, at the same time
                if fired && !self.push_row(w, t, outs) { return true; } // and after them
            }
            while w.ctl.multi && q < nt && times[q] < t { // CV_ONE_STEP: state from the step past the output
                self.record(w, times[q], outs, &mut out[q * no..(q + 1) * no]);
                q += 1;
            }
            if let Some((time, next_h)) = step_out {
                *time = t;
                // A positive value also tells the next call to keep BDF history.
                *next_h = if bdf { tn - times[0] } else { h };
                return true;
            }
        }
        false
    }

    /// With ctl.var, a row of `outs` at t into w.rows (as libroadrunner's variable-step simulate);
    /// false when it has ctl.max_rows rows (0: no limit).
    fn push_row(&self, w: &mut Work, t: f64, outs: &[u32]) -> bool {
        let Work { hist, ra, s, rows: (ts, ys), ctl, .. } = w;
        if !ctl.var { return true; }
        ts.push(t);
        let n = ys.len();
        ys.resize(n + outs.len(), 0.0);
        self.record_y(hist, ra, s, t, outs, &mut ys[n..]);
        ctl.max_rows == 0 || ts.len() < ctl.max_rows
    }

    /// Time course of a DAE model (algebraic states: M = diag(1 | 0)) with diffsol's BDF. Each
    /// segment starts from a consistent state (diffsol solves the algebraic part first) and is
    /// stepped until a trigger (or a comparison inside one) changes, located by bisection on
    /// diffsol's interpolant, a pending delayed event comes due, or the end time. Events then fire
    /// as in `simulate_cont` and diffsol restarts from the post-event state. delay() is not
    /// supported here (returns false).
    #[allow(clippy::too_many_arguments)]
    fn simulate_dae(&self, w: &mut Work, times: &[f64], rtol: f64, atol: f64, outs: &[u32], out: &mut [f64],
                    mut reached: Option<&mut f64>) -> bool {
        use diffsol::{NalgebraLU, NalgebraMat, NalgebraVec, OdeBuilder, OdeSolverMethod, OdeSolverStopReason, Vector};
        if !self.dly.is_empty() { return false; }
        let (nx, no, nt) = (self.nx, outs.len(), times.len());
        let has_ev = self.has_events();
        let tend = times[nt - 1];
        let (v, ev) = (RefCell::new(vec![0.0; self.nv]), RefCell::new(vec![0.0; self.jac.len() + self.dt.len()]));
        let hs = Hist::default(); // no delay() on this path
        let atols: Vec<f64> = (0..nx).map(|i| atol * w.atol_scale[i]).collect();
        let get = |y: &NalgebraVec<f64>, o: &mut [f64]| for (i, x) in o.iter_mut().enumerate() { *x = y.get_index(i) };
        let ctl = w.ctl;
        let (mut t, mut q, mut steps, mut n, mut row) = (times[0], 0, 0, 0, true); // row: a variable-step row due after the events
        'seg: loop {
            // a fresh problem per segment: events may have changed parameters and states
            let (rr, rf) = (RefCell::new(w.rr.clone()), RefCell::new(w.rf.clone()));
            let y0 = w.s.clone();
            let rhs = |x: &[f64], _p: &[f64], t: f64, y: &mut [f64]| {
                let mut v = v.borrow_mut();
                self.rates_at(&hs, &mut rr.borrow_mut(), &mut v, x, t);
                self.rhs(&v, y);
            };
            let jac = |x: &[f64], _p: &[f64], t: f64, u: &[f64], y: &mut [f64]| {
                let (mut v, mut ev) = (v.borrow_mut(), ev.borrow_mut());
                self.full_at(&hs, &mut rf.borrow_mut(), &mut v, &mut ev, x, t);
                y.fill(0.0);
                for &(d, e, c) in &self.jfull {
                    let (i, k) = (d as usize / nx, d as usize % nx);
                    y[i] += c * ev[e as usize] * u[k];
                }
            };
            let mass = |u: &[f64], _p: &[f64], _t: f64, beta: f64, y: &mut [f64]| {
                for i in 0..nx { y[i] = if self.alg[i] { 0.0 } else { u[i] } + beta * y[i]; }
            };
            let init = |_p: &[f64], _t: f64, y: &mut [f64]| y.copy_from_slice(&y0);
            let Ok(mut problem) = OdeBuilder::<NalgebraMat<f64>>::new().t0(t).rtol(rtol).atol(atols.clone())
                .rhs_implicit(rhs, jac).mass(mass).init(init, nx).build() else { return false };
            // Algebraic equations always use BDF, even with stiff=false.
            problem.ode_options.min_timestep = ctl.hmin.max(f64::EPSILON);
            if ctl.max_steps > 0 { problem.ode_options.max_number_timesteps = MAX_STEPS * 10; }
            let Ok(mut solver) = problem.bdf::<NalgebraLU<f64>>() else { return false };
            let h0 = ctl.lim(if ctl.h0 > 0.0 { ctl.h0 } else { solver.state().h });
            if h0 != solver.state().h { *solver.state_mut().h = h0; }
            get(solver.state().y, &mut w.s); // consistent algebraic states
            steps += 1;
            if steps > MAX_STEPS * 10 { return false; }
            if has_ev { // events see the consistent state; a firing restarts (algebraic states follow)
                match self.events_at(w, t) { None => return false, Some(true) => { row = true; continue 'seg } Some(false) => {} }
            }
            if std::mem::replace(&mut row, false) && !self.push_row(w, t, outs) { return true; }
            while q < nt && times[q] <= t {
                self.record(w, times[q], outs, &mut out[q * no..(q + 1) * no]);
                q += 1;
                n = 0;
            }
            if q == nt { return true; }
            let segment_end = w.pend.iter().fold(tend, |m, p| m.min(p.0));
            let mut ta = t; // start of the current step
            loop {
                steps += 1;
                n += 1;
                if steps > MAX_STEPS * 10 || n > ctl.step_limit() { return false; }
                // diffsol exposes its proposed order, but has no adaptive order cap.
                // Marking the state modified makes its next step restart at order one.
                if ctl.max_order > 0 && solver.order() > ctl.max_order { let _ = solver.state_mut(); }
                let cap = if ctl.hmax > 0.0 { segment_end.min(ta + ctl.hmax) } else { segment_end };
                if solver.set_stop_time(cap).is_err() { return false; }
                let Ok(reason) = solver.step() else { return false };
                let tn = solver.state().t;
                get(solver.state().y, &mut w.yn);
                let mut te = tn;
                if has_ev {
                    let Work { ra, yn, prev, prevc, hist, .. } = w;
                    if self.changed(hist, ra, prev, prevc, yn, tn) {
                        let at = |x: f64, yn: &mut [f64], ra: &mut [f64]| match solver.interpolate(x) {
                            Ok(y) => { get(&y, yn); Some(self.changed(hist, ra, prev, prevc, yn, x)) }
                            Err(_) => None,
                        };
                        let (mut lo, mut hi) = (ta, tn);
                        while hi - lo > 1e-12 * tn.abs().max(tn - ta) {
                            let mid = 0.5 * (lo + hi);
                            let Some(c) = at(mid, yn, ra) else { return false };
                            if c { hi = mid } else { lo = mid }
                        }
                        // a trigger like time >= 1 holds exactly at 1: prefer the round time if it already flips
                        let r = round12(hi);
                        if r > lo && r < hi && at(r, yn, ra) == Some(true) { hi = r; }
                        te = hi;
                    }
                }
                while q < nt && times[q] < te { // output points inside the step from the interpolant
                    if !ctl.multi || te < tn {
                        let Ok(y) = solver.interpolate(if ctl.multi { te } else { times[q] }) else { return false };
                        get(&y, &mut w.yn);
                    }
                    let Work { ra, yn, hist, .. } = w;
                    self.record_y(hist, ra, yn, times[q], outs, &mut out[q * no..(q + 1) * no]);
                    q += 1;
                    n = 0;
                }
                if let Some(time) = reached.take() {
                    if te < tn {
                        let Ok(y) = solver.interpolate(te) else { return false };
                        get(&y, &mut w.s);
                    } else {
                        get(solver.state().y, &mut w.s);
                    }
                    if has_ev && self.events_at(w, te).is_none() { return false; }
                    *time = te;
                    return true;
                }
                if te < tn || (cap == segment_end && matches!(reason, OdeSolverStopReason::TstopReached)) {
                    if te < tn {
                        let Ok(y) = solver.interpolate(te) else { return false };
                        get(&y, &mut w.s);
                    } else {
                        get(solver.state().y, &mut w.s);
                    }
                    t = te;
                    n = 0;
                    if !self.push_row(w, t, outs) { return true; }
                    continue 'seg;
                }
                if ctl.var || (q < nt && times[q] <= tn) {
                    get(solver.state().y, &mut w.s);
                    if !self.push_row(w, tn, outs) { return true; }
                    while q < nt && times[q] <= tn {
                        self.record(w, times[q], outs, &mut out[q * no..(q + 1) * no]);
                        q += 1;
                        n = 0;
                    }
                }
                ta = tn;
            }
        }
    }

    /// Gillespie direct method from the loaded state (concentrations; amounts = conc * volume,
    /// propensity = rate law at concentrations). Events are checked after every reaction and at
    /// their delayed execution times. Outputs as in `simulate`.
    /// The firing limit applies to each reporting interval (or maximum-time-step subdivision).
    pub fn gillespie(&self, w: &mut Work, times: &[f64], seed: u64, outs: &[u32], out: &mut [f64]) -> bool {
        let ctl = w.ctl;
        let hmax = if ctl.hmax > ctl.hmin { ctl.hmax } else { 0.0 };
        let (no, mut q, mut n, mut over) = (outs.len(), 0, 0, false);
        let mut window = times[0];
        self.ssa(w, times[0], seed, times, |w, _, stop, firing| {
            while q < times.len() && times[q] < stop {
                self.record(w, times[q], outs, &mut out[q * no..(q + 1) * no]);
                window = times[q];
                q += 1;
                n = 0;
            }
            if q == times.len() { return true; }
            if hmax > 0.0 && stop > window + hmax {
                // Constant-rate waiting times are memoryless. Keep the sampled firing, but
                // reset the budget at any reporting subdivisions that it passes.
                window += (((stop - window) / hmax).ceil() - 1.0) * hmax;
                n = 0;
            }
            n += usize::from(firing);
            let target = if hmax > 0.0 { times[q].min(window + hmax) } else { times[q] };
            over = firing && ctl.max_steps > 0 && n >= ctl.max_steps && stop < target;
            over
        }) && !over
    }

    /// Variable-step Gillespie from t0 to t1: an output row (time, then `outs`) at t0, after every
    /// reaction firing or event, and at t1; at most max_rows rows. None on failure.
    pub fn gillespie_steps(&self, w: &mut Work, t0: f64, t1: f64, seed: u64, outs: &[u32], max_rows: usize) -> Option<Vec<f64>> {
        let no = outs.len() + 1;
        let max_rows = if max_rows == 0 { usize::MAX } else { max_rows };
        let mut rows: Vec<f64> = vec![];
        let mut row = vec![0.0; outs.len()];
        let mut push = |w: &mut Work, t: f64, rows: &mut Vec<f64>| {
            self.record(w, t, outs, &mut row);
            rows.push(t);
            rows.extend_from_slice(&row);
        };
        let hmax = w.ctl.hmax;
        let ok = self.ssa(w, t0, seed, &[t1], |w, t, stop, _| {
            push(w, t.min(t1), &mut rows);
            if rows.len() / no == max_rows || t >= t1 { return true; }
            if stop > t1 {
                // A finite waiting time may exceed the cap in single-step mode. With
                // no more reactions, libroadrunner advances by the capped request.
                let mut k = 1.0;
                while stop.is_infinite() && hmax > 0.0 && t + k * hmax < t1 && rows.len() / no < max_rows {
                    push(w, t + k * hmax, &mut rows);
                    k += 1.0;
                }
                if rows.len() / no < max_rows { push(w, t1, &mut rows); }
                return true;
            }
            false
        });
        ok.then_some(rows)
    }

    pub(crate) fn gillespie_one(&self, w: &mut Work, t0: f64, end: f64,
                                uniform: impl FnMut() -> f64) -> Option<(f64, f64)> {
        let mut fired = f64::NAN;
        let mut last = t0;
        let ok = self.ssa_run(w, t0, &[end], uniform, Some(&mut fired), |_, t, stop, _| {
            last = t;
            stop.is_infinite() || (t >= end && (self.nv != self.nr || !self.dt.is_empty()))
        });
        ok.then(|| if fired.is_nan() { (last, end) } else { (fired, fired) })
    }

    /// SSA core: the callback sees the state valid on [t, stop) and whether the segment ends
    /// in a reaction attempt. It returns true to finish. With rate rules (`ssa_hybrid`) the state at t
    /// holds at t only, and stop never passes the next `grid` time.
    fn ssa(&self, w: &mut Work, t0: f64, seed: u64, grid: &[f64], at: impl FnMut(&mut Work, f64, f64, bool) -> bool) -> bool {
        let mut rng = Rng::new(seed);
        crate::rand::seed(seed ^ 0xD157_D157_D157_D157);
        self.events_reset(w);
        self.ssa_run(w, t0, grid, || rng.unif(), None, at)
    }

    fn ssa_run(&self, w: &mut Work, t0: f64, grid: &[f64], mut uniform: impl FnMut() -> f64,
               mut one: Option<&mut f64>, mut at: impl FnMut(&mut Work, f64, f64, bool) -> bool) -> bool {
        let (ns, nr) = (self.ns, self.nr);
        let mut t = t0;
        let resync = |w: &mut Work| for i in 0..ns { w.amt[i] = w.s[i] * self.vol[i]; };
        resync(w);
        if self.has_events() && self.events_at(w, t).is_none() { return false; }
        resync(w);
        if self.nv != nr || !self.dt.is_empty() { return self.ssa_hybrid(w, t, &mut uniform, grid, one, at); }
        loop {
            let Work { rr, v, s, hist, .. } = w;
            self.rates_at(hist, rr, v, s, t);
            let a0: f64 = v[..nr].iter().map(|&a| a.max(0.0)).sum();
            if !a0.is_finite() { return false; }
            let tn = if a0 > 0.0 { t - uniform().ln() / a0 } else { f64::INFINITY };
            let mut te = w.pend.iter().fold(f64::INFINITY, |m, p| m.min(p.0));
            if self.has_events() { // a trigger that turns with time alone while the state waits (t >= 25)
                let tb = tn.min(te).min(t + 1e6 * (1.0 + t.abs()));
                let Work { ra, s, prev, prevc, hist, .. } = w;
                if tb > t && self.changed(hist, ra, prev, prevc, s, tb) {
                    let (mut lo, mut hi) = (t, tb);
                    while hi - lo > 1e-12 * hi.abs().max(1.0) {
                        let mid = 0.5 * (lo + hi);
                        if self.changed(hist, ra, prev, prevc, s, mid) { hi = mid } else { lo = mid }
                    }
                    // a trigger like t >= 25 holds exactly at 25: prefer the round time if it already flips
                    let r = { let sc = 10f64.powi(12 - hi.abs().max(1e-300).log10().ceil() as i32); (hi * sc).round() / sc };
                    if r > t && r < hi && self.changed(hist, ra, prev, prevc, s, r) { hi = r; }
                    te = te.min(hi);
                }
            }
            if at(w, t, tn.min(te), tn < te) { return true; }
            let firing = tn < te;
            if te <= tn { // an event time: no reaction fires; the next waiting time is drawn afresh (memoryless)
                t = te;
            } else {
                t = tn;
                let (mut r, mut j) = (uniform() * a0, nr);
                for jj in 0..nr { // ponytail: linear search + full propensity recompute; dependency graph if nr is large
                    let a = w.v[jj].max(0.0);
                    if a > 0.0 { j = jj; if r <= a { break; } r -= a; }
                }
                self.fire(w, j);
            }
            if self.has_events() {
                match self.events_at(w, t) { None => return false, Some(true) => resync(w), _ => {} }
            }
            if firing {
                if let Some(time) = one.as_deref_mut() {
                    *time = t;
                    return true;
                }
            }
        }
    }

    /// Reaction j fires (w.ctl.nonneg: not if an amount would go negative, as libroadrunner's
    /// nonnegative).
    fn fire(&self, w: &mut Work, j: usize) {
        if w.ctl.nonneg && self.col[j].iter().any(|&(i, n)| w.amt[i as usize] + n < 0.0) { return; }
        for &(i, n) in &self.col[j] {
            let i = i as usize;
            w.amt[i] += n;
            w.s[i] = w.amt[i] / self.vol[i];
        }
    }

    /// SSA with rate rules (hybrid SSA, integrated-propensity method): the rate-rule states (columns
    /// nr..nv) and A(t) = int a0 dt are integrated deterministically (RK4 with step-doubling error
    /// control) from each segment start; a reaction fires when A reaches r = -ln U. A segment ends
    /// at a firing, a trigger change (bisection), a pending event or the next `grid` time; `at`
    /// then sees the state at the segment start. Algebraic states are not supported (false).
    fn ssa_hybrid(&self, w: &mut Work, mut t: f64, uniform: &mut impl FnMut() -> f64, grid: &[f64],
                  mut one: Option<&mut f64>, mut at: impl FnMut(&mut Work, f64, f64, bool) -> bool) -> bool {
        if self.alg.contains(&true) { return false; }
        let (ns, nr, ev) = (self.ns, self.nr, self.has_events());
        let (rtol, atol) = (1e-8, 1e-10);
        let rrc: Vec<(usize, usize)> = self.nz.iter().filter(|e| e.1 as usize >= nr).map(|e| (e.0 as usize, e.1 as usize)).collect();
        let nz = rrc.len() + 1; // rate-rule states, then A
        let load = |w: &mut Work, z: &[f64]| for (&(i, _), &x) in rrc.iter().zip(z) { w.s[i] = x };
        let f = |w: &mut Work, t: f64, z: &[f64], dz: &mut [f64]| {
            load(w, z);
            let Work { rr, v, s, hist, .. } = w;
            self.rates_at(hist, rr, v, s, t);
            for (d, &(_, j)) in dz.iter_mut().zip(&rrc) { *d = v[j]; }
            dz[nz - 1] = v[..nr].iter().map(|&a| a.max(0.0)).sum();
        };
        let rk4 = |w: &mut Work, t: f64, h: f64, z: &[f64], o: &mut [f64], b: &mut [f64]| {
            let (k, zt) = b.split_at_mut(4 * nz);
            let (k1, k) = k.split_at_mut(nz);
            let (k2, k) = k.split_at_mut(nz);
            let (k3, k4) = k.split_at_mut(nz);
            f(w, t, z, k1);
            for i in 0..nz { zt[i] = z[i] + 0.5 * h * k1[i]; }
            f(w, t + 0.5 * h, zt, k2);
            for i in 0..nz { zt[i] = z[i] + 0.5 * h * k2[i]; }
            f(w, t + 0.5 * h, zt, k3);
            for i in 0..nz { zt[i] = z[i] + h * k3[i]; }
            f(w, t + h, zt, k4);
            for i in 0..nz { o[i] = z[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]); }
        };
        // segment end condition at (t, z): A reached r or a trigger (or a comparison inside one) changed
        let hit = |w: &mut Work, t: f64, z: &[f64], r: f64| z[nz - 1] >= r || ev && {
            load(w, z);
            let Work { ra, s, prev, prevc, hist, .. } = w;
            self.changed(hist, ra, prev, prevc, s, t)
        };
        let (mut z, mut z1, mut z2, mut zh) = (vec![0.0; nz], vec![0.0; nz], vec![0.0; nz], vec![0.0; nz]);
        let mut b = vec![0.0; 5 * nz];
        let (mut h, mut steps) = (0.0f64, 0usize);
        loop {
            let tg = grid.get(grid.partition_point(|&g| g <= t)).copied().unwrap_or(f64::INFINITY);
            let cap = w.pend.iter().fold(tg, |m, p| m.min(p.0));
            if cap == f64::INFINITY { return at(w, t, cap, false); } // past the last grid time: `at` finishes
            let y0 = w.s.clone();
            let r = -uniform().ln();
            for (x, &(i, _)) in z.iter_mut().zip(&rrc) { *x = w.s[i]; }
            z[nz - 1] = 0.0;
            let mut tc = t;
            if h <= 0.0 { h = 1e-3 * (cap - t); }
            while tc < cap { // adaptive RK4 to the first hit or cap
                steps += 1;
                if steps > MAX_STEPS * 100 { return false; }
                if w.ctl.hmax > w.ctl.hmin { h = h.min(w.ctl.hmax); }
                let last = h >= cap - tc;
                if last { h = cap - tc; }
                rk4(w, tc, h, &z, &mut z1, &mut b);
                rk4(w, tc, 0.5 * h, &z, &mut zh, &mut b);
                rk4(w, tc + 0.5 * h, 0.5 * h, &zh, &mut z2, &mut b);
                let err = (0..nz).map(|i| (z2[i] - z1[i]).abs() / 15.0 / (atol + rtol * z[i].abs().max(z2[i].abs()))).fold(0.0, f64::max);
                if !err.is_finite() { return false; }
                if err > 1.0 {
                    h *= (0.9 * err.powf(-0.2)).max(0.2);
                    if h < 1e-14 * tc.abs().max(1.0) { return false; }
                    continue;
                }
                let tn = if last { cap } else { tc + h };
                if hit(w, tn, &z2, r) { // bisect on the step length from (tc, z)
                    let (mut lo, mut hi) = (0.0, tn - tc);
                    while hi - lo > 1e-12 * tn.abs().max(tn - tc) {
                        let mid = 0.5 * (lo + hi);
                        rk4(w, tc, mid, &z, &mut z1, &mut b);
                        if hit(w, tc + mid, &z1, r) { hi = mid } else { lo = mid }
                    }
                    // a trigger like k > 0.5 or time >= 25 holds at a round time: prefer it if it already flips
                    let rt = round12(tc + hi) - tc;
                    if rt > lo && rt < hi { rk4(w, tc, rt, &z, &mut z1, &mut b); if hit(w, tc + rt, &z1, r) { hi = rt } }
                    if hi < tn - tc { rk4(w, tc, hi, &z, &mut z2, &mut b); }
                    tc += hi;
                    z.copy_from_slice(&z2);
                    break;
                }
                tc = tn;
                z.copy_from_slice(&z2);
                h *= (0.9 * err.max(1e-10).powf(-0.2)).min(4.0);
            }
            // report the segment start, then move to its end
            let mut y1 = std::mem::replace(&mut w.s, y0);
            if at(w, t, tc, z[nz - 1] >= r) { return true; }
            std::mem::swap(&mut w.s, &mut y1);
            load(w, &z);
            t = tc;
            let firing = z[nz - 1] >= r;
            if firing { // a reaction fires at t, chosen by the propensities at t
                let Work { rr, v, s, hist, .. } = w;
                self.rates_at(hist, rr, v, s, t);
                let a0: f64 = v[..nr].iter().map(|&a| a.max(0.0)).sum();
                let (mut u, mut j) = (uniform() * a0, nr);
                for jj in 0..nr {
                    let a = w.v[jj].max(0.0);
                    if a > 0.0 { j = jj; if u <= a { break; } u -= a; }
                }
                if j < nr { self.fire(w, j); }
            }
            if ev {
                match self.events_at(w, t) { None => return false, Some(true) => for i in 0..ns { w.amt[i] = w.s[i] * self.vol[i] }, _ => {} }
            }
            if firing {
                if let Some(time) = one.as_deref_mut() {
                    *time = t;
                    return true;
                }
            }
        }
    }
}

/// Raw output buffer shared across threads; each sample writes its own disjoint row.
pub struct Rows { p: *mut f64, stride: usize }
unsafe impl Sync for Rows {}
unsafe impl Send for Rows {}
impl Rows {
    pub fn new(v: &mut [f64], stride: usize) -> Rows { Rows { p: v.as_mut_ptr(), stride } }
    /// SAFETY: callers use each row index from exactly one thread, within the buffer's bounds.
    #[allow(clippy::mut_from_ref)]
    pub unsafe fn row(&self, i: usize) -> &mut [f64] {
        std::slice::from_raw_parts_mut(self.p.add(i * self.stride), self.stride)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> { v.iter().map(|x| x.to_string()).collect() }

    // S1 -> S2 (k1*S1), S2 -> S1 (k2*S2): one moiety S1+S2.  Plus E+S <-> ES -> E+P style moiety.
    fn moiety() -> Model {
        // species: S1 S2 E ES P ; reactions: R1 S1->S2, R2 S2->S1, R3 E+S1->ES, R4 ES->E+S2, R5 $X->S1, R6 S2->
        let st = vec![
            -1., 1., -1., 0., 1., 0.,  // S1
             1., -1., 0., 1., 0., -1., // S2
             0., 0., -1., 1., 0., 0.,  // E
             0., 0., 1., -1., 0., 0.,  // ES
             0., 0., 0., 0., 0., 0.,   // P (no reactions)
        ];
        Model::new(s(&["S1", "S2", "E", "ES", "P"]), vec![1., 2., 1., 0.5, 3.], vec![1.; 5],
                   s(&["R1", "R2", "R3", "R4", "R5", "R6"]),
                   s(&["k1*S1", "k2*S2", "k3*E*S1", "k4*ES", "k5", "k6*S2/(1+S2)"]), st,
                   s(&["k1", "k2", "k3", "k4", "k5", "k6"]), vec![1., 0.5, 2., 1.5, 1., 2.], Rules::default()).unwrap()
    }

    #[test]
    fn conservation() {
        let m = moiety();
        // E + ES conserved, P has no reactions -> 2 dependent species
        assert_eq!(m.ni(), 3);
        assert_eq!(m.dep.len(), 2);
        // every dependent row of N equals l0 . NR
        for (d, &i) in m.dep.iter().enumerate() {
            for j in 0..m.nr {
                let r: f64 = (0..m.ni()).map(|b| m.l0[d * m.ni() + b] * m.stoich[m.indep[b] * m.nr + j]).sum();
                assert!((r - m.stoich[i * m.nr + j]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn steady_state_and_simulation() {
        let m = moiety();
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        assert!(m.steady_state(&mut w));
        let mut f = vec![0.0; m.ns];
        m.rhs(&w.v, &mut f);
        assert!(inf_norm(&f) < 1e-12, "{f:?}");
        assert!((w.s[2] + w.s[3] - 1.5).abs() < 1e-12); // moiety preserved
        // long simulation reaches the same steady state
        let ss = w.s.clone();
        m.load(&mut w, &[], &[]);
        let mut h = 0.0;
        assert!(m.integrate(&mut w, 0.0, 200.0, &mut h, 1e-10, 1e-12));
        for i in 0..m.ns { assert!((w.s[i] - ss[i]).abs() < 1e-7, "{i} {} {}", w.s[i], ss[i]); }
    }

    #[test]
    fn rodas_order() {
        // y' = -k y (1 + sin t)-free: decay with time-dependent rate law k*S*(1+time)
        let m = Model::new(s(&["S"]), vec![1.0], vec![1.0], s(&["R"]), s(&["k*S*(1+time)"]),
                           vec![-1.0], s(&["k"]), vec![1.0], Rules::default()).unwrap();
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let mut h = 0.0;
        assert!(m.integrate(&mut w, 0.0, 2.0, &mut h, 1e-10, 1e-14));
        let exact = (-4.0f64).exp(); // y = exp(-(t + t^2/2))
        assert!((w.s[0] - exact).abs() < 1e-9 * exact.max(1e-3), "{} {}", w.s[0], exact);
    }

    fn decay(x: Rules) -> Model {
        Model::new(s(&["S"]), vec![10.0], vec![1.0], s(&["R"]), s(&["k*S"]), vec![-1.0], s(&["k", "p"]), vec![1.0, 0.0], x).unwrap()
    }

    #[test]
    fn dae_nonlinear_algebraic_state() {
        // S' = -k S (S(0) = 10); algebraic p: p^3 + p - S = 0, started inconsistent at p = 0
        let m = decay(Rules { rate: vec![("p".into(), "p^3 + p - S".into())], algebraic: vec!["p".into()], ..Default::default() });
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let times = [0.0, 0.5, 1.0, 2.0];
        let outs = m.output_regs(&["S".into(), "p".into()]).unwrap();
        let mut out = vec![0.0; times.len() * 2];
        assert!(m.simulate(&mut w, &times, 1e-8, 1e-10, &outs, &mut out));
        for (q, &t) in times.iter().enumerate() {
            let (sv, p) = (out[2 * q], out[2 * q + 1]);
            assert!((sv - 10.0 * (-t as f64).exp()).abs() < 1e-5, "S({t}) = {sv}");
            assert!((p * p * p + p - sv).abs() < 1e-5, "p({t}) = {p} not on p^3 + p = {sv}");
        }
    }

    #[test]
    fn dae_with_event() {
        // S' = -k S, algebraic p^3 + p - S = 0; at time >= 1, S = 10 (p must follow on the constraint)
        let ev = EventDef { trigger: "time >= 1".into(), delay: None, priority: None, assign: vec![("S".into(), "10".into())],
            initial_value: true, persistent: true, use_trigger_values: true, div: vec![], rescale: vec![] };
        let m = decay(Rules { rate: vec![("p".into(), "p^3 + p - S".into())], algebraic: vec!["p".into()], events: vec![ev], ..Default::default() });
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let times = [0.0, 0.5, 0.999, 1.0, 1.5, 2.0];
        let outs = m.output_regs(&["S".into(), "p".into()]).unwrap();
        let mut out = vec![0.0; times.len() * 2];
        assert!(m.simulate(&mut w, &times, 1e-8, 1e-10, &outs, &mut out));
        for (q, &t) in times.iter().enumerate() {
            let ex = if t < 1.0 { 10.0 * (-t as f64).exp() } else { 10.0 * (-(t - 1.0) as f64).exp() };
            let (sv, p) = (out[2 * q], out[2 * q + 1]);
            assert!((sv - ex).abs() < 1e-5, "S({t}) = {sv} vs {ex}");
            assert!((p * p * p + p - ex).abs() < 1e-5, "p({t}) = {p} not on p^3 + p = {ex}");
        }
    }

    #[test]
    fn ssa_rate_rule() {
        // k' = 1 (k(0) = 0) drives S -> with propensity k S: E[S(t)] = S0 exp(-t^2/2)
        let check = |m: &Model, t: f64, s0: f64, pr: f64| {
            let mut w = m.work();
            let outs = m.output_regs(&s(&["S", "k"])).unwrap();
            let (n, mut acc) = (4000, 0.0);
            for i in 0..n {
                m.load(&mut w, &[], &[]);
                let mut out = [0.0; 4];
                assert!(m.gillespie(&mut w, &[0.0, t], i, &outs, &mut out));
                assert!((out[3] - t).abs() < 1e-9, "k({t}) = {}", out[3]);
                acc += out[2];
            }
            let mean = acc / n as f64;
            let se = (s0 * pr * (1.0 - pr) / n as f64).sqrt();
            assert!((mean - s0 * pr).abs() < 4.0 * se, "{mean} vs {}", s0 * pr);
        };
        let rk = || vec![("k".to_string(), "1".to_string())];
        let mut m = decay(Rules { rate: rk(), ..Default::default() });
        m.init[1] = 0.0; // k(0) = 0
        check(&m, 1.0, 10.0, (-0.5f64).exp());
        // variable-step rows: one per firing, ending at t1 with the integrated k(t1) = t1
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let rows = m.gillespie_steps(&mut w, 0.0, 1.5, 7, &m.output_regs(&s(&["S", "k"])).unwrap(), 1000).unwrap();
        let last = &rows[rows.len() - 3..];
        assert!(last[0] == 1.5 && (last[2] - 1.5).abs() < 1e-9, "{last:?}");
        for r in rows.chunks(3) { assert!((r[2] - r[0]).abs() < 1e-9, "k(t) = t: {r:?}"); }
        // a trigger on the rate-rule state: k > 0.5 sets S = 20, then E[S(1)] = 20 exp(-(1 - 0.25)/2)
        let ev = EventDef { trigger: "k > 0.5".into(), delay: None, priority: None, assign: vec![("S".into(), "20".into())],
            initial_value: true, persistent: true, use_trigger_values: true, div: vec![], rescale: vec![] };
        let mut m = decay(Rules { rate: rk(), events: vec![ev], ..Default::default() });
        m.init[1] = 0.0;
        check(&m, 1.0, 20.0, (-0.375f64).exp());
    }

    #[test]
    fn dense_output() {
        let m = decay(Rules::default());
        let mut w = m.work();
        for h0 in [0.2, 0.1] {
            m.load(&mut w, &[], &[]);
            let mut h = h0;
            assert_eq!(m.step(&mut w, 0.0, h0, &mut h, 1.0, 1.0), Some(h0));
            for i in 0..m.nx { w.d1[i] = (0..5).map(|j| H[0][j] * w.k[j][i]).sum(); w.d2[i] = (0..5).map(|j| H[1][j] * w.k[j][i]).sum(); }
            let th = 0.3;
            let y = (1.0 - th) * w.y0[0] + th * (w.s[0] + (1.0 - th) * (w.d1[0] + th * w.d2[0]));
            let ex = 10.0 * (-th * h0).exp();
            assert!((y - ex).abs() < 2e-5 * h0.powi(4) * 1e4, "h {h0}: {y} vs {ex}");
        }
    }

    #[test]
    fn events_and_rate_rules() {
        let ev = |trig: &str, asg: &[(&str, &str)]| EventDef { trigger: trig.into(), delay: None, priority: None,
            assign: asg.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect(),
            initial_value: true, persistent: true, use_trigger_values: true, div: vec![], rescale: vec![] };
        // S < 5 -> S = 10, and p counts firings; p' rate rule-free parameter
        let m = decay(Rules { events: vec![ev("S < 5", &[("S", "10"), ("p", "p + 1")])], ..Default::default() });
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let t = [0.0, 2f64.ln() + 0.1, 3.0 * 2f64.ln() + 0.1];
        let outs = m.output_regs(&s(&["S", "p", "R"])).unwrap();
        let mut out = vec![0.0; 9];
        assert!(m.simulate(&mut w, &t, 1e-10, 1e-12, &outs, &mut out));
        let ex = 10.0 * (-0.1f64).exp();
        for (a, b) in [(out[3], ex), (out[4], 1.0), (out[6], ex), (out[7], 3.0), (out[8], ex)] {
            assert!((a - b).abs() < 1e-7, "{out:?}");
        }
        // rate rule on parameter p' = 2*time + k*0, steady state rejects (no steady state), simulate exact
        let m = decay(Rules { rate: vec![("p".into(), "2*time".into())], assign: vec![("q".into(), "p*S".into())], ..Default::default() });
        assert_eq!(m.state_ids, s(&["S", "p"]));
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let outs = m.output_regs(&s(&["p", "q"])).unwrap();
        let mut out = vec![0.0; 4];
        assert!(m.simulate(&mut w, &[0.0, 1.5], 1e-10, 1e-12, &outs, &mut out));
        assert!((out[2] - 2.25).abs() < 1e-9 && (out[3] - 2.25 * 10.0 * (-1.5f64).exp()).abs() < 1e-7, "{out:?}");
        let mut j = vec![0.0; 4];
        m.jacobian(&mut w, 0.0, &mut j);
        assert_eq!(j, vec![-1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn ssa_mean() {
        // pure decay: E[S(t)] = S0 exp(-k t)
        let m = decay(Rules::default());
        let mut w = m.work();
        let outs = m.output_regs(&s(&["S"])).unwrap();
        let (n, mut acc) = (4000, 0.0);
        for i in 0..n {
            m.load(&mut w, &[], &[]);
            let mut out = [0.0; 2];
            assert!(m.gillespie(&mut w, &[0.0, 1.0], i, &outs, &mut out));
            acc += out[1];
        }
        let (mean, ex) = (acc / n as f64, 10.0 * (-1.0f64).exp());
        let se = (10.0 * (-1.0f64).exp() * (1.0 - (-1.0f64).exp()) / n as f64).sqrt();
        assert!((mean - ex).abs() < 4.0 * se, "{mean} vs {ex}");
    }
}

impl Model {
    /// Tape computing the rate vector v (outs, nv) from registers [state | params | time].
    pub fn rates_tape(&self) -> &Tape { &self.rates }
}

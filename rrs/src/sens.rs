//! Forward parameter sensitivities, exact rate derivatives and variable-step Gillespie on top of a
//! compiled `Model`. The rates tape is lifted back into a DAG in which parameters are variables
//! (Var(nx + k)), so `Dag::diff` yields exact d v / d parameter as well as d v / d state.
//!
//! Sensitivities: for each parameter (in parallel) an augmented model with states [x | s] and
//! ds/dt = J(x) s + df/dp_k is compiled and integrated by the same Rodas4 (its Jacobian, including
//! the d(J s)/dx block, is again exact).
use crate::expr::{Dag, Op, Tape, N};
use crate::model::{Key, Model, Rules};
use rayon::prelude::*;
use std::fmt::Write;

fn remap(op: Op, f: impl Fn(N) -> N) -> Op {
    use Op::*;
    match op {
        Const(_) | Var(_) | Par(_) | Time => op,
        Neg(a) => Neg(f(a)), Add(a, c) => Add(f(a), f(c)), Sub(a, c) => Sub(f(a), f(c)),
        Mul(a, c) => Mul(f(a), f(c)), Div(a, c) => Div(f(a), f(c)), Pow(a, c) => Pow(f(a), f(c)),
        PowI(a, k) => PowI(f(a), k), Exp(a) => Exp(f(a)), Ln(a) => Ln(f(a)), Log10(a) => Log10(f(a)),
        Sqrt(a) => Sqrt(f(a)), Abs(a) => Abs(f(a)), Sign(a) => Sign(f(a)), Floor(a) => Floor(f(a)),
        Ceil(a) => Ceil(f(a)), Sin(a) => Sin(f(a)), Cos(a) => Cos(f(a)), Tan(a) => Tan(f(a)),
        Lt(a, c) => Lt(f(a), f(c)), Le(a, c) => Le(f(a), f(c)), Eq(a, c) => Eq(f(a), f(c)),
        Ne(a, c) => Ne(f(a), f(c)), And(a, c) => And(f(a), f(c)), Or(a, c) => Or(f(a), f(c)),
        Xor(a, c) => Xor(f(a), f(c)), Not(a) => Not(f(a)), Sel(a, c, d) => Sel(f(a), f(c), f(d)),
        Min(a, c) => Min(f(a), f(c)), Max(a, c) => Max(f(a), f(c)),
    }
}

/// Rates tape -> (DAG, rate nodes). Registers [x (nx) | params (np)] become Var(0..nx+np).
fn lift(m: &Model) -> (Dag, Vec<N>) {
    let t = m.rates_tape();
    let nl = t.ns + t.np;
    let mut dag = Dag::default();
    let mut node = vec![0 as N; t.nregs];
    for (r, n) in node.iter_mut().enumerate().take(nl) { *n = dag.mk(Op::Var(r as u32)); }
    node[nl] = dag.mk(Op::Time);
    for &(r, x) in &t.consts { node[r as usize] = dag.cst(x); }
    for &(dst, op) in t.pre.iter().chain(&t.ins) { node[dst as usize] = dag.mk(remap(op, |r| node[r as usize])); }
    (dag, t.outs.iter().map(|&r| node[r as usize]).collect())
}

/// Base values with key overrides -> (x, params).
fn state(m: &Model, keys: &[Key], vals: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let (mut x, mut p) = (m.init.clone(), m.pvals.clone());
    for (k, &v) in keys.iter().zip(vals) {
        match *k { Key::Par(i) => p[i] = v, Key::Spec(i) => x[i] = v }
    }
    (x, p)
}

/// Exact rate derivatives, compiled once per model: one tape for v and every nonzero d v / d [x | params].
pub struct RateDerivs { tape: Tape, at: Vec<usize> }

impl RateDerivs {
    pub fn new(m: &Model) -> RateDerivs {
        let nz = m.nx + m.pnames.len();
        let (mut dag, v) = lift(m);
        let mut roots = v.clone();
        let mut at = vec![];
        for (j, &r) in v.iter().enumerate() {
            let mut has_t = false;
            for k in dag.deps(r, &mut has_t) {
                let d = dag.diff(r, Op::Var(k));
                if dag.val(d) != Some(0.0) { at.push(j * nz + k as usize); roots.push(d); }
            }
        }
        RateDerivs { tape: Tape::new(&dag, &roots, nz, 0), at }
    }

    /// (v [nv], d v / d [x | params] row-major nv x (nx + np)) at the keyed state, time t.
    pub fn eval(&self, m: &Model, keys: &[Key], vals: &[f64], t: f64) -> (Vec<f64>, Vec<f64>) {
        let (nx, nv, nz) = (m.nx, m.nv, m.nx + m.pnames.len());
        let (x, p) = state(m, keys, vals);
        let tape = &self.tape;
        let mut regs = tape.regs();
        regs[..nx].copy_from_slice(&x);
        regs[nx..nz].copy_from_slice(&p);
        regs[nz] = t;
        tape.run_all(&mut regs);
        let mut d = vec![0.0; nv * nz];
        for (q, &i) in self.at.iter().enumerate() { d[i] = regs[tape.outs[nv + q] as usize]; }
        (tape.outs[..nv].iter().map(|&r| regs[r as usize]).collect(), d)
    }
}

const MAX_FORMULA: usize = 1 << 24;

/// Infix text of node n (parsable by expr::parse). Stops early once `o` exceeds MAX_FORMULA.
fn show(d: &Dag, n: N, nm: &dyn Fn(u32) -> String, o: &mut String) {
    use Op::*;
    if o.len() > MAX_FORMULA { return; }
    let f = |o: &mut String, name: &str, args: &[N]| {
        o.push_str(name);
        o.push('(');
        for (i, &a) in args.iter().enumerate() {
            if i > 0 { o.push(','); }
            show(d, a, nm, o);
        }
        o.push(')');
    };
    let bin = |o: &mut String, a: N, s: &str, c: N| { o.push('('); show(d, a, nm, o); o.push_str(s); show(d, c, nm, o); o.push(')'); };
    match d.nodes[n as usize] {
        Const(b) => { let _ = write!(o, "({:?})", f64::from_bits(b)); }
        Var(r) => o.push_str(&nm(r)),
        Par(_) => unreachable!("lifted DAGs have no Par leaves"),
        Time => o.push_str("time"),
        Neg(a) => { o.push_str("(-"); show(d, a, nm, o); o.push(')'); }
        Add(a, c) => bin(o, a, "+", c), Sub(a, c) => bin(o, a, "-", c),
        Mul(a, c) => bin(o, a, "*", c), Div(a, c) => bin(o, a, "/", c),
        Pow(a, c) => f(o, "pow", &[a, c]),
        PowI(a, k) => { o.push_str("pow("); show(d, a, nm, o); let _ = write!(o, ",{k})"); }
        Exp(a) => f(o, "exp", &[a]), Ln(a) => f(o, "ln", &[a]), Log10(a) => f(o, "log10", &[a]),
        Sqrt(a) => f(o, "sqrt", &[a]), Abs(a) => f(o, "abs", &[a]), Floor(a) => f(o, "floor", &[a]),
        Ceil(a) => f(o, "ceil", &[a]), Sin(a) => f(o, "sin", &[a]), Cos(a) => f(o, "cos", &[a]),
        Tan(a) => f(o, "tan", &[a]),
        Sign(a) => { o.push_str("piecewise(1,"); show(d, a, nm, o); o.push_str(">0,-1,"); show(d, a, nm, o); o.push_str("<0,0)"); }
        Lt(a, c) => f(o, "lt", &[a, c]), Le(a, c) => f(o, "leq", &[a, c]), Eq(a, c) => f(o, "eq", &[a, c]),
        Ne(a, c) => f(o, "neq", &[a, c]), And(a, c) => f(o, "and", &[a, c]), Or(a, c) => f(o, "or", &[a, c]),
        Xor(a, c) => f(o, "xor", &[a, c]), Not(a) => f(o, "not", &[a]),
        Sel(c, x, y) => f(o, "piecewise", &[x, c, y]),
        Min(a, c) => f(o, "min", &[a, c]), Max(a, c) => f(o, "max", &[a, c]),
    }
}

/// Augmented model for parameter k (None: no parameter term, for initial-value sensitivities):
/// states [x | s], s = d x / d p_k, all driven by rate rules on parameters named as the original
/// states (then `__s{i}`), so the original events and `Key`s apply unchanged. Events act on x only
/// (no sensitivity jumps at events).
pub fn sens_model(m: &Model, k: Option<usize>) -> Result<Model, String> {
    let (nx, nv, np) = (m.nx, m.nv, m.pnames.len());
    let nl = nx + np;
    let (mut dag, v) = lift(m);
    let mut g = Vec::with_capacity(nv); // directional derivative of v_j along (s, e_k)
    for &r in &v {
        let mut has_t = false;
        let mut acc = match k { Some(k) => dag.diff(r, Op::Var((nx + k) as u32)), None => dag.cst(0.0) };
        for l in dag.deps(r, &mut has_t) {
            if l as usize >= nx { continue; }
            let dl = dag.diff(r, Op::Var(l));
            let sl = dag.mk(Op::Var((nl + l as usize) as u32));
            let t = dag.mk(Op::Mul(dl, sl));
            acc = dag.mk(Op::Add(acc, t));
        }
        g.push(acc);
    }
    let sname = |i: usize| format!("__s{i}");
    let nm = |r: u32| { let r = r as usize; if r < nx { m.state_ids[r].clone() } else if r < nl { m.pnames[r - nx].clone() } else { sname(r - nl) } };
    let mut rate = Vec::with_capacity(2 * nx);
    for (sens, terms) in [(false, &v), (true, &g)] {
        for i in 0..nx {
            let mut acc = dag.cst(0.0);
            for (j, &t) in terms.iter().enumerate() {
                let c = m.stoich[i * nv + j];
                if c == 0.0 { continue; }
                let c = dag.cst(if j < m.nr { c / m.vol[i] } else { c });
                let ct = dag.mk(Op::Mul(c, t));
                acc = dag.mk(Op::Add(acc, ct));
            }
            let mut s = String::new();
            show(&dag, acc, &nm, &mut s);
            if s.len() > MAX_FORMULA { return Err("sensitivity equations too large".into()); }
            rate.push((if sens { sname(i) } else { m.state_ids[i].clone() }, s));
        }
    }
    let (mut pn, mut pv) = (m.pnames.clone(), m.pvals.clone());
    for (i, x) in m.state_ids.iter().enumerate() {
        if !m.pnames.contains(x) { pn.push(x.clone()); pv.push(m.init[i]); }
    }
    pn.extend((0..nx).map(sname));
    pv.resize(pn.len(), 0.0);
    Model::new(vec![], vec![], vec![], vec![], vec![], vec![], pn, pv,
               Rules { rate, events: m.event_src.clone(), ..Default::default() })
}

/// Forward sensitivities d x(t) / d p for p in `params`: a parameter (x(0) independent of it) or a
/// state's initial value. Returns (x [nt x nx], s [nt x len(params) x nx]); concentrations.
/// order 1 gives d s / dt instead of s. Parameters run in parallel.
#[allow(clippy::too_many_arguments)]
pub fn time_sens(m: &Model, keys: &[Key], vals: &[f64], times: &[f64], params: &[Key], rtol: f64, atol: f64,
                 order: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    if order > 1 { return Err("sensitivities: only k = 0, 1".into()); }
    let (nx, nt, npar) = (m.nx, times.len(), params.len());
    let fail = || "sensitivities: integration failed".to_string();
    let mut x = vec![0.0; nt * nx];
    let mut w = m.work();
    m.load(&mut w, keys, vals);
    if !m.simulate(&mut w, times, rtol, atol, &(0..nx as u32).collect::<Vec<_>>(), &mut x) { return Err(fail()); }
    let runs: Vec<Vec<f64>> = params.par_iter().map(|&k| {
        let a = sens_model(m, match k { Key::Par(i) => Some(i), Key::Spec(_) => None })?;
        let mut w = a.work();
        a.load(&mut w, keys, vals);
        if let Key::Spec(i) = k { w.s[nx + i] = 1.0; }
        // order 1: every register up to time (states, parameters incl. event-changed ones, time)
        let outs: Vec<u32> = if order == 0 { (nx as u32..2 * nx as u32).collect() } else { (0..=(a.nx + a.pnames.len()) as u32).collect() };
        let no = outs.len();
        let mut out = vec![0.0; nt * no];
        if !a.simulate(&mut w, times, rtol, atol, &outs, &mut out) { return Err(fail()); }
        if order == 0 { return Ok(out); }
        let tape = a.rates_tape();
        let mut regs = tape.regs();
        Ok((0..nt).flat_map(|t| {
            regs[..no].copy_from_slice(&out[t * no..(t + 1) * no]);
            tape.run_all(&mut regs);
            tape.outs[nx..2 * nx].iter().map(|&r| regs[r as usize]).collect::<Vec<_>>()
        }).collect())
    }).collect::<Result<_, String>>()?;
    let mut s = vec![0.0; nt * npar * nx];
    for (q, r) in runs.iter().enumerate() {
        for t in 0..nt { s[(t * npar + q) * nx..(t * npar + q + 1) * nx].copy_from_slice(&r[t * nx..(t + 1) * nx]); }
    }
    Ok((x, s))
}

#[cfg(feature = "python")]
mod py {
    use super::*;
    use crate::py::PyModel;
    use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;

    fn arr<'py>(py: Python<'py>, v: Vec<f64>, shape: &[usize]) -> PyResult<Bound<'py, PyAny>> {
        Ok(PyArray1::from_vec(py, v).reshape(shape)?.into_any())
    }

    fn resolve(m: &Model, keys: &[String], values: &[f64]) -> PyResult<Vec<Key>> {
        let k = m.resolve(keys).map_err(PyValueError::new_err)?;
        if k.len() != values.len() { return Err(PyValueError::new_err("values and keys differ in length")); }
        Ok(k)
    }

    /// time_sensitivities(model, keys, values, times, params, rtol=1e-6, atol=1e-12, order=0)
    ///   -> (x[len(times), nState], sens[len(times), len(params), nState])
    /// Forward sensitivities d state / d parameter (concentrations, exact symbolic derivatives),
    /// state at keys/values. A state id in params gives d state / d its initial value; order=1 gives
    /// their time derivative. Events act on the state only. One stiff integration per parameter, in parallel.
    #[pyfunction]
    #[pyo3(signature = (model, keys, values, times, params, rtol=1e-6, atol=1e-12, order=0))]
    #[allow(clippy::too_many_arguments)]
    fn time_sensitivities<'py>(py: Python<'py>, model: &Bound<'py, PyModel>, keys: Vec<String>, values: Vec<f64>,
                               times: PyReadonlyArray1<f64>, params: Vec<String>, rtol: f64, atol: f64, order: usize)
                               -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
        let m = &model.get().m;
        let k = resolve(m, &keys, &values)?;
        let times: Vec<f64> = times.as_array().iter().copied().collect();
        if times.is_empty() || times.windows(2).any(|t| !(t[1] >= t[0])) {
            return Err(PyValueError::new_err("times must be non-empty and non-decreasing"));
        }
        let ps = m.resolve(&params).map_err(PyValueError::new_err)?;
        let (x, s) = py.detach(|| time_sens(m, &k, &values, &times, &ps, rtol, atol, order)).map_err(PyValueError::new_err)?;
        let (nt, nx) = (times.len(), m.nx);
        Ok((arr(py, x, &[nt, nx])?, arr(py, s, &[nt, ps.len(), nx])?))
    }

    /// RateDerivatives(model): exact rate derivatives compiled once.
    /// rd(keys, values, time=0.0) -> (v[nRates], dv_dx[nRates, nState], dv_dp[nRates, nParams], param_names)
    #[pyclass(frozen, name = "RateDerivatives")]
    struct PyRateDerivs { model: Py<PyModel>, rd: RateDerivs }

    #[pymethods]
    impl PyRateDerivs {
        #[new]
        fn new(model: Bound<'_, PyModel>) -> Self {
            let rd = RateDerivs::new(&model.get().m);
            PyRateDerivs { model: model.unbind(), rd }
        }

        #[pyo3(signature = (keys, values, time=0.0))]
        #[allow(clippy::type_complexity)]
        fn __call__<'py>(&self, py: Python<'py>, keys: Vec<String>, values: Vec<f64>, time: f64)
                         -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>, Bound<'py, PyAny>, Vec<String>)> {
            let m = &self.model.get().m;
            let k = resolve(m, &keys, &values)?;
            let (v, d) = self.rd.eval(m, &k, &values, time);
            let (nv, nx, nz) = (m.nv, m.nx, m.nx + m.pnames.len());
            let (dx, dp): (Vec<f64>, Vec<f64>) = (
                (0..nv).flat_map(|j| d[j * nz..j * nz + nx].to_vec()).collect(),
                (0..nv).flat_map(|j| d[j * nz + nx..(j + 1) * nz].to_vec()).collect());
            Ok((arr(py, v, &[nv])?, arr(py, dx, &[nv, nx])?, arr(py, dp, &[nv, nz - nx])?, m.pnames.clone()))
        }
    }

    /// gillespie_steps(model, keys, values, start, end, seed, outputs=None, max_rows=100000)
    ///   -> rows[n, 1 + len(outputs)] (time, outputs)
    /// Variable-step SSA: a row at start, after each reaction firing or event, and at end. Outputs
    /// as Model.output_ids (default: floating species, concentrations). No rate rules.
    #[pyfunction]
    #[pyo3(signature = (model, keys, values, start, end, seed, outputs=None, max_rows=100000))]
    #[allow(clippy::too_many_arguments)]
    fn gillespie_steps<'py>(py: Python<'py>, model: &Bound<'py, PyModel>, keys: Vec<String>, values: Vec<f64>,
                            start: f64, end: f64, seed: u64, outputs: Option<Vec<String>>, max_rows: usize)
                            -> PyResult<Bound<'py, PyAny>> {
        let m = &model.get().m;
        if m.nv != m.nr { return Err(PyValueError::new_err("gillespie: rate rules are not supported")); }
        let k = resolve(m, &keys, &values)?;
        let outs = match outputs {
            None => (0..m.ns as u32).collect(),
            Some(ids) => m.output_regs(&ids).map_err(PyValueError::new_err)?,
        };
        let rows = py.detach(|| {
            let mut w = m.work();
            m.load(&mut w, &k, &values);
            m.gillespie_steps(&mut w, start, end, seed, &outs, max_rows)
        }).ok_or_else(|| PyValueError::new_err("gillespie: failed"))?;
        let no = outs.len() + 1;
        arr(py, rows.clone(), &[rows.len() / no, no])
    }

    pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(time_sensitivities, m)?)?;
        m.add_class::<PyRateDerivs>()?;
        m.add_function(wrap_pyfunction!(gillespie_steps, m)?)
    }
}
#[cfg(feature = "python")]
pub use py::register;

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> { v.iter().map(|x| x.to_string()).collect() }

    // S1 -> S2 (k1*S1/(Km+S1)), S2 -> (k2*S2), volume 2
    fn mm() -> Model {
        Model::new(s(&["S1", "S2"]), vec![5.0, 0.0], vec![2.0, 2.0], s(&["J1", "J2"]),
                   s(&["k1*S1/(Km+S1)", "k2*S2*(1+0*time) + k1*Km"]), vec![-1.0, 0.0, 1.0, -1.0],
                   s(&["k1", "Km", "k2"]), vec![2.0, 0.5, 0.3], Rules::default()).unwrap()
    }

    #[test]
    fn rate_derivs_match_finite_difference() {
        let m = mm();
        let rd = RateDerivs::new(&m);
        let (v, d) = rd.eval(&m, &[], &[], 0.0);
        assert!((v[0] - 2.0 * 5.0 / 5.5).abs() < 1e-14);
        let nz = m.nx + m.pnames.len();
        for k in 0..nz {
            let key = if k < m.nx { Key::Spec(k) } else { Key::Par(k - m.nx) };
            let base = if k < m.nx { m.init[k] } else { m.pvals[k - m.nx] };
            let h = 1e-6 * base.abs().max(1e-3);
            let (vp, _) = rd.eval(&m, &[key], &[base + h], 0.0);
            let (vm, _) = rd.eval(&m, &[key], &[base - h], 0.0);
            for j in 0..m.nv {
                let fd = (vp[j] - vm[j]) / (2.0 * h);
                assert!((d[j * nz + k] - fd).abs() < 1e-6 * (1.0 + fd.abs()), "j {j} k {k}: {} vs {fd}", d[j * nz + k]);
            }
        }
    }

    #[test]
    fn sensitivities_match_finite_difference() {
        let m = mm();
        let times: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let (x, sv) = time_sens(&m, &[], &[], &times, &[Key::Par(0), Key::Par(1), Key::Par(2), Key::Spec(0)], 1e-10, 1e-14, 0).unwrap();
        for q in 0..4 {
            let (key, p0) = if q < 3 { (Key::Par(q), m.pvals[q]) } else { (Key::Spec(0), m.init[0]) };
            let h = 1e-5 * p0;
            let run = |p: f64| time_sens(&m, &[key], &[p], &times, &[], 1e-12, 1e-16, 0).unwrap().0;
            let (xp, xm) = (run(p0 + h), run(p0 - h));
            for i in 0..times.len() * m.nx {
                let fd = (xp[i] - xm[i]) / (2.0 * h);
                let (t, st) = (i / m.nx, i % m.nx);
                let a = sv[(t * 4 + q) * m.nx + st];
                assert!((a - fd).abs() < 1e-5 * (1.0 + fd.abs()), "p{q} t{t} x{st}: {a} vs {fd} (x {})", x[i]);
            }
        }
    }

    #[test]
    fn sensitivity_time_derivative() {
        let m = mm();
        let times: Vec<f64> = (0..41).map(|i| i as f64 * 0.1).collect();
        let (_, s0) = time_sens(&m, &[], &[], &times, &[Key::Par(0)], 1e-12, 1e-16, 0).unwrap();
        let (_, s1) = time_sens(&m, &[], &[], &times, &[Key::Par(0)], 1e-12, 1e-16, 1).unwrap();
        for t in 1..40 {
            for i in 0..m.nx {
                let fd = (s0[(t + 1) * m.nx + i] - s0[(t - 1) * m.nx + i]) / 0.2;
                assert!((s1[t * m.nx + i] - fd).abs() < 1e-3 * (1.0 + fd.abs()), "t{t} x{i}: {} vs {fd}", s1[t * m.nx + i]);
            }
        }
    }

    #[test]
    fn ssa_steps_mean() {
        // S1 -> S2 -> : E[S1(t)] = S1_0 exp(-k t) for linear decay
        let m = Model::new(s(&["S"]), vec![10.0], vec![1.0], s(&["R"]), s(&["k*S"]), vec![-1.0],
                           s(&["k"]), vec![1.0], Rules::default()).unwrap();
        let (n, mut acc) = (4000, 0.0);
        for i in 0..n {
            let mut w = m.work();
            m.load(&mut w, &[], &[]);
            let r = m.gillespie_steps(&mut w, 0.0, 1.0, i, &[0], 100_000).unwrap();
            let (t, x): (Vec<f64>, Vec<f64>) = (r.iter().step_by(2).copied().collect(), r.iter().skip(1).step_by(2).copied().collect());
            assert!(x.windows(2).take(x.len() - 2).all(|w| w[1] == w[0] - 1.0)); // one firing per row
            assert_eq!(t[0], 0.0);
            assert_eq!(*t.last().unwrap(), 1.0);
            assert!(t.windows(2).all(|w| w[1] >= w[0]));
            acc += x[x.len() - 1];
        }
        let (mean, ex) = (acc / n as f64, 10.0 * (-1.0f64).exp());
        let se = (ex * (1.0 - (-1.0f64).exp()) / n as f64).sqrt();
        assert!((mean - ex).abs() < 4.0 * se, "{mean} vs {ex}");
    }
}

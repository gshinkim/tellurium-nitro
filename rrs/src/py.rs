//! Python bindings (module `roadrunner_nitro`).
use crate::model::{EventDef, Key, Model, Rows, Rules, Work};
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString};
use rayon::prelude::*;

/// Compiled reaction network (Rust engine).
///
/// Model(species, init, volumes, reactions, rates, stoich, param_names, param_values,
///       rate_rules=None, assignment_rules=None, events=None)
///   species/init/volumes: floating species, initial concentrations, compartment volumes.
///   rates: kinetic-law formulas (amount/time) over species concentrations, parameters, `time`.
///   stoich: (nSpecies, nReactions).
///   rate_rules: [(target, formula)]; target is a reaction-less species (formula: conc/time) or a
///     parameter, which then becomes ODE state (appended to state_ids).
///   assignment_rules: [(id, formula)]: extra output variables.
///   events: [(trigger, [(target, formula)], delay|None, priority|None, initial_value, persistent,
///     use_values_from_trigger_time)]; targets are species (concentration) or parameters.
///
/// state_ids = species + promoted rate-rule parameters (columns of states, rows of `jacobian`).
/// output_ids = state_ids + assignment-rule ids + reaction ids (rates) + parameters + 'time'.
/// Formulas must not reference assignment-rule ids (inline them).
#[pyclass(frozen, name = "Model")]
pub(crate) struct PyModel { pub(crate) m: Model, pool: std::sync::Mutex<Vec<Work>>, names: std::sync::OnceLock<Vec<Py<PyString>>> }

fn verr(e: String) -> PyErr { PyValueError::new_err(e) }

/// Pre-resolved key or output ids (Model.resolve_keys / Model.resolve_outputs): skips per-call
/// string conversion and lookup in simulate_one.
#[pyclass(frozen)]
struct Resolved { keys: Vec<Key>, outs: Vec<u32> }

fn keys_of(m: &Model, k: &Bound<'_, PyAny>) -> PyResult<Vec<Key>> {
    match k.cast::<Resolved>() {
        Ok(r) => Ok(r.get().keys.clone()),
        Err(_) => m.resolve(&k.extract::<Vec<String>>()?).map_err(verr),
    }
}

type PyEvent = (String, Vec<(String, String)>, Option<String>, Option<String>, bool, bool, bool, Vec<Option<String>>, Vec<(String, String)>);

/// Run f(work, sample) for 0..n: serially for threads=1, else on rayon (all cores for None).
fn run_batch(m: &Model, n: usize, threads: Option<usize>, f: impl Fn(&mut Work, usize) + Sync + Send) -> PyResult<()> {
    match threads {
        Some(0) => Err(verr("threads must be >= 1".into())),
        Some(1) => { let mut w = m.work(); (0..n).for_each(|i| f(&mut w, i)); Ok(()) }
        None => { m.par(n, f); Ok(()) }
        Some(t) => {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(t).build().map_err(|e| verr(e.to_string()))?;
            pool.install(|| (0..n).into_par_iter().with_min_len(4).for_each_init(|| m.work(), |w, i| f(w, i)));
            Ok(())
        }
    }
}

fn arr<'py>(py: Python<'py>, v: Vec<f64>, shape: &[usize]) -> PyResult<Bound<'py, PyAny>> {
    Ok(PyArray1::from_vec(py, v).reshape(shape)?.into_any())
}

struct Batch { keys: Vec<Key>, vals: Vec<f64>, n: usize }

fn check_times(times: PyReadonlyArray1<f64>) -> PyResult<Vec<f64>> {
    let times: Vec<f64> = times.as_array().iter().copied().collect();
    if times.is_empty() || times.windows(2).any(|t| !(t[1] >= t[0])) {
        return Err(verr("times must be non-empty and non-decreasing".into()));
    }
    Ok(times)
}

impl PyModel {
    fn batch(&self, keys: Vec<String>, values: PyReadonlyArray2<f64>) -> PyResult<Batch> {
        let keys = self.m.resolve(&keys).map_err(verr)?;
        let a = values.as_array();
        let (n, k) = a.dim();
        if k != keys.len() { return Err(verr(format!("values has {k} columns but {} keys", keys.len()))); }
        Ok(Batch { keys, vals: a.iter().copied().collect(), n })
    }

    /// Output registers for `outputs` (default: floating species).
    fn outs(&self, outputs: Option<Vec<String>>) -> PyResult<Vec<u32>> {
        match outputs {
            None => Ok((0..self.m.ns as u32).collect()),
            Some(ids) => self.m.output_regs(&ids).map_err(verr),
        }
    }

    /// Time courses (ODE, or SSA with seed) for every sample -> (n, T, len(outs)), NaN rows on failure.
    #[allow(clippy::too_many_arguments)]
    fn time_courses<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                         times: PyReadonlyArray1<f64>, threads: Option<usize>, outputs: Option<Vec<String>>,
                         ode: Option<(f64, f64)>, seed: u64, atol_scale: Option<Vec<f64>>) -> PyResult<Bound<'py, PyAny>> {
        let b = self.batch(keys, values)?;
        if atol_scale.as_ref().is_some_and(|a| a.len() != self.m.nx) { return Err(verr("atol_scale must have one entry per state".into())); }
        let times = check_times(times)?;
        let outs = self.outs(outputs)?;
        let m = &self.m;
        let (no, nt, n, nk) = (outs.len(), times.len(), b.n, b.keys.len());
        let mut out = vec![f64::NAN; n * nt * no];
        py.detach(|| {
            let rows = Rows::new(&mut out, nt * no);
            run_batch(m, n, threads, |w, i| {
                let r = unsafe { rows.row(i) }; // SAFETY: one row per sample
                m.load(w, &b.keys, &b.vals[i * nk..(i + 1) * nk]);
                match &atol_scale { Some(a) => w.atol_scale.copy_from_slice(a), None => w.atol_scale.fill(1.0) }
                let ok = match ode {
                    Some((rtol, atol)) => { m.reseed(w, seed.wrapping_add(i as u64)); m.simulate(w, &times, rtol, atol, &outs, r) }
                    None => m.gillespie(w, &times, seed.wrapping_add(i as u64), &outs, r),
                };
                if !ok { r.fill(f64::NAN); }
            })
        })?;
        arr(py, out, &[n, nt, no])
    }

    /// Steady states (+ optional MCA) for every sample. Returns flat buffers.
    fn solve(&self, py: Python<'_>, b: &Batch, threads: Option<usize>, mca: bool) -> PyResult<Vec<Vec<f64>>> {
        let m = &self.m;
        let (nx, ns, nr, n, nk) = (m.nx, m.ns, m.nr, b.n, b.keys.len());
        let sizes = [nx, nr, nr * nx, 1, ns * nr, nr * nr, ns * nr, nr * nr, nr * ns];
        let nbuf = if mca { 9 } else { 4 };
        let mut bufs: Vec<Vec<f64>> = sizes[..nbuf].iter().map(|&s| vec![f64::NAN; s * n]).collect();
        py.detach(|| {
            let rows: Vec<Rows> = bufs.iter_mut().zip(sizes).map(|(v, s)| Rows::new(v, s)).collect();
            run_batch(m, n, threads, |w, i| {
                // SAFETY: sample i writes only row i of each buffer.
                let r = |q: usize| unsafe { rows[q].row(i) };
                m.load(w, &b.keys, &b.vals[i * nk..(i + 1) * nk]);
                let mut ok = m.steady_state(w);
                if ok && mca {
                    ok = m.mca(w, r(2), r(4), r(5), r(6), r(7), r(8));
                } else if ok {
                    m.elasticity_into(w, r(2));
                }
                if ok {
                    r(0).copy_from_slice(&w.s);
                    r(1).copy_from_slice(&w.v[..nr]);
                    r(3)[0] = 1.0;
                } else {
                    for q in [2, 4, 5, 6, 7, 8] { if q < rows.len() { r(q).fill(f64::NAN); } }
                }
            })
        })?;
        Ok(bufs)
    }
}

#[pymethods]
impl PyModel {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (species, init, volumes, reactions, rates, stoich, param_names, param_values,
                        rate_rules=None, assignment_rules=None, events=None, delays=None, algebraic=None))]
    fn new(species: Vec<String>, init: Vec<f64>, volumes: Vec<f64>, reactions: Vec<String>, rates: Vec<String>,
           stoich: PyReadonlyArray2<f64>, param_names: Vec<String>, param_values: Vec<f64>,
           rate_rules: Option<Vec<(String, String)>>, assignment_rules: Option<Vec<(String, String)>>,
           events: Option<Vec<PyEvent>>, delays: Option<Vec<(String, String, String, String)>>,
           algebraic: Option<Vec<String>>) -> PyResult<Self> {
        let a = stoich.as_array();
        if a.dim() != (species.len(), reactions.len()) {
            return Err(verr(format!("stoich shape {:?} != (nSpecies, nReactions)", a.dim())));
        }
        let st = a.iter().copied().collect();
        let x = Rules {
            rate: rate_rules.unwrap_or_default(),
            assign: assignment_rules.unwrap_or_default(),
            events: events.unwrap_or_default().into_iter().map(|(trigger, assign, delay, priority, initial_value, persistent, u, div, rescale)|
                EventDef { trigger, assign, delay, priority, initial_value, persistent, use_trigger_values: u, div, rescale }).collect(),
            delays: delays.unwrap_or_default(),
            algebraic: algebraic.unwrap_or_default(),
        };
        Model::new(species, init, volumes, reactions, rates, st, param_names, param_values, x)
            .map(|m| PyModel { m, pool: Default::default(), names: Default::default() }).map_err(verr)
    }

    /// Resolved keys for simulate_one (same as passing the id list, without per-call lookup).
    fn resolve_keys(&self, keys: Vec<String>) -> PyResult<Resolved> {
        Ok(Resolved { keys: self.m.resolve(&keys).map_err(verr)?, outs: vec![] })
    }

    /// Resolved outputs for simulate_one.
    fn resolve_outputs(&self, outputs: Vec<String>) -> PyResult<Resolved> {
        Ok(Resolved { keys: vec![], outs: self.m.output_regs(&outputs).map_err(verr)? })
    }

    #[getter]
    fn state_ids(&self) -> Vec<String> { self.m.state_ids.clone() }
    #[getter]
    fn output_ids(&self) -> Vec<String> { self.m.outputs.iter().map(|o| o.0.clone()).collect() }

    #[getter]
    fn species(&self) -> Vec<String> { self.m.species.clone() }
    #[getter]
    fn reactions(&self) -> Vec<String> { self.m.reactions.clone() }

    /// (NR, L, indep): reduced stoichiometry, link matrix (rows in state_ids order), independent states.
    /// With rate rules the columns are reactions + one unit column per rate rule.
    fn conservation<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>, Vec<usize>)> {
        let m = &self.m;
        let (ns, nr, ni) = (m.nx, m.nv, m.ni());
        let nrm: Vec<f64> = m.indep.iter().flat_map(|&i| m.stoich[i * nr..(i + 1) * nr].iter().copied()).collect();
        let mut l = vec![0.0; ns * ni];
        for (b, &i) in m.indep.iter().enumerate() { l[i * ni + b] = 1.0; }
        for (d, &i) in m.dep.iter().enumerate() { l[i * ni..(i + 1) * ni].copy_from_slice(&m.l0[d * ni..(d + 1) * ni]); }
        Ok((arr(py, nrm, &[ni, nr])?, arr(py, l, &[ns, ni])?, m.indep.clone()))
    }

    /// (S[n, nState], v[n, nReactions], E[n, nReactions, nState], ok[n]) at steady state for each
    /// sample row (events ignored, rate-rule states solved too; ok=False if singular/not found).
    #[pyo3(signature = (keys, values, threads=None))]
    fn steady_state_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                               threads: Option<usize>) -> PyResult<Bound<'py, PyAny>> {
        let b = self.batch(keys, values)?;
        let (ns, nr, n) = (self.m.nx, self.m.nr, b.n);
        let mut o = self.solve(py, &b, threads, false)?.into_iter();
        let mut nx = || o.next().unwrap();
        let (s, v, e, ok) = (nx(), nx(), nx(), nx());
        let ok: Vec<bool> = ok.iter().map(|&x| x == 1.0).collect();
        Ok((arr(py, s, &[n, ns])?, arr(py, v, &[n, nr])?, arr(py, e, &[n, nr, ns])?, PyArray1::from_vec(py, ok))
            .into_pyobject(py)?.into_any())
    }

    /// dict(ccc, fcc, elasticities, uccc, ufcc, S, v, ok). S is the full state (state_ids); the
    /// coefficients cover species x reactions. With rate rules, as roadrunner: steady state of the
    /// whole system, then MCA of the reaction network with rate-rule states held at their values.
    #[pyo3(signature = (keys, values, threads=None))]
    fn mca_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                      threads: Option<usize>) -> PyResult<Bound<'py, PyDict>> {
        let b = self.batch(keys, values)?;
        let (nx, ns, nr, n) = (self.m.nx, self.m.ns, self.m.nr, b.n);
        let o = self.solve(py, &b, threads, true)?;
        let [s, v, _e, ok, uc, uf, cc, fc, el]: [Vec<f64>; 9] = o.try_into().unwrap();
        let d = PyDict::new(py);
        d.set_item("ccc", arr(py, cc, &[n, ns, nr])?)?;
        d.set_item("fcc", arr(py, fc, &[n, nr, nr])?)?;
        d.set_item("elasticities", arr(py, el, &[n, nr, ns])?)?;
        d.set_item("uccc", arr(py, uc, &[n, ns, nr])?)?;
        d.set_item("ufcc", arr(py, uf, &[n, nr, nr])?)?;
        d.set_item("S", arr(py, s, &[n, nx])?)?;
        d.set_item("v", arr(py, v, &[n, nr])?)?;
        d.set_item("ok", PyArray1::from_vec(py, ok.iter().map(|&x| x == 1.0).collect::<Vec<bool>>()))?;
        Ok(d)
    }

    /// (n, len(times), len(outputs)) time courses (solver as in simulate, events root-located);
    /// outputs default to the floating species (concentrations). NaN rows where integration failed.
    /// atol_scale: per-state multipliers of atol, as in Session.advance.
    #[pyo3(signature = (keys, values, times, threads=None, rtol=1e-6, atol=1e-12, outputs=None, atol_scale=None, seed=0))]
    #[allow(clippy::too_many_arguments)]
    fn simulate_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                           times: PyReadonlyArray1<f64>, threads: Option<usize>, rtol: f64, atol: f64,
                           outputs: Option<Vec<String>>, atol_scale: Option<Vec<f64>>, seed: u64) -> PyResult<Bound<'py, PyAny>> {
        self.time_courses(py, keys, values, times, threads, outputs, Some((rtol, atol)), seed, atol_scale)
    }

    /// One time course -> (trajectory[len(times), len(outputs)], final) where final maps every
    /// output id (states, rules, rates, parameters incl. event-changed ones, time) to its value at
    /// times[-1]. Raises RuntimeError if integration fails.
    /// event_state: None (events start fresh, returns a 2-tuple) or the state returned by a previous
    /// call (an empty tuple () starts fresh) to continue as consecutive roadrunner simulate calls do;
    /// then returns (trajectory, final, event_state) with event_state = (trigger values,
    /// [(time, event index, values fixed at trigger time)]) for the next call.
    /// final_array=True gives final as an array in output_ids order instead of a dict.
    /// atol_scale: per-state (state_ids) multipliers of atol, as roadrunner's CVODE scales a scalar
    /// atol by the initial values.
    #[pyo3(signature = (keys, values, times, outputs=None, rtol=1e-6, atol=1e-12, event_state=None, final_array=false, atol_scale=None))]
    #[allow(clippy::too_many_arguments)]
    fn simulate_one<'py>(&self, py: Python<'py>, keys: Bound<'py, PyAny>, values: Vec<f64>, times: PyReadonlyArray1<f64>,
                         outputs: Option<Bound<'py, PyAny>>, rtol: f64, atol: f64, event_state: Option<Bound<'py, PyAny>>,
                         final_array: bool, atol_scale: Option<PyReadonlyArray1<f64>>) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.m;
        let keys = keys_of(m, &keys)?;
        if values.len() != keys.len() { return Err(verr("values and keys differ in length".into())); }
        let times = check_times(times)?;
        let outs = match outputs {
            Some(o) => match o.cast::<Resolved>() { Ok(r) => r.get().outs.clone(), Err(_) => self.outs(Some(o.extract()?))? },
            None => self.outs(None)?,
        };
        let all: Vec<u32> = m.outputs.iter().map(|o| o.1).collect();
        let (no, nt) = (outs.len(), times.len());
        let mut out = vec![0.0; nt * no];
        let mut fin = vec![0.0; all.len()];
        let mut w = self.pool.lock().unwrap().pop().unwrap_or_else(|| m.work()); // reuse workspaces across calls
        match &atol_scale {
            Some(a) if a.len()? != m.nx => return Err(verr("atol_scale must have one entry per state".into())),
            Some(a) => w.atol_scale.copy_from_slice(a.as_slice()?),
            None => w.atol_scale.fill(1.0),
        }
        let cont = match &event_state {
            None => false,
            Some(e) if e.len()? == 0 => { m.events_reset(&mut w); true }
            Some(e) => {
                let (prev, pend): (Vec<bool>, Vec<(f64, usize, Vec<f64>)>) = e.extract()?;
                m.set_event_state(&mut w, prev, pend).map_err(verr)?;
                true
            }
        };
        let ok = py.detach(|| {
            m.load(&mut w, &keys, &values);
            let ok = if cont { m.simulate_cont(&mut w, &times, rtol, atol, &outs, &mut out) }
                     else { m.simulate(&mut w, &times, rtol, atol, &outs, &mut out) };
            if ok { m.record(&mut w, times[nt - 1], &all, &mut fin); }
            ok
        });
        let es = m.event_state(&w);
        self.pool.lock().unwrap().push(w);
        if !ok { return Err(PyRuntimeError::new_err("roadrunner_nitro: integration failed")); }
        let d = if final_array { PyArray1::from_vec(py, fin).into_any() } else {
            let names = self.names.get_or_init(|| m.outputs.iter().map(|o| PyString::intern(py, &o.0).unbind()).collect());
            let d = PyDict::new(py);
            for (o, x) in names.iter().zip(fin) { d.set_item(o.bind(py), x)?; }
            d.into_any()
        };
        let y = arr(py, out, &[nt, no])?;
        Ok(if cont { (y, d, es).into_pyobject(py)?.into_any() } else { (y, d).into_pyobject(py)?.into_any() })
    }

    /// Outputs (all, in output_ids order, or `outputs`) at the keyed state and `time`, without integrating.
    #[pyo3(signature = (keys, values, time=0.0, outputs=None))]
    fn outputs_at<'py>(&self, py: Python<'py>, keys: Bound<'py, PyAny>, values: Vec<f64>, time: f64,
                       outputs: Option<Bound<'py, PyAny>>) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.m;
        let keys = keys_of(m, &keys)?;
        if values.len() != keys.len() { return Err(verr("values and keys differ in length".into())); }
        let outs = match outputs {
            Some(o) => match o.cast::<Resolved>() { Ok(r) => r.get().outs.clone(), Err(_) => self.outs(Some(o.extract()?))? },
            None => m.outputs.iter().map(|o| o.1).collect(),
        };
        let mut w = self.pool.lock().unwrap().pop().unwrap_or_else(|| m.work());
        let mut row = vec![0.0; outs.len()];
        m.load(&mut w, &keys, &values);
        m.record(&mut w, time, &outs, &mut row);
        self.pool.lock().unwrap().push(w);
        Ok(PyArray1::from_vec(py, row).into_any())
    }

    /// (n, len(times), len(outputs)) Gillespie direct-method trajectories (outputs as concentrations).
    /// Sample i uses an xoshiro256++ stream seeded from seed + i (seed=None: from the clock).
    #[pyo3(signature = (keys, values, times, seed=None, threads=None, outputs=None))]
    #[allow(clippy::too_many_arguments)]
    fn gillespie_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                            times: PyReadonlyArray1<f64>, seed: Option<u64>, threads: Option<usize>,
                            outputs: Option<Vec<String>>) -> PyResult<Bound<'py, PyAny>> {
        if self.m.has_algebraic() { return Err(verr("gillespie_batch: algebraic rules are not supported".into())); }
        let seed = seed.unwrap_or_else(|| std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64));
        self.time_courses(py, keys, values, times, threads, outputs, None, seed, None)
    }

    /// Full Jacobian d(dx/dt)/dx (nState x nState, concentrations, state_ids order) at the state
    /// given by the model values overridden by keys/values, at `time`. Exact (symbolic) derivatives.
    #[pyo3(signature = (keys, values, time=0.0))]
    fn jacobian<'py>(&self, py: Python<'py>, keys: Vec<String>, values: Vec<f64>, time: f64) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.m;
        let keys = m.resolve(&keys).map_err(verr)?;
        if values.len() != keys.len() { return Err(verr("values and keys differ in length".into())); }
        let mut w = m.work();
        m.load(&mut w, &keys, &values);
        let mut j = vec![0.0; m.nx * m.nx];
        m.jacobian(&mut w, time, &mut j);
        arr(py, j, &[m.nx, m.nx])
    }
}

/// A running simulation of one Model: keeps the event state (trigger values, pending delayed
/// events) between advance() calls, as consecutive roadrunner simulate/oneStep calls do.
///
/// Session(model); advance(keys, values, times, outputs=None, rtol=1e-6, atol=1e-12, atol_scale=None)
///   -> (trajectory[len(times), len(outputs)], final values of all outputs in output_ids order).
/// A failed run is retried once with atol <= 1e-12 (a loose atol can let tiny species go astray),
/// leaving the event state as it was before the call; RuntimeError if that fails too.
/// stats: solver counters of the last advance (solver, bdf_steps, rodas_steps, rejected, switched, atol).
#[pyclass(name = "Session")]
pub(crate) struct Session { model: Py<PyModel>, w: Work, fresh: bool, stats: ([usize; 4], f64, bool) }

#[pymethods]
impl Session {
    #[new]
    fn new(model: Py<PyModel>) -> Self {
        let w = model.get().m.work();
        Session { model, w, fresh: true, stats: ([0; 4], 0.0, false) }
    }

    /// Start the next advance with fresh event triggers (roadrunner's reset()).
    fn reset_events(&mut self) { self.fresh = true; }

    #[pyo3(signature = (keys, values, times, outputs=None, rtol=1e-6, atol=1e-12, atol_scale=None))]
    #[allow(clippy::too_many_arguments)]
    fn advance<'py>(&mut self, py: Python<'py>, keys: Bound<'py, PyAny>, values: Vec<f64>, times: PyReadonlyArray1<f64>,
                    outputs: Option<Bound<'py, PyAny>>, rtol: f64, atol: f64,
                    atol_scale: Option<PyReadonlyArray1<f64>>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
        let pm = self.model.get();
        let m = &pm.m;
        let keys = keys_of(m, &keys)?;
        if values.len() != keys.len() { return Err(verr("values and keys differ in length".into())); }
        let times = check_times(times)?;
        let outs = match outputs {
            Some(o) => match o.cast::<Resolved>() { Ok(r) => r.get().outs.clone(), Err(_) => pm.outs(Some(o.extract()?))? },
            None => pm.outs(None)?,
        };
        let all: Vec<u32> = m.outputs.iter().map(|o| o.1).collect();
        let (no, nt) = (outs.len(), times.len());
        let mut out = vec![0.0; nt * no];
        let mut fin = vec![0.0; all.len()];
        let w = &mut self.w;
        match &atol_scale {
            Some(a) if a.len()? != m.nx => return Err(verr("atol_scale must have one entry per state".into())),
            Some(a) => w.atol_scale.copy_from_slice(a.as_slice()?),
            None => w.atol_scale.fill(1.0),
        }
        if self.fresh { m.events_reset(w); }
        let saved = m.event_state(w);
        let fresh = &mut self.fresh;
        let stats = &mut self.stats;
        let ok = py.detach(|| {
            let mut tried = f64::NAN;
            for a in [atol, atol.min(1e-12)] {
                if a == tried { break; }
                tried = a;
                w.stats = [0; 4];
                let hist = w.hist.save(); // a continued run keeps its delay() history
                m.load(w, &keys, &values);
                if !*fresh { w.hist.restore(hist); }
                if m.simulate_cont(w, &times, rtol, a, &outs, &mut out) {
                    m.record(w, times[nt - 1], &all, &mut fin);
                    *fresh = false;
                    *stats = (w.stats, a, true);
                    return true;
                }
                let _ = m.set_event_state(w, saved.0.clone(), saved.1.clone()); // undo the failed attempt
            }
            *stats = (w.stats, tried, false);
            false
        });
        if !ok { return Err(PyRuntimeError::new_err("roadrunner_nitro: integration failed")); }
        Ok((arr(py, out, &[nt, no])?, PyArray1::from_vec(py, fin).into_any()))
    }

    /// Outputs at `time` from the session's own state (as the last advance left it, with its delay()
    /// history): (outputs) -> array.
    #[pyo3(signature = (time, outputs=None))]
    fn outputs_now<'py>(&mut self, py: Python<'py>, time: f64, outputs: Option<Bound<'py, PyAny>>) -> PyResult<Bound<'py, PyAny>> {
        let pm = self.model.get();
        let outs = match outputs {
            Some(o) => match o.cast::<Resolved>() { Ok(r) => r.get().outs.clone(), Err(_) => pm.outs(Some(o.extract()?))? },
            None => pm.m.outputs.iter().map(|o| o.1).collect(),
        };
        let mut row = vec![0.0; outs.len()];
        pm.m.record(&mut self.w, time, &outs, &mut row);
        Ok(PyArray1::from_vec(py, row).into_any())
    }

    /// Counters of the last advance: dict(solver='bdf'|'rodas4'|'bdf+rodas4', bdf_steps, rodas_steps,
    /// rejected, switched, atol (the one that succeeded), ok).
    #[getter]
    fn stats<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let ([b, r, rej, sw], a, ok) = self.stats;
        let d = PyDict::new(py);
        d.set_item("solver", if b > 0 && r > 0 { "bdf+rodas4" } else if b > 0 { "bdf" } else { "rodas4" })?;
        d.set_item("bdf_steps", b)?;
        d.set_item("rodas_steps", r)?;
        d.set_item("rejected", rej)?;
        d.set_item("switched", sw > 0)?;
        d.set_item("atol", a)?;
        d.set_item("ok", ok)?;
        Ok(d)
    }
}

/// Antimony source -> SBML Level 3 Version 2 text. ValueError lists the diagnostics on any error.
#[pyfunction]
fn antimony_to_sbml(source: &str) -> PyResult<String> {
    let loaded = antimony::parse(source);
    if loaded.has_errors() {
        return Err(pyo3::exceptions::PyValueError::new_err(loaded.diagnostics.to_string()));
    }
    Ok(antimony::sbml::write_sbml(&loaded.document))
}

/// Flatten a hierarchical (comp) SBML document into SBML core text. ``base``: directory external
/// model definitions are relative to. ValueError if it cannot be flattened.
#[pyfunction]
#[pyo3(signature = (sbml, base=None))]
fn flatten_comp(sbml: &str, base: Option<String>) -> PyResult<String> {
    sbml_io::comp::flatten(sbml, base.as_deref().map(std::path::Path::new)).map_err(verr)
}

/// SBML text -> (Model, info dict), without roadrunner or libsbml (see crate::sbml::from_sbml).
/// NotImplementedError: SBML the engine does not cover (comp: flatten first); ValueError: invalid SBML.
#[pyfunction]
fn from_sbml<'py>(py: Python<'py>, sbml: &str) -> PyResult<(PyModel, Bound<'py, PyDict>)> {
    use crate::sbml::Error;
    let (m, i) = crate::sbml::from_sbml(sbml).map_err(|e| match e {
        Error::NotImplemented(s) => pyo3::exceptions::PyNotImplementedError::new_err(s),
        Error::Invalid(s) => PyValueError::new_err(s),
    })?;
    let d = PyDict::new(py);
    let [f, b, g, c, r] = i.ids;
    d.set_item("ids", (f, b, g, c, r))?;
    let nr = i.n.len().checked_div(i.rows.len()).unwrap_or(0);
    let n = numpy::ndarray::Array2::from_shape_vec((i.rows.len(), nr), i.n).map_err(|e| verr(e.to_string()))?;
    d.set_item("N", numpy::PyArray2::from_owned_array(py, n))?;
    d.set_item("rows", i.rows)?;
    d.set_item("comp_of", i.comp_of)?;
    d.set_item("general", i.general)?;
    d.set_item("changing", i.changing)?;
    d.set_item("amount_params", i.amount_params)?;
    d.set_item("rate_rule_ids", i.rate_rule_ids)?;
    d.set_item("ia", i.ia.into_iter().collect::<std::collections::HashMap<_, _>>())?;
    d.set_item("fdefs", i.fdefs.into_iter().collect::<std::collections::HashMap<_, _>>())?;
    d.set_item("values", i.values)?;
    d.set_item("n_events", i.n_events)?;
    d.set_item("rule_vars", i.rule_vars)?;
    d.set_item("hosu", i.hosu)?;
    d.set_item("conversion_factors", i.conversion_factors)?;
    d.set_item("srefs", i.srefs)?;
    Ok((PyModel { m, pool: Default::default(), names: Default::default() }, d))
}

#[pymodule]
fn roadrunner_nitro(m: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::sens::register(m)?;
    m.add_function(wrap_pyfunction!(antimony_to_sbml, m)?)?;
    m.add_function(wrap_pyfunction!(flatten_comp, m)?)?;
    m.add_function(wrap_pyfunction!(from_sbml, m)?)?;
    m.add_class::<Session>()?;
    m.add_class::<PyRR>()?;
    m.add_class::<PyModel>()
}

// ---------------------------------------------------------------- RoadRunner (the Rust core)

use crate::rr::{Matrix, RoadRunner, Value};

type Mat<'py> = (Bound<'py, PyAny>, Vec<String>, Vec<String>);

/// libroadrunner-like instance on the Rust core (crate::rr): RoadRunner(sbml=None, base=None).
/// Matrices come back as (array, row ids, column ids).
#[pyclass(unsendable, name = "RoadRunner")]
struct PyRR { r: RoadRunner }

fn m3<'py>(py: Python<'py>, m: Result<Matrix, String>) -> PyResult<Mat<'py>> {
    let m = m.map_err(verr)?;
    Ok((arr(py, m.data, &[m.nr, m.nc])?, m.rows, m.cols))
}

fn to_value(v: &Bound<'_, PyAny>) -> PyResult<Value> {
    Ok(if let Ok(b) = v.cast::<pyo3::types::PyBool>() { Value::Bool(b.is_true()) }
       else if let Ok(i) = v.extract::<i64>() { Value::Int(i) }
       else if let Ok(u) = v.extract::<u64>() { Value::UInt(u) }
       else if let Ok(f) = v.extract::<f64>() { Value::Double(f) }
       else if let Ok(s) = v.extract::<String>() { Value::Str(s) }
       else { Value::Doubles(v.extract::<Vec<f64>>()?) })
}

fn from_value<'py>(py: Python<'py>, v: &Value) -> PyResult<Bound<'py, PyAny>> {
    Ok(match v {
        Value::Bool(b) => pyo3::types::PyBool::new(py, *b).to_owned().into_any(),
        Value::Int(i) => i.into_pyobject(py)?.into_any(),
        Value::UInt(u) => u.into_pyobject(py)?.into_any(),
        Value::Double(d) => d.into_pyobject(py)?.into_any(),
        Value::Str(s) => s.into_pyobject(py)?.into_any(),
        Value::Doubles(d) => d.clone().into_pyobject(py)?.into_any(),
    })
}

fn cplx(py: Python<'_>, e: Result<Vec<(f64, f64)>, String>) -> PyResult<Bound<'_, PyAny>> {
    let e = e.map_err(verr)?;
    Ok(PyArray1::from_vec(py, e.into_iter().map(|(re, im)| numpy::Complex64::new(re, im)).collect()).into_any())
}

#[allow(non_snake_case)]
#[pymethods]
impl PyRR {
    #[new]
    #[pyo3(signature = (sbml=None, base=None))]
    fn new(sbml: Option<&str>, base: Option<String>) -> PyResult<Self> {
        let mut r = RoadRunner::new();
        if let Some(s) = sbml { r.load_sbml(s, base.as_deref().map(std::path::Path::new)).map_err(verr)?; }
        Ok(PyRR { r })
    }
    #[pyo3(signature = (sbml, base=None))]
    fn load(&mut self, sbml: &str, base: Option<String>) -> PyResult<()> { self.r.load_sbml(sbml, base.as_deref().map(std::path::Path::new)).map_err(verr) }
    fn isModelLoaded(&self) -> bool { self.r.is_model_loaded() }
    fn clearModel(&mut self) { self.r.clear_model() }
    fn getSBML(&self) -> PyResult<String> { self.r.get_sbml().map_err(verr) }
    fn getCurrentSBML(&mut self) -> PyResult<String> { self.r.get_current_sbml().map_err(verr) }
    fn regenerate(&mut self) -> PyResult<()> { self.r.regenerate().map_err(verr) }
    fn saveState(&self, path: &str) -> PyResult<()> { self.r.save_state(std::path::Path::new(path)).map_err(verr) }
    fn loadState(&mut self, path: &str) -> PyResult<()> { self.r.load_state(std::path::Path::new(path)).map_err(verr) }

    /// simulate(start=None, end=None, points=None, selections=None) -> array (points x selections)
    #[pyo3(signature = (start=None, end=None, points=None, selections=None))]
    fn simulate<'py>(&mut self, py: Python<'py>, start: Option<f64>, end: Option<f64>, points: Option<usize>,
                     selections: Option<Vec<String>>) -> PyResult<Bound<'py, PyAny>> {
        if let Some(s) = start { self.r.start = s; }
        if let Some(e) = end { self.r.end = e; }
        if let Some(p) = points { self.r.points = p; }
        if let Some(s) = selections { self.r.set_selections(&s).map_err(verr)?; }
        let m = self.r.simulate().map_err(verr)?;
        arr(py, m.data.clone(), &[m.nr, m.nc])
    }
    fn oneStep(&mut self, t: f64, h: f64) -> PyResult<f64> { self.r.one_step(t, h).map_err(verr) }
    fn getCurrentTime(&self) -> f64 { self.r.time() }

    #[getter]
    fn get_timeCourseSelections(&self) -> Vec<String> { self.r.selections().to_vec() }
    #[setter]
    fn set_timeCourseSelections(&mut self, s: Vec<String>) -> PyResult<()> { self.r.set_selections(&s).map_err(verr) }
    #[getter]
    fn get_steadyStateSelections(&self) -> Vec<String> { self.r.steady_state_selections().to_vec() }
    #[setter]
    fn set_steadyStateSelections(&mut self, s: Vec<String>) -> PyResult<()> { self.r.set_steady_state_selections(&s).map_err(verr) }
    #[getter]
    fn get_conservedMoietyAnalysis(&self) -> bool { self.r.conservation_laws() }
    #[setter]
    fn set_conservedMoietyAnalysis(&mut self, on: bool) -> PyResult<()> { self.r.set_conservation_laws(on).map_err(verr) }
    #[getter]
    fn get_diffStepSize(&self) -> f64 { self.r.diff_step }
    #[setter]
    fn set_diffStepSize(&mut self, v: f64) { self.r.diff_step = v; }

    fn setIntegrator(&mut self, name: &str) -> PyResult<()> { self.r.set_integrator(name).map_err(verr) }
    fn getIntegratorName(&self) -> &'static str { self.r.integrator().name }
    fn setIntegratorSetting(&mut self, key: &str, v: &Bound<'_, PyAny>) -> PyResult<()> { self.r.integrator_mut().set(key, to_value(v)?).map_err(verr) }
    fn getIntegratorSetting<'py>(&self, py: Python<'py>, key: &str) -> PyResult<Bound<'py, PyAny>> { from_value(py, self.r.integrator().get(key).map_err(verr)?) }
    fn getIntegratorSettings(&self) -> Vec<&'static str> { self.r.integrator().params.iter().map(|p| p.name).collect() }
    fn setSteadyStateSolver(&mut self, name: &str) -> PyResult<()> { self.r.set_steady_state_solver(name).map_err(verr) }
    fn getSteadyStateSolverName(&self) -> &'static str { self.r.steady_state_solver().name }
    fn setSteadyStateSolverSetting(&mut self, key: &str, v: &Bound<'_, PyAny>) -> PyResult<()> { self.r.steady_state_solver_mut().set(key, to_value(v)?).map_err(verr) }
    fn getSteadyStateSolverSetting<'py>(&self, py: Python<'py>, key: &str) -> PyResult<Bound<'py, PyAny>> { from_value(py, self.r.steady_state_solver().get(key).map_err(verr)?) }
    fn getSeed(&mut self) -> u64 { self.r.seed() }
    fn setSeed(&mut self, s: u64) -> PyResult<()> { self.r.set_seed(s).map_err(verr) }

    fn getValue(&mut self, sel: &str) -> PyResult<f64> { self.r.get_value(sel).map_err(verr) }
    fn setValue(&mut self, sel: &str, v: f64) -> PyResult<()> { self.r.set_value(sel, v).map_err(verr) }
    fn __getitem__(&mut self, sel: &str) -> PyResult<f64> { self.getValue(sel) }
    fn __setitem__(&mut self, sel: &str, v: f64) -> PyResult<()> { self.setValue(sel, v) }
    fn getValues(&mut self, sels: Vec<String>) -> PyResult<Vec<f64>> { self.r.values(&sels).map_err(verr) }

    fn getFloatingSpeciesIds(&self) -> Vec<String> { self.r.floating_ids() }
    fn getBoundarySpeciesIds(&self) -> Vec<String> { self.r.boundary_ids() }
    fn getGlobalParameterIds(&self) -> Vec<String> { self.r.global_parameter_ids() }
    fn getCompartmentIds(&self) -> Vec<String> { self.r.compartment_ids() }
    fn getReactionIds(&self) -> Vec<String> { self.r.reaction_ids() }
    fn getRateRuleIds(&self) -> Vec<String> { self.r.rate_rule_ids() }
    fn getIndependentFloatingSpeciesIds(&self) -> Vec<String> { self.r.independent_ids() }
    fn getDependentFloatingSpeciesIds(&self) -> Vec<String> { self.r.dependent_ids() }
    fn getConservedMoietyIds(&self) -> Vec<String> { self.r.conserved_moiety_ids() }
    fn getEigenValueIds(&self) -> Vec<String> { self.r.eigenvalue_ids() }
    fn getRatesOfChangeIds(&self) -> Vec<String> { self.r.rates_of_change_ids() }
    fn getAvailableTimeCourseSymbols(&self) -> Vec<(&'static str, Vec<String>)> { self.r.available_time_course_symbols() }
    fn getFloatingSpeciesConcentrations(&mut self) -> PyResult<Vec<f64>> { self.r.floating_concentrations().map_err(verr) }
    fn getFloatingSpeciesAmounts(&mut self) -> PyResult<Vec<f64>> { self.r.floating_amounts().map_err(verr) }
    fn getFloatingSpeciesInitConcentrations(&mut self) -> PyResult<Vec<f64>> { self.r.floating_init_concentrations().map_err(verr) }
    fn getBoundarySpeciesConcentrations(&mut self) -> PyResult<Vec<f64>> { self.r.boundary_concentrations().map_err(verr) }
    fn getGlobalParameterValues(&mut self) -> PyResult<Vec<f64>> { self.r.global_parameter_values().map_err(verr) }
    fn getCompartmentVolumes(&mut self) -> PyResult<Vec<f64>> { self.r.compartment_volumes().map_err(verr) }
    fn getReactionRates(&mut self) -> PyResult<Vec<f64>> { self.r.reaction_rates().map_err(verr) }
    fn getRatesOfChange(&mut self) -> PyResult<Vec<f64>> { self.r.rates_of_change().map_err(verr) }

    fn reset(&mut self) -> PyResult<()> { self.r.reset().map_err(verr) }
    fn resetAll(&mut self) -> PyResult<()> { self.r.reset_all().map_err(verr) }
    fn resetToOrigin(&mut self) -> PyResult<()> { self.r.reset_to_origin().map_err(verr) }
    fn resetParameter(&mut self) -> PyResult<()> { self.r.reset_parameter().map_err(verr) }

    fn steadyState(&mut self) -> PyResult<f64> { self.r.steady_state().map_err(verr) }
    fn getSteadyStateValues(&mut self) -> PyResult<Vec<f64>> { self.r.steady_state_values().map_err(verr) }

    fn getFullStoichiometryMatrix<'py>(&self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.stoichiometry_matrix()) }
    fn getNrMatrix<'py>(&self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.reduced_stoichiometry_matrix()) }
    fn getLinkMatrix<'py>(&self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.link_matrix()) }
    fn getL0Matrix<'py>(&self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.l0_matrix()) }
    fn getConservationMatrix<'py>(&self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.conservation_matrix()) }
    fn getFullJacobian<'py>(&mut self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.full_jacobian()) }
    fn getReducedJacobian<'py>(&mut self, py: Python<'py>) -> PyResult<Mat<'py>> { m3(py, self.r.reduced_jacobian()) }
    fn getFullEigenValues<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> { cplx(py, self.r.full_eigenvalues()) }
    fn getReducedEigenValues<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> { cplx(py, self.r.reduced_eigenvalues()) }
    #[pyo3(signature = (scaled=false))]
    fn getElasticityMatrix<'py>(&mut self, py: Python<'py>, scaled: bool) -> PyResult<Mat<'py>> { m3(py, self.r.elasticity_matrix(scaled)) }
    #[pyo3(signature = (scaled=false))]
    fn getConcentrationControlCoefficientMatrix<'py>(&mut self, py: Python<'py>, scaled: bool) -> PyResult<Mat<'py>> { m3(py, self.r.cc_matrix(scaled)) }
    #[pyo3(signature = (scaled=false))]
    fn getFluxControlCoefficientMatrix<'py>(&mut self, py: Python<'py>, scaled: bool) -> PyResult<Mat<'py>> { m3(py, self.r.fcc_matrix(scaled)) }
    fn getCC(&mut self, var: &str, par: &str) -> PyResult<f64> { self.r.get_cc(var, par).map_err(verr) }
    fn getuCC(&mut self, var: &str, par: &str) -> PyResult<f64> { self.r.get_ucc(var, par).map_err(verr) }
    #[pyo3(signature = (reaction, x, steady=true))]
    fn getEE(&mut self, reaction: &str, x: &str, steady: bool) -> PyResult<f64> { self.r.get_ee(reaction, x, steady).map_err(verr) }
    #[pyo3(signature = (reaction, x, steady=true))]
    fn getuEE(&mut self, reaction: &str, x: &str, steady: bool) -> PyResult<f64> { self.r.get_uee(reaction, x, steady).map_err(verr) }

    /// gillespie(start, end, points=None) -> array: on a grid of `points`, or a row per firing.
    #[pyo3(signature = (start, end, points=None))]
    fn gillespie<'py>(&mut self, py: Python<'py>, start: f64, end: f64, points: Option<usize>) -> PyResult<Bound<'py, PyAny>> {
        let m = self.r.gillespie(start, end, points).map_err(verr)?;
        arr(py, m.data.clone(), &[m.nr, m.nc])
    }
    /// gillespieRuns(n, start, end, points, seed=None) -> array (n, points, selections), all cores.
    #[pyo3(signature = (n, start, end, points, seed=None))]
    fn gillespieRuns<'py>(&mut self, py: Python<'py>, n: usize, start: f64, end: f64, points: usize, seed: Option<u64>) -> PyResult<Bound<'py, PyAny>> {
        let nc = self.r.selections().len();
        let y = self.r.gillespie_runs(n, start, end, points, seed).map_err(verr)?;
        arr(py, y, &[n, points, nc])
    }
    fn gillespieMeanSD<'py>(&mut self, py: Python<'py>, n: usize, start: f64, end: f64, points: usize) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
        let (a, b) = self.r.gillespie_mean_sd(n, start, end, points).map_err(verr)?;
        Ok((arr(py, a.data, &[a.nr, a.nc])?, arr(py, b.data, &[b.nr, b.nc])?))
    }
}

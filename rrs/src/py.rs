//! Python bindings (module `roadrunner_rs`).
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

type PyEvent = (String, Vec<(String, String)>, Option<String>, Option<String>, bool, bool, bool);

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
                         ode: Option<(f64, f64)>, seed: u64) -> PyResult<Bound<'py, PyAny>> {
        let b = self.batch(keys, values)?;
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
                let ok = match ode {
                    Some((rtol, atol)) => m.simulate(w, &times, rtol, atol, &outs, r),
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
                        rate_rules=None, assignment_rules=None, events=None))]
    fn new(species: Vec<String>, init: Vec<f64>, volumes: Vec<f64>, reactions: Vec<String>, rates: Vec<String>,
           stoich: PyReadonlyArray2<f64>, param_names: Vec<String>, param_values: Vec<f64>,
           rate_rules: Option<Vec<(String, String)>>, assignment_rules: Option<Vec<(String, String)>>,
           events: Option<Vec<PyEvent>>) -> PyResult<Self> {
        let a = stoich.as_array();
        if a.dim() != (species.len(), reactions.len()) {
            return Err(verr(format!("stoich shape {:?} != (nSpecies, nReactions)", a.dim())));
        }
        let st = a.iter().copied().collect();
        let x = Rules {
            rate: rate_rules.unwrap_or_default(),
            assign: assignment_rules.unwrap_or_default(),
            events: events.unwrap_or_default().into_iter().map(|(trigger, assign, delay, priority, initial_value, persistent, u)|
                EventDef { trigger, assign, delay, priority, initial_value, persistent, use_trigger_values: u }).collect(),
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

    /// (n, len(times), len(outputs)) time courses (Rodas4, events root-located); outputs default to
    /// the floating species (concentrations). NaN rows where integration failed.
    #[pyo3(signature = (keys, values, times, threads=None, rtol=1e-6, atol=1e-12, outputs=None))]
    #[allow(clippy::too_many_arguments)]
    fn simulate_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                           times: PyReadonlyArray1<f64>, threads: Option<usize>, rtol: f64, atol: f64,
                           outputs: Option<Vec<String>>) -> PyResult<Bound<'py, PyAny>> {
        self.time_courses(py, keys, values, times, threads, outputs, Some((rtol, atol)), 0)
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
        if !ok { return Err(PyRuntimeError::new_err("roadrunner_rs: integration failed")); }
        let d = if final_array { PyArray1::from_vec(py, fin).into_any() } else {
            let names = self.names.get_or_init(|| m.outputs.iter().map(|o| PyString::intern(py, &o.0).unbind()).collect());
            let d = PyDict::new(py);
            for (o, x) in names.iter().zip(fin) { d.set_item(o.bind(py), x)?; }
            d.into_any()
        };
        let y = arr(py, out, &[nt, no])?;
        Ok(if cont { (y, d, es).into_pyobject(py)?.into_any() } else { (y, d).into_pyobject(py)?.into_any() })
    }

    /// (n, len(times), len(outputs)) Gillespie direct-method trajectories (outputs as concentrations).
    /// Sample i uses an xoshiro256++ stream seeded from seed + i (seed=None: from the clock).
    #[pyo3(signature = (keys, values, times, seed=None, threads=None, outputs=None))]
    #[allow(clippy::too_many_arguments)]
    fn gillespie_batch<'py>(&self, py: Python<'py>, keys: Vec<String>, values: PyReadonlyArray2<f64>,
                            times: PyReadonlyArray1<f64>, seed: Option<u64>, threads: Option<usize>,
                            outputs: Option<Vec<String>>) -> PyResult<Bound<'py, PyAny>> {
        if self.m.nv != self.m.nr { return Err(verr("gillespie_batch: rate rules are not supported".into())); }
        let seed = seed.unwrap_or_else(|| std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64));
        self.time_courses(py, keys, values, times, threads, outputs, None, seed)
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

#[pymodule]
fn roadrunner_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::sens::register(m)?;
    m.add_class::<PyModel>()
}

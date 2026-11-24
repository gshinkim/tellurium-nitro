//! More of libroadrunner's members: K and extended / reactant / product stoichiometry matrices,
//! stoichiometry ids and values, frequency response, internalOneStep, getIds by SelectionRecord type,
//! SBML level/version conversion, event roots and piecewise triggers, the model's random numbers.
use super::*;
use crate::sbml::math::{self, E, Fdefs};
use sbml_io::math::Math;

/// libroadrunner's SelectionRecord flags of the getIds groups (see `ids_of_type`).
const IND: i64 = 0x80000;
const DEP: i64 = 0x40000;

/// A species reference of the SBML: (species, reaction index, reactant, position in its list, id, stoichiometry).
struct SRef {
    species: String,
    rx: usize,
    reactant: bool,
    pos: usize,
    id: Option<String>,
    value: f64,
    math: Option<Math>,
}

/// Runtime coefficient changes have their own compiled document. Structural analysis and
/// getSBML continue to describe the loaded document, and resetAll restores that document.
pub(super) struct StoichEdit {
    pub(super) source: String,
    pub(super) runtime: String,
    initial: HashMap<(usize, bool, usize), f64>,
}

impl Core {
    pub(super) fn carry_stoich_edit(&mut self, old: &mut Core) {
        self.stoich_edit = old.stoich_edit.take();
        std::mem::swap(&mut self.ls_cell, &mut old.ls_cell);
        std::mem::swap(&mut self.ls_rows, &mut old.ls_rows);
        std::mem::swap(&mut self.ls_n, &mut old.ls_n);
        std::mem::swap(&mut self.ls_nr, &mut old.ls_nr);
    }

    fn integrate_one(&mut self, t: f64, reset: bool, name: &'static str,
                     step: impl FnOnce(&Model, &mut crate::model::Work, &mut f64) -> Result<f64, String>)
                     -> Result<f64, String> {
        let continuing = !reset && self.live && self.t == t;
        let mut next = self.internal_step.filter(|&(n, time, _)| continuing && n == name && time == t)
            .map_or(0.0, |(_, _, h)| h);
        if !continuing {
            self.m.load(&mut self.w, &self.kres, &self.x);
            self.m.events_reset(&mut self.w);
        }
        let reached = step(&self.m, &mut self.w, &mut next)?;
        self.m.record(&mut self.w, reached, &self.xo, &mut self.x);
        (self.t, self.fresh, self.qok, self.live) = (reached, false, false, true);
        self.internal_step = Some((name, reached, next));
        Ok(reached)
    }
}

impl RoadRunner {
    /// The runtime SBML document with comp models flattened.
    fn doc(&self) -> Result<sbml_io::Document, String> {
        let text = self.runtime_sbml()?;
        if text.contains("/comp/") {
            return read_doc(&sbml_io::comp::flatten(text, self.base.as_deref())?);
        }
        read_doc(text)
    }

    pub(super) fn runtime_sbml(&self) -> Result<&str, String> {
        Ok(self.core()?.stoich_edit.as_ref().filter(|e| e.source == self.sbml)
            .map_or(self.sbml.as_str(), |e| e.runtime.as_str()))
    }

    fn srefs(&self) -> Result<Vec<SRef>, String> {
        let doc = self.doc()?;
        let Some(m) = doc.model else { return Ok(vec![]) };
        Ok(m.reactions.iter().enumerate().flat_map(|(rx, r)| {
            let side = |l: &Vec<sbml_io::SpeciesReference>, reactant: bool| l.iter().enumerate().map(move |(pos, s)|
                SRef { species: s.species.clone(), rx, reactant, pos, id: s.id.clone(), value: s.stoichiometry.unwrap_or(1.0), math: s.stoichiometry_math.clone() }).collect::<Vec<_>>();
            side(&r.reactants, true).into_iter().chain(side(&r.products, false))
        }).collect())
    }

    // ------------------------------------------------------------ stoichiometry

    /// Species by reaction matrix with signed reactants and positive products.
    fn sref_matrix(&mut self, part: Option<bool>, boundary: bool) -> Result<Matrix, String> {
        let mut rows = self.floating_ids();
        if boundary { rows.extend(self.boundary_ids()); }
        let cols = self.reaction_ids();
        let mut data = vec![0.0; rows.len() * cols.len()];
        for s in self.srefs()? {
            if part.is_some_and(|p| p != s.reactant) { continue; }
            if let Some(i) = rows.iter().position(|r| *r == s.species) {
                let v = self.sref_value(&s, false)?;
                data[i * cols.len() + s.rx] += if s.reactant { -v } else { v };
            }
        }
        Ok(Matrix::new(rows, cols, data))
    }

    pub fn extended_stoichiometry_matrix(&mut self) -> Result<Matrix, String> {
        self.sref_matrix(None, true)
    }

    pub fn reactants_products_matrix(&mut self, reactants: bool, boundary: bool) -> Result<Matrix, String> {
        self.sref_matrix(Some(reactants), boundary)
    }

    /// A right nullspace basis in structural reaction order, with the identity block first.
    pub fn k_matrix(&self) -> Result<Matrix, String> {
        let c = self.core()?;
        let nr = c.info.ids[4].len();
        let ns = c.ls_rows.len();
        let n: Vec<_> = c.ls().order.iter().flat_map(|&i| c.ls_n[i * nr..(i + 1) * nr].iter().copied()).collect();
        let (rank, order, coeffs) = if n.iter().all(|&x| x == 0.0) {
            (0, (0..nr).collect::<Vec<_>>(), vec![])
        } else if ns == 1 && n[0].abs() > 1e-9 {
            (1, (0..nr).collect(), n[1..].iter().map(|v| v / n[0]).collect())
        } else {
            let nt: Vec<_> = (0..nr).flat_map(|r| (0..ns).map(move |s| (s, r)))
                .map(|(s, r)| n[s * nr + r]).collect();
            let qr = cons::analyze(&nt, nr, ns);
            (qr.ni, qr.order, qr.l0)
        };
        let (pivots, free) = order.split_at(rank);
        let nf = free.len();
        let mut data = vec![0.0; nr * nf];
        for f in 0..nf {
            data[f * nf + f] = 1.0;
            for b in 0..rank { data[(nf + b) * nf + f] = -coeffs[f * rank + b]; }
        }
        let name = |i: &usize| c.info.ids[4][*i].clone();
        Ok(Matrix::new(free.iter().chain(pivots).map(name).collect(), free.iter().map(name).collect(), data))
    }

    /// Each floating species/reaction pair occupies one executable stoichiometry index.
    fn stoich_entries(&self) -> Result<Vec<(Vec<SRef>, String)>, String> {
        let c = self.core()?;
        let deps = if self.cma { self.dependent_ids() } else { vec![] };
        let rids = self.reaction_ids();
        let mut entries: Vec<(Vec<SRef>, String)> = vec![];
        let mut positions = HashMap::new();
        for s in self.srefs()? {
            if !c.info.rows.contains(&s.species) || deps.contains(&s.species) { continue; }
            let key = (s.rx, s.species.clone());
            let i = *positions.entry(key).or_insert_with(|| {
                let id = format!("stoich({}, {})", s.species, rids[s.rx]);
                entries.push((vec![], id));
                entries.len() - 1
            });
            entries[i].0.push(s);
        }
        for (refs, id) in &mut entries {
            if let [s] = refs.as_slice() {
                if let Some(name) = &s.id { *id = name.clone(); }
            }
        }
        Ok(entries)
    }

    pub fn stoichiometry_ids(&self) -> Result<Vec<String>, String> {
        Ok(self.stoich_entries()?.into_iter().map(|e| e.1).collect())
    }

    pub fn stoichiometry_id(&self, index: usize) -> Result<String, String> {
        Ok(self.stoich_entries()?.into_iter().nth(index).ok_or("invalid stoichiometry index")?.1)
    }

    pub fn stoichiometry_index(&self, id: &str, reaction: Option<&str>) -> Result<isize, String> {
        let rids = self.reaction_ids();
        Ok(self.stoich_entries()?.iter().position(|(refs, _)| match reaction {
            Some(r) => refs[0].species == id && rids[refs[0].rx] == r,
            None => refs.len() == 1 && refs[0].id.as_deref() == Some(id),
        }).map_or(-1, |i| i as isize))
    }

    fn sref_value(&mut self, s: &SRef, initial: bool) -> Result<f64, String> {
        let c = self.core_mut()?;
        if initial {
            let key = (s.rx, s.reactant, s.pos);
            if let Some(&v) = c.stoich_edit.as_ref().and_then(|e| e.initial.get(&key)) {
                return Ok(v);
            }
        }
        if let Some(id) = &s.id {
            if initial {
                if let Some(v) = c.kpos.get(id).map(|&i| c.xi[i]).or_else(|| c.info.values.get(id).copied()) {
                    return Ok(v);
                }
            } else if let Ok(r) = c.reg(id) {
                return Ok(c.out(&[r])[0]);
            }
        }
        if let Some(expr) = &s.math {
            let doc = self.doc()?;
            let fd = function_definitions(doc.model.as_ref().ok_or("no model")?)?;
            let t = if initial { 0.0 } else { self.time() };
            let mut work = if initial {
                let c = self.core()?;
                let mut w = c.m.work();
                c.m.load(&mut w, &c.kres, &c.xi);
                w
            } else { self.root_work(t, None)? };
            let expr = math::conv(expr).map_err(|e| e.to_string())?;
            return value_at(self.core()?, &mut work, t, &expr, &fd);
        }
        Ok(s.value)
    }

    pub fn stoichiometry(&mut self, index: usize) -> Result<f64, String> {
        let (refs, _) = self.stoich_entries()?.into_iter().nth(index).ok_or("invalid stoichiometry index")?;
        let [s] = refs.as_slice() else { return Err("multiple species references share this index".into()); };
        self.sref_value(s, false)
    }

    pub fn init_stoichiometry(&mut self, index: usize) -> Result<f64, String> {
        let Some((refs, _)) = self.stoich_entries()?.into_iter().nth(index) else { return Ok(f64::NAN); };
        let [s] = refs.as_slice() else { return Err("multiple species references share this index".into()); };
        self.sref_value(s, true)
    }

    pub fn stoichiometry_cell(&mut self, si: usize, ri: usize, initial: bool) -> Result<f64, String> {
        let species = self.floating_ids();
        let Some(s) = species.get(si).filter(|_| ri < self.reaction_ids().len()) else {
            return Ok(f64::NAN);
        };
        if self.cma && self.dependent_ids().contains(s) { return Ok(f64::NAN); }
        self.srefs()?.into_iter().filter(|r| r.species == *s && r.rx == ri).try_fold(0.0, |sum, r| {
            let v = self.sref_value(&r, initial)?;
            Ok(sum + if r.reactant { -v } else { v })
        })
    }

    pub fn set_stoichiometry(&mut self, index: usize, value: f64) -> Result<i32, String> {
        if !value.is_finite() || value.is_sign_negative() { return Err("invalid stoichiometry value".into()); }
        let (refs, _) = self.stoich_entries()?.into_iter().nth(index).ok_or("invalid stoichiometry index")?;
        let [s] = refs.as_slice() else { return Err("multiple species references share this index".into()); };
        self.write_sref(s, value)
    }

    pub fn set_stoichiometry_cell(&mut self, si: usize, ri: usize, value: f64) -> Result<i32, String> {
        if self.cma { return Err("cannot set stoichiometry with conserved moiety analysis enabled".into()); }
        let species = self.floating_ids();
        let Some(s) = species.get(si).filter(|_| ri < self.reaction_ids().len()) else { return Ok(0); };
        let refs: Vec<_> = self.srefs()?.into_iter().filter(|r| r.species == *s && r.rx == ri).collect();
        match refs.as_slice() {
            [] => Ok(0),
            [s] => self.write_sref(s, if s.reactant { -value } else { value }),
            _ => Err("multiple species references share this cell".into()),
        }
    }

    fn write_sref(&mut self, s: &SRef, value: f64) -> Result<i32, String> {
        if self.cma { return Err("cannot set stoichiometry with conserved moiety analysis enabled".into()); }
        if !value.is_finite() { return Err("invalid stoichiometry value".into()); }
        if let Some(id) = &s.id {
            let c = self.core()?;
            if c.info.rate_rule_ids.contains(id) {
                self.set_value(id, value)?;
                return Ok(1);
            }
            if c.info.rule_vars.contains(id) { return Err(format!("'{id}' is governed by an assignment rule")); }
        }
        let initial = self.sref_value(s, true)?;
        let mut doc = self.doc()?;
        let reaction = &mut doc.model.as_mut().ok_or("no model")?.reactions[s.rx];
        let sr = &mut (if s.reactant { &mut reaction.reactants } else { &mut reaction.products })[s.pos];
        if sr.stoichiometry_math.is_some() { return Err("stoichiometry is governed by stoichiometryMath".into()); }
        sr.stoichiometry = Some(value);
        if let Some(id) = &s.id {
            doc.model.as_mut().unwrap().initial_assignments.retain(|a| &a.symbol != id);
        }
        let text = sbml_io::write(&doc);
        let mut core = self.recompile_stoich(&text)?;
        let old = self.core_mut()?;
        let (prev, pending) = old.m.event_state(&old.w);
        core.m.set_event_state(&mut core.w, prev, pending)?;
        core.w.hist.restore(old.w.hist.save());
        core.carry_stoich_edit(old);
        let edit = core.stoich_edit.get_or_insert_with(|| StoichEdit {
            source: self.sbml.clone(), runtime: String::new(), initial: HashMap::new(),
        });
        edit.runtime = text;
        edit.initial.entry((s.rx, s.reactant, s.pos)).or_insert(initial);
        if let Some(i) = s.id.as_ref().and_then(|id| core.kpos.get(id)) { core.x[*i] = value; }
        self.c = Some(core);
        Ok(1)
    }

    fn recompile_stoich(&self, text: &str) -> Result<Core, String> {
        let old = self.core()?;
        let sizes: Vec<_> = old.cconst.iter().filter_map(|id| {
            let i = *old.kpos.get(id)?;
            (old.x[i] != old.x0[i]).then(|| (id.clone(), old.x[i]))
        }).collect();
        let (model, info, compiled) = compile(text, self.base.as_deref(), &sizes)?;
        let mut core = Core::new(model, info, old.nss, &compiled)?;
        for (j, key) in core.keys.iter().enumerate() {
            if let Some(&i) = old.kpos.get(key) {
                (core.x0[j], core.xi[j], core.x[j]) = (old.x0[i], old.xi[i], old.x[i]);
            }
        }
        core.t = old.t;
        core.fresh = old.fresh;
        Ok(core)
    }

    pub(super) fn reset_stoichiometry(&mut self) -> Result<(), String> {
        if self.core()?.stoich_edit.is_some() {
            self.c = Some(self.recompile_stoich(&self.sbml)?);
        }
        Ok(())
    }

    // ------------------------------------------------------------ frequency response

    /// getFrequencyResponse: at the steady state (the model is left there), H(iw) = (iwI - J)^-1 df/dp
    /// of floating species `var` on parameter (or boundary species) `par`, at `points` frequencies w
    /// from `start` over `decades` decades. Rows (w, or w / 2pi with `hz`; |H|, in dB with `db`;
    /// phase in degrees in [0, 360)).
    #[allow(clippy::too_many_arguments)]
    pub fn frequency_response(&mut self, start: f64, decades: f64, points: usize, par: &str, var: &str, db: bool, hz: bool) -> Result<Vec<f64>, String> {
        use nalgebra::{Complex, DMatrix, DVector};
        if start <= 0.0 || !start.is_finite() || !decades.is_finite() || points == 0 {
            return Err("frequency response requires a positive start, finite decades and at least one point".into());
        }
        if !self.global_parameter_ids().iter().chain(&self.boundary_ids()).any(|p| p == par) {
            return Err(format!("frequency response: '{par}' is not a global parameter or boundary species"));
        }
        if !self.floating_ids().iter().any(|s| s == var) {
            return Err(format!("frequency response: '{var}' is not a floating species"));
        }
        if self.steady_state()? > self.steady_state_threshold() {
            return Err("steady-state residual exceeds steadyStateThresh".into());
        }
        let c = self.core_mut()?;
        if c.m.has_algebraic() || c.m.has_delays() {
            return Err("frequency response requires an ordinary differential equation".into());
        }
        let nx = c.m.nx;
        let k = c.m.state_ids.iter().position(|s| s == var)
            .ok_or_else(|| format!("frequency response: '{var}' has no differential equation"))?;
        let p = c.m.pnames.iter().position(|s| s == par)
            .ok_or_else(|| format!("frequency response: '{par}' is not an independent parameter"))?;
        let (nv, nz) = (c.m.nv, nx + c.m.pnames.len());
        let (_, d) = c.rate_derivs();
        let scale = if c.info.amount_params.contains(par) { c.species_raw(par)?.1 } else { 1.0 };
        let partial = |i: usize, q: usize| (0..nv)
            .map(|r| c.m.stoich[i * nv + r] * d[r * nz + q] / if r < c.m.nr { c.m.vol[i] } else { 1.0 })
            .sum::<f64>();
        let j = DMatrix::from_fn(nx, nx, |r, c| partial(r, c));
        let b = DVector::from_iterator(nx, (0..nx).map(|r| Complex::new(partial(r, nx + p) * scale, 0.0)));
        let size = points.checked_mul(3).ok_or("too many frequency points")?;
        let mut out = Vec::new();
        out.try_reserve_exact(size).map_err(|_| "too many frequency points")?;
        for p in 0..points {
            let w = 10f64.powf(start.log10() + if points > 1 { decades * p as f64 / (points - 1) as f64 } else { 0.0 });
            if !w.is_finite() || w <= 0.0 { return Err("frequency grid overflows or underflows".into()); }
            let a = DMatrix::from_fn(nx, nx, |r, c| Complex::new(-j[(r, c)], if r == c { w } else { 0.0 }));
            let h = a.lu().solve(&b).ok_or("frequency response: singular iwI - J")?[k];
            let phase = h.im.atan2(h.re).to_degrees();
            out.extend([if hz { w / (2.0 * std::f64::consts::PI) } else { w },
                        if db { 20.0 * h.norm().log10() } else { h.norm() }, if phase < 0.0 { phase + 360.0 } else { phase }]);
        }
        Ok(out)
    }

    // ------------------------------------------------------------ internalOneStep

    /// Advance one accepted step of the selected integrator.
    pub fn internal_one_step(&mut self, t: f64, h: f64, reset: bool) -> Result<f64, String> {
        self.core()?;
        let end = t + h;
        if !t.is_finite() || !end.is_finite() || end <= t {
            return Err("internalOneStep requires a positive, finite time interval".into());
        }
        let name = self.integ.name;
        let now = match name {
            "rk45" => {
                let (min, max, eps) = (self.integ.num("minimum_time_step"), self.integ.num("maximum_time_step"), self.integ.num("epsilon"));
                self.core_mut()?.integrate_one(t, reset, name, |m, w, next| {
                    crate::fixed::rk45_step(m, w, t, end, min, max, eps, next)
                })?
            }
            "cvode" => {
                let (rtol, tol) = (self.integ.num("relative_tolerance"), self.cvode_tolerances()?);
                let ctl = crate::model::Ctl { var: false, multi: true, ..self.cvode_ctl() };
                self.core_mut()?.integrate_one(t, reset, name, |m, w, next| {
                    w.atol_scale.copy_from_slice(&tol);
                    let old_ctl = std::mem::replace(&mut w.ctl, ctl);
                    let reached = m.integrate_one(w, t, end, rtol, next);
                    w.ctl = old_ctl;
                    reached.ok_or_else(|| "internal integration step failed".into())
                })?
            }
            "euler" | "rk4" => {
                use crate::fixed::Method;
                let (method, subdivision) = if name == "euler" {
                    (Method::Euler, h / self.integ.num("subdivision_steps").max(1.0))
                } else { (Method::Rk4, 0.0) };
                self.core_mut()?.integrate_one(t, reset, name, |m, w, _| {
                    crate::fixed::integrate_fixed(m, w, method, &[t, end], subdivision, 0.0, 0.0, &[], &mut [])?;
                    Ok(end)
                })?
            }
            "gillespie" => {
                let seed = self.seed();
                let ctl = self.integ.ctl();
                let (old_seed, mut rng) = self.step_rng.take()
                    .unwrap_or_else(|| (seed, crate::rand::Mt19937::new(seed as u32)));
                if old_seed != seed { rng = crate::rand::Mt19937::new(seed as u32); }
                let mut returned = end;
                let result = self.core_mut()?.integrate_one(t, reset, name, |m, w, _| {
                    let old_ctl = std::mem::replace(&mut w.ctl, ctl);
                    let times = m.gillespie_one(w, t, end, || rng.uniform());
                    w.ctl = old_ctl;
                    let (time, stop) = times.ok_or("Gillespie integration step failed")?;
                    returned = stop;
                    Ok(time)
                });
                self.step_rng = Some((seed, rng));
                result?;
                returned
            }
            name => return Err(format!("internalOneStep: unknown integrator {name}")),
        };
        Ok(now)
    }

    // ------------------------------------------------------------ getIds

    /// getIds(types): the selection ids whose SelectionRecord flags all lie in `types`, in
    /// libroadrunner's order; `eigen`: RoadRunner's list (with eigenvalues), else ExecutableModel's.
    pub fn ids_of_type(&self, types: i64, eigen: bool) -> Vec<String> {
        let Some(c) = &self.c else { return vec![] };
        let (fs, bs, gp, rx, cs) = (self.floating_ids(), self.boundary_ids(), self.global_parameter_ids(), self.reaction_ids(), self.compartment_ids());
        let deps = if self.cma { self.dependent_ids() } else { vec![] };
        let rr = self.rate_rule_ids();
        let assigned = |k: &String| c.info.rule_vars.contains(k) && !rr.contains(k);
        let ind = |ruled: bool| if ruled { DEP } else { IND };
        let csum = self.conserved_moiety_ids();
        let mut v: Vec<(i64, String)> = vec![];
        let mut add = |f: i64, ids: &mut dyn Iterator<Item = String>| v.extend(ids.map(|s| (f, s)));
        let dp = |s: &String| ind(deps.contains(s) || c.info.rule_vars.contains(s));
        for s in &fs { add(0x424 | dp(s), &mut std::iter::once(s.clone())); }
        for s in &bs { add(0x414 | dp(s), &mut std::iter::once(s.clone())); }
        for s in &fs { add(0x422 | dp(s), &mut std::iter::once(format!("[{s}]"))); }
        for s in &bs { add(0x412 | dp(s), &mut std::iter::once(format!("[{s}]"))); }
        for s in &cs { add(0x440 | ind(c.info.rule_vars.contains(s)), &mut std::iter::once(s.clone())); }
        for s in &gp { add(0x480 | ind(c.info.rule_vars.contains(s)), &mut std::iter::once(s.clone())); }
        add(0x40108, &mut rx.iter().cloned());
        add(0x222, &mut fs.iter().map(|s| format!("init([{s}])")));
        add(0x224, &mut fs.iter().map(|s| format!("init({s})")));
        add(0x212, &mut bs.iter().map(|s| format!("init([{s}])")));
        add(0x214, &mut bs.iter().map(|s| format!("init({s})")));
        add(0x240, &mut cs.iter().map(|s| format!("init({s})")));
        add(0xc0680, &mut gp.iter().map(|s| format!("init({s})")));
        for s in &rr {
            let flags = if fs.contains(s) { 0x40028 } else if bs.contains(s) { 0x40018 }
                else if cs.contains(s) { 0xc0448 } else if gp.contains(s) { 0x88 } else { continue };
            add(if types == 8 { 8 } else { flags }, &mut std::iter::once(format!("{s}'")));
        }
        add(0x40028, &mut fs.iter().filter(|s| !deps.contains(s) && !c.info.rule_vars.contains(*s)).map(|s| format!("{s}'")));
        add(0x4002a, &mut fs.iter().filter(|s| !deps.contains(s) && !assigned(s)).map(|s| format!("[{s}]'")));
        add(0x4001a, &mut bs.iter().filter(|s| !assigned(s)).map(|s| format!("[{s}]'")));
        if types & 0x1010000 != 0 {
            let entries = self.stoich_entries().unwrap_or_default();
            add(0x10000, &mut entries.iter().map(|(_, id)| id.clone()));
            add(0x1000000, &mut entries.iter().filter(|(refs, _)| refs.len() > 1)
                .flat_map(|(refs, _)| refs.iter().filter_map(|r| r.id.clone())));
        }
        let xs: Vec<&String> = fs.iter().chain(bs.iter().filter(|p| !assigned(p))).chain(gp.iter().filter(|p| !assigned(p) && !csum.contains(p))).collect();
        let ps: Vec<&String> = bs.iter().chain(&gp).filter(|p| !assigned(p)).collect();
        for (f, l) in [(0x1800, "uec"), (0x1000, "ec")] { add(f, &mut rx.iter().flat_map(|r| xs.iter().map(move |x| format!("{l}({r}, {x})")))); }
        for (f, l) in [(0x2800, "ucc"), (0x2000, "cc")] { add(f, &mut rx.iter().chain(&fs).flat_map(|r| ps.iter().map(move |x| format!("{l}({r}, {x})")))); }
        if eigen { add(0x4000, &mut self.eigenvalue_ids().into_iter()); }
        v.into_iter().filter(|(f, _)| f & types == *f).map(|x| x.1).collect()
    }

    pub fn supported_id_types(&self, road_runner: bool) -> i64 {
        let model = if self.c.is_some() { 0x10d07ff } else { 0 };
        model | if road_runner { 0x7800 } else { 0 }
    }

    pub fn steady_state_threshold(&self) -> f64 { self.steady_state_threshold }
    pub fn set_steady_state_threshold(&mut self, value: f64) { self.steady_state_threshold = value; }

    pub fn reset_seed(&mut self) -> Result<(), String> {
        let default = Solver::integrator("gillespie").unwrap();
        let Value::UInt(seed) = *default.get("seed")? else { unreachable!("Gillespie seed type"); };
        self.set_seed(seed)?;
        crate::rand::seed(seed);
        if self.c.is_some() {
            if self.sbml.contains("/distrib/") {
                self.rebuild(None)?;
            } else {
                let c = self.core_mut()?;
                c.xi.clone_from(&c.x0);
                self.reset_all()?;
            }
            let c = self.core_mut()?;
            c.m.reseed(&mut c.w, seed);
            c.m.reseed(&mut c.q, seed);
        }
        Ok(())
    }

    // ------------------------------------------------------------ SBML level / version

    /// Convert the loaded SBML to a supported level and version.
    pub fn sbml_converted(&self, level: u32, version: u32) -> Result<String, String> {
        self.convert_sbml(self.get_sbml()?, level, version)
    }

    /// Convert the SBML with the current state as its initial conditions.
    pub fn current_sbml_converted(&mut self, level: u32, version: u32) -> Result<String, String> {
        if level == 0 { return self.get_current_sbml(); }
        let text = self.current_sbml(true)?;
        self.convert_sbml(text, level, version)
    }

    fn convert_sbml(&self, text: String, level: u32, version: u32) -> Result<String, String> {
        if level == 0 { return Ok(text); }
        let mut doc = read_doc(&text)?;
        let version = if version == 0 { doc.version } else { version };
        let comp = text.contains("/comp/");
        let same_version = (level, version) == (doc.level, doc.version);
        let external = comp && text.contains("externalModelDefinition");
        if same_version && !external { return Ok(text); }
        if level == 2 && crate::sbml::convert::has_level_three_packages(&doc) {
            return Err(format!("Error converting SBML to level {level} version {version}: Level 3 packages"));
        }
        if comp {
            if !same_version && external {
                return Err(format!("Error converting SBML to level {level} version {version}: external model definitions"));
            }
            let flattened = sbml_io::comp::flatten(&text, self.base.as_deref())?;
            if same_version { return Ok(flattened); }
            doc = read_doc(&flattened)?;
        }
        crate::sbml::convert::convert(&mut doc, level, version)?;
        Ok(sbml_io::write(&doc))
    }

    // ------------------------------------------------------------ events, piecewise triggers

    pub fn event_ids(&self) -> Vec<String> {
        self.doc().ok().and_then(|d| d.model).map_or_else(Vec::new, |m| {
            m.events.iter().filter(|e| e.trigger.as_ref().is_some_and(|t| t.math.is_some()))
                .map(|e| e.id.clone().unwrap_or_default()).collect()
        })
    }

    pub fn num_events(&self) -> usize {
        self.c.as_ref().map_or(0, |c| c.m.event_src.len())
    }

    /// ExecutableModel vectors store rate-rule values first, then independent species amounts.
    pub fn state_vector_ids(&self) -> Vec<String> {
        let rr = self.rate_rule_ids();
        let deps = if self.cma { self.dependent_ids() } else { vec![] };
        let mut ids = rr.clone();
        if let Some(c) = &self.c {
            ids.extend(self.floating_ids().into_iter().filter(|s| {
                !rr.contains(s) && !deps.contains(s) && !c.info.rule_vars.contains(s)
            }));
        }
        ids
    }

    fn root_work(&mut self, t: f64, state: Option<&[f64]>) -> Result<crate::model::Work, String> {
        let ids = self.state_vector_ids();
        let c = self.core_mut()?;
        let mut work = c.m.work();
        let mut vals = c.x.clone();
        let mut history = c.live.then(|| c.w.hist.snapshot());
        if let Some(state) = state {
            if ids.len() != state.len() { return Err(format!("expected {} state values", ids.len())); }
            // Use the supplied compartment values when converting species amounts.
            for (s, &v) in ids.iter().zip(state) {
                if !c.is_species(s) {
                    if let Some(&k) = c.kpos.get(s) { vals[k] = v; }
                }
            }
            c.m.load(&mut work, &c.kres, &vals);
            if let Some(saved) = history.take() { work.hist.restore(saved); }
            for (s, &v) in ids.iter().zip(state) {
                if !c.is_species(s) { continue; }
                let Some(&k) = c.kpos.get(s) else { continue };
                let mut vol = [1.0];
                c.m.record(&mut work, t, &[c.reg(c.comp_of(s))?], &mut vol);
                vals[k] = if c.info.amount_params.contains(s) ||
                    (c.info.rate_rule_ids.contains(s) && !c.info.hosu.contains(s)) { v } else { v / vol[0] };
            }
            if c.live { history = Some(work.hist.save()); }
        }
        c.m.load(&mut work, &c.kres, &vals);
        if let Some(saved) = history { work.hist.restore(saved); }
        Ok(work)
    }

    pub fn event_roots(&mut self, time: Option<f64>, state: Option<&[f64]>) -> Result<Vec<f64>, String> {
        let t = time.unwrap_or(self.time());
        let mut work = self.root_work(t, state)?;
        Ok(self.core()?.m.event_roots(&mut work, t))
    }

    pub fn num_piecewise_triggers(&self) -> Result<usize, String> {
        Ok(self.core()?.info.piecewise_roots.len())
    }

    pub fn piecewise_trigger_roots(&mut self, time: Option<f64>, state: Option<&[f64]>) -> Result<Vec<f64>, String> {
        let t = time.unwrap_or(self.time());
        let mut w = self.root_work(t, state)?;
        let c = self.core()?;
        let mut roots = vec![0.0; c.info.piecewise_roots.len()];
        c.m.record(&mut w, t, &c.info.piecewise_roots, &mut roots);
        for root in &mut roots { *root = if *root != 0.0 { 1.0 } else { -1.0 }; }
        Ok(roots)
    }

    /// ExecutableModel::getRandom: a uniform number in (0, 1) from the model's seeded stream.
    pub fn random(&mut self) -> f64 { crate::model::Rng::new(self.run_seed()).unif() }
}

fn function_definitions(model: &sbml_io::Model) -> Result<Fdefs, String> {
    let mut fd = Fdefs::new();
    for f in &model.function_definitions {
        let mut expr = f.math.as_ref();
        while let Some(Math::Semantics { math, .. }) = expr { expr = Some(math); }
        if let Some(Math::Lambda { args, body }) = expr {
            fd.insert(f.id.clone(), (args.clone(), math::conv(body).map_err(|e| e.to_string())?));
        }
    }
    Ok(fd)
}

/// Evaluate document math at a scratch state, reading symbols in their SBML units.
fn value_at(c: &Core, w: &mut crate::model::Work, t: f64, expr: &E, fd: &Fdefs) -> Result<f64, String> {
    fn at_time(e: &E, t: f64, fd: &Fdefs, depth: usize) -> Result<E, String> {
        if depth > 50 { return Err("recursive SBML function".into()); }
        Ok(match e {
            E::Time => E::Num(t),
            E::Call(f, a) => {
                let (args, body) = fd.get(f).ok_or_else(|| format!("unknown SBML function '{f}'"))?;
                at_time(&body.subst(&args.iter().cloned().zip(a.iter().cloned()).collect()), t, fd, depth + 1)?
            }
            E::Op(f, _) if f == "delay" => return Err("delay() in a piecewise or stoichiometry expression needs compiled history".into()),
            E::Op(f, a) => E::Op(f.clone(), a.iter().map(|x| at_time(x, t, fd, depth)).collect::<Result<_, _>>()?),
            e => e.clone(),
        })
    }
    let expr = at_time(expr, t, fd, 0)?;
    math::eval(&expr, &mut |n| {
        let error = |s: String| crate::sbml::Error::Invalid(s);
        let mut v = [0.0; 2];
        let reg = c.reg(n).map_err(error)?;
        if c.is_species(n) && c.info.hosu.contains(n) != c.info.amount_params.contains(n) {
            let vol = c.reg(c.comp_of(n)).map_err(error)?;
            c.m.record(w, t, &[reg, vol], &mut v);
            return Ok(if c.info.hosu.contains(n) { v[0] * v[1] } else { v[0] / v[1] });
        }
        c.m.record(w, t, &[reg], &mut v[..1]);
        Ok(v[0])
    }, fd).map_err(|e| e.to_string())
}

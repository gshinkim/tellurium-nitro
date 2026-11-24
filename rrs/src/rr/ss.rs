//! Steady state with libroadrunner's solver settings (rrRoadRunner.cpp steadyState: the solver, then
//! PresimulationProgramDecorator, then ApproxSteadyStateDecorator):
//! 1. the solver's Newton from the current state (crate::model::NewtonOpts from the settings);
//! 2. allow_presimulation: for each t in presimulation_times, reset() (species back, time 0),
//!    integrate to t (CVODE defaults: rtol 1e-6, atol 1e-12), Newton again; the model is left at t;
//! 3. allow_approx: reset(), integrate to approx_time - h and approx_time (h = approx_time /
//!    approx_maximum_steps); accepted when sum |state(end - h) - state(end)| / h <= approx_tolerance
//!    (the model is left at approx_time; the value returned is that sum).
//!
//! Deliberate differences (Nitro is exact where libroadrunner is loose or buggy):
//! - newton / newton_linesearch converge to machine precision unless func_norm_tol / scaled_step_tol
//!   are set (KINSOL stops at |f| <= ~6e-6 and calls a state with |f| < 6e-8 a steady state);
//! - allow_approx applies to every solver (libroadrunner's NLEQ exceptions escape the approximation
//!   decorator when presimulation is off);
//! - an initial state with zero concentrations is a valid start for newton (KINSOL's positivity
//!   constraint rejects it as illegal input);
//! - auto_moiety_analysis: the engine always solves on the moiety-reduced system, the answer
//!   libroadrunner gets with it on.
// ponytail: settings that shape only the iteration path keep their values but cannot change the
// converged state: presimulation_time and presimulation_maximum_steps (libroadrunner 2.11 runs
// presimulation_times only), broyden_method, linearity, and KINSOL's print_level, no_init_setup,
// no_res_monitoring, max_setup_calls, max_subsetup_calls, eta_*, res_mon_*, no_min_eps,
// max_newton_step, max_beta_fails, rel_err_func.

use super::settings::Solver;
use super::Core;
use crate::model::NewtonOpts;

pub(super) struct SsCfg {
    pub newton: NewtonOpts,
    pub presim: Option<Vec<f64>>,
    /// (approx_tolerance, approx_maximum_steps, approx_time)
    pub approx: Option<(f64, f64, f64)>,
    /// NLEQ evaluates d/dt at time 0 (getStateVectorRate(0, ...)), which leaves the model at t = 0
    pub time0: bool,
}

impl SsCfg {
    pub fn of(s: &Solver) -> Result<SsCfg, String> {
        let int = |k: &str| s.num(k).trunc();
        let newton = if s.name.starts_with("nleq") {
            let rtol = s.num("relative_tolerance");
            if !(rtol > 0.0) { return Err("Nonpositive value for RTOL supplied to NLEQ".into()); }
            let iters = int("maximum_iterations");
            if iters < 0.0 { return Err("One or more fields specified in IOPT are invalid (NLEQ)".into()); }
            let md = s.num("minimum_damping");
            // NLEQ's default FCMIN when RWK(22) is 0: 1e-4, 1e-8 for extremely nonlinear problems
            let md = if md > 0.0 { md } else if int("linearity") >= 4.0 { 1e-8 } else { 1e-4 };
            NewtonOpts { max_iter: if iters == 0.0 { 50 } else { iters as usize }, rtol, min_damping: md, linesearch: true,
                         rank_deficient: s.name == "nleq2", allow_negative: false, fnorm_tol: 0.0, step_tol: 0.0 }
        } else {
            let eta = s.get("eta_form")?;
            if !matches!(eta, super::settings::Value::Str(e) if ["eta_choice1", "eta_choice2", "eta_constant"].contains(&e.as_str())) {
                return Err(format!("\"{eta:?}\". Valid options are \"eta_choice1\", \"eta_choice2\", \"eta_constant\", "));
            }
            let iters = int("num_max_iters");
            NewtonOpts { max_iter: if iters <= 0.0 { 200 } else { iters as usize }, rtol: 1e-12, min_damping: 1e-10,
                         linesearch: s.name == "newton_linesearch", rank_deficient: false,
                         allow_negative: s.get("allow_negative")?.bool(), fnorm_tol: s.num("func_norm_tol").max(0.0),
                         step_tol: s.num("scaled_step_tol").max(0.0) }
        };
        let presim = match s.get("allow_presimulation")?.bool() {
            true => match s.get("presimulation_times")? { super::settings::Value::Doubles(v) => Some(v.clone()), v => Some(vec![v.f64()]) },
            false => None,
        };
        let approx = s.get("allow_approx")?.bool().then(|| (s.num("approx_tolerance"), int("approx_maximum_steps"), s.num("approx_time")));
        Ok(SsCfg { newton, presim, approx, time0: s.name.starts_with("nleq") })
    }
}

impl Core {
    /// Newton from the current state into x (the solver's work is left in q).
    fn newton_here(&mut self, o: &NewtonOpts, time0: bool) -> bool {
        if time0 { self.t = 0.0; self.touch(); }
        self.m.load(&mut self.q, &self.kres, &self.x);
        let ok = self.m.newton_solve(&mut self.q, o);
        self.qok = false;
        if ok {
            let s = self.q.s.clone();
            self.take_states(&s);
        }
        ok
    }

    /// Integrate over times with a fresh CVODEIntegrator's defaults (rtol 1e-6, atol 1e-12 scaled
    /// per state as roadrunner's).
    fn advance_default(&mut self, times: &[f64]) -> Result<(), String> {
        let tol: Vec<f64> = self.atol_scale()?.iter().map(|s| 1e-12 * s).collect();
        self.advance(times, &[], 1e-6, &tol, crate::model::Ctl::default()).map(|_| ())
    }

    /// Floating species amounts (independent ones with moieties, as auto_moiety_analysis) and
    /// rate-rule states: roadrunner's state vector.
    fn state_vector(&mut self) -> Result<Vec<f64>, String> {
        let ids: Vec<String> = if self.ls().nd() > 0 { self.ls().order[..self.ls().ni].iter().map(|&i| self.ls_rows[i].clone()).collect() } else { self.fids().to_vec() };
        let mut out = vec![];
        let ids: Vec<String> = ids.into_iter().filter(|s| self.m.state_ids.contains(s)).collect();
        for s in &ids { out.push(self.species(s, false)?); }
        for r in self.info.rate_rule_ids.clone() {
            if !self.is_species(&r) { if let Some(&j) = self.kpos.get(&r) { out.push(self.x[j]); } }
        }
        Ok(out)
    }

    /// Steady state from the current state as the settings say (see the module doc). Returns
    /// libroadrunner's value: |d[S]/dt| (Newton) or the approximation's tolerance; q holds the
    /// solution (states, rates).
    pub(super) fn solve_with(&mut self, cfg: &SsCfg) -> Result<Option<f64>, String> {
        if self.nss { return Err("steady state: the model has events or time-dependent rates".into()); }
        let start = (self.x.clone(), self.t);
        match self.solve_settings(cfg) {
            Err(e) if self.robust_from(start) => { let _ = e; Ok(None) }
            r => r,
        }
    }

    /// Where libroadrunner's settings give up but a steady state exists (e.g. Newton diverges from a
    /// transient state), the engine's own solver finds it: damped Newton, then integration toward it
    /// with Newton after each horizon (model::steady_state). A model with no steady state still fails.
    fn robust_from(&mut self, (x, t): (Vec<f64>, f64)) -> bool {
        (self.x, self.t) = (x, t);
        self.touch();
        self.m.load(&mut self.q, &self.kres, &self.x);
        let ok = self.m.steady_state(&mut self.q);
        self.qok = false;
        if ok { let s = self.q.s.clone(); self.take_states(&s); }
        ok
    }

    fn solve_settings(&mut self, cfg: &SsCfg) -> Result<Option<f64>, String> {
        if self.newton_here(&cfg.newton, cfg.time0) { return Ok(None); }
        let mut err = "steady state not found (Newton failed from the current state)".to_string();
        if let Some(times) = &cfg.presim {
            for &t in times {
                self.reset()?;
                if t != 0.0 { self.advance_default(&[0.0, t])?; }
                if self.newton_here(&cfg.newton, cfg.time0) { return Ok(None); }
            }
            err = format!("SteadyStateSolver failed to solve presimulation program: {err}");
        }
        if let Some((thresh, steps, end)) = cfg.approx {
            let h = end / steps;
            self.reset()?;
            self.advance_default(&[0.0, end - h])?;
            let a = self.state_vector()?;
            self.advance_default(&[end - h, end])?;
            let b = self.state_vector()?;
            let tol: f64 = a.iter().zip(&b).map(|(x, y)| ((x - y) / h).abs()).sum();
            self.m.load(&mut self.q, &self.kres, &self.x);
            let mut r = vec![0.0; self.m.nx];
            self.m.state_rates(&mut self.q, self.t, &mut r);
            if !(tol <= thresh) {
                return Err(format!("Failed to converge while running steady state approximation. Tolerance {tol} is not \
                                    greater than threshold {thresh}. Try increasing the time point at which the approximation \
                                    is conducted (with the \"approx_time\" argument) or increasing the number of steps \
                                    parameter (with the \"approx_maximum_steps\" argument). Be aware that your model might \
                                    not have a steady state"));
            }
            return Ok(Some(tol));
        }
        Err(err)
    }
}

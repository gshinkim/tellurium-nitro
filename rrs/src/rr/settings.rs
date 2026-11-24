//! Integrator and steady-state solver registries with libroadrunner's names, setting types and
//! defaults.

#[derive(Debug, Clone, PartialEq)]
pub enum Value { Bool(bool), Int(i64), UInt(u64), Double(f64), Str(String), Doubles(Vec<f64>) }

impl Value {
    /// `v` converted to this value's type (numbers and bools convert when exact).
    fn coerce(&self, v: Value) -> Result<Value, String> {
        use Value::*;
        let num = match &v { Bool(b) => Some(*b as i64 as f64), Int(i) => Some(*i as f64), UInt(u) => Some(*u as f64), Double(d) => Some(*d), _ => None };
        let int = num.filter(|x| x.fract() == 0.0 && x.is_finite());
        Ok(match (self, v) {
            (Bool(_), v @ Bool(_)) | (Int(_), v @ Int(_)) | (UInt(_), v @ UInt(_)) | (Double(_), v @ Double(_))
            | (Str(_), v @ Str(_)) | (Doubles(_), v @ Doubles(_)) => v,
            (Bool(_), _) if int == Some(0.0) || int == Some(1.0) => Bool(int == Some(1.0)),
            (Int(_), _) if int.is_some() => Int(int.unwrap() as i64),
            (UInt(_), Int(i)) if i >= 0 => UInt(i as u64),
            (UInt(_), Double(d)) if int.is_some() && d >= 0.0 => UInt(d as u64),
            (Double(_), _) if num.is_some() => Double(num.unwrap()),
            (want, got) => return Err(format!("expected a value like {want:?}, got {got:?}")),
        })
    }
    pub fn f64(&self) -> f64 {
        match self { Value::Bool(b) => *b as i64 as f64, Value::Int(i) => *i as f64, Value::UInt(u) => *u as f64, Value::Double(d) => *d, _ => f64::NAN }
    }
    pub fn bool(&self) -> bool { self.f64() != 0.0 }
}

#[derive(Debug, Clone)]
pub struct Param { pub name: &'static str, pub display: &'static str, pub hint: &'static str, pub value: Value, default: Value }

#[derive(Debug, Clone)]
pub struct Solver { pub name: &'static str, pub hint: &'static str, pub description: &'static str, pub params: Vec<Param> }

pub const INTEGRATORS: [&str; 5] = ["cvode", "gillespie", "rk4", "rk45", "euler"];
pub const STEADY_STATE_SOLVERS: [&str; 4] = ["nleq1", "nleq2", "newton", "newton_linesearch"];

use Value::*;

fn p(name: &'static str, display: &'static str, hint: &'static str, v: Value) -> Param {
    Param { name, display, hint, default: v.clone(), value: v }
}

fn clock_seed() -> u64 {
    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_nanos() as u64) % (1 << 52)
}

impl Solver {
    pub fn integrator(name: &str) -> Option<Solver> {
        let (name, hint, description, params) = match name {
            "cvode" => ("cvode", "Deterministic ODE solver",
                "Deterministic ODE solver for stiff and non-stiff systems (Nitro: variable-order BDF / Rodas4 with root-located events, CVODE's settings).",
                vec![
                    p("relative_tolerance", "Relative Tolerance", "Specifies the scalar relative tolerance.", Double(1e-6)),
                    p("absolute_tolerance", "Absolute Tolerance", "Specifies the scalar absolute tolerance (scaled per species by its initial value).", Double(1e-12)),
                    p("stiff", "Stiff", "Specifies whether the integrator attempts to solve stiff equations.", Bool(true)),
                    p("maximum_bdf_order", "Maximum BDF Order", "Specifies the maximum order for BDF.", Int(5)),
                    p("maximum_adams_order", "Maximum Adams Order", "Specifies the maximum order for Adams-Moulton.", Int(12)),
                    p("maximum_num_steps", "Maximum Number of Steps", "Specifies the maximum number of steps per time step.", Int(20000)),
                    p("maximum_time_step", "Maximum Time Step", "Specifies the maximum absolute step size (0: no limit).", Double(0.0)),
                    p("minimum_time_step", "Minimum Time Step", "Specifies the minimum absolute step size (0: no limit).", Double(0.0)),
                    p("initial_time_step", "Initial Time Step", "Specifies the initial time step size (0: estimated).", Double(0.0)),
                    p("multiple_steps", "Multiple Steps", "Perform a multiple time step simulation.", Bool(false)),
                    p("variable_step_size", "Variable Step Size", "Perform a variable time step simulation.", Bool(false)),
                    p("max_output_rows", "Maximum Output Rows", "For variable step size simulations, the maximum number of rows in the result.", Int(100000)),
                ]),
            "gillespie" => ("gillespie", "Gillespie Direct Method SSA",
                "Gillespie direct method stochastic simulation; with variable_step_size, one output row per reaction firing.",
                vec![
                    p("seed", "Seed", "Seed of the random number generator.", UInt(clock_seed())),
                    p("variable_step_size", "Variable Step Size", "Perform a variable time step simulation.", Bool(true)),
                    p("minimum_time_step", "Minimum Time Step", "Specifies the minimum absolute step size.", Double(0.0)),
                    p("maximum_time_step", "Maximum Time Step", "Specifies the maximum absolute step size.", Double(0.0)),
                    p("nonnegative", "Non-negative species only", "Prevents species amounts from going negative.", Bool(false)),
                    p("max_output_rows", "Maximum Output Rows", "For variable step size simulations, the maximum number of rows in the result.", Int(100000)),
                    p("maximum_num_steps", "Maximum Number of Steps", "Maximum number of steps (0: no limit).", Int(0)),
                ]),
            "rk4" => ("rk4", "Internal RK4 ODE solver", "Fixed-step fourth-order Runge-Kutta.", vec![]),
            "rk45" => ("rk45", "Internal RK45 ODE solver", "Runge-Kutta-Fehlberg 4(5) with adaptive steps.", vec![
                p("variable_step_size", "Variable Step Size", "Perform a variable time step simulation.", Bool(true)),
                p("minimum_time_step", "Minimum Time Step", "Specifies the minimum absolute step size.", Double(1e-12)),
                p("maximum_time_step", "Maximum Time Step", "Specifies the maximum absolute step size.", Double(1.0)),
                p("epsilon", "Maximum error tolerance", "Maximum error tolerance.", Double(1e-12)),
                p("max_output_rows", "Maximum Output Rows", "For variable step size simulations, the maximum number of rows in the result.", Int(100000)),
            ]),
            "euler" => ("euler", "A simple Euler integrator", "Forward Euler with a fixed step.", vec![
                p("subdivision_steps", "Subdivision Steps", "Euler steps per output interval.", Int(1)),
            ]),
            _ => return None,
        };
        Some(Solver { name, hint, description, params })
    }

    pub fn steady_state(name: &str) -> Option<Solver> {
        let name = *STEADY_STATE_SOLVERS.iter().find(|&&n| n == name)?;
        let newton = name.starts_with("newton");
        let mut params = vec![
            p("auto_moiety_analysis", "Auto Moiety Analysis", "Use conserved moiety analysis during the steady state solve.", Bool(true)),
            p("allow_presimulation", "Allow Presimulation", "Integrate before solving when Newton fails.", Bool(newton)),
            p("presimulation_time", "Presimulation Time", "End time of presimulation.", Double(if newton { 5.0 } else { 100.0 })),
            p("presimulation_times", "Presimulation Times", "Successive presimulation end times.", Doubles(vec![0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0])),
            p("presimulation_maximum_steps", "Presimulation Maximum Steps", "Maximum integrator steps in presimulation.", if newton { Double(500.0) } else { Int(100) }),
            p("allow_approx", "Allow Approximation", "Accept a long simulation's end point as the steady state.", Bool(newton)),
            p("approx_tolerance", "Approximation Tolerance", "Tolerance of the approximate steady state.", Double(1e-6)),
            p("approx_maximum_steps", "Approximation Maximum Steps", "Maximum integrator steps when approximating.", Int(10000)),
            p("approx_time", "Approximation Time", "End time when approximating.", if newton { Int(10000) } else { Double(10000.0) }),
        ];
        if newton {
            params.push(p("num_max_iters", "Maximum Iterations", "Maximum Newton iterations.", Int(200)));
            params.push(p("allow_negative", "Allow Negative", "Accept negative concentrations.", Bool(false)));
        } else {
            params.push(p("relative_tolerance", "Relative Tolerance", "Relative tolerance of the solution.", Double(1e-12)));
            params.push(p("maximum_iterations", "Maximum Iterations", "Maximum Newton iterations.", Int(100)));
            params.push(p("minimum_damping", "Minimum Damping", "Smallest damping factor.", Double(1e-20)));
            params.push(p("broyden_method", "Broyden Method", "Use Broyden updates (0: off).", Int(0)));
            params.push(p("linearity", "Linearity", "Problem linearity hint (1-4).", Int(3)));
        }
        let hint = if newton { if name == "newton" { "Newton iteration" } else { "Newton iteration with linesearch strategy" } }
                   else { "Steady-state nonlinear systems of equations solver" };
        // ponytail: every solver runs the engine's damped Newton with presimulation fallback; the
        // settings are kept for the API, not honoured individually
        Some(Solver { name, hint, description: "Damped Newton on the reduced system, integrating towards the steady state when Newton fails.", params })
    }

    fn find(&self, key: &str) -> Result<usize, String> {
        self.params.iter().position(|p| p.name == key).ok_or_else(|| format!("{}: no setting '{key}'", self.name))
    }
    pub fn get(&self, key: &str) -> Result<&Value, String> { Ok(&self.params[self.find(key)?].value) }
    pub fn set(&mut self, key: &str, v: Value) -> Result<(), String> {
        let i = self.find(key)?;
        let p = &mut self.params[i];
        // cvode's absolute_tolerance: a scalar, or one value per state (as libroadrunner's)
        if let (true, Value::Doubles(_)) = (self.name == "cvode" && key == "absolute_tolerance", &v) { p.value = v; return Ok(()); }
        p.value = p.default.coerce(v).map_err(|e| format!("{}.{key}: {e}", self.name))?;
        Ok(())
    }
    pub fn reset(&mut self) { for p in &mut self.params { p.value = p.default.clone(); } }
    /// Numeric value of a setting (NaN if missing or not a number).
    pub fn num(&self, key: &str) -> f64 { self.get(key).map_or(f64::NAN, Value::f64) }

    /// The settings the engine's cvode / gillespie solvers take (model::Ctl; missing or <= 0: unset).
    /// Nonstiff ODE runs use Rodas4, a one-step method with a fixed order and the same
    /// error tolerances. maximum_adams_order is result-neutral because no Adams history is used.
    pub fn ctl(&self) -> crate::model::Ctl {
        let n = |k: &str| { let x = self.num(k); if x > 0.0 { x } else { 0.0 } };
        let b = |k: &str| self.get(k).is_ok_and(Value::bool);
        crate::model::Ctl {
            h0: n("initial_time_step"), hmin: n("minimum_time_step"), hmax: n("maximum_time_step"),
            max_steps: n("maximum_num_steps") as usize,
            max_order: if b("stiff") { n("maximum_bdf_order") as usize } else { 0 },
            nonstiff: self.name == "cvode" && !b("stiff"),
            var: b("variable_step_size"), max_rows: n("max_output_rows") as usize, multi: b("multiple_steps"), nonneg: b("nonnegative"),
        }
    }
}

//! libroadrunner's "euler", "rk4" and "rk45" (Runge-Kutta-Fehlberg) integrators, step for step as
//! roadrunner 2.11 (EulerIntegrator.h, RK4Integrator.cpp, RK45Integrator.cpp) driven by
//! RoadRunner::simulate's fixed output grid (rk45 with variable_step_size = false).
//!
//! Events, as libroadrunner: all three fire the events due at the start time. After that only euler
//! looks at events, at the end of each step (and only when some trigger is true there); rk4 and rk45
//! ignore them for the rest of the run. Algebraic rules and delay() are rejected (libroadrunner cannot
//! load algebraic rules; delays would need the history the adaptive solvers keep).
use crate::model::{Model, Work};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method { Euler, Rk4, Rk45 }

/// Time course on the output grid `times` from the loaded state into out (len(times) x len(outs),
/// row-major), as `Model::simulate`. Each output interval is integrated separately:
/// - Euler/Rk4: `step` is the step size; the interval is split into ceil(interval / step) equal
///   steps (step <= 0 or >= interval: one step per interval, libroadrunner's default; euler's
///   subdivision_steps = k is step = interval / k).
/// - Rk45: `step` is libroadrunner's maximum_time_step, which is also its first step (<= 0: 1.0);
///   `atol` is its epsilon (bound on the 2-norm of the per-step error estimate, in amounts);
///   minimum_time_step is 1e-12. `rtol` is unused; epsilon is clamped by `rk45_tol` where it is
///   meaningless for the state scale.
///
/// On success w.s / params hold the final state.
#[allow(clippy::too_many_arguments, clippy::neg_cmp_op_on_partial_ord)] // !(a <= b): NaN fails
pub fn integrate_fixed(m: &Model, w: &mut Work, method: Method, times: &[f64], step: f64, rtol: f64, atol: f64,
                       outs: &[u32], out: &mut [f64]) -> Result<(), String> {
    let _ = rtol;
    if m.has_algebraic() { return Err("euler/rk4/rk45: algebraic rules are not supported (libroadrunner cannot load them either)".into()); }
    if m.has_delays() { return Err("euler/rk4/rk45: delay() is not supported".into()); }
    let (nx, no, nt) = (m.nx, outs.len(), times.len());
    if nt == 0 || times.windows(2).any(|t| !(t[1] >= t[0])) { return Err("times must be non-empty and non-decreasing".into()); }
    if out.len() != nt * no { return Err(format!("out has {} entries, expected {}", out.len(), nt * no)); }
    let ev = m.has_events();
    m.events_reset(w);
    if ev && m.apply_events(w, times[0]).is_none() { return Err("event cascade does not terminate".into()); }
    let ne = m.event_state(w).0.len();
    if method == Method::Euler && ev { // EulerIntegrator's trigger memory starts all false
        let pend = m.event_state(w).1;
        m.set_event_state(w, vec![false; ne], pend)?;
    }
    m.record(w, times[0], outs, &mut out[..no]);
    let z = || vec![0.0; nx];
    let (mut y, mut yt, mut k1, mut k2, mut k3, mut k4, mut k5, mut k6) = (w.s.clone(), z(), z(), z(), z(), z(), z(), z());
    // rk45 settings and the step size it carries from one output interval to the next
    let hmax = if step > 0.0 { step } else { 1.0 };
    let (hmin, eps) = (1e-12f64.min(hmax), atol);
    let mut hcur = 0.0f64;
    let mut budget = 100_000_000usize;
    let axpy = |y: &mut [f64], a: f64, x: &[f64]| for (o, &v) in y.iter_mut().zip(x) { *o += a * v };
    for i in 1..nt {
        let (t0, t1) = (times[i - 1], times[i]);
        let dt = t1 - t0;
        if dt > 0.0 {
            let n = if step > 0.0 && step < dt { ((dt / step) * (1.0 - 1e-12)).ceil().max(1.0) as usize } else { 1 };
            match method {
                Method::Euler => {
                    let h = dt / n as f64;
                    for _ in 0..n {
                        // as EulerIntegrator::integrate: every sub-step evaluates the rates at the
                        // interval start t0 and ends at "t0 + h" (it never advances t0)
                        m.deriv(w, t0, &y, &mut k1);
                        axpy(&mut y, h, &k1);
                        if ev {
                            let te = t0 + h;
                            w.s.copy_from_slice(&y);
                            if m.any_trigger(w, te) {
                                if m.apply_events(w, te).is_none() { return Err("event cascade does not terminate".into()); }
                                y.copy_from_slice(&w.s);
                            } else { // no trigger true: libroadrunner skips applyEvents, remembers all-false
                                let pend = m.event_state(w).1;
                                m.set_event_state(w, vec![false; ne], pend)?;
                            }
                        }
                    }
                }
                Method::Rk4 => {
                    let h = dt / n as f64;
                    for s in 0..n {
                        let t = if s == 0 { t0 } else { t0 + s as f64 * h };
                        m.deriv(w, t, &y, &mut k1);
                        let a = h / 2.;
                        yt.copy_from_slice(&y); axpy(&mut yt, a, &k1);
                        m.deriv(w, t + a, &yt, &mut k2);
                        yt.copy_from_slice(&y); axpy(&mut yt, a, &k2);
                        m.deriv(w, t + a, &yt, &mut k3);
                        yt.copy_from_slice(&y); axpy(&mut yt, h, &k3);
                        m.deriv(w, t + h, &yt, &mut k4);
                        axpy(&mut k1, 2., &k2);
                        axpy(&mut k1, 2., &k3);
                        axpy(&mut k1, 1., &k4);
                        axpy(&mut y, h / 6., &k1);
                    }
                }
                Method::Rk45 => {
                    let mut hp = (if hcur > 0.0 { hcur } else { hmax }).max(hmin).min(hmax);
                    let mut tc = t0;
                    loop {
                        budget = budget.checked_sub(1).ok_or("rk45: too many steps")?;
                        let htry = hp;
                        let (mut h, mut reach) = (htry, false);
                        if h >= t1 - tc { h = t1 - tc; reach = true; }
                        m.deriv(w, tc, &y, &mut k1);
                        yt.copy_from_slice(&y); axpy(&mut yt, h / 4., &k1);
                        m.deriv(w, tc + h / 4., &yt, &mut k2);
                        yt.copy_from_slice(&y); axpy(&mut yt, 3. * h / 32., &k1); axpy(&mut yt, 9. * h / 32., &k2);
                        m.deriv(w, tc + 3. * h / 8., &yt, &mut k3);
                        yt.copy_from_slice(&y);
                        axpy(&mut yt, 1932. * h / 2197., &k1); axpy(&mut yt, -7200. * h / 2197., &k2); axpy(&mut yt, 7296. * h / 2197., &k3);
                        m.deriv(w, tc + 12. * h / 13., &yt, &mut k4);
                        yt.copy_from_slice(&y);
                        axpy(&mut yt, 439. * h / 216., &k1); axpy(&mut yt, -8. * h, &k2);
                        axpy(&mut yt, 3680. * h / 513., &k3); axpy(&mut yt, -845. * h / 4104., &k4);
                        m.deriv(w, tc + h, &yt, &mut k5);
                        yt.copy_from_slice(&y);
                        axpy(&mut yt, -8. * h / 27., &k1); axpy(&mut yt, 2. * h, &k2); axpy(&mut yt, -3544. * h / 2565., &k3);
                        axpy(&mut yt, 1859. * h / 4104., &k4); axpy(&mut yt, -11. * h / 40., &k5);
                        m.deriv(w, tc + h / 2., &yt, &mut k6);
                        // libroadrunner's error: |k1/360 - 128/4275 k3 - 2197/75240 k4 + k5/50 + 2/55 k6|_2
                        // (not scaled by h), over its state vector of amounts
                        let error = (0..nx).map(|j| {
                            let e = k1[j] / 360. - 128. / 4275. * k3[j] - 2197. / 75240. * k4[j] + k5[j] / 50. + 2. / 55. * k6[j];
                            (e * m.vol[j]).powi(2)
                        }).sum::<f64>().sqrt();
                        let tol = rk45_tol(eps, (0..nx).map(|j| (y[j] * m.vol[j]).powi(2)).sum::<f64>().sqrt());
                        let q = if !error.is_finite() { 0.1 } else if error == 0. { 4. } else { (0.84 * (tol / error).powf(0.25)).clamp(0.1, 4.) };
                        if !(error <= tol) && htry > hmin { hp = (q * h).max(hmin); continue; }
                        axpy(&mut y, 1408. * h / 2565., &k3);
                        axpy(&mut y, 2197. / 4104. * h, &k4);
                        axpy(&mut y, -1. / 5. * h, &k5);
                        axpy(&mut y, 25. / 216. * h, &k1);
                        tc = if reach { t1 } else { tc + h };
                        if tc >= t1 { tc = t1; reach = true; }
                        hp = q * h;
                        if h < htry && hp < htry { hp = htry; } // a step cut to end at t1 does not lower it
                        hp = hp.max(hmin).min(hmax);
                        hcur = hp;
                        if reach { break; }
                    }
                }
            }
            w.s.copy_from_slice(&y);
        }
        m.record(w, t1, outs, &mut out[i * no..(i + 1) * no]);
    }
    Ok(())
}

/// The per-step error bound for libroadrunner's absolute epsilon on a state of 2-norm `ynorm`
/// (amounts): epsilon itself where it means something (identical steps to libroadrunner), tightened
/// to epsilon relative where the state is tiny (epsilon would allow 100% errors) and loosened to what
/// rounding can reach where it is huge (else steps shrink to minimum_time_step and never end).
fn rk45_tol(eps: f64, ynorm: f64) -> f64 {
    let mut tol = eps;
    if ynorm > 0.0 { tol = tol.min(eps * ynorm); }
    tol.max(64.0 * f64::EPSILON * ynorm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Rules;

    /// S -> ; k*S, S(0) = 1, k = 1: S(t) = exp(-t).
    fn decay() -> Model {
        Model::new(vec!["S".into()], vec![1.0], vec![1.0], vec!["J".into()], vec!["k*S".into()], vec![-1.0],
                   vec!["k".into()], vec![1.0], Rules::default()).unwrap()
    }

    fn err_at_1(method: Method, step: f64, atol: f64) -> f64 {
        let m = decay();
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let mut out = [0.0; 2];
        integrate_fixed(&m, &mut w, method, &[0.0, 1.0], step, 0.0, atol, &[0], &mut out).unwrap();
        assert_eq!(out[0], 1.0);
        (out[1] - (-1f64).exp()).abs()
    }

    #[test]
    fn convergence_orders() {
        let (e1, e2) = (err_at_1(Method::Euler, 0.01, 0.0), err_at_1(Method::Euler, 0.005, 0.0));
        assert!((e1 / e2 - 2.0).abs() < 0.05, "euler ratio {}", e1 / e2);
        let (e1, e2) = (err_at_1(Method::Rk4, 0.1, 0.0), err_at_1(Method::Rk4, 0.05, 0.0));
        assert!((e1 / e2 - 16.0).abs() < 1.0, "rk4 ratio {}", e1 / e2);
        // one Euler step of 1: 1 - 1 = 0; one RK4 step: 1 - 1 + 1/2 - 1/6 + 1/24
        assert!((err_at_1(Method::Euler, 0.0, 0.0) - (-1f64).exp()).abs() < 1e-15);
        assert!((err_at_1(Method::Rk4, 0.0, 0.0) - (0.375 - (-1f64).exp()).abs()).abs() < 1e-15);
    }

    #[test]
    fn rk45_tolerance() {
        let (loose, tight) = (err_at_1(Method::Rk45, 0.0, 1e-6), err_at_1(Method::Rk45, 0.0, 1e-12));
        assert!(loose < 1e-5 && tight < 1e-10 && tight < loose, "{loose} {tight}");
    }
}

//! libroadrunner's "euler", "rk4" and "rk45" (Runge-Kutta-Fehlberg) integrators, step for step as
//! roadrunner 2.11 (EulerIntegrator.h, RK4Integrator.cpp, RK45Integrator.cpp) driven by
//! RoadRunner::simulate's fixed output grid (rk45 with variable_step_size = false).
//!
//! Event transitions are located inside steps. Algebraic rules and delay() are unsupported.
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
    if ev && m.apply_events(w, times[0]).is_none() { return Err("event cascade does not terminate".into()); }
    m.record(w, times[0], outs, &mut out[..no]);
    let z = || vec![0.0; nx];
    let (mut y, mut yt, mut k1, mut k2, mut k3, mut k4) = (w.s.clone(), z(), z(), z(), z(), z());
    // rk45: carries its step size from one output interval to the next
    let hmax = if step > 0.0 { step } else { 1.0 };
    let mut rk = Rk45::new(nx, 1e-12, hmax, atol);
    let axpy = |y: &mut [f64], a: f64, x: &[f64]| for (o, &v) in y.iter_mut().zip(x) { *o += a * v };
    for i in 1..nt {
        let (t0, t1) = (times[i - 1], times[i]);
        let dt = t1 - t0;
        if dt > 0.0 {
            let n = if step > 0.0 && step < dt { ((dt / step) * (1.0 - 1e-12)).ceil().max(1.0) as usize } else { 1 };
            let h = dt / n as f64;
            let mut t = t0;
            for j in 1..=if method == Method::Rk45 { 1 } else { n } {
                let end = if j == n || method == Method::Rk45 { t1 } else { t0 + j as f64 * h };
                while t < end {
                    let next = if ev { m.next_event_time(w, end) } else { end };
                    let hs = if ev { next - t } else { h };
                    let tn = match method {
                        Method::Euler => {
                            // Preserve EulerIntegrator's time convention for event-free intervals.
                            m.deriv(w, if ev { t } else { t0 }, &y, &mut k1);
                            axpy(&mut y, hs, &k1);
                            next
                        }
                        Method::Rk4 => {
                            m.deriv(w, t, &y, &mut k1);
                            let a = hs / 2.;
                            yt.copy_from_slice(&y); axpy(&mut yt, a, &k1);
                            m.deriv(w, t + a, &yt, &mut k2);
                            yt.copy_from_slice(&y); axpy(&mut yt, a, &k2);
                            m.deriv(w, t + a, &yt, &mut k3);
                            yt.copy_from_slice(&y); axpy(&mut yt, hs, &k3);
                            m.deriv(w, t + hs, &yt, &mut k4);
                            axpy(&mut k1, 2., &k2);
                            axpy(&mut k1, 2., &k3);
                            axpy(&mut k1, 1., &k4);
                            axpy(&mut y, hs / 6., &k1);
                            next
                        }
                        Method::Rk45 => rk.integrate(m, w, &mut y, t, next, ev)?,
                    };
                    if tn <= t { return Err("fixed integrator cannot advance time".into()); }
                    t = finish_step(m, w, method, t, tn, &mut y)?;
                }
            }
        }
        m.record(w, t1, outs, &mut out[i * no..(i + 1) * no]);
    }
    Ok(())
}

/// Locate a trigger transition on the accepted step's interpolant, then process pending
/// assignments with the same event queue as CVODE. The caller resumes from the post-event state.
fn finish_step(m: &Model, w: &mut Work, method: Method, t: f64, end: f64, y: &mut [f64]) -> Result<f64, String> {
    let mut te = end;
    if m.has_events() && m.trigger_changed(w, y, end) {
        let h = end - t;
        let y0 = w.s.clone();
        let (mut d0, mut d1, mut mid) = (vec![0.0; y.len()], vec![0.0; y.len()], y.to_vec());
        if method != Method::Euler {
            m.deriv(w, t, &y0, &mut d0);
            m.deriv(w, end, y, &mut d1);
        }
        let (mut lo, mut hi) = (t, end);
        while hi - lo > 1e-12 * end.abs().max(h) {
            let tm = lo + (hi - lo) / 2.0;
            if tm == lo || tm == hi { break; }
            let u = (tm - t) / h;
            for i in 0..y.len() {
                let delta = y[i] - y0[i];
                mid[i] = if method == Method::Euler { y0[i] + u * delta }
                    else { y0[i] + u * (delta + (1.0 - u) * ((1.0 - u) * (h * d0[i] - delta) + u * (delta - h * d1[i]))) };
            }
            if m.trigger_changed(w, &mid, tm) {
                hi = tm;
                w.s.copy_from_slice(&mid);
            } else { lo = tm; }
        }
        te = hi;
    }
    if te == end { w.s.copy_from_slice(y); }
    if m.has_events() {
        m.apply_events(w, te).ok_or("event cascade does not terminate")?;
        y.copy_from_slice(&w.s);
    }
    Ok(te)
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

/// libroadrunner's RK45Integrator::integrate(t, t_end - t): error-controlled Runge-Kutta-Fehlberg
/// steps from (t, y), the last one cut to end exactly at t_end; `one`: return after one accepted
/// step (variable_step_size), else after reaching t_end. The step size proposal carries over calls.
struct Rk45 { hmin: f64, hmax: f64, eps: f64, hcur: f64, budget: usize, k: [Vec<f64>; 6], yt: Vec<f64> }

impl Rk45 {
    fn new(nx: usize, hmin: f64, hmax: f64, eps: f64) -> Rk45 {
        let z = || vec![0.0; nx];
        Rk45 { hmin: hmin.min(hmax), hmax, eps, hcur: 0.0, budget: 100_000_000, k: [z(), z(), z(), z(), z(), z()], yt: z() }
    }

    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn integrate(&mut self, m: &Model, w: &mut Work, y: &mut [f64], t: f64, t1: f64, one: bool) -> Result<f64, String> {
        if !(t1 > t) { return Ok(t); }
        // fused multiply-add, as the BLAS daxpy libroadrunner calls (the error sum cancels: its rounding
        // steers the step size)
        let axpy = |y: &mut [f64], a: f64, x: &[f64]| for (o, &v) in y.iter_mut().zip(x) { *o = a.mul_add(v, *o) };
        let (hmin, hmax, eps) = (self.hmin, self.hmax, self.eps);
        let mut hp = (if self.hcur > 0.0 { self.hcur } else { hmax }).max(hmin).min(hmax);
        let mut tc = t;
        loop {
            self.budget = self.budget.checked_sub(1).ok_or("rk45: too many steps")?;
            let [k1, k2, k3, k4, k5, k6] = &mut self.k;
            let yt = &mut self.yt;
            let htry = hp;
            let (mut h, mut reach) = (htry, false);
            if h >= t1 - tc { h = t1 - tc; reach = true; }
            m.deriv(w, tc, y, k1);
            yt.copy_from_slice(y); axpy(yt, h / 4., k1);
            m.deriv(w, tc + h / 4., yt, k2);
            yt.copy_from_slice(y); axpy(yt, 3. * h / 32., k1); axpy(yt, 9. * h / 32., k2);
            m.deriv(w, tc + 3. * h / 8., yt, k3);
            yt.copy_from_slice(y);
            axpy(yt, 1932. * h / 2197., k1); axpy(yt, -7200. * h / 2197., k2); axpy(yt, 7296. * h / 2197., k3);
            m.deriv(w, tc + 12. * h / 13., yt, k4);
            yt.copy_from_slice(y);
            axpy(yt, 439. * h / 216., k1); axpy(yt, -8. * h, k2);
            axpy(yt, 3680. * h / 513., k3); axpy(yt, -845. * h / 4104., k4);
            m.deriv(w, tc + h, yt, k5);
            yt.copy_from_slice(y);
            axpy(yt, -8. * h / 27., k1); axpy(yt, 2. * h, k2); axpy(yt, -3544. * h / 2565., k3);
            axpy(yt, 1859. * h / 4104., k4); axpy(yt, -11. * h / 40., k5);
            m.deriv(w, tc + h / 2., yt, k6);
            // libroadrunner's error: |k1/360 - 128/4275 k3 - 2197/75240 k4 + k5/50 + 2/55 k6|_2
            // (not scaled by h), over its state vector of amounts
            let error = (0..y.len()).map(|j| {
                // daxpy order and fma, as libroadrunner
                let mut e = 0f64;
                for (a, kk) in [(1f64 / 360., &k1), (-128. / 4275., &k3), (-2197. / 75240., &k4), (1. / 50., &k5), (2. / 55., &k6)] { e = a.mul_add(kk[j] * m.vol[j], e); }
                e * e
            }).sum::<f64>().sqrt();
            let tol = rk45_tol(eps, (0..y.len()).map(|j| (y[j] * m.vol[j]).powi(2)).sum::<f64>().sqrt());
            let q = if !error.is_finite() { 0.1 } else if error == 0. { 4. } else { (0.84 * (tol / error).powf(0.25)).clamp(0.1, 4.) };
            if !(error <= tol) && htry > hmin { hp = (q * h).max(hmin); continue; }
            axpy(y, 1408. * h / 2565., k3);
            axpy(y, 2197. / 4104. * h, k4);
            axpy(y, -1. / 5. * h, k5);
            axpy(y, 25. / 216. * h, k1);
            tc = if reach { t1 } else { tc + h };
            if tc >= t1 { tc = t1; reach = true; }
            hp = q * h;
            if h < htry && hp < htry { hp = htry; } // a step cut to end at t1 does not lower it
            hp = hp.max(hmin).min(hmax);
            self.hcur = hp;
            if one || reach { return Ok(tc); }
        }
    }
}

/// libroadrunner's rk45 with variable_step_size = true (RoadRunner::simulate's variable step branch):
/// a row at t0, then one row per accepted step, the last step cut to end exactly at t1; at most
/// `max_rows` rows (0: no limit). `hmax` <= 0: 1. Events as `integrate_fixed`. Returns (times, rows
/// of outs, row-major). On success w.s / params hold the final state.
#[allow(clippy::too_many_arguments)]
pub fn integrate_rk45_variable(m: &Model, w: &mut Work, t0: f64, t1: f64, hmin: f64, hmax: f64, eps: f64,
                               outs: &[u32], max_rows: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    if m.has_algebraic() { return Err("rk45: algebraic rules are not supported (libroadrunner cannot load them either)".into()); }
    if m.has_delays() { return Err("rk45: delay() is not supported".into()); }
    if !(t1 > t0) { return Err("rk45 (variable step): need end > start".into()); }
    if m.has_events() && m.apply_events(w, t0).is_none() { return Err("event cascade does not terminate".into()); }
    let mut row = vec![0.0; outs.len()];
    m.record(w, t0, outs, &mut row);
    let (mut ts, mut rows) = (vec![t0], row.clone());
    let mut rk = Rk45::new(m.nx, hmin, if hmax > 0.0 { hmax } else { 1.0 }, eps);
    let mut y = w.s.clone();
    let mut t = t0;
    while max_rows == 0 || ts.len() < max_rows {
        let end = rk.integrate(m, w, &mut y, t, m.next_event_time(w, t1), true)?;
        if end <= t { return Err("rk45 cannot advance time".into()); }
        t = finish_step(m, w, Method::Rk45, t, end, &mut y)?;
        m.record(w, t, outs, &mut row);
        ts.push(t);
        rows.extend_from_slice(&row);
        if t >= t1 { break; }
    }
    Ok((ts, rows))
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

    #[test]
    fn rk45_variable_steps() {
        let m = decay();
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        let (t, y) = integrate_rk45_variable(&m, &mut w, 0.0, 3.0, 1e-12, 1.0, 1e-12, &[0], 0).unwrap();
        assert!(t.len() > 4 && t[0] == 0.0 && *t.last().unwrap() == 3.0 && t.windows(2).all(|p| p[1] > p[0]));
        for (ti, yi) in t.iter().zip(&y) { assert!((yi - (-ti).exp()).abs() < 1e-9, "{ti} {yi}"); }
        let mut w = m.work();
        m.load(&mut w, &[], &[]);
        assert_eq!(integrate_rk45_variable(&m, &mut w, 0.0, 3.0, 1e-12, 1.0, 1e-12, &[0], 3).unwrap().0.len(), 3);
    }
}

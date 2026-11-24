//! Random draws for SBML distrib functions (normal(...), poisson(...), ...): each evaluation draws
//! afresh from a per-thread stream that a run seeds (`seed`), so batches are reproducible per sample.
use crate::model::Rng;
use std::cell::RefCell;

thread_local!(static RNG: RefCell<Rng> = RefCell::new(Rng::new(0x5eed_d157)));

/// Start this thread's stream (at the start of each run).
pub fn seed(s: u64) { RNG.with(|r| *r.borrow_mut() = Rng::new(s)) }

/// Distribution kinds, in the order of `NAMES`, with their number of parameters.
pub const NAMES: [(&str, usize); 12] = [("normal", 2), ("uniform", 2), ("exponential", 1), ("gamma", 2), ("poisson", 1),
    ("binomial", 2), ("bernoulli", 1), ("cauchy", 2), ("chisquare", 1), ("laplace", 2), ("lognormal", 2), ("rayleigh", 1)];

fn normal(u: &mut impl FnMut() -> f64) -> f64 { // Box-Muller
    (-2.0 * u().ln()).sqrt() * (2.0 * std::f64::consts::PI * u()).cos()
}

fn gamma(k: f64, u: &mut impl FnMut() -> f64) -> f64 { // Marsaglia-Tsang, unit scale
    if k < 1.0 { return gamma(k + 1.0, u) * u().powf(1.0 / k); }
    let d = k - 1.0 / 3.0;
    let c = 1.0 / (9.0 * d).sqrt();
    loop {
        let x = normal(u);
        let v = (1.0 + c * x).powi(3);
        if v > 0.0 && u().ln() < 0.5 * x * x + d - d * v + d * v.ln() { return d * v; }
    }
}

fn poisson(l: f64, u: &mut impl FnMut() -> f64) -> f64 { // Knuth on chunks of <= 30 (sums of Poissons are Poisson)
    let (mut l, mut n) = (l, 0.0);
    while l > 0.0 {
        let c = l.min(30.0);
        l -= c;
        let (lim, mut p) = ((-c).exp(), u());
        while p > lim { n += 1.0; p *= u(); }
    }
    n
}

fn once(k: u8, a: f64, b: f64, u: &mut impl FnMut() -> f64) -> f64 {
    match k {
        0 => a + b * normal(u),
        1 => a + (b - a) * u(),
        2 => -u().ln() / a,
        3 => gamma(a, u) * b,
        4 => poisson(a, u),
        5 => (0..a.round().max(0.0) as u64).filter(|_| u() <= b).count() as f64,
        6 => if u() <= a { 1.0 } else { 0.0 },
        7 => a + b * (std::f64::consts::PI * (u() - 0.5)).tan(),
        8 => 2.0 * gamma(0.5 * a, u),
        9 => { let v = u() - 0.5; a - b * v.signum() * (1.0 - 2.0 * v.abs()).max(f64::MIN_POSITIVE).ln() }
        10 => (a + b * normal(u)).exp(),
        11 => a * (-2.0 * u().ln()).sqrt(),
        _ => f64::NAN,
    }
}

/// One draw of kind k with parameters (a, b), truncated to [lo, hi] by rejection (NaN if the
/// interval is (practically) never hit).
pub fn draw(k: u8, a: f64, b: f64, lo: f64, hi: f64) -> f64 {
    RNG.with(|r| {
        let mut r = r.borrow_mut();
        let mut u = || r.unif();
        for _ in 0..1_000_000 {
            let x = once(k, a, b, &mut u);
            if !(x < lo || x > hi) { return x; }
        }
        f64::NAN
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stats(k: u8, a: f64, b: f64, lo: f64, hi: f64) -> (f64, f64) {
        seed(7);
        let n = 200_000;
        let xs: Vec<f64> = (0..n).map(|_| draw(k, a, b, lo, hi)).collect();
        let m = xs.iter().sum::<f64>() / n as f64;
        (m, (xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt())
    }

    #[test]
    fn moments() {
        let inf = f64::INFINITY;
        // (kind, a, b, mean, sd)
        for (k, a, b, mu, sd) in [(0, 1.0, 2.0, 1.0, 2.0), (1, 0.0, 1.0, 0.5, (1.0f64 / 12.0).sqrt()), (2, 2.0, 0.0, 0.5, 0.5),
                                  (3, 2.0, 3.0, 6.0, 18f64.sqrt()), (4, 1.5, 0.0, 1.5, 1.5f64.sqrt()), (4, 45.0, 0.0, 45.0, 45f64.sqrt()),
                                  (5, 100.0, 0.2, 20.0, 16f64.sqrt()), (6, 0.23, 0.0, 0.23, (0.23f64 * 0.77).sqrt()),
                                  (8, 3.0, 0.0, 3.0, 6f64.sqrt()), (9, 1.0, 2.0, 1.0, 8f64.sqrt()),
                                  (11, 2.0, 0.0, 2.0 * (std::f64::consts::PI / 2.0).sqrt(), 2.0 * ((4.0 - std::f64::consts::PI) / 2.0).sqrt())] {
            let (m, s) = stats(k, a, b, -inf, inf);
            assert!((m - mu).abs() < 5.0 * sd / 447.0, "kind {k}: mean {m} vs {mu}");
            assert!((s - sd).abs() < 0.02 * sd, "kind {k}: sd {s} vs {sd}");
        }
        let (m, _) = stats(0, 0.0, 1.0, 0.0, inf); // half-normal: mean sqrt(2/pi)
        assert!((m - (2.0 / std::f64::consts::PI).sqrt()).abs() < 0.01);
    }
}

/// The 32-bit MT19937 stream used by libroadrunner's Gillespie integrator.
pub(crate) struct Mt19937 {
    state: Box<[u32; 624]>,
    next: usize,
}

impl Mt19937 {
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = Box::new([0u32; 624]);
        state[0] = seed;
        for i in 1..state.len() {
            let prev = state[i - 1];
            state[i] = 1812433253u32.wrapping_mul(prev ^ (prev >> 30)).wrapping_add(i as u32);
        }
        Self { next: state.len(), state }
    }

    pub(crate) fn uniform(&mut self) -> f64 {
        let n = self.state.len();
        if self.next == n {
            for i in 0..n {
                let bits = (self.state[i] & 0x80000000) | (self.state[(i + 1) % n] & 0x7fffffff);
                self.state[i] = self.state[(i + 397) % n] ^ (bits >> 1)
                    ^ if bits & 1 == 0 { 0 } else { 0x9908b0df };
            }
            self.next = 0;
        }
        let mut value = self.state[self.next];
        self.next += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c5680;
        value ^= (value << 15) & 0xefc60000;
        value ^= value >> 18;
        f64::from(value.max(1)) / f64::from(u32::MAX)
    }
}

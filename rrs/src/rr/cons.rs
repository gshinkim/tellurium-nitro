//! Conserved moieties as libroadrunner computes them: a port of libstruct's analyzeWithQR
//! (rr-libstruct lsLibStructural.cpp): QR with column pivoting (LAPACK dgeqp3, unblocked dlaqp2) of
//! N^T, R rounded to 1e-12, Gauss-Jordan at 1e-9, rank from the zero rows. The pivot order is the
//! species order roadrunner reports (independent species first), and L0 comes from the reduced R.

/// Structural analysis of a stoichiometry matrix (species x reactions).
#[derive(Debug, Clone, PartialEq)]
pub struct Ls {
    /// Species (row indices of N) in libstruct's order: independent first, then dependent.
    pub order: Vec<usize>,
    /// Number of independent species.
    pub ni: usize,
    /// L0: dependent x independent (row-major), dependent amounts = L0 * independent amounts + totals.
    pub l0: Vec<f64>,
}

impl Ls {
    pub fn nd(&self) -> usize { self.order.len() - self.ni }
}

pub(crate) fn round_tol(x: f64, tol: f64) -> f64 {
    if x.abs() < tol { 0.0 } else if (x.ceil() - x).abs() < tol { x.ceil() } else if (x - x.floor()).abs() < tol { x.floor() } else { x }
}

/// a * b + c as roadrunner's LAPACK computes it: its macOS arm64 build (clang, fp-contract=on) fuses
/// these into FMAs, which decides near-ties between pivot columns.
// ponytail: fused on aarch64 only; matches roadrunner's x86 builds only where those don't contract
#[inline]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 { if cfg!(target_arch = "aarch64") { a.mul_add(b, c) } else { a * b + c } }

/// Reference BLAS dnrm2 (scaled sum of squares).
fn dnrm2(x: impl Iterator<Item = f64>) -> f64 {
    let (mut scale, mut ssq) = (0.0f64, 1.0f64);
    for v in x {
        if v != 0.0 {
            let a = v.abs();
            if scale < a { let d = scale / a; ssq = fma(ssq, d * d, 1.0); scale = a; } else { let d = a / scale; ssq = fma(d, d, ssq); }
        }
    }
    scale * ssq.sqrt()
}

fn dlapy2(x: f64, y: f64) -> f64 {
    let (w, z) = (x.abs().max(y.abs()), x.abs().min(y.abs()));
    if z == 0.0 { w } else { w * fma(z / w, z / w, 1.0).sqrt() }
}

/// N: rows x cols (species x reactions), row-major.
pub fn analyze(n: &[f64], rows: usize, cols: usize) -> Ls {
    if rows == 0 || cols == 0 { // libstruct: "no floating species / no reactions": nothing to reduce
        return Ls { order: (0..rows).collect(), ni: rows, l0: vec![] };
    }
    // A = N^T (m reactions x k species), row-major
    let (m, k) = (cols, rows);
    let mut a: Vec<f64> = (0..m * k).map(|t| n[(t % k) * cols + t / k]).collect();
    let at = |a: &[f64], i: usize, j: usize| a[i * k + j];
    let norm = |a: &[f64], j: usize, from: usize| dnrm2((from..m).map(|i| at(a, i, j)));
    let mut jpvt: Vec<usize> = (0..k).collect();
    let mut vn1: Vec<f64> = (0..k).map(|j| norm(&a, j, 0)).collect();
    let mut vn2 = vn1.clone();
    let tol3z = (f64::EPSILON * 0.5).sqrt();
    for i in 0..m.min(k) {
        let mut pvt = i;
        for j in i + 1..k { if vn1[j] > vn1[pvt] { pvt = j; } }
        if pvt != i {
            for r in 0..m { a.swap(r * k + pvt, r * k + i); }
            jpvt.swap(pvt, i);
            vn1[pvt] = vn1[i];
            vn2[pvt] = vn2[i];
        }
        // Householder reflector for A[i.., i] with a non-negative diagonal (dlarfp)
        let mut alpha = at(&a, i, i);
        let xnorm = if i + 1 < m { norm(&a, i, i + 1) } else { 0.0 };
        let tau;
        if xnorm == 0.0 {
            tau = if alpha >= 0.0 { 0.0 } else { 2.0 };
            if alpha < 0.0 { for r in i + 1..m { a[r * k + i] = 0.0; } a[i * k + i] = -alpha; }
        } else {
            let d = dlapy2(alpha, xnorm);
            let mut beta = if alpha >= 0.0 { d } else { -d };
            alpha += beta;
            if beta < 0.0 {
                beta = -beta;
                tau = -alpha / beta;
            } else {
                alpha = xnorm * (xnorm / alpha);
                tau = alpha / beta;
                alpha = -alpha;
            }
            let sc = 1.0 / alpha;
            for r in i + 1..m { a[r * k + i] *= sc; }
            a[i * k + i] = beta;
        }
        // apply H = I - tau v v^T (v = [1; A[i+1.., i]]) to A[i.., i+1..] (dlarf: dgemv, then dger)
        if tau != 0.0 {
            for j in i + 1..k {
                let mut w = 0.0;
                w += at(&a, i, j);
                for r in i + 1..m { w = fma(at(&a, r, j), at(&a, r, i), w); }
                let t = -tau * w;
                a[i * k + j] += t;
                for r in i + 1..m { a[r * k + j] = fma(at(&a, r, i), t, a[r * k + j]); }
            }
        }
        // downdate the partial column norms
        for j in i + 1..k {
            if vn1[j] == 0.0 { continue; }
            let d = at(&a, i, j).abs() / vn1[j];
            let temp = fma(-d, d, 1.0).max(0.0);
            let q = vn1[j] / vn2[j];
            let temp2 = temp * (q * q);
            if temp2 <= tol3z {
                vn1[j] = if i + 1 < m { norm(&a, j, i + 1) } else { 0.0 };
                vn2[j] = vn1[j];
            } else {
                vn1[j] *= temp.sqrt();
            }
        }
    }
    // R: upper trapezoid, rounded as libstruct
    let mut r: Vec<f64> = (0..m * k).map(|t| if t % k < t / k { 0.0 } else { round_tol(a[t], 1e-12) }).collect();
    gauss_jordan(&mut r, m, k, 1e-9);
    // rank: zero rows at the bottom (row 0 never counts as zero, as ls::findRank)
    let mut ni = m;
    let mut i = m - 1;
    while i != 0 && (0..k).map(|j| r[i * k + j].abs()).sum::<f64>() < 1e-9 { ni -= 1; i -= 1; }
    let ni = ni.min(k);
    let nd = k - ni;
    let mut l0 = vec![0.0; nd * ni];
    for d in 0..nd { for b in 0..ni { l0[d * ni + b] = r[b * k + ni + d]; } }
    Ls { order: jpvt, ni, l0 }
}

/// ls::gaussJordan: partial pivoting, eliminates above the pivots only, then rounds.
fn gauss_jordan(a: &mut [f64], rows: usize, cols: usize, tol: f64) {
    let (mut pr, mut pc) = (0, 0);
    while pr < rows && pc < cols {
        let mut x = pr;
        for c in pr..rows { if a[c * cols + pc].abs() > a[x * cols + pc].abs() { x = c; } }
        if x != pr { for b in 0..cols { a.swap(pr * cols + b, x * cols + b); } }
        let p = a[pr * cols + pc];
        if p.abs() > tol {
            for b in 0..cols { a[pr * cols + b] /= p; }
            for b in 0..pr {
                let f = a[b * cols + pc];
                for d in pr..cols { a[b * cols + d] = fma(-a[pr * cols + d], f, a[b * cols + d]); }
            }
            pr += 1;
        }
        pc += 1;
    }
    for x in a.iter_mut() { *x = round_tol(*x, tol); }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn libstruct_order_and_l0() {
        // S1, S2, S3, A, B x J1..J5 (roadrunner: S1, S2, A independent; B, S3 dependent)
        let n = [-1., 0., 1., -1., 1., 1., -1., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 0., -1., 1., 0., 0., 0., 1., -1.];
        let ls = analyze(&n, 5, 5);
        assert_eq!(ls.order, vec![0, 1, 3, 4, 2]);
        assert_eq!(ls.ni, 3);
        assert_eq!(ls.l0, vec![0., 0., -1., -1., -1., 1.]);
    }
}

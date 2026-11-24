//! Eigenvalues as libroadrunner computes them: rr-libstruct ls::getEigenValues (lsLibla.cpp) copies
//! the (real) Jacobian M into a complex column-major array as A[j + n*i] = M(i, j), i.e. A = M^T,
//! calls CLAPACK zgeev('N', 'N') with lwork = 2n and rounds each part to 1e-12 (ls::RoundToTolerance).
//! zgeev = zlascl (extreme norms only) + zgebal('B') + zgehrd (its unblocked zgehd2 at this lwork) +
//! zhseqr('E', 'N') (zlahqr up to n = 75, the multishift zlaqr0 above: eig_aed.rs). Complex
//! arithmetic as f2c writes it, contracted into FMAs like roadrunner's arm64 build (see cons::fma).
//!
//! Matrices live in one column-major buffer `h: &mut [C]`; a view `V` is LAPACK's (pointer, ld) pair
//! and indexes 1-based like the Fortran. Vectors (w, tau, work) are separate slices, 0-based.

use super::cons::{fma, round_tol};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub(crate) struct C { pub r: f64, pub i: f64 }

pub(crate) const Z: C = C { r: 0.0, i: 0.0 };
pub(crate) const ONE: C = C { r: 1.0, i: 0.0 };
pub(crate) fn c(r: f64, i: f64) -> C { C { r, i } }
/// f2c's a*b: (ar*br - ai*bi, ar*bi + ai*br), the left product fused.
#[inline] pub(crate) fn mul(a: C, b: C) -> C { c(fma(a.r, b.r, -(a.i * b.i)), fma(a.r, b.i, a.i * b.r)) }
#[inline] pub(crate) fn add(a: C, b: C) -> C { c(a.r + b.r, a.i + b.i) }
#[inline] pub(crate) fn sub(a: C, b: C) -> C { c(a.r - b.r, a.i - b.i) }
#[inline] pub(crate) fn neg(a: C) -> C { c(-a.r, -a.i) }
#[inline] pub(crate) fn cj(a: C) -> C { c(a.r, -a.i) }
#[inline] pub(crate) fn scl(s: f64, a: C) -> C { c(s * a.r, s * a.i) }
#[inline] pub(crate) fn div_r(a: C, s: f64) -> C { c(a.r / s, a.i / s) }
#[inline] pub(crate) fn abs1(a: C) -> f64 { a.r.abs() + a.i.abs() }
#[inline] pub(crate) fn nz(a: C) -> bool { a.r != 0.0 || a.i != 0.0 }

/// A column-major matrix inside a buffer: element (i, j) (1-based) at o + (i-1) + (j-1)*ld.
#[derive(Debug, Clone, Copy)]
pub(crate) struct V { pub o: usize, pub ld: usize }
impl V {
    #[inline] pub(crate) fn i(self, r: usize, c: usize) -> usize { self.o + (r - 1) + (c - 1) * self.ld }
    #[inline] pub(crate) fn at(self, r: usize, c: usize) -> V { V { o: self.i(r, c), ld: self.ld } }
}

/// libf2c f__cabs.
pub(crate) fn cabs(re: f64, im: f64) -> f64 {
    let (mut r, mut i) = (re.abs(), im.abs());
    if i > r { std::mem::swap(&mut r, &mut i); }
    if r + i == r { return r; }
    let t = i / r;
    r * fma(t, t, 1.0).sqrt()
}
pub(crate) fn zabs(a: C) -> f64 { cabs(a.r, a.i) }

/// libf2c z_sqrt.
pub(crate) fn zsqrt(z: C) -> C {
    let mag = cabs(z.r, z.i);
    if mag == 0.0 { Z } else if z.r > 0.0 {
        let r = (0.5 * (mag + z.r)).sqrt();
        c(r, z.i / r / 2.0)
    } else {
        let mut i = (0.5 * (mag - z.r)).sqrt();
        if z.i < 0.0 { i = -i; }
        c(z.i / i / 2.0, i)
    }
}

/// libf2c pow_zi(x, 2).
pub(crate) fn sq(x: C) -> C { c(fma(x.r, x.r, -(x.i * x.i)), 2.0 * x.r * x.i) }

/// LAPACK 3.2 dladiv via zladiv: x / y.
pub(crate) fn zladiv(x: C, y: C) -> C {
    let (a, b, cc, d) = (x.r, x.i, y.r, y.i);
    if d.abs() < cc.abs() {
        let e = d / cc;
        let f = fma(d, e, cc);
        c(fma(b, e, a) / f, fma(-a, e, b) / f)
    } else {
        let e = cc / d;
        let f = fma(cc, e, d);
        c(fma(a, e, b) / f, fma(b, e, -a) / f)
    }
}

pub(crate) fn dlapy2(x: f64, y: f64) -> f64 {
    let (w, z) = (x.abs().max(y.abs()), x.abs().min(y.abs()));
    if z == 0.0 { w } else { w * fma(z / w, z / w, 1.0).sqrt() }
}

fn dlapy3(x: f64, y: f64, z: f64) -> f64 {
    let (xa, ya, za) = (x.abs(), y.abs(), z.abs());
    let w = xa.max(ya).max(za);
    if w == 0.0 { return xa + ya + za; }
    let (a, b, d) = (xa / w, ya / w, za / w);
    w * fma(d, d, fma(a, a, b * b)).sqrt()
}

/// Reference BLAS dznrm2.
pub(crate) fn dznrm2(x: impl Iterator<Item = C>) -> f64 {
    let (mut scale, mut ssq) = (0.0f64, 1.0f64);
    let mut acc = |v: f64| {
        if v != 0.0 {
            let t = v.abs();
            if scale < t { let d = scale / t; ssq = fma(ssq, d * d, 1.0); scale = t; } else { let d = t / scale; ssq = fma(d, d, ssq); }
        }
    };
    for z in x { acc(z.r); acc(z.i); }
    scale * ssq.sqrt()
}

pub(crate) const SAFMIN: f64 = f64::MIN_POSITIVE; // dlamch('S')
pub(crate) const EPS: f64 = f64::EPSILON * 0.5; // dlamch('E')
pub(crate) const ULP: f64 = f64::EPSILON; // dlamch('P')

/// zlarfg: alpha, x (n-1 values) -> tau; alpha becomes beta, x becomes v(2:n).
pub(crate) fn zlarfg(alpha: &mut C, x: &mut [C]) -> C {
    let mut xnorm = dznrm2(x.iter().copied());
    let (mut alphr, mut alphi) = (alpha.r, alpha.i);
    if xnorm == 0.0 && alphi == 0.0 { return Z; }
    let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };
    let mut beta = -sign(dlapy3(alphr, alphi, xnorm), alphr);
    let safmin = SAFMIN / EPS;
    let rsafmn = 1.0 / safmin;
    let mut knt = 0;
    if beta.abs() < safmin {
        loop {
            knt += 1;
            for v in x.iter_mut() { *v = scl(rsafmn, *v); }
            beta *= rsafmn; alphi *= rsafmn; alphr *= rsafmn;
            if beta.abs() >= safmin { break; }
        }
        xnorm = dznrm2(x.iter().copied());
        *alpha = c(alphr, alphi);
        beta = -sign(dlapy3(alphr, alphi, xnorm), alphr);
    }
    let tau = c((beta - alphr) / beta, -alphi / beta);
    let s = zladiv(ONE, c(alpha.r - beta, alpha.i));
    for v in x.iter_mut() { *v = mul(s, *v); }
    for _ in 0..knt { beta *= safmin; }
    *alpha = c(beta, 0.0);
    tau
}

/// zlarf on the m x n block cv: H = I - tau v v^H from the left (`left`) or the right.
pub(crate) fn zlarf(left: bool, m: usize, n: usize, v: &[C], tau: C, h: &mut [C], cv: V) {
    if !nz(tau) { return; }
    let mut lastv = if left { m } else { n };
    while lastv > 0 && !nz(v[lastv - 1]) { lastv -= 1; }
    let ix = |i: usize, j: usize| cv.i(i, j);
    let alpha = neg(tau);
    if left {
        // ilazlc(lastv, n)
        let lastc = if n == 0 { 0 } else if nz(h[ix(1, n)]) || nz(h[ix(lastv, n)]) { n } else {
            let mut k = n;
            while k >= 1 && !(1..=lastv).any(|i| nz(h[ix(i, k)])) { k -= 1; }
            k
        };
        if lastc == 0 { return; }
        let mut w = vec![Z; lastc];
        for j in 1..=lastc {
            let mut t = Z;
            for i in 1..=lastv { t = add(t, mul(cj(h[ix(i, j)]), v[i - 1])); }
            w[j - 1] = add(Z, mul(ONE, t));
        }
        for j in 1..=lastc {
            if nz(w[j - 1]) {
                let t = mul(alpha, cj(w[j - 1]));
                for i in 1..=lastv { let k = ix(i, j); h[k] = add(h[k], mul(v[i - 1], t)); }
            }
        }
    } else {
        // ilazlr(m, lastv)
        let lastc = if m == 0 { 0 } else if nz(h[ix(m, 1)]) || nz(h[ix(m, lastv)]) { m } else {
            (1..=lastv).map(|j| { let mut i = m; while i >= 1 && !nz(h[ix(i, j)]) { i -= 1; } i }).max().unwrap_or(0)
        };
        if lastc == 0 { return; }
        let mut w = vec![Z; lastc];
        for j in 1..=lastv {
            if nz(v[j - 1]) {
                let t = mul(ONE, v[j - 1]);
                for i in 1..=lastc { w[i - 1] = add(w[i - 1], mul(t, h[ix(i, j)])); }
            }
        }
        for j in 1..=lastv {
            if nz(v[j - 1]) {
                let t = mul(alpha, cj(v[j - 1]));
                for i in 1..=lastc { let k = ix(i, j); h[k] = add(h[k], mul(w[i - 1], t)); }
            }
        }
    }
}

/// izamax over an iterator (1-based index of the first largest |re| + |im|).
fn izamax(x: impl Iterator<Item = C>) -> usize {
    let (mut k, mut best) = (0, -1.0);
    for (i, z) in x.enumerate() { if i == 0 || abs1(z) > best { k = i + 1; best = abs1(z); } }
    k
}

/// zgebal('B'): permute and scale the n x n matrix at a; returns (ilo, ihi).
fn zgebal(n: usize, h: &mut [C], a: V) -> (usize, usize) {
    let mut scale = vec![1.0f64; n + 1];
    let (mut k, mut l) = (1usize, n);
    let ix = |i: usize, j: usize| a.i(i, j);
    let exc = |h: &mut [C], scale: &mut [f64], j: usize, m: usize, k: usize, l: usize| {
        scale[m] = j as f64;
        if j != m {
            for i in 1..=l { h.swap(ix(i, j), ix(i, m)); }
            for i in k..=n { h.swap(ix(j, i), ix(m, i)); }
        }
    };
    // rows isolating an eigenvalue: push down
    'rows: loop {
        for j in (1..=l).rev() {
            if (1..=l).all(|i| i == j || !nz(h[ix(j, i)])) {
                exc(h, &mut scale, j, l, k, l);
                if l == 1 { return (k, l); }
                l -= 1;
                continue 'rows;
            }
        }
        break;
    }
    // columns isolating an eigenvalue: push left
    'cols: loop {
        for j in k..=l {
            if (k..=l).all(|i| i == j || !nz(h[ix(i, j)])) {
                exc(h, &mut scale, j, k, k, l);
                k += 1;
                continue 'cols;
            }
        }
        break;
    }
    for s in &mut scale[k..=l] { *s = 1.0; }
    let sfmin1 = SAFMIN / ULP;
    let sfmax1 = 1.0 / sfmin1;
    let sfmin2 = sfmin1 * 2.0;
    let sfmax2 = 1.0 / sfmin2;
    loop {
        let mut noconv = false;
        for i in k..=l {
            let (mut cc, mut r) = (0.0f64, 0.0f64);
            for j in k..=l {
                if j == i { continue; }
                cc += h[ix(j, i)].r.abs() + h[ix(j, i)].i.abs();
                r += h[ix(i, j)].r.abs() + h[ix(i, j)].i.abs();
            }
            let ica = izamax((1..=l).map(|q| h[ix(q, i)]));
            let mut ca = zabs(h[ix(ica, i)]);
            let ira = izamax((k..=n).map(|q| h[ix(i, q)]));
            let mut ra = zabs(h[ix(i, ira + k - 1)]);
            if cc == 0.0 || r == 0.0 { continue; }
            let mut g = r / 2.0;
            let mut f = 1.0f64;
            let s = cc + r;
            while !(cc >= g || f.max(cc).max(ca) >= sfmax2 || r.min(g).min(ra) <= sfmin2) {
                f *= 2.0; cc *= 2.0; ca *= 2.0; r /= 2.0; g /= 2.0; ra /= 2.0;
            }
            g = cc / 2.0;
            while !(g < r || r.max(ra) >= sfmax2 || f.min(cc).min(g).min(ca) <= sfmin2) {
                f /= 2.0; cc /= 2.0; g /= 2.0; ca /= 2.0; r *= 2.0; ra *= 2.0;
            }
            if cc + r >= s * 0.95 { continue; }
            if f < 1.0 && scale[i] < 1.0 && f * scale[i] <= sfmin1 { continue; }
            if f > 1.0 && scale[i] > 1.0 && scale[i] >= sfmax1 / f { continue; }
            let g = 1.0 / f;
            scale[i] *= f;
            noconv = true;
            for q in k..=n { h[ix(i, q)] = scl(g, h[ix(i, q)]); }
            for q in 1..=l { h[ix(q, i)] = scl(f, h[ix(q, i)]); }
        }
        if !noconv { break; }
    }
    (k, l)
}

/// zgehd2: Hessenberg reduction of rows/columns ilo..ihi of the n x n matrix a; tau[0..n-1].
pub(crate) fn zgehd2(n: usize, ilo: usize, ihi: usize, h: &mut [C], a: V, tau: &mut [C]) {
    for i in ilo..ihi {
        let mut alpha = h[a.i(i + 1, i)];
        let mut x: Vec<C> = (i + 2..=ihi).map(|q| h[a.i(q, i)]).collect();
        tau[i - 1] = zlarfg(&mut alpha, &mut x);
        for (q, v) in (i + 2..=ihi).zip(&x) { h[a.i(q, i)] = *v; }
        let v: Vec<C> = std::iter::once(ONE).chain(x).collect();
        zlarf(false, ihi, ihi - i, &v, tau[i - 1], h, a.at(1, i + 1));
        zlarf(true, ihi - i, n - i, &v, cj(tau[i - 1]), h, a.at(i + 1, i + 1));
        h[a.i(i + 1, i)] = alpha;
    }
}

/// zlahqr: eigenvalues (and Schur form with `wantt`, Schur vectors into z with `wantz`) of the
/// Hessenberg block ilo..ihi of the n x n matrix hv; w[i-1] = eigenvalue i. Returns LAPACK's info.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zlahqr(wantt: bool, wantz: bool, n: usize, ilo: usize, ihi: usize, h: &mut [C], hv: V, w: &mut [C],
                     iloz: usize, ihiz: usize, zv: V) -> usize {
    if n == 0 { return 0; }
    let ix = |i: usize, j: usize| hv.i(i, j);
    let zx = |i: usize, j: usize| zv.i(i, j);
    if ilo == ihi { w[ilo - 1] = h[ix(ilo, ilo)]; return 0; }
    let mut j = ilo;
    while j + 3 <= ihi { h[ix(j + 2, j)] = Z; h[ix(j + 3, j)] = Z; j += 1; }
    if ilo + 2 <= ihi { h[ix(ihi, ihi - 2)] = Z; }
    let (jlo, jhi) = if wantt { (1, n) } else { (ilo, ihi) };
    for i in ilo + 1..=ihi {
        if h[ix(i, i - 1)].i != 0.0 {
            let d = abs1(h[ix(i, i - 1)]);
            let mut sc = div_r(h[ix(i, i - 1)], d);
            sc = div_r(cj(sc), zabs(sc));
            h[ix(i, i - 1)] = c(zabs(h[ix(i, i - 1)]), 0.0);
            for q in i..=jhi { h[ix(i, q)] = mul(sc, h[ix(i, q)]); }
            for q in jlo..=jhi.min(i + 1) { h[ix(q, i)] = mul(cj(sc), h[ix(q, i)]); }
            if wantz { for q in iloz..=ihiz { h[zx(q, i)] = mul(cj(sc), h[zx(q, i)]); } }
        }
    }
    let nh = ihi - ilo + 1;
    let nz_ = ihiz as isize - iloz as isize + 1;
    let smlnum = SAFMIN * (nh as f64 / ULP);
    let (mut i1, mut i2) = (1, n);
    let mut i = ihi;
    while i >= ilo {
        let mut l = ilo;
        let mut converged = false;
        for its in 0..=30 {
            // look for a single small subdiagonal element
            let mut k = i;
            while k > l {
                if abs1(h[ix(k, k - 1)]) <= smlnum { break; }
                let mut tst = abs1(h[ix(k - 1, k - 1)]) + abs1(h[ix(k, k)]);
                if tst == 0.0 {
                    if k >= ilo + 2 { tst += h[ix(k - 1, k - 2)].r.abs(); }
                    if k < ihi { tst += h[ix(k + 1, k)].r.abs(); }
                }
                if h[ix(k, k - 1)].r.abs() <= ULP * tst {
                    let (p, q) = (abs1(h[ix(k, k - 1)]), abs1(h[ix(k - 1, k)]));
                    let (ab, ba) = (p.max(q), p.min(q));
                    let zz = sub(h[ix(k - 1, k - 1)], h[ix(k, k)]);
                    let (p, q) = (abs1(h[ix(k, k)]), abs1(zz));
                    let (aa, bb) = (p.max(q), p.min(q));
                    let s = aa + ab;
                    if ba * (ab / s) <= smlnum.max(ULP * (bb * (aa / s))) { break; }
                }
                k -= 1;
            }
            l = k;
            if l > ilo { h[ix(l, l - 1)] = Z; }
            if l >= i { converged = true; break; }
            if !wantt { i1 = l; i2 = i; }
            let mut t;
            if its == 10 {
                let s = h[ix(l + 1, l)].r.abs() * 0.75;
                t = c(s + h[ix(l, l)].r, h[ix(l, l)].i);
            } else if its == 20 {
                let s = h[ix(i, i - 1)].r.abs() * 0.75;
                t = c(s + h[ix(i, i)].r, h[ix(i, i)].i);
            } else {
                t = h[ix(i, i)];
                let u = mul(zsqrt(h[ix(i - 1, i)]), zsqrt(h[ix(i, i - 1)]));
                let mut s = abs1(u);
                if s != 0.0 {
                    let x = scl(0.5, sub(h[ix(i - 1, i - 1)], t));
                    let sx = abs1(x);
                    s = s.max(abs1(x));
                    let mut y = scl(s, zsqrt(add(sq(div_r(x, s)), sq(div_r(u, s)))));
                    if sx > 0.0 {
                        let xs = div_r(x, sx);
                        if fma(xs.r, y.r, xs.i * y.i) < 0.0 { y = neg(y); }
                    }
                    t = sub(t, mul(u, zladiv(u, add(x, y))));
                }
            }
            // look for two consecutive small subdiagonal elements
            let mut m = i - 1;
            let mut v = loop {
                let mm = m.max(l);
                let (h11, h22) = (h[ix(mm, mm)], h[ix(mm + 1, mm + 1)]);
                let h11s = sub(h11, t);
                let h21 = h[ix(mm + 1, mm)].r;
                let s = abs1(h11s) + h21.abs();
                let (h11s, h21) = (div_r(h11s, s), h21 / s);
                if m <= l { m = l; break [h11s, c(h21, 0.0)]; }
                let h10 = h[ix(m, m - 1)].r;
                if h10.abs() * h21.abs() <= ULP * (abs1(h11s) * (abs1(h11) + abs1(h22))) { break [h11s, c(h21, 0.0)]; }
                m -= 1;
            };
            // single-shift QR step
            for k in m..i {
                if k > m { v = [h[ix(k, k - 1)], h[ix(k + 1, k - 1)]]; }
                let (mut a0, mut x1) = (v[0], [v[1]]);
                let t1 = zlarfg(&mut a0, &mut x1);
                v = [a0, x1[0]];
                if k > m { h[ix(k, k - 1)] = v[0]; h[ix(k + 1, k - 1)] = Z; }
                let v2 = v[1];
                let t2 = fma(t1.r, v2.r, -(t1.i * v2.i));
                for j in k..=i2 {
                    let sum = add(mul(cj(t1), h[ix(k, j)]), scl(t2, h[ix(k + 1, j)]));
                    h[ix(k, j)] = sub(h[ix(k, j)], sum);
                    h[ix(k + 1, j)] = sub(h[ix(k + 1, j)], mul(sum, v2));
                }
                for j in i1..=(k + 2).min(i) {
                    let sum = add(mul(t1, h[ix(j, k)]), scl(t2, h[ix(j, k + 1)]));
                    h[ix(j, k)] = sub(h[ix(j, k)], sum);
                    h[ix(j, k + 1)] = sub(h[ix(j, k + 1)], mul(sum, cj(v2)));
                }
                if wantz {
                    for j in iloz..=ihiz {
                        let sum = add(mul(t1, h[zx(j, k)]), scl(t2, h[zx(j, k + 1)]));
                        h[zx(j, k)] = sub(h[zx(j, k)], sum);
                        h[zx(j, k + 1)] = sub(h[zx(j, k + 1)], mul(sum, cj(v2)));
                    }
                }
                if k == m && m > l {
                    let mut temp = c(1.0 - t1.r, 0.0 - t1.i);
                    temp = div_r(temp, zabs(temp));
                    h[ix(m + 1, m)] = mul(h[ix(m + 1, m)], cj(temp));
                    if m + 2 <= i { h[ix(m + 2, m + 1)] = mul(h[ix(m + 2, m + 1)], temp); }
                    for j in m..=i {
                        if j != m + 1 {
                            if i2 > j { for q in j + 1..=i2 { h[ix(j, q)] = mul(temp, h[ix(j, q)]); } }
                            for q in i1..j { h[ix(q, j)] = mul(cj(temp), h[ix(q, j)]); }
                            if wantz && nz_ > 0 { for q in iloz..=ihiz { h[zx(q, j)] = mul(cj(temp), h[zx(q, j)]); } }
                        }
                    }
                }
            }
            let temp = h[ix(i, i - 1)];
            if temp.i != 0.0 {
                let rtemp = zabs(temp);
                h[ix(i, i - 1)] = c(rtemp, 0.0);
                let temp = div_r(temp, rtemp);
                if i2 > i { for q in i + 1..=i2 { h[ix(i, q)] = mul(cj(temp), h[ix(i, q)]); } }
                for q in i1..i { h[ix(q, i)] = mul(temp, h[ix(q, i)]); }
                if wantz && nz_ > 0 { for q in iloz..=ihiz { h[zx(q, i)] = mul(temp, h[zx(q, i)]); } }
            }
        }
        if !converged { return i; }
        w[i - 1] = h[ix(i, i)];
        if l <= ilo { break; }
        i = l - 1;
    }
    0
}

/// zhseqr('E', 'N') on the balanced Hessenberg n x n matrix (ld n at offset 0).
fn zhseqr(n: usize, ilo: usize, ihi: usize, h: &mut [C], w: &mut [C]) -> Result<(), String> {
    let hv = V { o: 0, ld: n };
    for i in 1..ilo { w[i - 1] = h[hv.i(i, i)]; }
    for i in ihi + 1..=n { w[i - 1] = h[hv.i(i, i)]; }
    if ilo == ihi { w[ilo - 1] = h[hv.i(ilo, ilo)]; return Ok(()); }
    let mut work = vec![Z; 2 * n];
    let info = if n > 75 {
        super::eig_aed::zlaqr0(false, false, n, ilo, ihi, h, hv, w, ilo, ihi, hv, &mut work)
    } else {
        let info = zlahqr(false, false, n, ilo, ihi, h, hv, w, ilo, ihi, hv);
        if info == 0 { 0 } else if n >= 49 {
            super::eig_aed::zlaqr0(false, false, n, ilo, info, h, hv, w, ilo, ihi, hv, &mut work)
        } else {
            // zlaqr0 on a 49 x 49 copy padded with zeros
            let mut hl = vec![Z; 49 * 49];
            let lv = V { o: 0, ld: 49 };
            for j in 1..=n { for i in 1..=n { hl[lv.i(i, j)] = h[hv.i(i, j)]; } }
            let mut wl = w.to_vec();
            wl.resize(49, Z);
            let mut workl = vec![Z; 49];
            let info = super::eig_aed::zlaqr0(false, false, 49, ilo, info, &mut hl, lv, &mut wl, ilo, ihi, lv, &mut workl);
            w.copy_from_slice(&wl[..n]);
            info
        }
    };
    if info > 0 { return Err(format!("eigenvalues: LAPACK zhseqr failed to converge (info {info})")); }
    Ok(())
}

/// Eigenvalues of the row-major n x n matrix `m` in libroadrunner's order, rounded as libstruct.
pub fn eigenvalues(m: &[f64], n: usize) -> Result<Vec<(f64, f64)>, String> {
    if n == 0 { return Ok(vec![]); }
    if m.iter().any(|x| !x.is_finite()) { return Err("eigenvalues: the Jacobian has non-finite entries".into()); }
    // A = M^T column-major = M row-major read as is
    let mut a: Vec<C> = m.iter().map(|&x| c(x, 0.0)).collect();
    let anrm = m.iter().fold(0.0f64, |s, x| s.max(x.abs())); // zlange('M')
    let smlnum = SAFMIN.sqrt() / ULP;
    let cscale = if anrm > 0.0 && anrm < smlnum { smlnum } else if anrm > 1.0 / smlnum { 1.0 / smlnum } else { 0.0 };
    if cscale != 0.0 { zlascl(anrm, cscale, &mut a); }
    let av = V { o: 0, ld: n };
    let (ilo, ihi) = zgebal(n, &mut a, av);
    let mut tau = vec![Z; n];
    zgehd2(n, ilo, ihi, &mut a, av, &mut tau);
    let mut w = vec![Z; n];
    zhseqr(n, ilo, ihi, &mut a, &mut w)?;
    if cscale != 0.0 { zlascl(cscale, anrm, &mut w); }
    Ok(w.iter().map(|z| (round_tol(z.r, 1e-12), round_tol(z.i, 1e-12))).collect())
}

/// zlascl('G'): x *= cto / cfrom without over/underflow.
fn zlascl(cfrom: f64, cto: f64, x: &mut [C]) {
    let (smlnum, bignum) = (SAFMIN, 1.0 / SAFMIN);
    let (mut cf, mut ct) = (cfrom, cto);
    loop {
        let cf1 = cf * smlnum;
        let (mul, done);
        if cf1 == cf { mul = ct / cf; done = true; } else {
            let ct1 = ct / bignum;
            if ct1 == ct { mul = ct; done = true; cf = 1.0; }
            else if cf1.abs() > ct.abs() && ct != 0.0 { mul = smlnum; done = false; cf = cf1; }
            else if ct1.abs() > cf.abs() { mul = bignum; done = false; ct = ct1; }
            else { mul = ct / cf; done = true; }
        }
        for z in x.iter_mut() { *z = scl(mul, *z); }
        if done { break; }
    }
}

#[cfg(test)]
mod tests {
    /// Values and order from CLAPACK zgeev compiled from libroadrunner-deps (as libroadrunner).
    #[test]
    fn lapack_order() {
        let e = super::eigenvalues(&[0., 1., 0., 0., 0., 1., -6., -11., -6.], 3).unwrap();
        assert_eq!(e, vec![(-1., 0.), (-3., 0.), (-2., 0.)]);
        let e = super::eigenvalues(&[1., 2., 0., 0., 3., 4., 0., 0., 0., 0., 0., -1., 0., 0., 1., 0.], 4).unwrap();
        assert_eq!(e, vec![(-0.37228132326901453, 0.), (5.372281323269014, 0.), (0., 1.), (0., -1.)]);
    }
}

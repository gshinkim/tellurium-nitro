//! CLAPACK's multishift QR with aggressive early deflation, as zhseqr runs it above n = 75: zlaqr0 /
//! zlaqr4, zlaqr3 / zlaqr2, zlaqr5, zlaqr1, ztrexc, zlartg, zunmhr (zunmqr, zunm2r), zgehrd (blocked:
//! zlahr2) and the reference BLAS / LAPACK kernels they call (zgemm, ztrmm, ztrmv, zgemv, zlarft,
//! zlarfb). Workspace splits and block sizes (ilaenv / iparmq) follow the C code, so the rounding and
//! hence the eigenvalue order match libroadrunner. Conventions: see eig.rs.

use super::eig::*;
use super::cons::fma;

const NMIN: usize = 75; // iparmq ISPEC 12
const NIBBLE: usize = 14; // iparmq ISPEC 14

/// iparmq's shift count for an active block of size nh (ISPEC 13 / 15 / 16 share it).
fn iparmq_ns(nh: usize) -> usize {
    let mut ns = 2;
    if nh >= 30 { ns = 4; }
    if nh >= 60 { ns = 10; }
    if nh >= 150 {
        let r = (((nh as f64).ln() / 2f64.ln()) as f32).round() as usize; // i_nint of a real
        ns = 10.max(nh / r);
    }
    if nh >= 590 { ns = 64; }
    if nh >= 3000 { ns = 128; }
    if nh >= 6000 { ns = 256; }
    2.max(ns - ns % 2)
}

/// A copy of the m x n block at v (ld m), for operands that alias the output buffer.
fn take(h: &[C], v: V, m: usize, n: usize) -> (Vec<C>, V) {
    let ld = m.max(1);
    let mut out = vec![Z; ld * n.max(1)];
    for j in 1..=n { for i in 1..=m { out[(i - 1) + (j - 1) * ld] = h[v.i(i, j)]; } }
    (out, V { o: 0, ld })
}
fn lacpy(h: &mut [C], from: V, to: V, m: usize, n: usize) {
    for j in 1..=n { for i in 1..=m { h[to.i(i, j)] = h[from.i(i, j)]; } }
}
fn laset(h: &mut [C], v: V, m: usize, n: usize, off: C, diag: C) {
    for j in 1..=n { for i in 1..=m { h[v.i(i, j)] = if i == j { diag } else { off }; } }
}

/// Reference zgemm: C = alpha op(A) op(B) + beta C; t: b'N', b'T' or b'C'.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zgemm(ta: u8, tb: u8, m: usize, n: usize, k: usize, alpha: C, a: &[C], av: V, b: &[C], bv: V, beta: C, cb: &mut [C], cv: V) {
    if m == 0 || n == 0 || ((!nz(alpha) || k == 0) && beta == ONE) { return; }
    let ci = |i: usize, j: usize| cv.i(i, j);
    let opa = |x: C| if ta == b'C' { cj(x) } else { x };
    let opb = |x: C| if tb == b'C' { cj(x) } else { x };
    if !nz(alpha) {
        for j in 1..=n { for i in 1..=m { let q = ci(i, j); cb[q] = if nz(beta) { mul(beta, cb[q]) } else { Z }; } }
        return;
    }
    let scale_col = |cb: &mut [C], j: usize| {
        if !nz(beta) { for i in 1..=m { cb[ci(i, j)] = Z; } } else if beta != ONE { for i in 1..=m { let q = ci(i, j); cb[q] = mul(beta, cb[q]); } }
    };
    if ta == b'N' {
        // axpy form: C(:,j) += (alpha op(B)(l,j)) A(:,l)
        for j in 1..=n {
            scale_col(cb, j);
            for l in 1..=k {
                let bl = if tb == b'N' { b[bv.i(l, j)] } else { b[bv.i(j, l)] };
                if nz(bl) {
                    let temp = mul(alpha, opb(bl));
                    for i in 1..=m { let q = ci(i, j); cb[q] = add(cb[q], mul(temp, a[av.i(i, l)])); }
                }
            }
        }
    } else {
        // dot form
        for j in 1..=n {
            for i in 1..=m {
                let mut temp = Z;
                for l in 1..=k {
                    let bl = if tb == b'N' { b[bv.i(l, j)] } else { opb(b[bv.i(j, l)]) };
                    temp = add(temp, mul(opa(a[av.i(l, i)]), bl));
                }
                let q = ci(i, j);
                cb[q] = if !nz(beta) { mul(alpha, temp) } else { add(mul(alpha, temp), mul(beta, cb[q])) };
            }
        }
    }
}

/// Reference ztrmm: B = alpha op(A) B (left) or alpha B op(A) (right), A triangular.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ztrmm(left: bool, upper: bool, tr: u8, unit: bool, m: usize, n: usize, alpha: C, a: &[C], av: V, b: &mut [C], bv: V) {
    if m == 0 || n == 0 { return; }
    let ai = |i: usize, j: usize| a[av.i(i, j)];
    let bi = |i: usize, j: usize| bv.i(i, j);
    let op = |x: C| if tr == b'C' { cj(x) } else { x };
    if !nz(alpha) { for j in 1..=n { for i in 1..=m { b[bi(i, j)] = Z; } } return; }
    if left {
        if tr == b'N' {
            if upper {
                for j in 1..=n { for k in 1..=m {
                    if nz(b[bi(k, j)]) {
                        let mut temp = mul(alpha, b[bi(k, j)]);
                        for i in 1..k { let q = bi(i, j); b[q] = add(b[q], mul(temp, ai(i, k))); }
                        if !unit { temp = mul(temp, ai(k, k)); }
                        b[bi(k, j)] = temp;
                    }
                } }
            } else {
                for j in 1..=n { for k in (1..=m).rev() {
                    if nz(b[bi(k, j)]) {
                        let temp = mul(alpha, b[bi(k, j)]);
                        b[bi(k, j)] = temp;
                        if !unit { b[bi(k, j)] = mul(b[bi(k, j)], ai(k, k)); }
                        for i in k + 1..=m { let q = bi(i, j); b[q] = add(b[q], mul(temp, ai(i, k))); }
                    }
                } }
            }
        } else if upper {
            for j in 1..=n { for i in (1..=m).rev() {
                let mut temp = b[bi(i, j)];
                if !unit { temp = mul(temp, op(ai(i, i))); }
                for k in 1..i { temp = add(temp, mul(op(ai(k, i)), b[bi(k, j)])); }
                b[bi(i, j)] = mul(alpha, temp);
            } }
        } else {
            for j in 1..=n { for i in 1..=m {
                let mut temp = b[bi(i, j)];
                if !unit { temp = mul(temp, op(ai(i, i))); }
                for k in i + 1..=m { temp = add(temp, mul(op(ai(k, i)), b[bi(k, j)])); }
                b[bi(i, j)] = mul(alpha, temp);
            } }
        }
    } else if tr == b'N' {
        if upper {
            for j in (1..=n).rev() {
                let mut temp = alpha;
                if !unit { temp = mul(temp, ai(j, j)); }
                for i in 1..=m { let q = bi(i, j); b[q] = mul(temp, b[q]); }
                for k in 1..j {
                    if nz(ai(k, j)) {
                        let temp = mul(alpha, ai(k, j));
                        for i in 1..=m { let q = bi(i, j); b[q] = add(b[q], mul(temp, b[bi(i, k)])); }
                    }
                }
            }
        } else {
            for j in 1..=n {
                let mut temp = alpha;
                if !unit { temp = mul(temp, ai(j, j)); }
                for i in 1..=m { let q = bi(i, j); b[q] = mul(temp, b[q]); }
                for k in j + 1..=n {
                    if nz(ai(k, j)) {
                        let temp = mul(alpha, ai(k, j));
                        for i in 1..=m { let q = bi(i, j); b[q] = add(b[q], mul(temp, b[bi(i, k)])); }
                    }
                }
            }
        }
    } else if upper {
        for k in 1..=n {
            for j in 1..k {
                if nz(ai(j, k)) {
                    let temp = mul(alpha, op(ai(j, k)));
                    for i in 1..=m { let q = bi(i, j); b[q] = add(b[q], mul(temp, b[bi(i, k)])); }
                }
            }
            let mut temp = alpha;
            if !unit { temp = mul(temp, op(ai(k, k))); }
            if temp != ONE { for i in 1..=m { let q = bi(i, k); b[q] = mul(temp, b[q]); } }
        }
    } else {
        for k in (1..=n).rev() {
            for j in k + 1..=n {
                if nz(ai(j, k)) {
                    let temp = mul(alpha, op(ai(j, k)));
                    for i in 1..=m { let q = bi(i, j); b[q] = add(b[q], mul(temp, b[bi(i, k)])); }
                }
            }
            let mut temp = alpha;
            if !unit { temp = mul(temp, op(ai(k, k))); }
            if temp != ONE { for i in 1..=m { let q = bi(i, k); b[q] = mul(temp, b[q]); } }
        }
    }
}

/// Reference ztrmv (incx = 1): x = op(A) x.
fn ztrmv(upper: bool, tr: u8, unit: bool, n: usize, a: &[C], av: V, x: &mut [C]) {
    let ai = |i: usize, j: usize| a[av.i(i, j)];
    let op = |z: C| if tr == b'C' { cj(z) } else { z };
    if tr == b'N' {
        if upper {
            for j in 1..=n {
                if nz(x[j - 1]) {
                    let temp = x[j - 1];
                    for i in 1..j { x[i - 1] = add(x[i - 1], mul(temp, ai(i, j))); }
                    if !unit { x[j - 1] = mul(x[j - 1], ai(j, j)); }
                }
            }
        } else {
            for j in (1..=n).rev() {
                if nz(x[j - 1]) {
                    let temp = x[j - 1];
                    for i in (j + 1..=n).rev() { x[i - 1] = add(x[i - 1], mul(temp, ai(i, j))); }
                    if !unit { x[j - 1] = mul(x[j - 1], ai(j, j)); }
                }
            }
        }
    } else if upper {
        for j in (1..=n).rev() {
            let mut temp = x[j - 1];
            if !unit { temp = mul(temp, op(ai(j, j))); }
            for i in (1..j).rev() { temp = add(temp, mul(op(ai(i, j)), x[i - 1])); }
            x[j - 1] = temp;
        }
    } else {
        for j in 1..=n {
            let mut temp = x[j - 1];
            if !unit { temp = mul(temp, op(ai(j, j))); }
            for i in j + 1..=n { temp = add(temp, mul(op(ai(i, j)), x[i - 1])); }
            x[j - 1] = temp;
        }
    }
}

/// Reference zgemv: y = alpha op(A) x + beta y; y[(i-1)*incy].
#[allow(clippy::too_many_arguments)]
fn zgemv(tr: u8, m: usize, n: usize, alpha: C, a: &[C], av: V, x: &[C], beta: C, y: &mut [C], incy: usize) {
    if m == 0 || n == 0 || (!nz(alpha) && beta == ONE) { return; }
    let leny = if tr == b'N' { m } else { n };
    if beta != ONE {
        for i in 0..leny { y[i * incy] = if nz(beta) { mul(beta, y[i * incy]) } else { Z }; }
    }
    if !nz(alpha) { return; }
    if tr == b'N' {
        for j in 1..=n {
            if nz(x[j - 1]) {
                let temp = mul(alpha, x[j - 1]);
                for i in 1..=m { let q = (i - 1) * incy; y[q] = add(y[q], mul(temp, a[av.i(i, j)])); }
            }
        }
    } else {
        for j in 1..=n {
            let mut temp = Z;
            for i in 1..=m {
                let aij = if tr == b'C' { cj(a[av.i(i, j)]) } else { a[av.i(i, j)] };
                temp = add(temp, mul(aij, x[i - 1]));
            }
            let q = (j - 1) * incy;
            y[q] = add(y[q], mul(alpha, temp));
        }
    }
}

/// ilazlr: last non-zero row of the m x n block.
fn ilazlr(h: &[C], v: V, m: usize, n: usize) -> usize {
    if m == 0 { return 0; }
    if nz(h[v.i(m, 1)]) || nz(h[v.i(m, n)]) { return m; }
    (1..=n).map(|j| { let mut i = m; while i >= 1 && !nz(h[v.i(i, j)]) { i -= 1; } i }).max().unwrap_or(0)
}
/// ilazlc: last non-zero column of the m x n block.
fn ilazlc(h: &[C], v: V, m: usize, n: usize) -> usize {
    if n == 0 { return 0; }
    if nz(h[v.i(1, n)]) || nz(h[v.i(m, n)]) { return n; }
    let mut k = n;
    while k >= 1 && !(1..=m).any(|i| nz(h[v.i(i, k)])) { k -= 1; }
    k
}

/// zlarft('F', 'C'): T (k x k, upper) of the block reflector with vectors in vv (n x k).
fn zlarft(n: usize, k: usize, h: &mut [C], vv: V, tau: &[C], t: &mut [C], tv: V) {
    if n == 0 { return; }
    let mut prevlastv = n;

    for i in 1..=k {
        prevlastv = prevlastv.max(i);
        if !nz(tau[i - 1]) {
            for j in 1..=i { t[tv.i(j, i)] = Z; }
        } else {
            let vii = h[vv.i(i, i)];
            h[vv.i(i, i)] = ONE;
            let mut lastv = n;
            while lastv > i && !nz(h[vv.i(lastv, i)]) { lastv -= 1; }
            let j = lastv.min(prevlastv);
            let (a, av) = take(h, vv.at(i, 1), j - i + 1, i - 1);
            let x: Vec<C> = (i..=j).map(|q| h[vv.i(q, i)]).collect();
            let o = tv.i(1, i);
            zgemv(b'C', j - i + 1, i - 1, neg(tau[i - 1]), &a, av, &x, Z, &mut t[o..], 1);
            h[vv.i(i, i)] = vii;
            let (tt, ttv) = take(t, tv, i - 1, i - 1);
            ztrmv(true, b'N', false, i - 1, &tt, ttv, &mut t[o..]);
            t[tv.i(i, i)] = tau[i - 1];
            prevlastv = if i > 1 { prevlastv.max(lastv) } else { lastv };
        }
    }
}

/// zlarfb(side, trans, 'F', 'C'): apply H or H^H (block reflector V, T) to the m x n block cv.
#[allow(clippy::too_many_arguments)]
fn zlarfb(left: bool, trans: u8, m: usize, n: usize, k: usize, h: &mut [C], vv: V, t: &[C], tv: V, cv: V, work: &mut [C], ldw: usize) {
    if m == 0 || n == 0 { return; }
    let transt = if trans == b'N' { b'C' } else { b'N' };
    let wv = V { o: 0, ld: ldw };
    let m1 = neg(ONE);
    if left {
        let lastv = k.max(ilazlr(h, vv, m, k));
        let lastc = ilazlc(h, cv, lastv, n);
        for j in 1..=k { for i in 1..=lastc { work[wv.i(i, j)] = cj(h[cv.i(j, i)]); } }
        let (v1, v1v) = take(h, vv, k, k);
        ztrmm(false, false, b'N', true, lastc, k, ONE, &v1, v1v, work, wv);
        if lastv > k {
            let (cc, ccv) = take(h, cv.at(k + 1, 1), lastv - k, lastc);
            let (v2, v2v) = take(h, vv.at(k + 1, 1), lastv - k, k);
            zgemm(b'C', b'N', lastc, k, lastv - k, ONE, &cc, ccv, &v2, v2v, ONE, work, wv);
        }
        ztrmm(false, true, transt, false, lastc, k, ONE, t, tv, work, wv);
        if m > k {
            let (v2, v2v) = take(h, vv.at(k + 1, 1), lastv - k, k);
            let (w2, w2v) = take(work, wv, lastc, k);
            zgemm(b'N', b'C', lastv - k, lastc, k, m1, &v2, v2v, &w2, w2v, ONE, h, cv.at(k + 1, 1));
        }
        ztrmm(false, false, b'C', true, lastc, k, ONE, &v1, v1v, work, wv);
        for j in 1..=k { for i in 1..=lastc { let q = cv.i(j, i); h[q] = sub(h[q], cj(work[wv.i(i, j)])); } }
    } else {
        let lastv = k.max(ilazlr(h, vv, n, k));
        let lastc = ilazlr(h, cv, m, lastv);
        for j in 1..=k { for i in 1..=lastc { work[wv.i(i, j)] = h[cv.i(i, j)]; } }
        let (v1, v1v) = take(h, vv, k, k);
        ztrmm(false, false, b'N', true, lastc, k, ONE, &v1, v1v, work, wv);
        if lastv > k {
            let (cc, ccv) = take(h, cv.at(1, k + 1), lastc, lastv - k);
            let (v2, v2v) = take(h, vv.at(k + 1, 1), lastv - k, k);
            zgemm(b'N', b'N', lastc, k, lastv - k, ONE, &cc, ccv, &v2, v2v, ONE, work, wv);
        }
        ztrmm(false, true, trans, false, lastc, k, ONE, t, tv, work, wv);
        if lastv > k {
            let (v2, v2v) = take(h, vv.at(k + 1, 1), lastv - k, k);
            let (w2, w2v) = take(work, wv, lastc, k);
            zgemm(b'N', b'C', lastc, lastv - k, k, m1, &w2, w2v, &v2, v2v, ONE, h, cv.at(1, k + 1));
        }
        ztrmm(false, false, b'C', true, lastc, k, ONE, &v1, v1v, work, wv);
        for j in 1..=k { for i in 1..=lastc { let q = cv.i(i, j); h[q] = sub(h[q], work[wv.i(i, j)]); } }
    }
}

/// zunmhr('R', 'N', m, n, ilo = 1, ihi): C (m x n) = C Q, Q from zgehrd's reflectors in a / tau.
#[allow(clippy::too_many_arguments)]
fn zunmhr_rn(m: usize, n: usize, ihi: usize, h: &mut [C], a: V, tau: &[C], cv: V, work: &mut [C]) {
    let nh = ihi - 1;
    if m == 0 || n == 0 || nh == 0 { return; }
    // zunmqr('R', 'N', mi = m, ni = nh, k = nh, A(2, 1), tau, C(1, 2))
    let (ni, k, av, cv) = (nh, nh, a.at(2, 1), cv.at(1, 2));
    let nq = ni;
    let nw = m;
    let lwork = work.len();
    let mut nb = 32usize.min(64);
    let mut nbmin = 2;
    if nb > 1 && nb < k && lwork < nw * nb { nb = lwork / nw; nbmin = 2; }
    if nb < nbmin || nb >= k {
        // zunm2r, forward
        for i in 1..=k {
            let ni_i = ni - i + 1;
            let aii = h[av.i(i, i)];
            h[av.i(i, i)] = ONE;
            let v: Vec<C> = (0..ni_i).map(|q| h[av.i(i + q, i)]).collect();
            zlarf(false, m, ni_i, &v, tau[i - 1], h, cv.at(1, i));
            h[av.i(i, i)] = aii;
        }
    } else {
        let mut t = vec![Z; 65 * 64];
        let tv = V { o: 0, ld: 65 };
        let mut i = 1;
        while i <= k {
            let ib = nb.min(k - i + 1);
            zlarft(nq - i + 1, ib, h, av.at(i, i), &tau[i - 1..], &mut t, tv);
            zlarfb(false, b'N', m, ni - i + 1, ib, h, av.at(i, i), &t, tv, cv.at(1, i), work, nw);
            i += nb;
        }
    }
}

/// zgehrd(n, 1, ihi): blocked Hessenberg reduction (zlahr2 + zlarfb) when the workspace allows.
fn zgehrd(n: usize, ihi: usize, h: &mut [C], a: V, tau: &mut [C], work: &mut [C]) {
    let ilo = 1;
    for t in tau.iter_mut().take(n.saturating_sub(1)).skip(ihi.max(1) - 1) { *t = Z; }
    let nh = ihi - ilo + 1;
    if nh <= 1 { return; }
    let lwork = work.len();
    let mut nb = 32usize;
    let nbmin = 2;
    let mut nx = 0;
    if nb > 1 && nb < nh {
        nx = nb.max(128);
        if nx < nh && lwork < n * nb { nb = if lwork >= n * nbmin { lwork / n } else { 1 }; }
    }
    let ldwork = n;
    let mut i = ilo;
    if !(nb < nbmin || nb >= nh) {
        let mut t = vec![Z; 65 * 64];
        let tv = V { o: 0, ld: 65 };
        let yv = V { o: 0, ld: ldwork };
        while i + nx < ihi {
            let ib = nb.min(ihi - i);
            zlahr2(ihi, i, ib, h, a.at(1, i), &mut tau[i - 1..], &mut t, tv, work, yv);
            let q = a.i(i + ib, i + ib - 1);
            let ei = h[q];
            h[q] = ONE;
            let (y, yy) = take(work, yv, ihi, ib);
            let (vv, vvv) = take(h, a.at(i + ib, i), ihi - i - ib + 1, ib);
            zgemm(b'N', b'C', ihi, ihi - i - ib + 1, ib, neg(ONE), &y, yy, &vv, vvv, ONE, h, a.at(1, i + ib));
            h[q] = ei;
            let (l, lv) = take(h, a.at(i + 1, i), ib - 1, ib - 1);
            ztrmm(false, false, b'C', true, i, ib - 1, ONE, &l, lv, work, yv);
            for j in 0..ib - 1 {
                for r in 1..=i { let p = a.i(r, i + j + 1); h[p] = add(h[p], mul(neg(ONE), work[ldwork * j + r - 1])); }
            }
            let tt = t.clone();
            zlarfb(true, b'C', ihi - i, n - i - ib + 1, ib, h, a.at(i + 1, i), &tt, tv, a.at(i + 1, i + ib), work, ldwork);
            i += nb;
        }
    }
    zgehd2(n, i, ihi, h, a, tau);
}

/// zlahr2: reduce nb columns of the n-column block a (k rows above the reduced part).
#[allow(clippy::too_many_arguments)]
fn zlahr2(n: usize, k: usize, nb: usize, h: &mut [C], a: V, tau: &mut [C], t: &mut [C], tv: V, y: &mut [C], yv: V) {
    if n <= 1 { return; }
    let m1 = neg(ONE);
    let mut ei = Z;
    for i in 1..=nb {
        if i > 1 {
            // A(k+1:n, i) -= Y(k+1:n, 1:i-1) A(k+i-1, 1:i-1)^H
            let x: Vec<C> = (1..i).map(|q| cj(h[a.i(k + i - 1, q)])).collect();
            let (yy, yyv) = take(y, yv.at(k + 1, 1), n - k, i - 1);
            let o = a.i(k + 1, i);
            zgemv(b'N', n - k, i - 1, m1, &yy, yyv, &x, ONE, &mut h[o..], 1);
            // apply I - V T^H V^H from the left, with T(1:i-1, nb) as workspace
            let w0 = tv.i(1, nb);
            for q in 1..i { t[w0 + q - 1] = h[a.i(k + q, i)]; }
            let (l, lv) = take(h, a.at(k + 1, 1), i - 1, i - 1);
            ztrmv(false, b'C', true, i - 1, &l, lv, &mut t[w0..]);
            let (aa, aav) = take(h, a.at(k + i, 1), n - k - i + 1, i - 1);
            let x: Vec<C> = (k + i..=n).map(|q| h[a.i(q, i)]).collect();
            zgemv(b'C', n - k - i + 1, i - 1, ONE, &aa, aav, &x, ONE, &mut t[w0..], 1);
            let (tt, ttv) = take(t, tv, i - 1, i - 1);
            ztrmv(true, b'C', false, i - 1, &tt, ttv, &mut t[w0..]);
            let x: Vec<C> = t[w0..w0 + i - 1].to_vec();
            let o = a.i(k + i, i);
            zgemv(b'N', n - k - i + 1, i - 1, m1, &aa, aav, &x, ONE, &mut h[o..], 1);
            ztrmv(false, b'N', true, i - 1, &l, lv, &mut t[w0..]);
            for q in 1..i { let p = a.i(k + q, i); h[p] = add(h[p], mul(m1, t[w0 + q - 1])); }
            h[a.i(k + i - 1, i - 1)] = ei;
        }
        let mut alpha = h[a.i(k + i, i)];
        let mut x: Vec<C> = (k + i + 1..=n).map(|q| h[a.i(q, i)]).collect();
        tau[i - 1] = zlarfg(&mut alpha, &mut x);
        h[a.i(k + i, i)] = alpha;
        for (q, v) in (k + i + 1..=n).zip(&x) { h[a.i(q, i)] = *v; }
        ei = h[a.i(k + i, i)];
        h[a.i(k + i, i)] = ONE;
        let vcol: Vec<C> = (k + i..=n).map(|q| h[a.i(q, i)]).collect();
        let (a2, a2v) = take(h, a.at(k + 1, i + 1), n - k, n - k - i + 1);
        let o = yv.i(k + 1, i);
        zgemv(b'N', n - k, n - k - i + 1, ONE, &a2, a2v, &vcol, Z, &mut y[o..], 1);
        let (aa, aav) = take(h, a.at(k + i, 1), n - k - i + 1, i - 1);
        let o = tv.i(1, i);
        zgemv(b'C', n - k - i + 1, i - 1, ONE, &aa, aav, &vcol, Z, &mut t[o..], 1);
        let (yy, yyv) = take(y, yv.at(k + 1, 1), n - k, i - 1);
        let x: Vec<C> = t[o..o + i - 1].to_vec();
        let oy = yv.i(k + 1, i);
        zgemv(b'N', n - k, i - 1, m1, &yy, yyv, &x, ONE, &mut y[oy..], 1);
        for q in 0..n - k { y[oy + q] = mul(tau[i - 1], y[oy + q]); }
        for q in 0..i - 1 { t[o + q] = mul(neg(tau[i - 1]), t[o + q]); }
        let (tt, ttv) = take(t, tv, i - 1, i - 1);
        ztrmv(true, b'N', false, i - 1, &tt, ttv, &mut t[o..]);
        t[tv.i(i, i)] = tau[i - 1];
    }
    h[a.i(k + nb, nb)] = ei;
    // Y(1:k, 1:nb)
    for j in 1..=nb { for r in 1..=k { y[yv.i(r, j)] = h[a.i(r, j + 1)]; } }
    let (l, lv) = take(h, a.at(k + 1, 1), nb, nb);
    ztrmm(false, false, b'N', true, k, nb, ONE, &l, lv, y, yv);
    if n > k + nb {
        let (a1, a1v) = take(h, a.at(1, nb + 2), k, n - k - nb);
        let (a2, a2v) = take(h, a.at(k + 1 + nb, 1), n - k - nb, nb);
        zgemm(b'N', b'N', k, nb, n - k - nb, ONE, &a1, a1v, &a2, a2v, ONE, y, yv);
    }
    let (tt, ttv) = take(t, tv, nb, nb);
    ztrmm(false, true, b'N', false, k, nb, ONE, &tt, ttv, y, yv);
}

/// zlartg: plane rotation [cs sn; -conj(sn) cs] [f; g] = [r; 0].
fn zlartg(f: C, g: C) -> (f64, C, C) {
    let safmn2 = 2f64.powi(-484); // dlamch('B')^int(log(safmin/eps)/log(2)/2)
    let safmx2 = 1.0 / safmn2;
    let mut scale = f.r.abs().max(f.i.abs()).max(g.r.abs().max(g.i.abs()));
    let (mut fs, mut gs) = (f, g);
    let mut count = 0i32;
    if scale >= safmx2 {
        loop { count += 1; fs = scl(safmn2, fs); gs = scl(safmn2, gs); scale *= safmn2; if scale < safmx2 { break; } }
    } else if scale <= safmn2 {
        if !nz(g) { return (1.0, Z, f); }
        loop { count -= 1; fs = scl(safmx2, fs); gs = scl(safmx2, gs); scale *= safmx2; if scale > safmn2 { break; } }
    }
    let f2 = fma(fs.r, fs.r, fs.i * fs.i);
    let g2 = fma(gs.r, gs.r, gs.i * gs.i);
    if f2 <= g2.max(1.0) * SAFMIN {
        if !nz(f) {
            let r = c(dlapy2(g.r, g.i), 0.0);
            let d = dlapy2(gs.r, gs.i);
            return (0.0, c(gs.r / d, -gs.i / d), r);
        }
        let f2s = dlapy2(fs.r, fs.i);
        let g2s = g2.sqrt();
        let cs = f2s / g2s;
        let ff = if f.r.abs().max(f.i.abs()) > 1.0 {
            let d = dlapy2(f.r, f.i);
            c(f.r / d, f.i / d)
        } else {
            let (dr, di) = (safmx2 * f.r, safmx2 * f.i);
            let d = dlapy2(dr, di);
            c(dr / d, di / d)
        };
        let sn = mul(ff, c(gs.r / g2s, -gs.i / g2s));
        let r = add(scl(cs, f), mul(sn, g));
        (cs, sn, r)
    } else {
        let f2s = (g2 / f2 + 1.0).sqrt();
        let mut r = c(f2s * fs.r, f2s * fs.i);
        let cs = 1.0 / f2s;
        let d = f2 + g2;
        let sn = mul(c(r.r / d, r.i / d), cj(gs));
        if count > 0 { for _ in 0..count { r = scl(safmx2, r); } }
        if count < 0 { for _ in 0..-count { r = scl(safmn2, r); } }
        (cs, sn, r)
    }
}

/// zrot on the element index pairs (x, y).
fn zrot(h: &mut [C], pairs: impl Iterator<Item = (usize, usize)>, cs: f64, s: C) {
    for (x, y) in pairs {
        let stemp = add(scl(cs, h[x]), mul(s, h[y]));
        h[y] = sub(scl(cs, h[y]), mul(cj(s), h[x]));
        h[x] = stemp;
    }
}

/// ztrexc('V'): move the diagonal entry ifst of the upper triangular n x n t to ilst; Q updated.
fn ztrexc(n: usize, h: &mut [C], t: V, q: V, ifst: usize, ilst: usize) {
    if n == 1 || ifst == ilst { return; }
    let ks: Vec<usize> = if ifst < ilst { (ifst..ilst).collect() } else { (ilst..ifst).rev().collect() };
    for k in ks {
        let (t11, t22) = (h[t.i(k, k)], h[t.i(k + 1, k + 1)]);
        let (cs, sn, _) = zlartg(h[t.i(k, k + 1)], sub(t22, t11)); // (T(k, k+1) keeps its value, as ztrexc)
        if k + 2 <= n { zrot(h, (k + 2..=n).map(|j| (t.i(k, j), t.i(k + 1, j))), cs, sn); }
        zrot(h, (1..k).map(|i| (t.i(i, k), t.i(i, k + 1))), cs, cj(sn));
        h[t.i(k, k)] = t22;
        h[t.i(k + 1, k + 1)] = t11;
        zrot(h, (1..=n).map(|i| (q.i(i, k), q.i(i, k + 1))), cs, cj(sn));
    }
}

/// zlaqr1: first column of (H - s1)(H - s2) for the n = 2 or 3 block at hv, scaled.
fn zlaqr1(n: usize, h: &[C], hv: V, s1: C, s2: C) -> [C; 3] {
    let x = |i: usize, j: usize| h[hv.i(i, j)];
    if n == 2 {
        let s = abs1(sub(x(1, 1), s2)) + abs1(x(2, 1));
        if s == 0.0 { return [Z; 3]; }
        let h21s = div_r(x(2, 1), s);
        let v1 = add(mul(h21s, x(1, 2)), mul(sub(x(1, 1), s1), div_r(sub(x(1, 1), s2), s)));
        let v2 = mul(h21s, sub(sub(add(x(1, 1), x(2, 2)), s1), s2));
        [v1, v2, Z]
    } else {
        let s = abs1(sub(x(1, 1), s2)) + abs1(x(2, 1)) + abs1(x(3, 1));
        if s == 0.0 { return [Z; 3]; }
        let h21s = div_r(x(2, 1), s);
        let h31s = div_r(x(3, 1), s);
        let v1 = add(add(mul(sub(x(1, 1), s1), div_r(sub(x(1, 1), s2), s)), mul(x(1, 2), h21s)), mul(x(1, 3), h31s));
        let v2 = add(mul(h21s, sub(sub(add(x(1, 1), x(2, 2)), s1), s2)), mul(x(2, 3), h31s));
        let v3 = add(mul(h31s, sub(sub(add(x(1, 1), x(3, 3)), s1), s2)), mul(h21s, x(3, 2)));
        [v1, v2, v3]
    }
}

/// zlaqr0 (`aed3`: aggressive early deflation by zlaqr3) or zlaqr4 (zlaqr2). Returns info.
#[allow(clippy::too_many_arguments)]
pub(crate) fn zlaqr0(wantt: bool, wantz: bool, n: usize, ilo: usize, ihi: usize, h: &mut [C], hv: V, w: &mut [C],
                     iloz: usize, ihiz: usize, zv: V, work: &mut [C]) -> usize {
    zlaqr04(true, wantt, wantz, n, ilo, ihi, h, hv, w, iloz, ihiz, zv, work)
}

#[allow(clippy::too_many_arguments)]
fn zlaqr04(aed3: bool, wantt: bool, wantz: bool, n: usize, ilo: usize, ihi: usize, h: &mut [C], hv: V, w: &mut [C],
           iloz: usize, ihiz: usize, zv: V, work: &mut [C]) -> usize {
    if n == 0 { return 0; }
    if n <= 11 { return zlahqr(wantt, wantz, n, ilo, ihi, h, hv, w, iloz, ihiz, zv); }
    let lwork = work.len();
    let ix = |i: usize, j: usize| hv.i(i, j);
    let nh0 = ihi - ilo + 1;
    let ns0 = iparmq_ns(nh0);
    let mut nwr = if nh0 <= 500 { ns0 } else { ns0 * 3 / 2 };
    nwr = 2.max(nwr);
    nwr = nwr.min(nh0.min((n - 1) / 3));
    let mut nsr = ns0;
    nsr = nsr.min((n + 6) / 9).min(ihi - ilo);
    nsr = 2.max(nsr - nsr % 2);
    let nmin = NMIN.max(11);
    let nibble = NIBBLE;
    let kacc22 = { let mut k = 0; if ns0 >= 14 { k = 1; } if ns0 >= 14 { k = 2; } k };
    let nwmax = ((n - 1) / 3).min(lwork / 2);
    let mut nw = nwmax;
    let mut nsmax = ((n + 6) / 9).min(2 * lwork / 3);
    nsmax -= nsmax % 2;
    let mut ndfl = 1usize;
    let itmax = 10.max(ihi - ilo + 1) * 30;
    let mut kbot = ihi as isize;
    let mut ndec: isize = -1;
    for _it in 1..=itmax {
        if kbot < ilo as isize { return 0; }
        let kb = kbot as usize;
        let mut k = kb;
        while k > ilo && nz(h[ix(k, k - 1)]) { k -= 1; }
        let ktop = k;
        let nh = kb - ktop + 1;
        let nwupbd = nh.min(nwmax);
        if ndfl < 5 { nw = nwupbd.min(nwr); } else { nw = nwupbd.min(2 * nw); }
        if nw < nwmax {
            if nw + 1 >= nh { nw = nh; } else {
                let kwtop = kb - nw + 1;
                if abs1(h[ix(kwtop, kwtop - 1)]) > abs1(h[ix(kwtop - 1, kwtop - 2)]) { nw += 1; }
            }
        }
        if ndfl < 5 { ndec = -1; } else if ndec >= 0 || nw >= nwupbd {
            ndec += 1;
            if (nw as isize) - ndec < 2 { ndec = 0; }
            nw = (nw as isize - ndec) as usize;
        }
        let kv = n - nw + 1;
        let kt = nw + 1;
        let nho = (n as isize - nw as isize - 1 - kt as isize + 1).max(0) as usize;
        let kwv = nw + 2;
        let nve = (n as isize - nw as isize - kwv as isize + 1).max(0) as usize;
        let (ls, ld) = zlaqr23(aed3, wantt, wantz, n, ktop, kb, nw, h, hv, iloz, ihiz, zv, w,
                               hv.at(kv, 1), nho, hv.at(kv, kt), nve, hv.at(kwv, 1), work);
        kbot -= ld as isize;
        let kb = kbot as usize;
        let mut ks = kb as isize - ls as isize + 1;
        if ld == 0 || (ld * 100 <= nw * nibble && (kbot - ktop as isize + 1) > nmin.min(nwmax) as isize) {
            let mut ns = nsmax.min(nsr).min(2.max(kb.saturating_sub(ktop)));
            ns -= ns % 2;
            if ndfl % 6 == 0 {
                ks = kb as isize - ns as isize + 1;
                let mut i = kb as isize;
                while i >= ks + 1 {
                    let iu = i as usize;
                    let d = abs1(h[ix(iu, iu - 1)]) * 0.75;
                    w[iu - 1] = c(h[ix(iu, iu)].r + d, h[ix(iu, iu)].i);
                    w[iu - 2] = w[iu - 1];
                    i -= 2;
                }
            } else {
                if (kb as isize - ks + 1) <= (ns / 2) as isize {
                    ks = kb as isize - ns as isize + 1;
                    let ksu = ks as usize;
                    let kt = n - ns + 1;
                    lacpy(h, hv.at(ksu, ksu), hv.at(kt, 1), ns, ns);
                    let inf = if aed3 && ns > nmin {
                        zlaqr04(false, false, false, ns, 1, ns, h, hv.at(kt, 1), &mut w[ksu - 1..], 1, 1, hv.at(kt, 1), work)
                    } else {
                        zlahqr(false, false, ns, 1, ns, h, hv.at(kt, 1), &mut w[ksu - 1..], 1, 1, hv.at(kt, 1))
                    };
                    ks += inf as isize;
                    if ks >= kb as isize {
                        let s = abs1(h[ix(kb - 1, kb - 1)]) + abs1(h[ix(kb, kb - 1)]) + abs1(h[ix(kb - 1, kb)]) + abs1(h[ix(kb, kb)]);
                        let aa = div_r(h[ix(kb - 1, kb - 1)], s);
                        let cc = div_r(h[ix(kb, kb - 1)], s);
                        let bb = div_r(h[ix(kb - 1, kb)], s);
                        let dd = div_r(h[ix(kb, kb)], s);
                        let tr2 = div_r(add(aa, dd), 2.0);
                        let det = sub(mul(sub(aa, tr2), sub(dd, tr2)), mul(bb, cc));
                        let rtdisc = zsqrt(neg(det));
                        w[kb - 2] = scl(s, add(tr2, rtdisc));
                        w[kb - 1] = scl(s, sub(tr2, rtdisc));
                        ks = kb as isize - 1;
                    }
                }
                if (kb as isize - ks + 1) > ns as isize {
                    // bubble sort by decreasing |re| + |im|
                    let ksu = ks as usize;
                    let mut sorted = false;
                    let mut k = kb;
                    while k > ksu {
                        if sorted { break; }
                        sorted = true;
                        for i in ksu..k {
                            if abs1(w[i - 1]) < abs1(w[i]) { sorted = false; w.swap(i - 1, i); }
                        }
                        k -= 1;
                    }
                }
            }
            if (kb as isize - ks + 1) == 2 {
                if abs1(sub(w[kb - 1], h[ix(kb, kb)])) < abs1(sub(w[kb - 2], h[ix(kb, kb)])) { w[kb - 2] = w[kb - 1]; } else { w[kb - 1] = w[kb - 2]; }
            }
            ns = ns.min((kb as isize - ks + 1) as usize);
            ns -= ns % 2;
            let ks = kb - ns + 1;
            let kdu = 3 * ns - 3;
            let ku = n - kdu + 1;
            let kwh = kdu + 1;
            let nho = (n as isize - kdu as isize - 3 - (kdu as isize + 1) + 1).max(0) as usize;
            let kwv = kdu + 4;
            let nve = (n as isize - kdu as isize - kwv as isize + 1).max(0) as usize;
            let s: Vec<C> = w[ks - 1..ks - 1 + ns].to_vec();
            zlaqr5(wantt, wantz, kacc22, n, ktop, kb, ns, &s, h, hv, iloz, ihiz, zv, work,
                   hv.at(ku, 1), nve, hv.at(kwv, 1), nho, hv.at(ku, kwh));
        }
        if ld > 0 { ndfl = 1; } else { ndfl += 1; }
    }
    kbot.max(0) as usize
}

/// zlaqr3 (`aed3`) / zlaqr2: aggressive early deflation of the nw-window ending at kbot.
/// Returns (ns, nd); shifts into w[kwtop-1..].
#[allow(clippy::too_many_arguments)]
fn zlaqr23(aed3: bool, wantt: bool, wantz: bool, n: usize, ktop: usize, kbot: usize, nw: usize, h: &mut [C], hv: V,
           iloz: usize, ihiz: usize, zv: V, sh: &mut [C], vv: V, nh: usize, tv: V, nv: usize, wvv: V, work: &mut [C]) -> (usize, usize) {
    let ix = |i: usize, j: usize| hv.i(i, j);
    if ktop > kbot || nw < 1 { return (0, 0); }
    let smlnum = SAFMIN * (n as f64 / ULP);
    let jw = nw.min(kbot - ktop + 1);
    let kwtop = kbot - jw + 1;
    let mut s = if kwtop == ktop { Z } else { h[ix(kwtop, kwtop - 1)] };
    if kbot == kwtop {
        sh[kwtop - 1] = h[ix(kwtop, kwtop)];
        let (mut ns, mut nd) = (1, 0);
        if abs1(s) <= smlnum.max(ULP * abs1(h[ix(kwtop, kwtop)])) {
            ns = 0; nd = 1;
            if kwtop > ktop { h[ix(kwtop, kwtop - 1)] = Z; }
        }
        return (ns, nd);
    }
    // T = upper part of the window plus its subdiagonal; V = I
    for j in 1..=jw { for i in 1..=j.min(jw) { h[tv.i(i, j)] = h[ix(kwtop + i - 1, kwtop + j - 1)]; } }
    for q in 1..jw { h[tv.i(q + 1, q)] = h[ix(kwtop + q, kwtop + q - 1)]; }
    laset(h, vv, jw, jw, Z, ONE);
    let infqr = if aed3 && jw > NMIN {
        zlaqr04(false, true, true, jw, 1, jw, h, tv, &mut sh[kwtop - 1..], 1, jw, vv, work)
    } else {
        zlahqr(true, true, jw, 1, jw, h, tv, &mut sh[kwtop - 1..], 1, jw, vv)
    };
    // deflation detection
    let mut ns = jw;
    let mut ilst = infqr + 1;
    for _knt in infqr + 1..=jw {
        let mut foo = abs1(h[tv.i(ns, ns)]);
        if foo == 0.0 { foo = abs1(s); }
        if abs1(s) * abs1(h[vv.i(1, ns)]) <= smlnum.max(ULP * foo) {
            ns -= 1;
        } else {
            ztrexc(jw, h, tv, vv, ns, ilst);
            ilst += 1;
        }
    }
    if ns == 0 { s = Z; }
    if ns < jw {
        // sort the undeflatable eigenvalues by decreasing magnitude (selection with ztrexc)
        for i in infqr + 1..=ns {
            let mut ifst = i;
            for j in i + 1..=ns { if abs1(h[tv.i(j, j)]) > abs1(h[tv.i(ifst, ifst)]) { ifst = j; } }
            if ifst != i { ztrexc(jw, h, tv, vv, ifst, i); }
        }
    }
    for i in infqr + 1..=jw { sh[kwtop + i - 2] = h[tv.i(i, i)]; }
    if ns < jw || !nz(s) {
        if ns > 1 && nz(s) {
            // reflect the spike back into upper Hessenberg form
            let mut v: Vec<C> = (1..=ns).map(|q| cj(h[vv.i(1, q)])).collect();
            let mut beta = v[0];
            let tau = zlarfg(&mut beta, &mut v[1..]);
            v[0] = ONE;
            // zlaset('L'): the lower triangle (diagonal included) of T(3:, 1:jw-2)
            for j in 1..=jw.saturating_sub(2) { for i in j..=jw - 2 { h[tv.i(i + 2, j)] = Z; } }
            zlarf(true, ns, jw, &v, cj(tau), h, tv);
            zlarf(false, ns, ns, &v, tau, h, tv);
            zlarf(false, jw, ns, &v, tau, h, vv);
            let (tauv, wk) = work.split_at_mut(jw);
            zgehrd(jw, ns, h, tv, tauv, wk);
        }
        if kwtop > 1 { h[ix(kwtop, kwtop - 1)] = mul(s, cj(h[vv.i(1, 1)])); }
        for j in 1..=jw { for i in 1..=j { h[ix(kwtop + i - 1, kwtop + j - 1)] = h[tv.i(i, j)]; } }
        for q in 1..jw { h[ix(kwtop + q, kwtop + q - 1)] = h[tv.i(q + 1, q)]; }
        if ns > 1 && nz(s) {
            let (tauv, wk) = work.split_at_mut(jw);
            zunmhr_rn(jw, ns, ns, h, tv, tauv, vv, wk);
        }
        // update the rest of H and Z
        let ltop = if wantt { 1 } else { ktop };
        let mut krow = ltop;
        while krow + 1 <= kwtop {
            let kln = nv.min(kwtop - krow);
            let (a, av) = take(h, hv.at(krow, kwtop), kln, jw);
            let (b, bv) = take(h, vv, jw, jw);
            zgemm(b'N', b'N', kln, jw, jw, ONE, &a, av, &b, bv, Z, h, wvv);
            lacpy(h, wvv, hv.at(krow, kwtop), kln, jw);
            krow += nv;
        }
        if wantt {
            let mut kcol = kbot + 1;
            while kcol <= n {
                let kln = nh.min(n - kcol + 1);
                let (a, av) = take(h, vv, jw, jw);
                let (b, bv) = take(h, hv.at(kwtop, kcol), jw, kln);
                zgemm(b'C', b'N', jw, kln, jw, ONE, &a, av, &b, bv, Z, h, tv);
                lacpy(h, tv, hv.at(kwtop, kcol), jw, kln);
                kcol += nh;
            }
        }
        if wantz {
            let mut krow = iloz;
            while krow <= ihiz {
                let kln = nv.min(ihiz - krow + 1);
                let (a, av) = take(h, zv.at(krow, kwtop), kln, jw);
                let (b, bv) = take(h, vv, jw, jw);
                zgemm(b'N', b'N', kln, jw, jw, ONE, &a, av, &b, bv, Z, h, wvv);
                lacpy(h, wvv, zv.at(krow, kwtop), kln, jw);
                krow += nv;
            }
        }
    }
    (ns - infqr, jw - ns)
}

/// zlaqr5: one multishift QR sweep with ns shifts s on the active block ktop..kbot.
#[allow(clippy::too_many_arguments)]
fn zlaqr5(wantt: bool, wantz: bool, kacc22: usize, n: usize, ktop: usize, kbot: usize, nshfts: usize, s: &[C], h: &mut [C], hv: V,
          iloz: usize, ihiz: usize, zv: V, vbuf: &mut [C], uv: V, nv: usize, wvv: V, nh: usize, whv: V) {
    if nshfts < 2 || ktop >= kbot { return; }
    let ix = |i: usize, j: usize| hv.i(i, j);
    let zx = |i: usize, j: usize| zv.i(i, j);
    let ns = nshfts - nshfts % 2;
    let smlnum = SAFMIN * (n as f64 / ULP);
    let accum = kacc22 == 1 || kacc22 == 2;
    let blk22 = ns > 2 && kacc22 == 2;
    if ktop + 2 <= kbot { h[ix(ktop + 2, ktop)] = Z; }
    let nbmps = ns / 2;
    let kdu = 6 * nbmps - 3;
    let vw = V { o: 0, ld: 3 };
    let vi = |r: usize, m: usize| vw.i(r, m);
    let ui = |i: usize, j: usize| uv.i(i, j);
    let (kt, kb) = (ktop as isize, kbot as isize);
    let mut incol = 3 * (1 - nbmps as isize) + kt - 1;
    while incol <= kb - 2 {
        let ndcol = incol + kdu as isize;
        if accum { laset(h, uv, kdu, kdu, Z, ONE); }
        for krcol in incol..=(incol + 3 * nbmps as isize - 3).min(kb - 2) {
            let mtop = 1.max((kt - 1 - krcol + 2) / 3 + 1) as usize;
            let mbot = (nbmps as isize).min((kb - krcol) / 3) as usize;
            let m22 = mbot + 1;
            let bmp22 = mbot < nbmps && krcol + 3 * (m22 as isize - 1) == kb - 2;
            // bulges: chase / introduce
            for m in mtop..=mbot {
                let k = (krcol + 3 * (m as isize - 1)) as usize; // k >= ktop - 1 >= 0
                if k + 1 == ktop {
                    let v = zlaqr1(3, h, hv.at(ktop, ktop), s[2 * m - 2], s[2 * m - 1]);
                    let mut alpha = v[0];
                    let mut x = [v[1], v[2]];
                    let tau = zlarfg(&mut alpha, &mut x);
                    vbuf[vi(1, m)] = tau; vbuf[vi(2, m)] = x[0]; vbuf[vi(3, m)] = x[1];
                } else {
                    let mut beta = h[ix(k + 1, k)];
                    let mut x = [h[ix(k + 2, k)], h[ix(k + 3, k)]];
                    let tau = zlarfg(&mut beta, &mut x);
                    vbuf[vi(1, m)] = tau; vbuf[vi(2, m)] = x[0]; vbuf[vi(3, m)] = x[1];
                    if nz(h[ix(k + 3, k)]) || nz(h[ix(k + 3, k + 1)]) || !nz(h[ix(k + 3, k + 2)]) {
                        h[ix(k + 1, k)] = beta; h[ix(k + 2, k)] = Z; h[ix(k + 3, k)] = Z;
                    } else {
                        // a bulge may collapse: try a fresh one
                        let vt0 = zlaqr1(3, h, hv.at(k + 1, k + 1), s[2 * m - 2], s[2 * m - 1]);
                        let mut alpha = vt0[0];
                        let mut x = [vt0[1], vt0[2]];
                        let t0 = zlarfg(&mut alpha, &mut x);
                        let vt = [t0, x[0], x[1]];
                        let refsum = mul(cj(vt[0]), add(h[ix(k + 1, k)], mul(cj(vt[1]), h[ix(k + 2, k)])));
                        let z1 = sub(h[ix(k + 2, k)], mul(refsum, vt[1]));
                        let z4 = mul(refsum, vt[2]);
                        if abs1(z1) + abs1(z4) > ULP * (abs1(h[ix(k, k)]) + abs1(h[ix(k + 1, k + 1)]) + abs1(h[ix(k + 2, k + 2)])) {
                            h[ix(k + 1, k)] = beta; h[ix(k + 2, k)] = Z; h[ix(k + 3, k)] = Z;
                        } else {
                            h[ix(k + 1, k)] = sub(h[ix(k + 1, k)], refsum);
                            h[ix(k + 2, k)] = Z; h[ix(k + 3, k)] = Z;
                            vbuf[vi(1, m)] = vt[0]; vbuf[vi(2, m)] = vt[1]; vbuf[vi(3, m)] = vt[2];
                        }
                    }
                }
            }
            let k22 = krcol + 3 * (m22 as isize - 1);
            if bmp22 {
                let k = k22 as usize;
                if k + 1 == ktop {
                    let v = zlaqr1(2, h, hv.at(k + 1, k + 1), s[2 * m22 - 2], s[2 * m22 - 1]);
                    let mut beta = v[0];
                    let mut x = [v[1]];
                    let tau = zlarfg(&mut beta, &mut x);
                    vbuf[vi(1, m22)] = tau; vbuf[vi(2, m22)] = x[0];
                } else {
                    let mut beta = h[ix(k + 1, k)];
                    let mut x = [h[ix(k + 2, k)]];
                    let tau = zlarfg(&mut beta, &mut x);
                    vbuf[vi(1, m22)] = tau; vbuf[vi(2, m22)] = x[0];
                    h[ix(k + 1, k)] = beta; h[ix(k + 2, k)] = Z;
                }
            }
            // apply reflections from the left
            let jbot = if accum { (ndcol.min(kb)) as usize } else if wantt { n } else { kbot };
            for j in (ktop as isize).max(krcol) as usize..=jbot {
                let mend = (mbot as isize).min((j as isize - krcol + 2) / 3);
                if mend < mtop as isize { continue; }
                for m in mtop..=mend as usize {
                    let k = (krcol + 3 * (m as isize - 1)) as usize;
                    let refsum = mul(cj(vbuf[vi(1, m)]), add(add(h[ix(k + 1, j)], mul(cj(vbuf[vi(2, m)]), h[ix(k + 2, j)])), mul(cj(vbuf[vi(3, m)]), h[ix(k + 3, j)])));
                    h[ix(k + 1, j)] = sub(h[ix(k + 1, j)], refsum);
                    h[ix(k + 2, j)] = sub(h[ix(k + 2, j)], mul(refsum, vbuf[vi(2, m)]));
                    h[ix(k + 3, j)] = sub(h[ix(k + 3, j)], mul(refsum, vbuf[vi(3, m)]));
                }
            }
            if bmp22 {
                let k = k22 as usize;
                for j in (k + 1).max(ktop)..=jbot {
                    let refsum = mul(cj(vbuf[vi(1, m22)]), add(h[ix(k + 1, j)], mul(cj(vbuf[vi(2, m22)]), h[ix(k + 2, j)])));
                    h[ix(k + 1, j)] = sub(h[ix(k + 1, j)], refsum);
                    h[ix(k + 2, j)] = sub(h[ix(k + 2, j)], mul(refsum, vbuf[vi(2, m22)]));
                }
            }
            // apply reflections from the right
            let jtop = if accum { (kt.max(incol)) as usize } else if wantt { 1 } else { ktop };
            for m in mtop..=mbot {
                if !nz(vbuf[vi(1, m)]) { continue; }
                let k = (krcol + 3 * (m as isize - 1)) as usize;
                let (v1, v2, v3) = (vbuf[vi(1, m)], vbuf[vi(2, m)], vbuf[vi(3, m)]);
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = mul(v1, add(add(h[ix(j, k + 1)], mul(v2, h[ix(j, k + 2)])), mul(v3, h[ix(j, k + 3)])));
                    h[ix(j, k + 1)] = sub(h[ix(j, k + 1)], refsum);
                    h[ix(j, k + 2)] = sub(h[ix(j, k + 2)], mul(refsum, cj(v2)));
                    h[ix(j, k + 3)] = sub(h[ix(j, k + 3)], mul(refsum, cj(v3)));
                }
                if accum {
                    let kms = (k as isize - incol) as usize;
                    for j in (1.max(kt - incol)) as usize..=kdu {
                        let refsum = mul(v1, add(add(h[ui(j, kms + 1)], mul(v2, h[ui(j, kms + 2)])), mul(v3, h[ui(j, kms + 3)])));
                        h[ui(j, kms + 1)] = sub(h[ui(j, kms + 1)], refsum);
                        h[ui(j, kms + 2)] = sub(h[ui(j, kms + 2)], mul(refsum, cj(v2)));
                        h[ui(j, kms + 3)] = sub(h[ui(j, kms + 3)], mul(refsum, cj(v3)));
                    }
                } else if wantz {
                    for j in iloz..=ihiz {
                        let refsum = mul(v1, add(add(h[zx(j, k + 1)], mul(v2, h[zx(j, k + 2)])), mul(v3, h[zx(j, k + 3)])));
                        h[zx(j, k + 1)] = sub(h[zx(j, k + 1)], refsum);
                        h[zx(j, k + 2)] = sub(h[zx(j, k + 2)], mul(refsum, cj(v2)));
                        h[zx(j, k + 3)] = sub(h[zx(j, k + 3)], mul(refsum, cj(v3)));
                    }
                }
            }
            if bmp22 && nz(vbuf[vi(1, m22)]) {
                let k = k22 as usize;
                let (v1, v2) = (vbuf[vi(1, m22)], vbuf[vi(2, m22)]);
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = mul(v1, add(h[ix(j, k + 1)], mul(v2, h[ix(j, k + 2)])));
                    h[ix(j, k + 1)] = sub(h[ix(j, k + 1)], refsum);
                    h[ix(j, k + 2)] = sub(h[ix(j, k + 2)], mul(refsum, cj(v2)));
                }
                if accum {
                    let kms = (k as isize - incol) as usize;
                    for j in (1.max(kt - incol)) as usize..=kdu {
                        let refsum = mul(v1, add(h[ui(j, kms + 1)], mul(v2, h[ui(j, kms + 2)])));
                        h[ui(j, kms + 1)] = sub(h[ui(j, kms + 1)], refsum);
                        h[ui(j, kms + 2)] = sub(h[ui(j, kms + 2)], mul(refsum, cj(v2)));
                    }
                } else if wantz {
                    for j in iloz..=ihiz {
                        let refsum = mul(v1, add(h[zx(j, k + 1)], mul(v2, h[zx(j, k + 2)])));
                        h[zx(j, k + 1)] = sub(h[zx(j, k + 1)], refsum);
                        h[zx(j, k + 2)] = sub(h[zx(j, k + 2)], mul(refsum, cj(v2)));
                    }
                }
            }
            // vigilant deflation check
            let mut mstart = mtop;
            if krcol + 3 * (mstart as isize - 1) < kt { mstart += 1; }
            let mut mend = mbot;
            if bmp22 { mend += 1; }
            if krcol == kb - 2 { mend += 1; }
            for m in mstart..=mend {
                let k = (kb - 1).min(krcol + 3 * (m as isize - 1)) as usize;
                if nz(h[ix(k + 1, k)]) {
                    let mut tst1 = abs1(h[ix(k, k)]) + abs1(h[ix(k + 1, k + 1)]);
                    if tst1 == 0.0 {
                        if k > ktop { tst1 += abs1(h[ix(k, k - 1)]); }
                        if k >= ktop + 2 { tst1 += abs1(h[ix(k, k - 2)]); }
                        if k >= ktop + 3 { tst1 += abs1(h[ix(k, k - 3)]); }
                        if k + 2 <= kbot { tst1 += abs1(h[ix(k + 2, k + 1)]); }
                        if k + 3 <= kbot { tst1 += abs1(h[ix(k + 3, k + 1)]); }
                        if k + 4 <= kbot { tst1 += abs1(h[ix(k + 4, k + 1)]); }
                    }
                    if abs1(h[ix(k + 1, k)]) <= smlnum.max(ULP * tst1) {
                        let (p, q) = (abs1(h[ix(k + 1, k)]), abs1(h[ix(k, k + 1)]));
                        let (h12, h21) = (p.max(q), p.min(q));
                        let zz = sub(h[ix(k, k)], h[ix(k + 1, k + 1)]);
                        let (p, q) = (abs1(h[ix(k + 1, k + 1)]), abs1(zz));
                        let (h11, h22) = (p.max(q), p.min(q));
                        let scl_ = h11 + h12;
                        let tst2 = h22 * (h11 / scl_);
                        if tst2 == 0.0 || h21 * (h12 / scl_) <= smlnum.max(ULP * tst2) { h[ix(k + 1, k)] = Z; }
                    }
                }
            }
            // fill in the last row of each bulge
            let mend = (nbmps as isize).min((kb - krcol - 1) / 3);
            if mend >= mtop as isize {
                for m in mtop..=mend as usize {
                    let k = (krcol + 3 * (m as isize - 1)) as usize;
                    let refsum = mul(mul(vbuf[vi(1, m)], vbuf[vi(3, m)]), h[ix(k + 4, k + 3)]);
                    h[ix(k + 4, k + 1)] = neg(refsum);
                    h[ix(k + 4, k + 2)] = mul(neg(refsum), cj(vbuf[vi(2, m)]));
                    h[ix(k + 4, k + 3)] = sub(h[ix(k + 4, k + 3)], mul(refsum, cj(vbuf[vi(3, m)])));
                }
            }
        }
        // use U to update far-from-diagonal entries
        if accum {
            let (jtop, jbot) = if wantt { (1, n) } else { (ktop, kbot) };
            if !blk22 || incol < kt || ndcol > kb || ns <= 2 {
                let k1 = 1.max(kt - incol) as usize;
                let nu = kdu - 0.max(ndcol - kb) as usize - k1 + 1;
                let mut jcol = (ndcol.min(kb) + 1) as usize;
                while jcol <= jbot {
                    let jlen = nh.min(jbot - jcol + 1);
                    let (a, av) = take(h, uv.at(k1, k1), nu, nu);
                    let r0 = (incol + k1 as isize) as usize;
                    let (b, bv) = take(h, hv.at(r0, jcol), nu, jlen);
                    zgemm(b'C', b'N', nu, jlen, nu, ONE, &a, av, &b, bv, Z, h, whv);
                    lacpy(h, whv, hv.at(r0, jcol), nu, jlen);
                    jcol += nh;
                }
                let top = (kt.max(incol) - 1) as usize;
                let mut jrow = jtop;
                while jrow <= top {
                    let jlen = nv.min(kt.max(incol) as usize - jrow);
                    let c0 = (incol + k1 as isize) as usize;
                    let (a, av) = take(h, hv.at(jrow, c0), jlen, nu);
                    let (b, bv) = take(h, uv.at(k1, k1), nu, nu);
                    zgemm(b'N', b'N', jlen, nu, nu, ONE, &a, av, &b, bv, Z, h, wvv);
                    lacpy(h, wvv, hv.at(jrow, c0), jlen, nu);
                    jrow += nv;
                }
                if wantz {
                    let mut jrow = iloz;
                    while jrow <= ihiz {
                        let jlen = nv.min(ihiz - jrow + 1);
                        let c0 = (incol + k1 as isize) as usize;
                        let (a, av) = take(h, zv.at(jrow, c0), jlen, nu);
                        let (b, bv) = take(h, uv.at(k1, k1), nu, nu);
                        zgemm(b'N', b'N', jlen, nu, nu, ONE, &a, av, &b, bv, Z, h, wvv);
                        lacpy(h, wvv, zv.at(jrow, c0), jlen, nu);
                        jrow += nv;
                    }
                }
            } else {
                // U's 2 x 2 block structure
                let incol = incol as usize;
                let i2 = (kdu + 1) / 2;
                let i4 = kdu;
                let j2 = i4 - i2;
                let j4 = kdu;
                let kzs = j4 - j2 - (ns + 1);
                let knz = ns + 1;
                let mut jcol = (ndcol.min(kb) + 1) as usize;
                while jcol <= jbot {
                    let jlen = nh.min(jbot - jcol + 1);
                    lacpy(h, hv.at(incol + 1 + j2, jcol), whv.at(kzs + 1, 1), knz, jlen);
                    laset(h, whv, kzs, jlen, Z, Z);
                    let (u, uvv) = take(h, uv.at(j2 + 1, kzs + 1), knz, knz);
                    ztrmm(true, true, b'C', false, knz, jlen, ONE, &u, uvv, h, whv.at(kzs + 1, 1));
                    let (a, av) = take(h, uv, j2, i2);
                    let (b, bv) = take(h, hv.at(incol + 1, jcol), j2, jlen);
                    zgemm(b'C', b'N', i2, jlen, j2, ONE, &a, av, &b, bv, ONE, h, whv);
                    lacpy(h, hv.at(incol + 1, jcol), whv.at(i2 + 1, 1), j2, jlen);
                    let (u, uvv) = take(h, uv.at(1, i2 + 1), j2, j2);
                    ztrmm(true, false, b'C', false, j2, jlen, ONE, &u, uvv, h, whv.at(i2 + 1, 1));
                    let (a, av) = take(h, uv.at(j2 + 1, i2 + 1), j4 - j2, i4 - i2);
                    let (b, bv) = take(h, hv.at(incol + 1 + j2, jcol), j4 - j2, jlen);
                    zgemm(b'C', b'N', i4 - i2, jlen, j4 - j2, ONE, &a, av, &b, bv, ONE, h, whv.at(i2 + 1, 1));
                    lacpy(h, whv, hv.at(incol + 1, jcol), kdu, jlen);
                    jcol += nh;
                }
                let rows = |h: &mut [C], jrow: usize, jlen: usize, mv: V| {
                    lacpy(h, mv.at(jrow, incol + 1 + j2), wvv.at(1, kzs + 1), jlen, knz);
                    laset(h, wvv, jlen, kzs, Z, Z);
                    let (u, uvv) = take(h, uv.at(j2 + 1, kzs + 1), knz, knz);
                    ztrmm(false, true, b'N', false, jlen, knz, ONE, &u, uvv, h, wvv.at(1, kzs + 1));
                    let (a, av) = take(h, mv.at(jrow, incol + 1), jlen, j2);
                    let (b, bv) = take(h, uv, j2, i2);
                    zgemm(b'N', b'N', jlen, i2, j2, ONE, &a, av, &b, bv, ONE, h, wvv);
                    lacpy(h, mv.at(jrow, incol + 1), wvv.at(1, i2 + 1), jlen, j2);
                    let (u, uvv) = take(h, uv.at(1, i2 + 1), i4 - i2, i4 - i2);
                    ztrmm(false, false, b'N', false, jlen, i4 - i2, ONE, &u, uvv, h, wvv.at(1, i2 + 1));
                    let (a, av) = take(h, mv.at(jrow, incol + 1 + j2), jlen, j4 - j2);
                    let (b, bv) = take(h, uv.at(j2 + 1, i2 + 1), j4 - j2, i4 - i2);
                    zgemm(b'N', b'N', jlen, i4 - i2, j4 - j2, ONE, &a, av, &b, bv, ONE, h, wvv.at(1, i2 + 1));
                    lacpy(h, wvv, mv.at(jrow, incol + 1), jlen, kdu);
                };
                let top = incol.max(ktop) - 1;
                let mut jrow = jtop;
                while jrow <= top {
                    let jlen = nv.min(incol.max(ktop) - jrow);
                    rows(h, jrow, jlen, hv);
                    jrow += nv;
                }
                if wantz {
                    let mut jrow = iloz;
                    while jrow <= ihiz {
                        let jlen = nv.min(ihiz - jrow + 1);
                        rows(h, jrow, jlen, zv);
                        jrow += nv;
                    }
                }
            }
        }
        incol += 3 * nbmps as isize - 2;
    }
}

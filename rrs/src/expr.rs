//! Rate-law expressions: infix parser (libsbml L3 syntax) -> hash-consed DAG with constant
//! folding -> symbolic derivatives -> one flat register tape that evaluates every rate and every
//! nonzero Jacobian entry in a single allocation-free pass.
use crate::jit::{compile, Jit};
use std::collections::HashMap;
use std::sync::{atomic::{AtomicU32, Ordering}, OnceLock};

pub type N = u32;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Op {
    Const(u64), // f64 bits
    Var(u32),   // species concentration
    Par(u32),   // parameter
    Time,
    Neg(N), Add(N, N), Sub(N, N), Mul(N, N), Div(N, N), Pow(N, N), PowI(N, i32),
    Exp(N), Ln(N), Log10(N), Sqrt(N), Abs(N), Sign(N), Floor(N), Ceil(N), Sin(N), Cos(N), Tan(N),
    Lt(N, N), Le(N, N), Eq(N, N), Ne(N, N), And(N, N), Or(N, N), Xor(N, N), Not(N),
    Sel(N, N, N), Min(N, N), Max(N, N),
}

#[inline(always)]
fn b(x: bool) -> f64 { if x { 1.0 } else { 0.0 } }

/// Apply a non-leaf op to already computed operand values.
#[inline(always)]
fn apply(op: &Op, r: impl Fn(N) -> f64) -> f64 {
    use Op::*;
    match *op {
        Neg(a) => -r(a), Add(a, c) => r(a) + r(c), Sub(a, c) => r(a) - r(c),
        Mul(a, c) => r(a) * r(c), Div(a, c) => r(a) / r(c), Pow(a, c) => r(a).powf(r(c)),
        PowI(a, k) => r(a).powi(k), Exp(a) => r(a).exp(), Ln(a) => r(a).ln(),
        Log10(a) => r(a).log10(), Sqrt(a) => r(a).sqrt(), Abs(a) => r(a).abs(),
        Sign(a) => { let x = r(a); if x > 0.0 { 1.0 } else if x < 0.0 { -1.0 } else { 0.0 } }
        Floor(a) => r(a).floor(), Ceil(a) => r(a).ceil(),
        Sin(a) => r(a).sin(), Cos(a) => r(a).cos(), Tan(a) => r(a).tan(),
        Lt(a, c) => b(r(a) < r(c)), Le(a, c) => b(r(a) <= r(c)),
        Eq(a, c) => b(r(a) == r(c)), Ne(a, c) => b(r(a) != r(c)),
        And(a, c) => b(r(a) != 0.0 && r(c) != 0.0), Or(a, c) => b(r(a) != 0.0 || r(c) != 0.0),
        Xor(a, c) => b((r(a) != 0.0) != (r(c) != 0.0)), Not(a) => b(r(a) == 0.0),
        Sel(c, x, y) => if r(c) != 0.0 { r(x) } else { r(y) },
        Min(a, c) => r(a).min(r(c)), Max(a, c) => r(a).max(r(c)),
        Const(_) | Var(_) | Par(_) | Time => unreachable!(),
    }
}

fn children(op: &Op) -> [Option<N>; 3] {
    use Op::*;
    match *op {
        Const(_) | Var(_) | Par(_) | Time => [None, None, None],
        Neg(a) | PowI(a, _) | Exp(a) | Ln(a) | Log10(a) | Sqrt(a) | Abs(a) | Sign(a) | Floor(a)
        | Ceil(a) | Sin(a) | Cos(a) | Tan(a) | Not(a) => [Some(a), None, None],
        Sel(a, c, d) => [Some(a), Some(c), Some(d)],
        Add(a, c) | Sub(a, c) | Mul(a, c) | Div(a, c) | Pow(a, c) | Lt(a, c) | Le(a, c) | Eq(a, c)
        | Ne(a, c) | And(a, c) | Or(a, c) | Xor(a, c) | Min(a, c) | Max(a, c) => [Some(a), Some(c), None],
    }
}

#[derive(Default)]
pub struct Dag {
    pub nodes: Vec<Op>,
    memo: HashMap<Op, N>,
    dmemo: HashMap<(N, Op), N>,
}

impl Dag {
    pub fn cst(&mut self, x: f64) -> N { self.raw(Op::Const(x.to_bits())) }
    fn raw(&mut self, op: Op) -> N {
        if let Some(&i) = self.memo.get(&op) { return i; }
        let i = self.nodes.len() as N;
        self.nodes.push(op);
        self.memo.insert(op, i);
        i
    }
    pub fn val(&self, n: N) -> Option<f64> {
        if let Op::Const(x) = self.nodes[n as usize] { Some(f64::from_bits(x)) } else { None }
    }
    fn is(&self, n: N, x: f64) -> bool { self.val(n) == Some(x) }

    /// Build a node with algebraic simplification and constant folding.
    pub fn mk(&mut self, op: Op) -> N {
        use Op::*;
        let ch = children(&op);
        if ch[0].is_some() && ch.iter().flatten().all(|&c| self.val(c).is_some()) {
            let v = apply(&op, |c| self.val(c).unwrap());
            return self.cst(v);
        }
        match op {
            Add(a, c) if self.is(a, 0.0) => return c,
            Add(a, c) if self.is(c, 0.0) => return a,
            Sub(a, c) if self.is(c, 0.0) => return a,
            Sub(a, c) if self.is(a, 0.0) => return self.mk(Neg(c)),
            Sub(a, c) if a == c => return self.cst(0.0),
            Mul(a, c) if self.is(a, 0.0) || self.is(c, 0.0) => return self.cst(0.0),
            Mul(a, c) if self.is(a, 1.0) => return c,
            Mul(a, c) if self.is(c, 1.0) => return a,
            Mul(a, c) if self.is(a, -1.0) => return self.mk(Neg(c)),
            Mul(a, c) if self.is(c, -1.0) => return self.mk(Neg(a)),
            Div(a, _) if self.is(a, 0.0) => return self.cst(0.0),
            Div(a, c) if self.is(c, 1.0) => return a,
            Neg(a) => if let Neg(x) = self.nodes[a as usize] { return x },
            Pow(a, c) => match self.val(c) {
                Some(e) if e == 1.0 => return a,
                Some(e) if e == 0.0 => return self.cst(1.0),
                Some(e) if e == 0.5 => return self.mk(Sqrt(a)),
                Some(e) if e.fract() == 0.0 && e.abs() <= 64.0 => return self.raw(PowI(a, e as i32)),
                _ => {}
            },
            PowI(a, 1) => return a,
            PowI(_, 0) => return self.cst(1.0),
            Sel(c, x, y) if x == y => { let _ = c; return x }
            Sel(c, x, y) => if let Some(cv) = self.val(c) { return if cv != 0.0 { x } else { y } },
            _ => {}
        }
        self.raw(op)
    }

    /// d node / d wrt, where wrt is Var(k) or Time. Memoized.
    pub fn diff(&mut self, n: N, wrt: Op) -> N {
        if let Some(&d) = self.dmemo.get(&(n, wrt)) { return d; }
        use Op::*;
        let op = self.nodes[n as usize];
        let zero = self.cst(0.0);
        let d = match op {
            Const(_) | Par(_) => zero,
            Var(_) | Time => self.cst(if op == wrt { 1.0 } else { 0.0 }),
            Neg(a) => { let da = self.diff(a, wrt); self.mk(Neg(da)) }
            Add(a, c) => { let (x, y) = (self.diff(a, wrt), self.diff(c, wrt)); self.mk(Add(x, y)) }
            Sub(a, c) => { let (x, y) = (self.diff(a, wrt), self.diff(c, wrt)); self.mk(Sub(x, y)) }
            Mul(a, c) => {
                let (x, y) = (self.diff(a, wrt), self.diff(c, wrt));
                let (p, q) = (self.mk(Mul(x, c)), self.mk(Mul(a, y)));
                self.mk(Add(p, q))
            }
            Div(a, c) => { // (da - (a/c) dc) / c
                let (x, y) = (self.diff(a, wrt), self.diff(c, wrt));
                let t = self.mk(Mul(n, y)); let u = self.mk(Sub(x, t)); self.mk(Div(u, c))
            }
            PowI(a, k) => { // k a^(k-1) da
                let da = self.diff(a, wrt);
                let p = self.mk(PowI(a, k - 1)); let kc = self.cst(k as f64);
                let t = self.mk(Mul(kc, p)); self.mk(Mul(t, da))
            }
            Pow(a, c) => {
                let (x, y) = (self.diff(a, wrt), self.diff(c, wrt));
                if self.is(y, 0.0) { // c a^(c-1) da
                    let one = self.cst(1.0); let cm = self.mk(Sub(c, one));
                    let p = self.mk(Pow(a, cm)); let t = self.mk(Mul(c, p)); self.mk(Mul(t, x))
                } else { // a^c (dc ln a + c da / a)
                    let la = self.mk(Ln(a)); let t1 = self.mk(Mul(y, la));
                    let q = self.mk(Div(x, a)); let t2 = self.mk(Mul(c, q));
                    let s = self.mk(Add(t1, t2)); self.mk(Mul(n, s))
                }
            }
            Exp(a) => { let da = self.diff(a, wrt); self.mk(Mul(n, da)) }
            Ln(a) => { let da = self.diff(a, wrt); self.mk(Div(da, a)) }
            Log10(a) => {
                let da = self.diff(a, wrt); let l = self.cst(std::f64::consts::LN_10);
                let den = self.mk(Mul(a, l)); self.mk(Div(da, den))
            }
            Sqrt(a) => {
                let da = self.diff(a, wrt); let two = self.cst(2.0);
                let den = self.mk(Mul(two, n)); self.mk(Div(da, den))
            }
            Abs(a) => { let da = self.diff(a, wrt); let s = self.mk(Sign(a)); self.mk(Mul(s, da)) }
            Sin(a) => { let da = self.diff(a, wrt); let c = self.mk(Cos(a)); self.mk(Mul(c, da)) }
            Cos(a) => {
                let da = self.diff(a, wrt); let s = self.mk(Sin(a));
                let t = self.mk(Mul(s, da)); self.mk(Neg(t))
            }
            Tan(a) => {
                let da = self.diff(a, wrt); let one = self.cst(1.0);
                let t2 = self.mk(PowI(n, 2)); let s = self.mk(Add(one, t2)); self.mk(Mul(s, da))
            }
            Sel(c, x, y) => {
                let (dx, dy) = (self.diff(x, wrt), self.diff(y, wrt)); self.mk(Sel(c, dx, dy))
            }
            Min(a, c) => {
                let (x, y) = (self.diff(a, wrt), self.diff(c, wrt));
                let cond = self.mk(Le(a, c)); self.mk(Sel(cond, x, y))
            }
            Max(a, c) => {
                let (x, y) = (self.diff(a, wrt), self.diff(c, wrt));
                let cond = self.mk(Le(c, a)); self.mk(Sel(cond, x, y))
            }
            Sign(_) | Floor(_) | Ceil(_) | Lt(..) | Le(..) | Eq(..) | Ne(..) | And(..) | Or(..)
            | Xor(..) | Not(_) => zero,
        };
        self.dmemo.insert((n, wrt), d);
        d
    }

    /// Species (Var indices) a node depends on, sorted.
    pub fn deps(&self, n: N, has_time: &mut bool) -> Vec<u32> {
        let mut seen = vec![false; n as usize + 1];
        let mut st = vec![n];
        let mut out = vec![];
        while let Some(i) = st.pop() {
            if std::mem::replace(&mut seen[i as usize], true) { continue; }
            match self.nodes[i as usize] {
                Op::Var(k) => out.push(k),
                Op::Time => *has_time = true,
                op => st.extend(children(&op).iter().flatten()),
            }
        }
        out.sort();
        out
    }
}

// ---------------------------------------------------------------- parser

#[derive(Clone, Debug, PartialEq)]
enum Tok { Num(f64), Id(String), Sym(&'static str) }

fn lex(s: &str) -> Result<Vec<Tok>, String> {
    let c: Vec<char> = s.chars().collect();
    let (mut i, mut out) = (0, vec![]);
    const SYMS: [&str; 17] = ["<=", ">=", "==", "!=", "&&", "||", "<", ">", "+", "-", "*", "/", "^", "(", ")", ",", "!"];
    while i < c.len() {
        let ch = c[i];
        if ch.is_whitespace() { i += 1; continue; }
        if ch.is_ascii_digit() || (ch == '.' && i + 1 < c.len() && c[i + 1].is_ascii_digit()) {
            let st = i;
            while i < c.len() && (c[i].is_ascii_digit() || c[i] == '.') { i += 1; }
            if i < c.len() && (c[i] == 'e' || c[i] == 'E') {
                let mut j = i + 1;
                if j < c.len() && (c[j] == '+' || c[j] == '-') { j += 1; }
                if j < c.len() && c[j].is_ascii_digit() {
                    i = j;
                    while i < c.len() && c[i].is_ascii_digit() { i += 1; }
                }
            }
            let t: String = c[st..i].iter().collect();
            out.push(Tok::Num(t.parse().map_err(|_| format!("bad number '{t}'"))?));
        } else if ch.is_alphabetic() || ch == '_' {
            let st = i;
            while i < c.len() && (c[i].is_alphanumeric() || c[i] == '_') { i += 1; }
            out.push(Tok::Id(c[st..i].iter().collect()));
        } else {
            let rest: String = c[i..(i + 2).min(c.len())].iter().collect();
            let sym = SYMS.iter().find(|s| rest.starts_with(**s)).ok_or(format!("unexpected '{ch}' in '{s}'"))?;
            out.push(Tok::Sym(sym));
            i += sym.len();
        }
    }
    Ok(out)
}

struct P<'a, F: FnMut(&str) -> Option<Op>> { t: Vec<Tok>, i: usize, d: &'a mut Dag, res: F }

impl<'a, F: FnMut(&str) -> Option<Op>> P<'a, F> {
    fn peek(&self, s: &str) -> bool { matches!(self.t.get(self.i), Some(Tok::Sym(x)) if *x == s) }
    fn eat(&mut self, s: &str) -> bool { let m = self.peek(s); if m { self.i += 1; } m }
    fn expect(&mut self, s: &str) -> Result<(), String> {
        if self.eat(s) { Ok(()) } else { Err(format!("expected '{s}' at token {}", self.i)) }
    }
    fn or(&mut self) -> Result<N, String> {
        let mut a = self.and()?;
        while self.eat("||") { let c = self.and()?; a = self.d.mk(Op::Or(a, c)); }
        Ok(a)
    }
    fn and(&mut self) -> Result<N, String> {
        let mut a = self.cmp()?;
        while self.eat("&&") { let c = self.cmp()?; a = self.d.mk(Op::And(a, c)); }
        Ok(a)
    }
    fn cmp(&mut self) -> Result<N, String> {
        let mut a = self.sum()?;
        loop {
            let op = ["<=", ">=", "==", "!=", "<", ">"].into_iter().find(|s| self.peek(s));
            let Some(op) = op else { return Ok(a) };
            self.i += 1;
            let c = self.sum()?;
            a = self.d.mk(rel(op, a, c));
        }
    }
    fn sum(&mut self) -> Result<N, String> {
        let mut a = self.prod()?;
        loop {
            if self.eat("+") { let c = self.prod()?; a = self.d.mk(Op::Add(a, c)); }
            else if self.eat("-") { let c = self.prod()?; a = self.d.mk(Op::Sub(a, c)); }
            else { return Ok(a); }
        }
    }
    fn prod(&mut self) -> Result<N, String> {
        let mut a = self.unary()?;
        loop {
            if self.eat("*") { let c = self.unary()?; a = self.d.mk(Op::Mul(a, c)); }
            else if self.eat("/") { let c = self.unary()?; a = self.d.mk(Op::Div(a, c)); }
            else { return Ok(a); }
        }
    }
    fn unary(&mut self) -> Result<N, String> {
        if self.eat("-") { let a = self.unary()?; return Ok(self.d.mk(Op::Neg(a))); }
        if self.eat("+") { return self.unary(); }
        if self.eat("!") { let a = self.unary()?; return Ok(self.d.mk(Op::Not(a))); }
        let a = self.atom()?;
        if self.eat("^") { let e = self.unary()?; return Ok(self.d.mk(Op::Pow(a, e))); }
        Ok(a)
    }
    fn atom(&mut self) -> Result<N, String> {
        let tok = self.t.get(self.i).cloned().ok_or("unexpected end of formula")?;
        self.i += 1;
        match tok {
            Tok::Num(x) => Ok(self.d.cst(x)),
            Tok::Sym("(") => { let a = self.or()?; self.expect(")")?; Ok(a) }
            Tok::Id(name) if self.peek("(") => {
                self.i += 1;
                let mut args = vec![];
                if !self.eat(")") {
                    loop { args.push(self.or()?); if self.eat(")") { break; } self.expect(",")?; }
                }
                self.call(&name, args)
            }
            Tok::Id(name) => {
                let c = match name.as_str() {
                    "true" => 1.0, "false" => 0.0, "pi" => std::f64::consts::PI,
                    "exponentiale" => std::f64::consts::E, "infinity" | "INF" | "inf" => f64::INFINITY,
                    "NaN" | "nan" => f64::NAN,
                    _ => {
                        let op = (self.res)(&name).ok_or(format!("unknown symbol '{name}'"))?;
                        return Ok(self.d.mk(op));
                    }
                };
                Ok(self.d.cst(c))
            }
            Tok::Sym(s) => Err(format!("unexpected '{s}'")),
        }
    }
    fn call(&mut self, f: &str, a: Vec<N>) -> Result<N, String> {
        use Op::*;
        let n = a.len();
        let need = |k: usize| if n == k { Ok(()) } else { Err(format!("{f}() takes {k} argument(s), got {n}")) };
        let d = &mut *self.d;
        let fold = |d: &mut Dag, mk: fn(N, N) -> Op| -> Result<N, String> {
            if a.is_empty() { return Err(format!("{f}() needs arguments")); }
            Ok(a[1..].iter().fold(a[0], |x, &y| d.mk(mk(x, y))))
        };
        Ok(match f {
            "pow" | "power" => { need(2)?; d.mk(Pow(a[0], a[1])) }
            "exp" => { need(1)?; d.mk(Exp(a[0])) }
            "ln" => { need(1)?; d.mk(Ln(a[0])) }
            "log10" => { need(1)?; d.mk(Log10(a[0])) }
            "log" if n == 1 => d.mk(Log10(a[0])),
            "log" => { need(2)?; let x = d.mk(Ln(a[1])); let y = d.mk(Ln(a[0])); d.mk(Div(x, y)) }
            "sqrt" => { need(1)?; d.mk(Sqrt(a[0])) }
            "root" if n == 1 => d.mk(Sqrt(a[0])),
            "root" => { need(2)?; let one = d.cst(1.0); let e = d.mk(Div(one, a[0])); d.mk(Pow(a[1], e)) }
            "abs" => { need(1)?; d.mk(Abs(a[0])) }
            "floor" => { need(1)?; d.mk(Floor(a[0])) }
            "ceil" | "ceiling" => { need(1)?; d.mk(Ceil(a[0])) }
            "sin" => { need(1)?; d.mk(Sin(a[0])) }
            "cos" => { need(1)?; d.mk(Cos(a[0])) }
            "tan" => { need(1)?; d.mk(Tan(a[0])) }
            "min" => fold(d, Min)?,
            "max" => fold(d, Max)?,
            "plus" => fold(d, Add)?,
            "times" => fold(d, Mul)?,
            "minus" => { need(2)?; d.mk(Sub(a[0], a[1])) }
            "divide" => { need(2)?; d.mk(Div(a[0], a[1])) }
            "and" => fold(d, And)?,
            "or" => fold(d, Or)?,
            "xor" => fold(d, Xor)?,
            "not" => { need(1)?; d.mk(Not(a[0])) }
            "lt" | "gt" | "leq" | "geq" | "eq" | "neq" => {
                need(2)?;
                let op = match f { "lt" => "<", "gt" => ">", "leq" => "<=", "geq" => ">=", "eq" => "==", _ => "!=" };
                d.mk(rel(op, a[0], a[1]))
            }
            "piecewise" => {
                if n == 0 { return Err("piecewise() needs arguments".into()); }
                let mut acc = if n % 2 == 1 { a[n - 1] } else { d.cst(f64::NAN) };
                for k in (0..n / 2).rev() { acc = d.mk(Sel(a[2 * k + 1], a[2 * k], acc)); }
                acc
            }
            _ => return Err(format!("unsupported function '{f}'")),
        })
    }
}

fn rel(op: &str, a: N, c: N) -> Op {
    match op {
        "<" => Op::Lt(a, c), ">" => Op::Lt(c, a), "<=" => Op::Le(a, c), ">=" => Op::Le(c, a),
        "==" => Op::Eq(a, c), _ => Op::Ne(a, c),
    }
}

/// Parse an infix formula into `dag`; `res` maps identifiers to Var/Par/Time.
pub fn parse(src: &str, dag: &mut Dag, res: impl FnMut(&str) -> Option<Op>) -> Result<N, String> {
    let mut p = P { t: lex(src)?, i: 0, d: dag, res };
    let r = p.or()?;
    if p.i != p.t.len() { return Err(format!("trailing input in '{src}'")); }
    Ok(r)
}

// ---------------------------------------------------------------- tape

/// Register layout: [species | params | time | constants | computed].
pub struct Tape {
    pub ns: usize,
    pub np: usize,
    pub nregs: usize,
    pub consts: Vec<(u32, f64)>,
    pub ins: Vec<(u32, Op)>, // destination register, op with operands remapped to registers
    pub pre: Vec<(u32, Op)>, // instructions depending on parameters and constants only (see `prep`)
    pub outs: Vec<u32>,
    runs: AtomicU32,
    jit: OnceLock<Option<Jit>>,
}

/// `run` interprets the first JIT_AFTER calls, then compiles the tape to native code. Compiling
/// costs about as much as 2000 interpreted runs save, so short runs (and the many small
/// sensitivity models) never pay for it and long ones lose at most 2x the optimum.
const JIT_AFTER: u32 = 2000;

impl Tape {
    pub fn new(dag: &Dag, roots: &[N], ns: usize, np: usize) -> Tape {
        let mut need = vec![false; dag.nodes.len()];
        let mut st: Vec<N> = roots.to_vec();
        while let Some(i) = st.pop() {
            if std::mem::replace(&mut need[i as usize], true) { continue; }
            st.extend(children(&dag.nodes[i as usize]).iter().flatten());
        }
        let time = (ns + np) as u32;
        let mut reg = vec![u32::MAX; dag.nodes.len()];
        let mut next = time + 1;
        let mut consts = vec![];
        let mut ins = vec![];
        for (i, op) in dag.nodes.iter().enumerate() {
            if !need[i] { continue; }
            reg[i] = match *op {
                Op::Var(k) => k,
                Op::Par(k) => ns as u32 + k,
                Op::Time => time,
                Op::Const(x) => { consts.push((next, f64::from_bits(x))); next += 1; next - 1 }
                _ => { next += 1; next - 1 }
            };
        }
        let m = |n: N| reg[n as usize];
        // parameter/constant-only nodes are hoisted out of `run` into `prep`
        let mut stat = vec![false; dag.nodes.len()];
        let mut pre = vec![];
        for (i, op) in dag.nodes.iter().enumerate() {
            if !need[i] { continue; }
            stat[i] = match op {
                Op::Const(_) | Op::Par(_) => true,
                Op::Var(_) | Op::Time => false,
                _ => children(op).iter().flatten().all(|&c| stat[c as usize]),
            };
        }
        for (i, op) in dag.nodes.iter().enumerate() {
            if !need[i] { continue; }
            use Op::*;
            let r = match *op {
                Const(_) | Var(_) | Par(_) | Time => continue,
                Neg(a) => Neg(m(a)), Add(a, c) => Add(m(a), m(c)), Sub(a, c) => Sub(m(a), m(c)),
                Mul(a, c) => Mul(m(a), m(c)), Div(a, c) => Div(m(a), m(c)), Pow(a, c) => Pow(m(a), m(c)),
                PowI(a, k) => PowI(m(a), k), Exp(a) => Exp(m(a)), Ln(a) => Ln(m(a)), Log10(a) => Log10(m(a)),
                Sqrt(a) => Sqrt(m(a)), Abs(a) => Abs(m(a)), Sign(a) => Sign(m(a)), Floor(a) => Floor(m(a)),
                Ceil(a) => Ceil(m(a)), Sin(a) => Sin(m(a)), Cos(a) => Cos(m(a)), Tan(a) => Tan(m(a)),
                Lt(a, c) => Lt(m(a), m(c)), Le(a, c) => Le(m(a), m(c)), Eq(a, c) => Eq(m(a), m(c)),
                Ne(a, c) => Ne(m(a), m(c)), And(a, c) => And(m(a), m(c)), Or(a, c) => Or(m(a), m(c)),
                Xor(a, c) => Xor(m(a), m(c)), Not(a) => Not(m(a)), Sel(a, c, d) => Sel(m(a), m(c), m(d)),
                Min(a, c) => Min(m(a), m(c)), Max(a, c) => Max(m(a), m(c)),
            };
            if stat[i] { pre.push((reg[i], r)) } else { ins.push((reg[i], r)) }
        }
        Tape { ns, np, nregs: next as usize, consts, ins, pre, outs: roots.iter().map(|&r| m(r)).collect(), runs: AtomicU32::new(0), jit: OnceLock::new() }
    }

    /// Fresh register file with constants loaded.
    pub fn regs(&self) -> Vec<f64> {
        let mut r = vec![0.0; self.nregs];
        for &(i, x) in &self.consts { r[i as usize] = x; }
        r
    }

    /// Evaluate the parameter-only part; call after any parameter register changes.
    pub fn prep(&self, regs: &mut [f64]) { Self::exec(&self.pre, regs, self.nregs) }

    /// `prep` then `run`, for one-off evaluations.
    pub fn run_all(&self, regs: &mut [f64]) { self.prep(regs); self.run(regs) }

    /// Evaluate; species values must already be in regs[..ns], params in regs[ns..ns+np], and
    /// `prep` must have run since the parameters last changed.
    #[inline]
    pub fn run(&self, regs: &mut [f64]) {
        assert!(regs.len() >= self.nregs);
        match self.jit.get() {
            // SAFETY: the compiled code touches registers < nregs only (checked above).
            Some(Some(j)) => unsafe { (j.f)(regs.as_mut_ptr()) },
            Some(None) => Self::exec(&self.ins, regs, self.nregs),
            None => {
                if self.runs.fetch_add(1, Ordering::Relaxed) >= JIT_AFTER {
                    self.jit.get_or_init(|| compile(&self.ins, &self.outs, self.nregs));
                }
                Self::exec(&self.ins, regs, self.nregs)
            }
        }
    }

    #[inline]
    fn exec(ins: &[(u32, Op)], regs: &mut [f64], nregs: usize) {
        assert!(regs.len() >= nregs);
        for (dst, op) in ins {
            // SAFETY: every register index is < nregs by construction (checked by the assert).
            let v = apply(op, |r| unsafe { *regs.get_unchecked(r as usize) });
            unsafe { *regs.get_unchecked_mut(*dst as usize) = v; }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(src: &str) -> (Dag, N) {
        let mut d = Dag::default();
        let r = parse(src, &mut d, |s| match s { "x" => Some(Op::Var(0)), "y" => Some(Op::Var(1)),
            "k" => Some(Op::Par(0)), "time" => Some(Op::Time), _ => None }).unwrap();
        (d, r)
    }

    fn eval(src: &str, x: f64, y: f64, k: f64) -> (f64, f64, f64) {
        let (mut d, r) = build(src);
        let dx = d.diff(r, Op::Var(0));
        let dy = d.diff(r, Op::Var(1));
        let t = Tape::new(&d, &[r, dx, dy], 2, 1);
        let mut regs = t.regs();
        regs[0] = x; regs[1] = y; regs[2] = k; regs[3] = 0.5;
        t.run_all(&mut regs);
        (regs[t.outs[0] as usize], regs[t.outs[1] as usize], regs[t.outs[2] as usize])
    }

    #[test]
    fn values() {
        let e = |s| eval(s, 2.0, 3.0, 4.0).0;
        assert_eq!(e("1 + 2 * 3"), 7.0);
        assert_eq!(e("-2^2"), -4.0);
        assert_eq!(e("2^3^2"), 512.0);
        assert_eq!(e("1e-3 * 1000"), 1.0);
        assert_eq!(e("x * y - k"), 2.0);
        assert_eq!(e("piecewise(1, x > y, 2, x < y, 3)"), 2.0);
        assert_eq!(e("piecewise(1, x > y)").is_nan(), true);
        assert_eq!(e("log(100) + log(2, 8) + root(3, 27)"), 2.0 + 3.0 + 3.0);
        assert_eq!(e("(x <= 2 && !(y == 4)) || false"), 1.0);
        assert_eq!(e("min(x, y, 1) + max(x, k)"), 5.0);
        assert!((e("exp(ln(x)) + time") - 2.5).abs() < 1e-15);
        assert!(parse("foo(1)", &mut Dag::default(), |_| None).is_err());
        assert!(parse("x + ", &mut Dag::default(), |_| Some(Op::Var(0))).is_err());
    }

    #[test]
    fn jit_matches_interpreter() {
        let fs = ["k*x*y/(1+x+y^2)", "pow(x, 2.5)*exp(-y)", "x^y", "sqrt(x)*ln(y)+log10(x*y)",
                  "abs(x-y)*sin(x)+cos(y)*tan(x/5)", "piecewise(x^2, x > y, y^3)", "min(x,y)*max(x,y)",
                  "x^-3 * (1+y)^7 + floor(x) - ceil(y) + time", "(x <= 2 && !(y == 4)) || (x != y) + xor(x > 1, y < 1)",
                  "min(x, k) + max(k, y) - abs(x - 3) * (x - y)^13", "piecewise(1, x > y, 2, x < y, 3) + exp(x)/(1+exp(y))"];
        for f in fs {
            let (mut d, r) = build(f);
            let (dx, dy) = (d.diff(r, Op::Var(0)), d.diff(r, Op::Var(1)));
            let t = Tape::new(&d, &[r, dx, dy], 2, 1);
            let j = compile(&t.ins, &t.outs, t.nregs).unwrap();
            for &(x, y, k) in &[(0.7, 1.3, 2.0), (2.0, 0.4, f64::NAN), (-1.5, 1.6, 0.0), (2.0, 2.0, 1.0), (0.0, -0.0, -3.0)] {
                let mut a = t.regs();
                a[0] = x; a[1] = y; a[2] = k; a[3] = 0.25;
                t.prep(&mut a);
                let mut b = a.clone();
                Tape::exec(&t.ins, &mut a, t.nregs);
                unsafe { (j.f)(b.as_mut_ptr()) };
                for &o in &t.outs { assert_eq!(a[o as usize].to_bits(), b[o as usize].to_bits(), "{f} at {x} {y} {k}"); }
            }
        }
    }

    #[test]
    fn derivatives_match_finite_difference() {
        let fs = ["k*x*y/(1+x+y^2)", "pow(x, 2.5)*exp(-y)", "x^y", "sqrt(x)*ln(y)+log10(x*y)",
                  "abs(x-y)*sin(x)+cos(y)*tan(x/5)", "piecewise(x^2, x > y, y^3)", "min(x,y)*max(x,y)",
                  "root(3, x*y)", "k/(x*x) - 1/y", "x^-2 * (1+y)^3"];
        for f in fs {
            for &(x, y) in &[(0.7, 1.3), (2.0, 0.4), (1.5, 1.6)] {
                let (_, dx, dy) = eval(f, x, y, 1.7);
                let h = 1e-6;
                let fdx = (eval(f, x + h, y, 1.7).0 - eval(f, x - h, y, 1.7).0) / (2.0 * h);
                let fdy = (eval(f, x, y + h, 1.7).0 - eval(f, x, y - h, 1.7).0) / (2.0 * h);
                assert!((dx - fdx).abs() < 1e-6 * (1.0 + fdx.abs()), "{f} d/dx {dx} vs {fdx}");
                assert!((dy - fdy).abs() < 1e-6 * (1.0 + fdy.abs()), "{f} d/dy {dy} vs {fdy}");
            }
        }
    }
}

//! Native code for a tape (Cranelift): one straight-line function `f(regs)` over the tape's
//! instructions. Temporaries stay in SSA values; only output registers are stored. Results are
//! bit-identical to the interpreter (same operations in the same order, math via Rust's libm).
use crate::expr::Op;
use cranelift_codegen::ir::{condcodes::FloatCC, types::F64, AbiParam, InstBuilder, MemFlagsData, Signature, Value};
use cranelift_codegen::isa::OwnedTargetIsa;
use cranelift_codegen::settings::{self, Configurable};
use cranelift_frontend::{FunctionBuilder, FunctionBuilderContext};
use cranelift_jit::{JITBuilder, JITModule};
use cranelift_module::{default_libcall_names, Module};
use std::sync::OnceLock;

pub struct Jit { module: Option<JITModule>, pub f: unsafe extern "C" fn(*mut f64) }

// SAFETY: the module only owns finalized, immutable code; `f` is a pure function of its argument.
unsafe impl Send for Jit {}
unsafe impl Sync for Jit {}

impl Drop for Jit {
    fn drop(&mut self) {
        // SAFETY: `f` is dropped with self; nothing else points into the code.
        if let Some(m) = self.module.take() { unsafe { m.free_memory() } }
    }
}

extern "C" fn powf(a: f64, b: f64) -> f64 { a.powf(b) }
extern "C" fn exp(a: f64) -> f64 { a.exp() }
extern "C" fn ln(a: f64) -> f64 { a.ln() }
extern "C" fn log10(a: f64) -> f64 { a.log10() }
extern "C" fn sin(a: f64) -> f64 { a.sin() }
extern "C" fn cos(a: f64) -> f64 { a.cos() }
extern "C" fn tan(a: f64) -> f64 { a.tan() }

/// Compile `ins` (as in `Tape::ins`); registers in `outs` that `ins` writes are stored back.
pub fn compile(ins: &[(u32, Op)], outs: &[u32], nregs: usize) -> Option<Jit> {
    static ISA: OnceLock<Option<OwnedTargetIsa>> = OnceLock::new();
    let isa = ISA.get_or_init(|| {
        let mut fb = settings::builder();
        fb.set("opt_level", "none").ok()?; // the tape is already CSE'd and folded; "speed" runs no faster
        fb.set("enable_verifier", "false").ok()?;
        cranelift_native::builder().ok()?.finish(settings::Flags::new(fb)).ok()
    }).clone()?;
    let mut module = JITModule::new(JITBuilder::with_isa(isa, default_libcall_names()));
    let ptr = module.target_config().pointer_type();
    let cc = module.target_config().default_call_conv;
    let mut ctx = module.make_context();
    ctx.func.signature.params.push(AbiParam::new(ptr));
    let sig = |n: usize| { let mut s = Signature::new(cc); s.params.extend((0..n).map(|_| AbiParam::new(F64))); s.returns.push(AbiParam::new(F64)); s };
    let mut fctx = FunctionBuilderContext::new();
    {
        let mut b = FunctionBuilder::new(&mut ctx.func, &mut fctx);
        let (s1, s2) = (b.import_signature(sig(1)), b.import_signature(sig(2)));
        let blk = b.create_block();
        b.append_block_params_for_function_params(blk);
        b.switch_to_block(blk);
        let base = b.block_params(blk)[0];
        let mf = MemFlagsData::trusted();
        let mut val: Vec<Option<Value>> = vec![None; nregs];
        let mut store = vec![false; nregs];
        for &o in outs { store[o as usize] = true; }
        let (zero, one) = (b.ins().f64const(0.0), b.ins().f64const(1.0));
        for &(dst, op) in ins {
            let mut r = |b: &mut FunctionBuilder, i: u32| *val[i as usize].get_or_insert_with(|| b.ins().load(F64, mf, base, 8 * i as i32));
            macro_rules! call {
                ($f:expr, $($a:expr),+) => {{
                    let args = [$(r(&mut b, $a)),+];
                    let p = b.ins().iconst(ptr, $f as usize as i64);
                    let c = b.ins().call_indirect(if args.len() == 1 { s1 } else { s2 }, p, &args);
                    b.inst_results(c)[0]
                }};
            }
            let cmp = |b: &mut FunctionBuilder, cc: FloatCC, x: Value, y: Value| { let c = b.ins().fcmp(cc, x, y); b.ins().select(c, one, zero) };
            let nz = |b: &mut FunctionBuilder, x: Value| b.ins().fcmp(FloatCC::NotEqual, x, zero);
            use Op::*;
            let v = match op {
                Neg(a) => { let a = r(&mut b, a); b.ins().fneg(a) }
                Add(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); b.ins().fadd(a, c) }
                Sub(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); b.ins().fsub(a, c) }
                Mul(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); b.ins().fmul(a, c) }
                Div(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); b.ins().fdiv(a, c) }
                Pow(a, c) => call!(powf as extern "C" fn(f64, f64) -> f64, a, c),
                PowI(a, k) => { // compiler-rt __powidf2, unrolled
                    let mut x = r(&mut b, a);
                    let (mut acc, mut e) = (None::<Value>, k.unsigned_abs());
                    loop {
                        if e & 1 == 1 { acc = Some(match acc { None => x, Some(y) => b.ins().fmul(y, x) }); }
                        e /= 2;
                        if e == 0 { break; }
                        x = b.ins().fmul(x, x);
                    }
                    let p = acc.unwrap_or(one);
                    if k < 0 { b.ins().fdiv(one, p) } else { p }
                }
                Exp(a) => call!(exp as extern "C" fn(f64) -> f64, a),
                Ln(a) => call!(ln as extern "C" fn(f64) -> f64, a),
                Log10(a) => call!(log10 as extern "C" fn(f64) -> f64, a),
                Sin(a) => call!(sin as extern "C" fn(f64) -> f64, a),
                Cos(a) => call!(cos as extern "C" fn(f64) -> f64, a),
                Tan(a) => call!(tan as extern "C" fn(f64) -> f64, a),
                Sqrt(a) => { let a = r(&mut b, a); b.ins().sqrt(a) }
                Abs(a) => { let a = r(&mut b, a); b.ins().fabs(a) }
                Floor(a) => { let a = r(&mut b, a); b.ins().floor(a) }
                Ceil(a) => { let a = r(&mut b, a); b.ins().ceil(a) }
                Sign(a) => {
                    let a = r(&mut b, a);
                    let m1 = b.ins().f64const(-1.0);
                    let (gt, lt) = (b.ins().fcmp(FloatCC::GreaterThan, a, zero), b.ins().fcmp(FloatCC::LessThan, a, zero));
                    let n = b.ins().select(lt, m1, zero);
                    b.ins().select(gt, one, n)
                }
                Lt(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); cmp(&mut b, FloatCC::LessThan, a, c) }
                Le(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); cmp(&mut b, FloatCC::LessThanOrEqual, a, c) }
                Eq(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); cmp(&mut b, FloatCC::Equal, a, c) }
                Ne(a, c) => { let (a, c) = (r(&mut b, a), r(&mut b, c)); cmp(&mut b, FloatCC::NotEqual, a, c) }
                And(a, c) | Or(a, c) | Xor(a, c) => {
                    let (a, c) = (r(&mut b, a), r(&mut b, c));
                    let (x, y) = (nz(&mut b, a), nz(&mut b, c));
                    let z = match op { And(..) => b.ins().band(x, y), Or(..) => b.ins().bor(x, y), _ => b.ins().bxor(x, y) };
                    b.ins().select(z, one, zero)
                }
                Not(a) => { let a = r(&mut b, a); cmp(&mut b, FloatCC::Equal, a, zero) }
                Sel(c, x, y) => {
                    let (c, x, y) = (r(&mut b, c), r(&mut b, x), r(&mut b, y));
                    let t = nz(&mut b, c);
                    b.ins().select(t, x, y)
                }
                Min(a, c) | Max(a, c) => { // Rust f64::min/max (fminnm): a NaN operand yields the other one
                    let (a, c) = (r(&mut b, a), r(&mut b, c));
                    let m = if matches!(op, Min(..)) { b.ins().fmin(a, c) } else { b.ins().fmax(a, c) };
                    let (na, nc) = (b.ins().fcmp(FloatCC::Unordered, a, a), b.ins().fcmp(FloatCC::Unordered, c, c));
                    let m = b.ins().select(nc, a, m);
                    b.ins().select(na, c, m)
                }
                Const(_) | Var(_) | Par(_) | Time => unreachable!(),
            };
            val[dst as usize] = Some(v);
            if store[dst as usize] { b.ins().store(mf, v, base, 8 * dst as i32); }
        }
        b.ins().return_(&[]);
        b.seal_all_blocks();
        b.finalize(module.target_config());
    }
    let id = module.declare_anonymous_function(&ctx.func.signature).ok()?;
    module.define_function(id, &mut ctx).ok()?;
    module.clear_context(&mut ctx);
    module.finalize_definitions().ok()?;
    let p = module.get_finalized_function(id);
    // SAFETY: the function was built with signature (ptr) -> () in the default call convention.
    let f = unsafe { std::mem::transmute::<*const u8, unsafe extern "C" fn(*mut f64)>(p) };
    Some(Jit { module: Some(module), f })
}

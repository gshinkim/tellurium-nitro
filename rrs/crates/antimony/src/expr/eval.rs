//! Evaluating an expression (Delphi `TAntExprNode.TryEvaluate`).
//!
//! This is constant folding, not simulation. `None` means "no value": a
//! symbol the caller cannot resolve, a call that is not a built-in we can
//! compute (a user function, `delay`, `rateOf`), or an explicit division by
//! zero. An unknown name is never silently a zero (Delphi §6.5).

use super::{BinaryOp, Expr, ExprKind, UnaryOp};

fn truth(value: bool) -> f64 {
    if value { 1.0 } else { 0.0 }
}

/// `1 / x`, or `None` for a zero denominator rather than an invented
/// infinity.
fn reciprocal(x: f64) -> Option<f64> {
    (x != 0.0).then(|| 1.0 / x)
}

impl Expr {
    /// The value, if it can be computed without knowing any symbol.
    pub fn evaluate(&self) -> Option<f64> {
        self.evaluate_with(&|_| None)
    }

    /// The value, asking `resolve` for each symbol's value. `resolve`
    /// returning `None` makes the whole evaluation `None`, except where
    /// `&&`, `||` or `piecewise` never needed that part.
    pub fn evaluate_with(&self, resolve: &dyn Fn(&str) -> Option<f64>) -> Option<f64> {
        match &self.kind {
            ExprKind::Number { value, .. } => Some(*value),
            ExprKind::Constant(c) => Some(c.value()),
            ExprKind::Symbol(name) => resolve(name),
            ExprKind::Unary { op, operand } => {
                let x = operand.evaluate_with(resolve)?;
                Some(match op {
                    UnaryOp::Plus => x,
                    UnaryOp::Minus => -x,
                    UnaryOp::Not => truth(x == 0.0),
                })
            }
            ExprKind::Conditional {
                condition,
                then,
                otherwise,
            } => {
                if condition.evaluate_with(resolve)? != 0.0 {
                    then.evaluate_with(resolve)
                } else {
                    otherwise.evaluate_with(resolve)
                }
            }
            ExprKind::Call { name, args } => evaluate_call(name, args, resolve),
            ExprKind::Binary { op, left, right } => {
                let l = left.evaluate_with(resolve)?;
                // `&&` and `||` short-circuit, so `0 && unknown` is still 0.
                match op {
                    BinaryOp::And if l == 0.0 => return Some(0.0),
                    BinaryOp::Or if l != 0.0 => return Some(1.0),
                    _ => {}
                }
                let r = right.evaluate_with(resolve)?;
                Some(match op {
                    BinaryOp::And | BinaryOp::Or => truth(r != 0.0),
                    BinaryOp::Add => l + r,
                    BinaryOp::Subtract => l - r,
                    BinaryOp::Multiply => l * r,
                    BinaryOp::Divide => return reciprocal(r).map(|_| l / r),
                    BinaryOp::Modulo => return reciprocal(r).map(|_| l - (l / r).trunc() * r),
                    BinaryOp::Power => l.powf(r),
                    BinaryOp::Equal => truth(l == r),
                    BinaryOp::NotEqual => truth(l != r),
                    BinaryOp::Less => truth(l < r),
                    BinaryOp::LessEqual => truth(l <= r),
                    BinaryOp::Greater => truth(l > r),
                    BinaryOp::GreaterEqual => truth(l >= r),
                })
            }
        }
    }
}

fn evaluate_call(name: &str, args: &[Expr], resolve: &dyn Fn(&str) -> Option<f64>) -> Option<f64> {
    // piecewise(v1, c1, v2, c2, …, otherwise), evaluated lazily: a branch that
    // is not selected may be impossible to evaluate without failing the call.
    if name == "piecewise" {
        let mut i = 0;
        while i + 1 < args.len() {
            if args[i + 1].evaluate_with(resolve)? != 0.0 {
                return args[i].evaluate_with(resolve);
            }
            i += 2;
        }
        return args.get(i)?.evaluate_with(resolve);
    }

    let values = args
        .iter()
        .map(|arg| arg.evaluate_with(resolve))
        .collect::<Option<Vec<f64>>>()?;

    // MathML's n-ary operators, callable by name. The corpus uses
    // `xor(a, b)` 21 times and has no infix xor at all.
    match name {
        "and" | "or" | "xor" if values.is_empty() => return None,
        "and" => return Some(truth(values.iter().all(|&v| v != 0.0))),
        "or" => return Some(truth(values.iter().any(|&v| v != 0.0))),
        // n-ary xor is parity.
        "xor" => return Some(truth(values.iter().filter(|&&v| v != 0.0).count() % 2 == 1)),
        "implies" if values.len() == 2 => {
            return Some(truth(values[0] == 0.0 || values[1] != 0.0));
        }
        "implies" => return None,
        "plus" => return Some(values.iter().sum()),
        "times" => return Some(values.iter().product()),
        "min" | "max" if !values.is_empty() => {
            let pick = |a: f64, b: f64| {
                if (name == "min") == (b < a) { b } else { a }
            };
            return values.iter().copied().reduce(pick);
        }
        _ => {}
    }

    match *values.as_slice() {
        [x] => evaluate_unary_call(name, x),
        [a, b] => evaluate_binary_call(name, a, b),
        _ => None,
    }
}

fn evaluate_unary_call(name: &str, x: f64) -> Option<f64> {
    Some(match name {
        "abs" => x.abs(),
        "exp" => x.exp(),
        "ln" => x.ln(),
        "log10" => x.log10(),
        "sqrt" => x.sqrt(),
        "ceil" | "ceiling" => x.ceil(),
        "floor" => x.floor(),
        // A bare `log` has no stated base. libSBML's formula parser reads it
        // as base 10, so interoperability says base 10. The corpus never
        // uses it; it writes log10.
        "log" => x.log10(),
        "sin" => x.sin(),
        "cos" => x.cos(),
        "tan" => x.tan(),
        "sinh" => x.sinh(),
        "cosh" => x.cosh(),
        "tanh" => x.tanh(),
        "sec" => return reciprocal(x.cos()),
        "csc" => return reciprocal(x.sin()),
        "cot" => return reciprocal(x.tan()),
        "sech" => return reciprocal(x.cosh()),
        "csch" => return reciprocal(x.sinh()),
        "coth" => return reciprocal(x.tanh()),
        // libantimony writes `asin` where the documentation says `arcsin`.
        "asin" | "arcsin" => x.asin(),
        "acos" | "arccos" => x.acos(),
        "atan" | "arctan" => x.atan(),
        "arcsinh" => x.asinh(),
        "arccosh" => x.acosh(),
        "arctanh" => x.atanh(),
        // arcsec(x) = arccos(1/x), and so on.
        "arcsec" => return reciprocal(x).map(f64::acos),
        "arccsc" => return reciprocal(x).map(f64::asin),
        "arccot" => return reciprocal(x).map(f64::atan),
        "arcsech" => return reciprocal(x).map(f64::acosh),
        "arccsch" => return reciprocal(x).map(f64::asinh),
        "arccoth" => return reciprocal(x).map(f64::atanh),
        "factorial" => {
            // 171! overflows a double.
            if x < 0.0 || x.trunc() != x || x > 170.0 {
                return None;
            }
            (2..=x as u32).map(f64::from).product()
        }
        "minus" => -x,
        "not" => truth(x == 0.0),
        _ => return None,
    })
}

fn evaluate_binary_call(name: &str, a: f64, b: f64) -> Option<f64> {
    Some(match name {
        "pow" | "power" => a.powf(b),
        "root" => return reciprocal(a).map(|inverse| b.powf(inverse)),
        // MathML's log takes the base first.
        "log" => b.ln() / a.ln(),
        "minus" => a - b,
        "divide" => return reciprocal(b).map(|_| a / b),
        "quotient" => return reciprocal(b).map(|_| (a / b).trunc()),
        "rem" => return reciprocal(b).map(|_| a - (a / b).trunc() * b),
        // The comparisons are callable by name as well as infix.
        "eq" => truth(a == b),
        "neq" => truth(a != b),
        "lt" => truth(a < b),
        "leq" => truth(a <= b),
        "gt" => truth(a > b),
        "geq" => truth(a >= b),
        _ => return None,
    })
}

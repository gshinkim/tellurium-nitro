//! Expressions (Delphi `uAntExprNode`, Delphi §5.3).
//!
//! An [`Expr`] is an owned tree: copying one is `.clone()`, and it can never
//! be shared by two owners, which is what Delphi §3.3's clone-on-assign rule
//! asked for by hand. Every node carries its [`Span`].
//!
//! Names are never resolved here. A symbol is its name, and binding names to
//! model items is the semantic pass's job (SPEC.md §3.1).

mod eval;
mod parse;

use std::collections::HashSet;
use std::fmt;

pub use parse::{ExprParser, MAX_DEPTH, MAX_NESTING, can_start_expression, parse_expression};

use crate::number;
use crate::span::Span;

/// One node of an expression, and where it came from.
#[derive(Clone, Debug, PartialEq)]
pub struct Expr {
    pub kind: ExprKind,
    pub span: Span,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ExprKind {
    /// `2.5`, `inf`, `nan`, or a number with a unit: `40 foo`. The unit is
    /// SBML's `<cn units="…">`, kept as written; it is one identifier.
    Number {
        value: f64,
        units: Option<String>,
    },
    /// `pi`, `exponentiale`, `true`, `false`.
    Constant(Constant),
    /// A name, possibly dotted: `S1`, `A.S1`. `time` is a symbol too.
    Symbol(String),
    Unary {
        op: UnaryOp,
        operand: Box<Expr>,
    },
    Binary {
        op: BinaryOp,
        left: Box<Expr>,
        right: Box<Expr>,
    },
    /// `condition ? then : otherwise`.
    Conditional {
        condition: Box<Expr>,
        then: Box<Expr>,
        otherwise: Box<Expr>,
    },
    /// `f(…)`, including `piecewise` and the MathML operators called by name.
    Call {
        name: String,
        args: Vec<Expr>,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum UnaryOp {
    Plus,
    Minus,
    Not,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BinaryOp {
    Add,
    Subtract,
    Multiply,
    Divide,
    Modulo,
    Power,
    Equal,
    NotEqual,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    And,
    Or,
}

/// The MathML named constants. They are their own kind of node, not
/// symbols, so that inlining a function can substitute them and no function
/// argument escapes as an invented parameter (Delphi §9.1.4 #4). `inf` and
/// `nan` are numbers instead.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Constant {
    Pi,
    ExponentialE,
    True,
    False,
}

/// Operator precedence (Delphi §5.3), loosest first. Two levels differ from
/// most languages, deliberately: `!` sits between `&&` and `==`, so `!a == b`
/// is `!(a == b)`; and unary `+`/`-` bind tighter than `^`, so `-2^2` is 4.
pub mod precedence {
    pub const CONDITIONAL: u8 = 1;
    pub const OR: u8 = 2;
    pub const AND: u8 = 3;
    pub const NOT: u8 = 4;
    pub const EQUALITY: u8 = 5;
    pub const RELATIONAL: u8 = 6;
    pub const ADDITIVE: u8 = 7;
    pub const MULTIPLICATIVE: u8 = 8;
    pub const POWER: u8 = 9;
    pub const UNARY: u8 = 10;
    pub const PRIMARY: u8 = 11;
}

impl UnaryOp {
    pub fn spelling(self) -> &'static str {
        match self {
            UnaryOp::Plus => "+",
            UnaryOp::Minus => "-",
            UnaryOp::Not => "!",
        }
    }

    pub fn precedence(self) -> u8 {
        match self {
            UnaryOp::Not => precedence::NOT,
            UnaryOp::Plus | UnaryOp::Minus => precedence::UNARY,
        }
    }
}

impl BinaryOp {
    pub fn spelling(self) -> &'static str {
        use BinaryOp::*;
        match self {
            Add => "+",
            Subtract => "-",
            Multiply => "*",
            Divide => "/",
            Modulo => "%",
            Power => "^",
            Equal => "==",
            NotEqual => "!=",
            Less => "<",
            LessEqual => "<=",
            Greater => ">",
            GreaterEqual => ">=",
            And => "&&",
            Or => "||",
        }
    }

    pub fn precedence(self) -> u8 {
        use BinaryOp::*;
        match self {
            Add | Subtract => precedence::ADDITIVE,
            Multiply | Divide | Modulo => precedence::MULTIPLICATIVE,
            Power => precedence::POWER,
            Equal | NotEqual => precedence::EQUALITY,
            Less | LessEqual | Greater | GreaterEqual => precedence::RELATIONAL,
            And => precedence::AND,
            Or => precedence::OR,
        }
    }
}

impl Constant {
    pub fn spelling(self) -> &'static str {
        match self {
            Constant::Pi => "pi",
            Constant::ExponentialE => "exponentiale",
            Constant::True => "true",
            Constant::False => "false",
        }
    }

    pub fn value(self) -> f64 {
        match self {
            Constant::Pi => std::f64::consts::PI,
            Constant::ExponentialE => std::f64::consts::E,
            Constant::True => 1.0,
            Constant::False => 0.0,
        }
    }

    /// The constant `word` spells, case-sensitively.
    pub fn from_word(word: &str) -> Option<Constant> {
        match word {
            "pi" => Some(Constant::Pi),
            "exponentiale" => Some(Constant::ExponentialE),
            "true" => Some(Constant::True),
            "false" => Some(Constant::False),
            _ => None,
        }
    }
}

impl Expr {
    pub fn new(kind: ExprKind, span: Span) -> Expr {
        Expr { kind, span }
    }

    pub fn number(value: f64, span: Span) -> Expr {
        Expr::new(ExprKind::Number { value, units: None }, span)
    }

    pub fn constant(constant: Constant, span: Span) -> Expr {
        Expr::new(ExprKind::Constant(constant), span)
    }

    pub fn symbol(name: impl Into<String>, span: Span) -> Expr {
        Expr::new(ExprKind::Symbol(name.into()), span)
    }

    pub fn unary(op: UnaryOp, operand: Expr, span: Span) -> Expr {
        Expr::new(
            ExprKind::Unary {
                op,
                operand: Box::new(operand),
            },
            span,
        )
    }

    pub fn binary(op: BinaryOp, left: Expr, right: Expr, span: Span) -> Expr {
        Expr::new(
            ExprKind::Binary {
                op,
                left: Box::new(left),
                right: Box::new(right),
            },
            span,
        )
    }

    pub fn conditional(condition: Expr, then: Expr, otherwise: Expr, span: Span) -> Expr {
        Expr::new(
            ExprKind::Conditional {
                condition: Box::new(condition),
                then: Box::new(then),
                otherwise: Box::new(otherwise),
            },
            span,
        )
    }

    pub fn call(name: impl Into<String>, args: Vec<Expr>, span: Span) -> Expr {
        Expr::new(
            ExprKind::Call {
                name: name.into(),
                args,
            },
            span,
        )
    }

    /// The node's children, in order.
    pub fn children(&self) -> Vec<&Expr> {
        match &self.kind {
            ExprKind::Number { .. } | ExprKind::Constant(_) | ExprKind::Symbol(_) => Vec::new(),
            ExprKind::Unary { operand, .. } => vec![operand],
            ExprKind::Binary { left, right, .. } => vec![left, right],
            ExprKind::Conditional {
                condition,
                then,
                otherwise,
            } => vec![condition, then, otherwise],
            ExprKind::Call { args, .. } => args.iter().collect(),
        }
    }

    /// The same tree, ignoring spans. Two `nan` literals are the same tree,
    /// although NaN is not equal to itself as a number.
    pub fn same_tree(&self, other: &Expr) -> bool {
        use ExprKind::*;
        match (&self.kind, &other.kind) {
            (
                Number {
                    value: a,
                    units: ua,
                },
                Number {
                    value: b,
                    units: ub,
                },
            ) => ua == ub && (a == b || (a.is_nan() && b.is_nan())),
            (Constant(a), Constant(b)) => a == b,
            (Symbol(a), Symbol(b)) => a == b,
            (Unary { op: a, operand: x }, Unary { op: b, operand: y }) => a == b && x.same_tree(y),
            (
                Binary {
                    op: a,
                    left: l1,
                    right: r1,
                },
                Binary {
                    op: b,
                    left: l2,
                    right: r2,
                },
            ) => a == b && l1.same_tree(l2) && r1.same_tree(r2),
            (
                Conditional {
                    condition: c1,
                    then: t1,
                    otherwise: o1,
                },
                Conditional {
                    condition: c2,
                    then: t2,
                    otherwise: o2,
                },
            ) => c1.same_tree(c2) && t1.same_tree(t2) && o1.same_tree(o2),
            (Call { name: a, args: x }, Call { name: b, args: y }) => {
                a == b && x.len() == y.len() && x.iter().zip(y).all(|(p, q)| p.same_tree(q))
            }
            _ => false,
        }
    }

    /// No symbols and no calls, so the value is fixed when the file is read.
    /// This is what tells `x = 3` (a value) from `x = 2*k` (an initial
    /// assignment). A call is not constant: it may not be pure.
    pub fn is_constant(&self) -> bool {
        match &self.kind {
            ExprKind::Number { .. } | ExprKind::Constant(_) => true,
            ExprKind::Symbol(_) | ExprKind::Call { .. } => false,
            _ => self.children().into_iter().all(Expr::is_constant),
        }
    }

    /// The symbols used, each once, in order of first appearance. Named
    /// constants are not symbols, and neither is a called function's name.
    pub fn symbols(&self) -> Vec<&str> {
        let mut seen = HashSet::new();
        let mut order = Vec::new();
        self.collect_symbols(&mut seen, &mut order);
        order
    }

    fn collect_symbols<'a>(&'a self, seen: &mut HashSet<&'a str>, order: &mut Vec<&'a str>) {
        if let ExprKind::Symbol(name) = &self.kind {
            if seen.insert(name) {
                order.push(name);
            }
            return;
        }
        for child in self.children() {
            child.collect_symbols(seen, order);
        }
    }

    pub fn node_count(&self) -> usize {
        1 + self
            .children()
            .into_iter()
            .map(Expr::node_count)
            .sum::<usize>()
    }

    /// 1 for a leaf.
    pub fn depth(&self) -> usize {
        1 + self
            .children()
            .into_iter()
            .map(Expr::depth)
            .max()
            .unwrap_or(0)
    }

    fn precedence(&self) -> u8 {
        match &self.kind {
            ExprKind::Number { .. }
            | ExprKind::Constant(_)
            | ExprKind::Symbol(_)
            | ExprKind::Call { .. } => precedence::PRIMARY,
            ExprKind::Unary { op, .. } => op.precedence(),
            ExprKind::Binary { op, .. } => op.precedence(),
            ExprKind::Conditional { .. } => precedence::CONDITIONAL,
        }
    }

    /// Writes the expression with as few brackets as its meaning needs, so
    /// writing, reading and writing again gives the same text.
    fn emit(&self, out: &mut fmt::Formatter<'_>, parent: u8, right_child: bool) -> fmt::Result {
        let own = self.precedence();
        // At equal precedence the wrong side needs brackets too: `a - (b - c)`
        // and `a ^ (b ^ c)` differ from `a - b - c` and `(a ^ b) ^ c`.
        let brackets = own < parent
            || (own == parent
                && match &self.kind {
                    ExprKind::Binary {
                        op: BinaryOp::Power,
                        ..
                    } => !right_child,
                    ExprKind::Binary { .. } => right_child,
                    _ => false,
                });
        if brackets {
            out.write_str("(")?;
        }
        match &self.kind {
            ExprKind::Number { value, units } => {
                out.write_str(&number::format(*value))?;
                if let Some(units) = units {
                    write!(out, " {units}")?;
                }
            }
            ExprKind::Constant(c) => out.write_str(c.spelling())?,
            ExprKind::Symbol(name) => out.write_str(name)?,
            ExprKind::Unary { op, operand } => {
                out.write_str(op.spelling())?;
                // A space before a unary operand of a unary operator, or
                // `- -x` would come out as `--x`, the DNA operator. Delphi
                // bracketed it, `-(-x)`, but every bracket is a level of
                // nesting, and what is written must never be nested deeper
                // than what was read, or it could not be read back.
                if matches!(operand.kind, ExprKind::Unary { .. }) {
                    out.write_str(" ")?;
                }
                operand.emit(out, own, false)?;
            }
            ExprKind::Binary { op, left, right } => {
                left.emit(out, own, false)?;
                write!(out, " {} ", op.spelling())?;
                right.emit(out, own, true)?;
            }
            ExprKind::Conditional {
                condition,
                then,
                otherwise,
            } => {
                condition.emit(out, precedence::CONDITIONAL + 1, false)?;
                out.write_str(" ? ")?;
                then.emit(out, precedence::CONDITIONAL, false)?;
                out.write_str(" : ")?;
                otherwise.emit(out, precedence::CONDITIONAL, true)?;
            }
            ExprKind::Call { name, args } => {
                write!(out, "{name}(")?;
                for (i, arg) in args.iter().enumerate() {
                    if i > 0 {
                        out.write_str(", ")?;
                    }
                    arg.emit(out, precedence::CONDITIONAL, false)?;
                }
                out.write_str(")")?;
            }
        }
        if brackets {
            out.write_str(")")?;
        }
        Ok(())
    }
}

impl fmt::Display for Expr {
    /// Antimony text, with only the brackets the meaning needs:
    /// `k1 * S1 ^ 2 * S2`, `(a + b) * c`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.emit(f, 0, false)
    }
}

//! MathML, the subset SBML allows, as a tree.

/// A MathML expression.
#[derive(Debug, Clone, PartialEq)]
pub enum Math {
    /// `<cn>`, with the units Level 3 allows on it.
    Number {
        value: Number,
        units: Option<String>,
    },
    /// `<ci>`: a reference to an id.
    Ident(String),
    /// `<csymbol>` on its own: time or Avogadro's number.
    Symbol(Symbol),
    /// `<true/>`, `<pi/>` and the rest.
    Constant(Constant),
    /// A built-in operator applied to its arguments.
    Apply(Op, Vec<Math>),
    /// A call to a function definition.
    Call(String, Vec<Math>),
    /// A `<csymbol>` applied to arguments: `delay` or `rateOf`.
    CallSymbol(Symbol, Vec<Math>),
    /// `<root>`, with its optional `<degree>`.
    Root {
        degree: Option<Box<Math>>,
        arg: Box<Math>,
    },
    /// `<log>`, with its optional `<logbase>`.
    Log {
        base: Option<Box<Math>>,
        arg: Box<Math>,
    },
    /// `<piecewise>`: (value, condition) pairs, then the value otherwise.
    Piecewise {
        pieces: Vec<(Math, Math)>,
        otherwise: Option<Box<Math>>,
    },
    /// `<lambda>`: the body of a function definition.
    Lambda { args: Vec<String>, body: Box<Math> },
    /// `<semantics>`: an expression with annotations, kept as raw XML.
    Semantics {
        math: Box<Math>,
        annotations: Vec<String>,
    },
    /// MathML that `sbml-io` doesn't understand, kept as raw XML.
    Raw(String),
}

/// The value of a `<cn>`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Number {
    Integer(i64),
    Real(f64),
    /// mantissa × 10^exponent.
    ENotation(f64, i64),
    /// numerator / denominator.
    Rational(i64, i64),
}

impl Number {
    pub fn value(self) -> f64 {
        match self {
            Number::Integer(i) => i as f64,
            Number::Real(r) => r,
            Number::ENotation(m, e) => format!("{m}e{e}").parse().unwrap_or(f64::NAN),
            Number::Rational(n, d) => n as f64 / d as f64,
        }
    }
}

/// A `<csymbol>`: its definition URL and the name written in it.
#[derive(Debug, Clone, PartialEq)]
pub struct Symbol {
    pub url: String,
    pub name: String,
}

impl Symbol {
    pub const TIME: &str = "http://www.sbml.org/sbml/symbols/time";
    pub const DELAY: &str = "http://www.sbml.org/sbml/symbols/delay";
    pub const AVOGADRO: &str = "http://www.sbml.org/sbml/symbols/avogadro";
    /// Level 3 Version 2 only.
    pub const RATE_OF: &str = "http://www.sbml.org/sbml/symbols/rateOf";

    pub fn new(url: &str, name: &str) -> Self {
        Symbol {
            url: url.to_string(),
            name: name.to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Constant {
    True,
    False,
    Pi,
    ExponentialE,
    NotANumber,
    Infinity,
}

impl Constant {
    pub(crate) const ALL: [(Constant, &str); 6] = [
        (Constant::True, "true"),
        (Constant::False, "false"),
        (Constant::Pi, "pi"),
        (Constant::ExponentialE, "exponentiale"),
        (Constant::NotANumber, "notanumber"),
        (Constant::Infinity, "infinity"),
    ];

    /// The MathML element name.
    pub fn tag(self) -> &'static str {
        Self::ALL.iter().find(|(c, _)| *c == self).unwrap().1
    }

    pub fn from_tag(tag: &str) -> Option<Self> {
        Self::ALL.iter().find(|(_, t)| *t == tag).map(|(c, _)| *c)
    }
}

macro_rules! ops {
    ($($op:ident = $tag:literal),* $(,)?) => {
        /// A built-in MathML operator. `root` and `log` have their own
        /// `Math` variants, because of their qualifiers.
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub enum Op { $($op),* }

        impl Op {
            /// The MathML element name.
            pub fn tag(self) -> &'static str {
                match self { $(Op::$op => $tag),* }
            }

            pub fn from_tag(tag: &str) -> Option<Self> {
                match tag { $($tag => Some(Op::$op),)* _ => None }
            }
        }
    };
}

ops! {
    Plus = "plus", Minus = "minus", Times = "times", Divide = "divide", Power = "power",
    Abs = "abs", Exp = "exp", Ln = "ln", Floor = "floor", Ceiling = "ceiling",
    Factorial = "factorial",
    Eq = "eq", Neq = "neq", Gt = "gt", Lt = "lt", Geq = "geq", Leq = "leq",
    And = "and", Or = "or", Xor = "xor", Not = "not", Implies = "implies",
    Quotient = "quotient", Rem = "rem", Max = "max", Min = "min",
    Sin = "sin", Cos = "cos", Tan = "tan", Sec = "sec", Csc = "csc", Cot = "cot",
    Sinh = "sinh", Cosh = "cosh", Tanh = "tanh", Sech = "sech", Csch = "csch", Coth = "coth",
    ArcSin = "arcsin", ArcCos = "arccos", ArcTan = "arctan",
    ArcSec = "arcsec", ArcCsc = "arccsc", ArcCot = "arccot",
    ArcSinh = "arcsinh", ArcCosh = "arccosh", ArcTanh = "arctanh",
    ArcSech = "arcsech", ArcCsch = "arccsch", ArcCoth = "arccoth",
}

impl Math {
    pub fn int(i: i64) -> Self {
        Math::Number {
            value: Number::Integer(i),
            units: None,
        }
    }

    pub fn real(r: f64) -> Self {
        Math::Number {
            value: Number::Real(r),
            units: None,
        }
    }

    pub fn ident(id: &str) -> Self {
        Math::Ident(id.to_string())
    }

    /// Calls `f` on every identifier the expression refers to: `<ci>`s and
    /// called function names. Lambda arguments are reported too, where used.
    pub fn visit_idents<'a>(&'a self, f: &mut impl FnMut(&'a str)) {
        match self {
            Math::Ident(id) => f(id),
            Math::Number { .. } | Math::Symbol(_) | Math::Constant(_) | Math::Raw(_) => {}
            Math::Apply(_, args) | Math::CallSymbol(_, args) => {
                args.iter().for_each(|a| a.visit_idents(f))
            }
            Math::Call(name, args) => {
                f(name);
                args.iter().for_each(|a| a.visit_idents(f));
            }
            Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
                if let Some(q) = q {
                    q.visit_idents(f);
                }
                arg.visit_idents(f);
            }
            Math::Piecewise { pieces, otherwise } => {
                for (v, c) in pieces {
                    v.visit_idents(f);
                    c.visit_idents(f);
                }
                if let Some(o) = otherwise {
                    o.visit_idents(f);
                }
            }
            Math::Lambda { body, .. } => body.visit_idents(f),
            Math::Semantics { math, .. } => math.visit_idents(f),
        }
    }
}

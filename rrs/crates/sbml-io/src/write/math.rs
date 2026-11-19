//! Writing `Math` as MathML.

use super::{W, escape, fmt_f64};
use crate::math::{Math, Number, Symbol};
use crate::ns;

impl W {
    /// A `<math>` element.
    pub(crate) fn math(&mut self, m: &Math) {
        // Level 3 allows units on numbers, in the core namespace.
        let units = self.level >= 3 && has_units(m);
        if units {
            self.line(&format!(
                "<math xmlns=\"{}\" xmlns:sbml=\"{}\">",
                ns::MATHML,
                self.core
            ));
        } else {
            self.line(&format!("<math xmlns=\"{}\">", ns::MATHML));
        }
        self.push_depth();
        self.expr(m);
        self.pop_depth();
        self.line("</math>");
    }

    fn expr(&mut self, m: &Math) {
        match m {
            Math::Number { value, units } => {
                let units = match units {
                    Some(u) if self.level >= 3 => format!(" sbml:units=\"{}\"", escape(u, true)),
                    _ => String::new(),
                };
                let s = match value {
                    Number::Integer(i) => format!("<cn type=\"integer\"{units}> {i} </cn>"),
                    Number::Real(r) => format!("<cn{units}> {} </cn>", fmt_f64(*r)),
                    Number::ENotation(m, e) => format!(
                        "<cn type=\"e-notation\"{units}> {} <sep/> {e} </cn>",
                        fmt_f64(*m)
                    ),
                    Number::Rational(n, d) => {
                        format!("<cn type=\"rational\"{units}> {n} <sep/> {d} </cn>")
                    }
                };
                self.line(&s);
            }
            Math::Ident(id) => self.line(&format!("<ci> {} </ci>", escape(id, false))),
            Math::Symbol(s) => self.csymbol(s),
            Math::Constant(c) => self.line(&format!("<{}/>", c.tag())),
            Math::Apply(op, args) => self.apply(&format!("<{}/>", op.tag()), None, args),
            Math::Call(f, args) => {
                self.apply(&format!("<ci> {} </ci>", escape(f, false)), None, args)
            }
            Math::CallSymbol(s, args) => {
                self.line("<apply>");
                self.push_depth();
                self.csymbol(s);
                for a in args {
                    self.expr(a);
                }
                self.pop_depth();
                self.line("</apply>");
            }
            Math::Root { degree, arg } => {
                let q = degree.as_deref().map(|d| ("degree", d));
                self.apply("<root/>", q, std::slice::from_ref(arg.as_ref()));
            }
            Math::Log { base, arg } => {
                let q = base.as_deref().map(|b| ("logbase", b));
                self.apply("<log/>", q, std::slice::from_ref(arg.as_ref()));
            }
            Math::Piecewise { pieces, otherwise } => {
                self.line("<piecewise>");
                self.push_depth();
                for (v, c) in pieces {
                    self.line("<piece>");
                    self.push_depth();
                    self.expr(v);
                    self.expr(c);
                    self.pop_depth();
                    self.line("</piece>");
                }
                if let Some(o) = otherwise {
                    self.line("<otherwise>");
                    self.push_depth();
                    self.expr(o);
                    self.pop_depth();
                    self.line("</otherwise>");
                }
                self.pop_depth();
                self.line("</piecewise>");
            }
            Math::Lambda { args, body } => {
                self.line("<lambda>");
                self.push_depth();
                for a in args {
                    self.line(&format!("<bvar><ci> {} </ci></bvar>", escape(a, false)));
                }
                self.expr(body);
                self.pop_depth();
                self.line("</lambda>");
            }
            Math::Semantics { math, annotations } => {
                self.line("<semantics>");
                self.push_depth();
                self.expr(math);
                for a in annotations {
                    self.line(a);
                }
                self.pop_depth();
                self.line("</semantics>");
            }
            Math::Raw(xml) => self.line(xml),
        }
    }

    fn csymbol(&mut self, s: &Symbol) {
        self.line(&format!(
            "<csymbol encoding=\"text\" definitionURL=\"{}\"> {} </csymbol>",
            escape(&s.url, true),
            escape(&s.name, false)
        ));
    }

    fn apply(&mut self, head: &str, qualifier: Option<(&str, &Math)>, args: &[Math]) {
        self.line("<apply>");
        self.push_depth();
        self.line(head);
        if let Some((tag, q)) = qualifier {
            self.line(&format!("<{tag}>"));
            self.push_depth();
            self.expr(q);
            self.pop_depth();
            self.line(&format!("</{tag}>"));
        }
        for a in args {
            self.expr(a);
        }
        self.pop_depth();
        self.line("</apply>");
    }
}

fn has_units(m: &Math) -> bool {
    match m {
        Math::Number { units, .. } => units.is_some(),
        Math::Ident(_) | Math::Symbol(_) | Math::Constant(_) | Math::Raw(_) => false,
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            args.iter().any(has_units)
        }
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            q.as_deref().is_some_and(has_units) || has_units(arg)
        }
        Math::Piecewise { pieces, otherwise } => {
            pieces.iter().any(|(v, c)| has_units(v) || has_units(c))
                || otherwise.as_deref().is_some_and(has_units)
        }
        Math::Lambda { body, .. } => has_units(body),
        Math::Semantics { math, .. } => has_units(math),
    }
}

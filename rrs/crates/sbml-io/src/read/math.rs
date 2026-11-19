//! Reading MathML into `Math`. Anything outside the subset `Math` models is
//! kept as `Math::Raw`.

use roxmltree::Node;

use super::{Reader, parse_double};
use crate::math::{Constant, Math, Number, Op, Symbol};
use crate::ns;

impl<'a, 'i: 'a> Reader<'a, 'i> {
    /// A `<math>` element. `None` if it's empty.
    pub(super) fn math(&mut self, n: Node<'a, 'i>) -> Option<Math> {
        let mut kids = n.children().filter(|c| c.is_element());
        let first = kids.next()?;
        if let Some(more) = kids.next() {
            self.warn(more, "only the first expression in <math> is read");
        }
        Some(self.expr(first))
    }

    fn expr(&mut self, n: Node<'a, 'i>) -> Math {
        self.try_expr(n).unwrap_or_else(|| Math::Raw(self.raw(n)))
    }

    /// `None` means "keep it as raw XML".
    fn try_expr(&mut self, n: Node<'a, 'i>) -> Option<Math> {
        let name = self.name_in(n, ns::MATHML)?;
        match name {
            "cn" => self.cn(n),
            "ci" => {
                only_attrs(n, &[])?;
                Some(Math::Ident(text(n)?.trim().to_string()))
            }
            "csymbol" => self.csymbol(n).map(Math::Symbol),
            "apply" => self.apply(n),
            "piecewise" => {
                only_attrs(n, &[])?;
                let mut pieces = Vec::new();
                let mut otherwise = None;
                for c in elements(n) {
                    match self.name_in(c, ns::MATHML)? {
                        "piece" => {
                            let [v, cond] = elements(c).collect::<Vec<_>>().try_into().ok()?;
                            pieces.push((self.expr(v), self.expr(cond)));
                        }
                        "otherwise" if otherwise.is_none() => {
                            let [v] = elements(c).collect::<Vec<_>>().try_into().ok()?;
                            otherwise = Some(Box::new(self.expr(v)));
                        }
                        _ => return None,
                    }
                }
                Some(Math::Piecewise { pieces, otherwise })
            }
            "lambda" => {
                only_attrs(n, &[])?;
                let mut args = Vec::new();
                let mut body = None;
                for c in elements(n) {
                    if body.is_some() {
                        return None;
                    }
                    if self.name_in(c, ns::MATHML) == Some("bvar") {
                        let [ci] = elements(c).collect::<Vec<_>>().try_into().ok()?;
                        if self.name_in(ci, ns::MATHML) != Some("ci") {
                            return None;
                        }
                        args.push(text(ci)?.trim().to_string());
                    } else {
                        body = Some(Box::new(self.expr(c)));
                    }
                }
                Some(Math::Lambda { args, body: body? })
            }
            "semantics" => {
                let mut kids = elements(n);
                let math = Box::new(self.expr(kids.next()?));
                let annotations = kids.map(|c| self.raw(c)).collect();
                Some(Math::Semantics { math, annotations })
            }
            _ => {
                only_attrs(n, &[])?;
                elements(n).next().is_none().then_some(())?;
                Constant::from_tag(name).map(Math::Constant)
            }
        }
    }

    fn cn(&mut self, n: Node<'a, 'i>) -> Option<Math> {
        let mut kind = None;
        let mut units = None;
        for a in n.attributes() {
            match (a.namespace(), a.name()) {
                (None, "type") => kind = Some(a.value()),
                (Some(u), "units") if u == self.core => units = Some(a.value().to_string()),
                _ => return None,
            }
        }
        let parts = sep_parts(n)?;
        let value = match (kind.unwrap_or("real"), parts.as_slice()) {
            ("integer", [v]) => match v.parse() {
                Ok(i) => Number::Integer(i),
                Err(_) => {
                    self.error(n, &format!("\"{v}\" isn't an integer"));
                    Number::Integer(0)
                }
            },
            ("real", [v]) => match parse_double(v) {
                Some(r) => Number::Real(r),
                None => {
                    self.error(n, &format!("\"{v}\" isn't a number"));
                    Number::Real(0.0)
                }
            },
            ("e-notation", [m, e]) => Number::ENotation(parse_double(m)?, e.parse().ok()?),
            ("rational", [num, den]) => Number::Rational(num.parse().ok()?, den.parse().ok()?),
            _ => return None,
        };
        Some(Math::Number { value, units })
    }

    fn csymbol(&mut self, n: Node<'a, 'i>) -> Option<Symbol> {
        let mut url = None;
        for a in n.attributes() {
            match (a.namespace(), a.name()) {
                (None, "definitionURL") => url = Some(a.value()),
                (None, "encoding") if a.value() == "text" => {}
                _ => return None,
            }
        }
        Some(Symbol {
            url: url?.to_string(),
            name: text(n)?.trim().to_string(),
        })
    }

    fn apply(&mut self, n: Node<'a, 'i>) -> Option<Math> {
        only_attrs(n, &[])?;
        let mut kids = elements(n);
        let head = kids.next()?;
        let head_name = self.name_in(head, ns::MATHML)?;
        let rest: Vec<_> = kids.collect();
        let qualifier = |r: &mut Self, q: &str| -> Option<(Option<Box<Math>>, Vec<Node<'a, 'i>>)> {
            let mut found = None;
            let mut args = Vec::new();
            for c in &rest {
                if r.name_in(*c, ns::MATHML) == Some(q) && found.is_none() {
                    let [v] = elements(*c).collect::<Vec<_>>().try_into().ok()?;
                    found = Some(Box::new(r.expr(v)));
                } else {
                    args.push(*c);
                }
            }
            Some((found, args))
        };
        // Qualifiers belong only to root and log.
        let qualifiers = [
            "degree",
            "logbase",
            "bvar",
            "lowlimit",
            "uplimit",
            "condition",
        ];
        let has_qualifier = |r: &Self, rest: &[Node<'a, 'i>]| {
            rest.iter().any(|c| {
                r.name_in(*c, ns::MATHML)
                    .is_some_and(|x| qualifiers.contains(&x))
            })
        };
        match head_name {
            "root" | "log" => {
                let q = if head_name == "root" {
                    "degree"
                } else {
                    "logbase"
                };
                let (found, args) = qualifier(self, q)?;
                if has_qualifier(self, &args) {
                    return None;
                }
                let [arg] = args.try_into().ok()?;
                let arg = Box::new(self.expr(arg));
                Some(if head_name == "root" {
                    Math::Root { degree: found, arg }
                } else {
                    Math::Log { base: found, arg }
                })
            }
            _ => {
                if has_qualifier(self, &rest) {
                    return None;
                }
                let args = rest.iter().map(|c| self.expr(*c)).collect();
                match head_name {
                    "ci" => {
                        only_attrs(head, &[])?;
                        Some(Math::Call(text(head)?.trim().to_string(), args))
                    }
                    "csymbol" => Some(Math::CallSymbol(self.csymbol(head)?, args)),
                    op => {
                        only_attrs(head, &[])?;
                        Some(Math::Apply(Op::from_tag(op)?, args))
                    }
                }
            }
        }
    }
}

fn elements<'a, 'i>(n: Node<'a, 'i>) -> impl Iterator<Item = Node<'a, 'i>> {
    n.children().filter(|c| c.is_element())
}

/// `Some(())` if the element has no attributes other than these.
fn only_attrs(n: Node, allowed: &[&str]) -> Option<()> {
    n.attributes()
        .all(|a| a.namespace().is_none() && allowed.contains(&a.name()))
        .then_some(())
}

/// The text of an element with no child elements.
fn text<'a>(n: Node<'a, '_>) -> Option<String> {
    if elements(n).next().is_some() {
        return None;
    }
    Some(
        n.children()
            .filter(|c| c.is_text())
            .filter_map(|c| c.text())
            .collect(),
    )
}

/// The text of a `<cn>`, split at its `<sep/>`s and trimmed.
fn sep_parts(n: Node) -> Option<Vec<String>> {
    let mut parts = vec![String::new()];
    for c in n.children() {
        if c.is_element() {
            if c.tag_name().name() != "sep" || c.tag_name().namespace() != Some(ns::MATHML) {
                return None;
            }
            parts.push(String::new());
        } else if let Some(t) = c.text().filter(|_| c.is_text()) {
            parts.last_mut().unwrap().push_str(t);
        }
    }
    Some(parts.into_iter().map(|p| p.trim().to_string()).collect())
}

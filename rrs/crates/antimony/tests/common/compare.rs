//! Deep structural comparison of two modules, ported from the Delphi
//! `Tests/uAntModelCompare.pas`, for the round-trip tests.
//!
//! Text comparison would be the wrong test twice over: it would fail on a
//! harmless reordering, and it would pass a writer and reader that agreed
//! about something wrong. Spans are ignored, since the second parse reads
//! different text, and so is the order of symbols, which the writer groups
//! by kind. Everything else the model holds is compared, in order where
//! order is meaning: a writer that dropped something would otherwise pass.
//!
//! Beyond Delphi's comparer: a reaction's `name_generated`, the module's SBO
//! term, and the raw sections (layout and unknown annotations), which
//! Delphi's did not look at. Comments are not compared here: a comment can
//! only come back attached to something that is written, which the writer
//! tests check separately.

use antimony::expr::Expr;
use antimony::model::{ConstraintKind, CvTerm, Metadata, Module, Participant};

/// The ways `b` differs from `a`; empty when they are the same model.
pub fn differences(a: &Module, b: &Module) -> Vec<String> {
    let mut c = Comparer(Vec::new());
    c.modules(a, b);
    c.0
}

struct Comparer(Vec<String>);

fn show(e: Option<&Expr>) -> String {
    e.map(|e| e.to_string()).unwrap_or_else(|| "(none)".into())
}

impl Comparer {
    fn note(&mut self, what: String) {
        self.0.push(what);
    }

    fn same<T: PartialEq + std::fmt::Debug>(&mut self, what: &str, a: &T, b: &T) {
        if a != b {
            self.note(format!("{what}: {a:?} became {b:?}"));
        }
    }

    fn expr(&mut self, what: &str, a: Option<&Expr>, b: Option<&Expr>) {
        let same = match (a, b) {
            (None, None) => true,
            (Some(x), Some(y)) => x.same_tree(y),
            _ => false,
        };
        if !same {
            self.note(format!("{what}: \"{}\" became \"{}\"", show(a), show(b)));
        }
    }

    fn modules(&mut self, a: &Module, b: &Module) {
        self.same("module name", &a.name, &b.name);
        self.same("main-module marker", &a.is_main, &b.is_main);
        self.same("module display name", &a.display_name, &b.display_name);
        self.same("module notes", &a.notes, &b.notes);
        self.same(
            "module conversionFactor",
            &a.conversion_factor,
            &b.conversion_factor,
        );
        self.same("module formal arguments", &a.formal_args, &b.formal_args);
        self.same("module sboTerm", &a.sbo_term, &b.sbo_term);
        self.annotations("module", &a.cv_terms, &a.metadata, &b.cv_terms, &b.metadata);
        self.symbols(a, b);
        self.reactions(a, b);
        self.interactions(a, b);
        self.events(a, b);
        self.functions(a, b);
        self.algebraic_rules(a, b);
        self.constraints(a, b);
        self.objectives(a, b);
        self.units(a, b);
        let raw = |m: &Module| -> Vec<(antimony::model::RawKind, String)> {
            m.unparsed_sections
                .iter()
                .map(|s| (s.kind, s.text.clone()))
                .collect()
        };
        self.same("raw sections", &raw(a), &raw(b));
    }

    fn symbols(&mut self, a: &Module, b: &Module) {
        self.same("symbol count", &a.symbols.len(), &b.symbols.len());
        for s in &a.symbols {
            let Some(o) = b.symbol(&s.name) else {
                self.note(format!("symbol \"{}\" was lost", s.name));
                continue;
            };
            let w = |field: &str| format!("symbol \"{}\" {field}", s.name);
            self.same(&w("kind"), &s.kind, &o.kind);
            self.same(&w("constness"), &s.constness, &o.constness);
            self.same(&w("boundary"), &s.is_boundary, &o.is_boundary);
            self.same(
                &w("substanceOnly"),
                &s.is_substance_only,
                &o.is_substance_only,
            );
            self.same(&w("compartment"), &s.compartment, &o.compartment);
            self.same(&w("units"), &s.units, &o.units);
            self.same(
                &w("conversionFactor"),
                &s.conversion_factor,
                &o.conversion_factor,
            );
            self.same(&w("display name"), &s.display_name, &o.display_name);
            self.same(&w("notes"), &s.notes, &o.notes);
            self.same(&w("sboTerm"), &s.sbo_term, &o.sbo_term);
            // An empty right-hand side is information: losing it is the
            // libantimony defect this must not repeat.
            self.same(&w("empty value"), &s.has_empty_value, &o.has_empty_value);
            self.same(
                &w("rule kind"),
                &s.rule.as_ref().map(|r| r.kind),
                &o.rule.as_ref().map(|r| r.kind),
            );
            self.expr(
                &w("value"),
                s.initial_value.as_ref(),
                o.initial_value.as_ref(),
            );
            self.expr(
                &w("rule"),
                s.rule.as_ref().and_then(|r| r.expr.as_ref()),
                o.rule.as_ref().and_then(|r| r.expr.as_ref()),
            );
            self.annotations(&w(""), &s.cv_terms, &s.metadata, &o.cv_terms, &o.metadata);
        }
    }

    fn side(&mut self, what: &str, a: &[Participant], b: &[Participant]) {
        if a.len() != b.len() {
            self.note(format!(
                "{what}: {} participant(s) became {}",
                a.len(),
                b.len()
            ));
            return;
        }
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            self.same(&format!("{what}[{i}] species"), &x.species, &y.species);
            self.same(
                &format!("{what}[{i}] boundary"),
                &x.is_boundary,
                &y.is_boundary,
            );
            self.expr(
                &format!("{what}[{i}] stoichiometry"),
                x.stoichiometry.as_ref(),
                y.stoichiometry.as_ref(),
            );
        }
    }

    fn reactions(&mut self, a: &Module, b: &Module) {
        if a.reactions.len() != b.reactions.len() {
            self.note(format!(
                "reaction count: {} became {}",
                a.reactions.len(),
                b.reactions.len()
            ));
            return;
        }
        for (x, y) in a.reactions.iter().zip(&b.reactions) {
            let w = |field: &str| format!("reaction \"{}\" {field}", x.name);
            self.same(&w("name"), &x.name, &y.name);
            self.same(&w("generated name"), &x.name_generated, &y.name_generated);
            // Getting this backwards inverts the reaction.
            self.same(&w("kind"), &x.kind, &y.kind);
            self.same(
                &w("empty law"),
                &x.has_empty_kinetic_law,
                &y.has_empty_kinetic_law,
            );
            self.side(&w("reactants"), &x.reactants, &y.reactants);
            self.side(&w("products"), &x.products, &y.products);
            self.side(&w("modifiers"), &x.modifiers, &y.modifiers);
            self.expr(&w("law"), x.kinetic_law.as_ref(), y.kinetic_law.as_ref());
        }
    }

    fn interactions(&mut self, a: &Module, b: &Module) {
        let flat = |m: &Module| -> Vec<String> {
            m.interactions
                .iter()
                .map(|i| format!("{:?} {} {:?} {}", i.name, i.species, i.kind, i.reaction))
                .collect()
        };
        self.same("interactions", &flat(a), &flat(b));
    }

    fn events(&mut self, a: &Module, b: &Module) {
        if a.events.len() != b.events.len() {
            self.note(format!(
                "event count: {} became {}",
                a.events.len(),
                b.events.len()
            ));
            return;
        }
        for (x, y) in a.events.iter().zip(&b.events) {
            let w = |field: &str| format!("event {:?} {field}", x.name);
            self.same(&w("name"), &x.name, &y.name);
            self.expr(&w("trigger"), Some(&x.trigger), Some(&y.trigger));
            self.expr(&w("delay"), x.delay.as_ref(), y.delay.as_ref());
            self.expr(&w("priority"), x.priority.as_ref(), y.priority.as_ref());
            self.expr(
                &w("persistent"),
                x.persistent.as_ref(),
                y.persistent.as_ref(),
            );
            self.expr(&w("t0"), x.initial_value.as_ref(), y.initial_value.as_ref());
            self.expr(
                &w("fromTrigger"),
                x.from_trigger.as_ref(),
                y.from_trigger.as_ref(),
            );
            if x.assignments.len() != y.assignments.len() {
                self.note(w("assignment count"));
                continue;
            }
            for (i, (p, q)) in x.assignments.iter().zip(&y.assignments).enumerate() {
                self.same(&w(&format!("assignment {i} target")), &p.target, &q.target);
                self.expr(
                    &w(&format!("assignment {i}")),
                    p.expr.as_ref(),
                    q.expr.as_ref(),
                );
            }
        }
    }

    fn functions(&mut self, a: &Module, b: &Module) {
        self.same("function count", &a.functions.len(), &b.functions.len());
        for f in &a.functions {
            let Some(o) = b.function(&f.name) else {
                self.note(format!("function \"{}\" was lost", f.name));
                continue;
            };
            self.same(&format!("function \"{}\" args", f.name), &f.args, &o.args);
            self.expr(
                &format!("function \"{}\" body", f.name),
                f.body.as_ref(),
                o.body.as_ref(),
            );
        }
    }

    fn algebraic_rules(&mut self, a: &Module, b: &Module) {
        if a.algebraic_rules.len() != b.algebraic_rules.len() {
            self.note("algebraic rule count changed".into());
            return;
        }
        for (x, y) in a.algebraic_rules.iter().zip(&b.algebraic_rules) {
            self.same("algebraic rule name", &x.name, &y.name);
            self.expr(
                &format!("algebraic rule \"{}\"", x.name),
                x.expr.as_ref(),
                y.expr.as_ref(),
            );
        }
    }

    fn constraints(&mut self, a: &Module, b: &Module) {
        if a.constraints.len() != b.constraints.len() {
            self.note("constraint count changed".into());
            return;
        }
        for (x, y) in a.constraints.iter().zip(&b.constraints) {
            let w = format!("constraint {:?}", x.name);
            self.same(&format!("{w} name"), &x.name, &y.name);
            self.same(&format!("{w} message"), &x.message, &y.message);
            match (&x.kind, &y.kind) {
                (ConstraintKind::Expr(p), ConstraintKind::Expr(q)) => {
                    self.expr(&w, Some(p), Some(q))
                }
                (
                    ConstraintKind::Range {
                        lower: l1,
                        subject: s1,
                        upper: u1,
                        lower_inclusive: li1,
                        upper_inclusive: ui1,
                    },
                    ConstraintKind::Range {
                        lower: l2,
                        subject: s2,
                        upper: u2,
                        lower_inclusive: li2,
                        upper_inclusive: ui2,
                    },
                ) => {
                    self.expr(&format!("{w} lower"), Some(l1), Some(l2));
                    self.expr(&format!("{w} subject"), Some(s1), Some(s2));
                    self.expr(&format!("{w} upper"), Some(u1), Some(u2));
                    self.same(&format!("{w} lower strictness"), li1, li2);
                    self.same(&format!("{w} upper strictness"), ui1, ui2);
                }
                _ => self.note(format!("{w}: range form changed")),
            }
        }
    }

    fn objectives(&mut self, a: &Module, b: &Module) {
        if a.objectives.len() != b.objectives.len() {
            self.note("objective count changed".into());
            return;
        }
        for (x, y) in a.objectives.iter().zip(&b.objectives) {
            let w = format!("objective {:?}", x.name);
            self.same(&format!("{w} name"), &x.name, &y.name);
            self.same(&format!("{w} direction"), &x.maximise, &y.maximise);
            self.expr(&w, Some(&x.expr), Some(&y.expr));
        }
    }

    fn units(&mut self, a: &Module, b: &Module) {
        self.same("unit count", &a.units.len(), &b.units.len());
        for u in &a.units {
            let Some(o) = b.unit(&u.name) else {
                self.note(format!("unit \"{}\" was lost", u.name));
                continue;
            };
            let w = |field: &str| format!("unit \"{}\" {field}", u.name);
            self.same(&w("definition"), &u.definition, &o.definition);
            self.same(&w("display name"), &u.display_name, &o.display_name);
            self.annotations(&w(""), &u.cv_terms, &u.metadata, &o.cv_terms, &o.metadata);
        }
    }

    /// In order: `S1 part "a"` then `S1 part "b"` is a different file from the
    /// reverse, and SBML's `rdf:Bag` keeps the order.
    fn annotations(
        &mut self,
        what: &str,
        terms_a: &[CvTerm],
        meta_a: &Metadata,
        terms_b: &[CvTerm],
        meta_b: &Metadata,
    ) {
        let flat = |terms: &[CvTerm]| -> Vec<String> {
            terms
                .iter()
                .map(|t| format!("{:?} {} {}", t.qualifier, t.spelling, t.resource))
                .collect()
        };
        self.same(&format!("{what} CV terms"), &flat(terms_a), &flat(terms_b));
        self.same(&format!("{what} created"), &meta_a.created, &meta_b.created);
        self.same(
            &format!("{what} modified"),
            &meta_a.modified,
            &meta_b.modified,
        );
        self.same(
            &format!("{what} metadata notes"),
            &meta_a.notes,
            &meta_b.notes,
        );
        self.same(
            &format!("{what} creators"),
            &meta_a.creators,
            &meta_b.creators,
        );
    }
}

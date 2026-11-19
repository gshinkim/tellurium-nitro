//! Answers about a model's quantities: the species, compartments and
//! parameters by kind, and their values at the start (Delphi `uAntFacade`'s
//! query surface, re-shaped; SPEC.md §6).
//!
//! The model's own fields say what the file wrote, with `Option` for what
//! it did not. These methods answer the question a caller usually means,
//! "what is it worth?", with Delphi's defaults where the model is silent: a
//! species with no value is 0, a compartment with no size 1 (0 would make
//! every concentration in it a division by zero), a parameter 0 (and the
//! semantic pass warned about it), a stoichiometry 1. [`Module::value`] is
//! the raw answer, `None` where nothing can be computed.
//!
//! A method asked about an id that is not of its kind answers `None`: that
//! is a mistake in the question, not a fact about the model.

use super::{Event, Module, Participant, RuleKind, Symbol, SymbolKind};
use crate::expr::Expr;

/// How deep [`Module::value`] follows one symbol's value into another's.
/// `a = b; b = a` is not a circular rule, so the semantic pass allows it,
/// but following it would never end. Ten is far past any real chain.
const MAX_RESOLVE_DEPTH: usize = 10;

impl Module {
    /// The species, in order.
    pub fn species(&self) -> impl Iterator<Item = &Symbol> {
        self.symbols_of(SymbolKind::Species)
    }

    /// The species that are not boundary species (`$S`): the ones a
    /// simulation changes. Index `i` here is not index `i` of
    /// [`Module::species`] once a boundary species comes first.
    pub fn floating_species(&self) -> impl Iterator<Item = &Symbol> {
        self.species().filter(|s| !s.is_boundary)
    }

    /// The boundary species (`$S`), held fixed or set by rules. `is_boundary`
    /// is the `$`, not constness: `const S` is a different fact.
    pub fn boundary_species(&self) -> impl Iterator<Item = &Symbol> {
        self.species().filter(|s| s.is_boundary)
    }

    pub fn compartments(&self) -> impl Iterator<Item = &Symbol> {
        self.symbols_of(SymbolKind::Compartment)
    }

    pub fn parameters(&self) -> impl Iterator<Item = &Symbol> {
        self.symbols_of(SymbolKind::Parameter)
    }

    /// The event with this name. An unnamed event has no name to look up;
    /// find it in `events` by position.
    pub fn event(&self, name: &str) -> Option<&Event> {
        self.events.iter().find(|e| e.name.as_deref() == Some(name))
    }

    /// A symbol's value at the start: its assignment rule if it has one
    /// (that holds at all times, the start included), otherwise its initial
    /// value, worked out through the values of the symbols it names. `None`
    /// when there is no such symbol, it has neither, or the expression
    /// cannot be computed (it names `time`, a function, or a symbol with no
    /// value).
    pub fn value(&self, name: &str) -> Option<f64> {
        self.value_at(name, 0)
    }

    fn value_at(&self, name: &str, depth: usize) -> Option<f64> {
        let symbol = self.symbol(name)?;
        let expr = match symbol.assignment_rule() {
            Some(rule) => rule,
            None => symbol.initial_value.as_ref()?,
        };
        self.fold(expr, depth)
    }

    /// An expression's value, with each symbol's value from this module.
    pub fn evaluate(&self, expr: &Expr) -> Option<f64> {
        self.fold(expr, 0)
    }

    fn fold(&self, expr: &Expr, depth: usize) -> Option<f64> {
        if depth >= MAX_RESOLVE_DEPTH {
            return None;
        }
        expr.evaluate_with(&|name| self.value_at(name, depth + 1))
    }

    /// A species' concentration at the start, or `None` if `id` is not a
    /// species. A species stores one number, an amount if it is
    /// `substanceOnly` and a concentration otherwise; the other is worked
    /// out through its compartment's size. No value is 0, and so is the
    /// concentration in a compartment of size 0.
    pub fn concentration(&self, id: &str) -> Option<f64> {
        let species = self.symbol_of(id, SymbolKind::Species)?;
        let value = self.value(id).unwrap_or(0.0);
        if !species.is_substance_only {
            return Some(value);
        }
        let size = self.size_of(species);
        Some(if size == 0.0 { 0.0 } else { value / size })
    }

    /// A species' amount at the start, or `None` if `id` is not a species:
    /// as [`Module::concentration`], the other way round.
    pub fn amount(&self, id: &str) -> Option<f64> {
        let species = self.symbol_of(id, SymbolKind::Species)?;
        let value = self.value(id).unwrap_or(0.0);
        Some(if species.is_substance_only {
            value
        } else {
            value * self.size_of(species)
        })
    }

    /// A compartment's size at the start, 1 when the model gives none, or
    /// `None` if `id` is not a compartment.
    pub fn compartment_size(&self, id: &str) -> Option<f64> {
        self.symbol_of(id, SymbolKind::Compartment)?;
        Some(self.value(id).unwrap_or(1.0))
    }

    /// A parameter's value at the start, 0 when the model gives none (the
    /// semantic pass warns about that), or `None` if `id` is not a
    /// parameter.
    pub fn parameter_value(&self, id: &str) -> Option<f64> {
        self.symbol_of(id, SymbolKind::Parameter)?;
        Some(self.value(id).unwrap_or(0.0))
    }

    /// A participant's stoichiometry: 1 when none is written, a number as
    /// written, a named one (`sr = 2; J1: sr S1 => S2`) through the symbol's
    /// value; 1 too when it cannot be computed. For what was written, see
    /// [`Participant::stoichiometry_text`].
    pub fn stoichiometry(&self, participant: &Participant) -> f64 {
        participant
            .stoichiometry
            .as_ref()
            .and_then(|s| self.evaluate(s))
            .unwrap_or(1.0)
    }

    fn symbol_of(&self, id: &str, kind: SymbolKind) -> Option<&Symbol> {
        self.symbol(id).filter(|s| s.kind == kind)
    }

    /// The size of a species' compartment: 1 when it names none.
    fn size_of(&self, species: &Symbol) -> f64 {
        species
            .compartment
            .as_deref()
            .and_then(|c| self.compartment_size(c))
            .unwrap_or(1.0)
    }
}

impl Symbol {
    /// `x := …`, if that is how this symbol is set.
    pub fn assignment_rule(&self) -> Option<&Expr> {
        self.rule_of(RuleKind::Assignment)
    }

    /// `x' = …`, if that is how this symbol changes.
    pub fn rate_rule(&self) -> Option<&Expr> {
        self.rule_of(RuleKind::Rate)
    }

    fn rule_of(&self, kind: RuleKind) -> Option<&Expr> {
        self.rule
            .as_ref()
            .filter(|r| r.kind == kind)
            .and_then(|r| r.expr.as_ref())
    }
}

impl Participant {
    /// The stoichiometry as written, `1` when none is: `sr` for a named one,
    /// where [`Module::stoichiometry`] gives its value.
    pub fn stoichiometry_text(&self) -> String {
        self.stoichiometry
            .as_ref()
            .map_or_else(|| "1".to_string(), Expr::to_string)
    }
}

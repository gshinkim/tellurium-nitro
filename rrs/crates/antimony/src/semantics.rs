//! The semantic pass (Delphi `uAntSemantics`, Delphi §6).
//!
//! The parser does Delphi §6.1's passes 1 and 2: it reads the statements and
//! declares symbols as it meets them, which is what lets a forward reference
//! work. This does the rest:
//!
//! - **resolve**: every name used is looked up, type inference is completed,
//!   and a name never declared becomes a parameter (Antimony's rule) with a
//!   warning, since a typo silently becoming a parameter is how a model comes
//!   to run and give nonsense (Delphi §6.5);
//! - **validate**: Delphi §6.3's checks.
//!
//! Nothing here fails: every finding is a diagnostic.

use std::collections::{HashMap, HashSet};

use crate::builtins::is_builtin_function;
use crate::diagnostics::{Code, Diagnostics, Severity};
use crate::expr::{Expr, ExprKind};
use crate::messages;
use crate::model::{Constness, ConstraintKind, Document, Module, RuleKind, Symbol, SymbolKind};
use crate::number;
use crate::span::Span;

/// How a source is read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Options {
    /// Report the uninitialised-symbol warnings (`dcUninitialisedSymbol`,
    /// `dcImplicitParameter`) as errors (Delphi `StrictMode`). By default
    /// they are warnings: valid published models must load.
    pub strict: bool,
    /// Give every symbol that has no value one: 0, or 1 for a compartment
    /// (Delphi §6.5.1, `--init-zero`). Off by default, because inventing a
    /// value is what §6.5 exists to stop; the warning is still given, saying
    /// what was filled in. A compartment gets 1 because a concentration is an
    /// amount over its compartment's size, and 0 would divide by zero.
    pub initialise_missing: bool,
}

/// Resolves and checks every module of a parsed document, changing it
/// where inference completes it, and returns what was found.
pub fn analyse(document: &mut Document, options: &Options) -> Diagnostics {
    let mut diagnostics = Diagnostics::new();
    for module in &mut document.modules {
        Analyser::new(module, &mut diagnostics, options).resolve();
    }
    for module in &mut document.modules {
        Analyser::new(module, &mut diagnostics, options).validate();
    }
    diagnostics
}

/// The symbol a name used in an expression refers to, after [`analyse`]:
/// `None` for a built-in function or csymbol (`sin`, `time`), and for a name
/// that is no symbol (a module's formal argument). A function's formal
/// arguments are local to its body, so resolve those against the function.
pub fn resolve<'m>(module: &'m Module, name: &str) -> Option<&'m Symbol> {
    if is_builtin_function(name) {
        return None;
    }
    module.symbol(name)
}

/// Names used in an expression, with where each is used.
type Names = Vec<(String, Span)>;

/// Every name used in `expr`, each occurrence, in order. A call's name is
/// not a reference; its arguments are.
fn references(expr: Option<&Expr>, into: &mut Names) {
    let Some(expr) = expr else {
        return;
    };
    if let ExprKind::Symbol(name) = &expr.kind {
        into.push((name.clone(), expr.span));
        return;
    }
    for child in expr.children() {
        references(Some(child), into);
    }
}

/// Every call in `expr`, in order: (name, number of arguments, span).
fn calls(expr: Option<&Expr>, into: &mut Vec<(String, usize, Span)>) {
    let Some(expr) = expr else {
        return;
    };
    if let ExprKind::Call { name, args } = &expr.kind {
        into.push((name.clone(), args.len(), expr.span));
    }
    for child in expr.children() {
        calls(Some(child), into);
    }
}

struct Analyser<'a> {
    module: &'a mut Module,
    /// Symbol name to position in `module.symbols`.
    index: HashMap<String, usize>,
    diagnostics: &'a mut Diagnostics,
    options: &'a Options,
    /// Names in scope only in the construct being resolved: a function's
    /// formal arguments. Without them `function MM(S, Vm, Km)` would invent
    /// three parameters (the `bvar` leak, Delphi §9.1.4 #4).
    locals: HashSet<String>,
}

impl<'a> Analyser<'a> {
    fn new(
        module: &'a mut Module,
        diagnostics: &'a mut Diagnostics,
        options: &'a Options,
    ) -> Analyser<'a> {
        let index = module
            .symbols
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.clone(), i))
            .collect();
        Analyser {
            module,
            index,
            diagnostics,
            options,
            locals: HashSet::new(),
        }
    }

    fn find(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    fn declare(&mut self, name: &str, span: Span) -> usize {
        if let Some(at) = self.find(name) {
            return at;
        }
        self.module.symbols.push(Symbol::new(name, span));
        let at = self.module.symbols.len() - 1;
        self.index.insert(name.to_string(), at);
        at
    }

    /// A warning that strict mode makes an error.
    fn warn(&mut self, code: Code, span: Span, message: String) {
        let severity = if self.options.strict {
            Severity::Error
        } else {
            code.default_severity()
        };
        self.diagnostics.add_with(code, severity, span, message);
    }

    // --- resolve -------------------------------------------------------------

    /// The symbol a used name refers to, inventing a parameter for a name
    /// never declared (Antimony's rule), marked so §6.5 can warn about it.
    fn resolve_reference(&mut self, name: &str, span: Span) -> Option<usize> {
        if is_builtin_function(name) || self.locals.contains(name) {
            return None;
        }
        if let Some(at) = self.find(name) {
            return Some(at);
        }
        if self.module.formal_args.iter().any(|a| a == name) {
            return None;
        }
        let at = self.declare(name, span);
        let symbol = &mut self.module.symbols[at];
        symbol.is_implicit = true;
        symbol.kind = SymbolKind::Parameter;
        Some(at)
    }

    fn resolve_names(&mut self, names: Names) {
        for (name, span) in names {
            self.resolve_reference(&name, span);
        }
    }

    /// The values and rules of every symbol, including any added meanwhile.
    fn resolve_symbol_values(&mut self) {
        let mut i = 0;
        while i < self.module.symbols.len() {
            let symbol = &self.module.symbols[i];
            let mut names = Vec::new();
            references(symbol.initial_value.as_ref(), &mut names);
            references(
                symbol.rule.as_ref().and_then(|r| r.expr.as_ref()),
                &mut names,
            );
            self.resolve_names(names);
            i += 1;
        }
    }

    fn resolve(&mut self) {
        // Values and rules first, so a symbol used only in a rule is present
        // before the reactions are resolved.
        self.resolve_symbol_values();

        for r in 0..self.module.reactions.len() {
            let reaction = &self.module.reactions[r];
            // Each participant's species, then its stoichiometry, in order.
            let participants: Vec<(String, Span, Names)> = reaction
                .reactants
                .iter()
                .chain(&reaction.products)
                .chain(&reaction.modifiers)
                .map(|p| {
                    let mut names = Vec::new();
                    references(p.stoichiometry.as_ref(), &mut names);
                    (p.species.clone(), p.span, names)
                })
                .collect();
            let mut law_names = Vec::new();
            references(reaction.kinetic_law.as_ref(), &mut law_names);

            for (species, span, stoichiometry_names) in participants {
                if self.find(&species).is_none() {
                    let at = self.declare(&species, span);
                    self.module.symbols[at].kind = SymbolKind::Species;
                }
                self.resolve_names(stoichiometry_names);
            }
            self.resolve_names(law_names);
        }

        for e in 0..self.module.events.len() {
            let event = &self.module.events[e];
            let mut names = Vec::new();
            references(Some(&event.trigger), &mut names);
            for modifier in [
                &event.delay,
                &event.priority,
                &event.persistent,
                &event.initial_value,
                &event.from_trigger,
            ] {
                references(modifier.as_ref(), &mut names);
            }
            self.resolve_names(names);
            // An event assignment's target is something the model contains,
            // as an ordinary assignment's is.
            let assignments: Vec<(String, Span, Names)> = self.module.events[e]
                .assignments
                .iter()
                .map(|a| {
                    let mut names = Vec::new();
                    references(a.expr.as_ref(), &mut names);
                    (a.target.clone(), a.span, names)
                })
                .collect();
            for (target, span, names) in assignments {
                self.resolve_reference(&target, span);
                self.resolve_names(names);
            }
        }

        for f in 0..self.module.functions.len() {
            let function = &self.module.functions[f];
            self.locals = function.args.iter().cloned().collect();
            let mut names = Vec::new();
            references(function.body.as_ref(), &mut names);
            self.resolve_names(names);
            self.locals.clear();
        }

        let mut names = Vec::new();
        for rule in &self.module.algebraic_rules {
            references(rule.expr.as_ref(), &mut names);
        }
        for constraint in &self.module.constraints {
            match &constraint.kind {
                ConstraintKind::Expr(e) => references(Some(e), &mut names),
                ConstraintKind::Range {
                    lower,
                    subject,
                    upper,
                    ..
                } => {
                    references(Some(lower), &mut names);
                    references(Some(subject), &mut names);
                    references(Some(upper), &mut names);
                }
            }
        }
        for objective in &self.module.objectives {
            references(Some(&objective.expr), &mut names);
        }
        self.resolve_names(names);

        // Again: resolving reactions and events may have added symbols.
        self.resolve_symbol_values();
        self.resolve_attachments();

        // Anything still unclassified is a parameter.
        for symbol in &mut self.module.symbols {
            if symbol.kind == SymbolKind::Unknown {
                symbol.kind = SymbolKind::Parameter;
            }
        }
    }

    fn resolve_attachments(&mut self) {
        let mut i = 0;
        while i < self.module.symbols.len() {
            let symbol = &self.module.symbols[i];
            let (span, compartment, factor) = (
                symbol.span,
                symbol.compartment.clone(),
                symbol.conversion_factor.clone(),
            );
            if let Some(name) = compartment {
                // `S1 in cell` says `cell` is a compartment: inference, not
                // invention, so it is not an implicit parameter.
                match self.find(&name) {
                    None => {
                        let at = self.declare(&name, span);
                        self.module.symbols[at].kind = SymbolKind::Compartment;
                    }
                    Some(at) => {
                        let target = &mut self.module.symbols[at];
                        if matches!(target.kind, SymbolKind::Unknown | SymbolKind::Parameter) {
                            target.kind = SymbolKind::Compartment;
                        }
                    }
                }
            }
            // A unit names a unit definition or a base unit, never a symbol,
            // so `has` resolves nothing here.
            if let Some(name) = factor {
                self.resolve_reference(&name, span);
            }
            i += 1;
        }
    }

    // --- validate --------------------------------------------------------------

    fn validate(&mut self) {
        for i in 0..self.module.symbols.len() {
            self.check_rule_conflicts(i);
            self.check_compartment(i);
            self.check_conversion_factor(i);
            self.check_initialisation(i);
        }
        self.check_circular_rules();
        self.check_functions();
    }

    /// A constant cannot be driven by a rule.
    fn check_rule_conflicts(&mut self, i: usize) {
        let symbol = &self.module.symbols[i];
        if let (Some(rule), Constness::Const) = (&symbol.rule, symbol.constness) {
            let message = messages::const_with_rule(&symbol.name, rule.kind.name());
            let span = symbol.span;
            self.diagnostics.add(Code::ConstWithRule, span, message);
        }
    }

    fn check_compartment(&mut self, i: usize) {
        let symbol = &self.module.symbols[i];
        let Some(target) = symbol.compartment.as_deref().and_then(|c| self.find(c)) else {
            return;
        };
        let target = &self.module.symbols[target];
        if target.kind != SymbolKind::Compartment {
            let message = messages::not_a_compartment(&target.name, target.kind.name());
            let span = symbol.span;
            self.diagnostics.add(Code::NotACompartment, span, message);
        }
    }

    /// Delphi §5.8: a conversion factor must be a parameter.
    fn check_conversion_factor(&mut self, i: usize) {
        let symbol = &self.module.symbols[i];
        let Some(name) = symbol.conversion_factor.as_deref() else {
            return;
        };
        if is_builtin_function(name) {
            return;
        }
        let Some(target) = self.find(name) else {
            return;
        };
        let target = &self.module.symbols[target];
        if target.kind != SymbolKind::Parameter {
            let message = messages::bad_conversion_factor(&target.name, target.kind.name());
            let span = symbol.span;
            self.diagnostics
                .add(Code::BadConversionFactor, span, message);
        }
    }

    /// Delphi §6.5: a symbol with no value and no rule. A name that was never
    /// declared (the typo case) and one declared but given nothing are
    /// reported separately. With `initialise_missing` the value is filled in
    /// and the report says so.
    fn check_initialisation(&mut self, i: usize) {
        let symbol = &self.module.symbols[i];
        // Names, not quantities: asking their value is a category error.
        if matches!(
            symbol.kind,
            SymbolKind::Reaction
                | SymbolKind::Event
                | SymbolKind::Function
                | SymbolKind::Submodule
                | SymbolKind::UnitDef
                | SymbolKind::Constraint
                | SymbolKind::Deleted
                | SymbolKind::Interaction
        ) || symbol.is_initialised()
        {
            return;
        }
        let (name, span, implicit) = (symbol.name.clone(), symbol.span, symbol.is_implicit);
        let is_compartment = symbol.kind == SymbolKind::Compartment;
        let fill = self.options.initialise_missing;
        if implicit {
            // A likely typo, reported the same way with or without filling in.
            self.warn(
                Code::ImplicitParameter,
                span,
                messages::implicit_parameter(&name),
            );
        } else if !fill {
            self.warn(
                Code::UninitialisedSymbol,
                span,
                messages::uninitialised_symbol(&name),
            );
        }
        if !fill {
            return;
        }
        let value = if is_compartment { 1.0 } else { 0.0 };
        self.module.symbols[i].initial_value = Some(Expr::number(value, span));
        if !implicit {
            self.warn(
                Code::UninitialisedSymbol,
                span,
                messages::initialised_to_default(&name, &number::format(value)),
            );
        }
    }

    /// The symbols with an assignment rule that symbol `i`'s rule uses, in
    /// order.
    fn rule_dependencies(&self, i: usize) -> Vec<usize> {
        let mut names = Vec::new();
        let rule = self.module.symbols[i].rule.as_ref();
        references(rule.and_then(|r| r.expr.as_ref()), &mut names);
        names
            .iter()
            .filter(|(name, _)| !is_builtin_function(name))
            .filter_map(|(name, _)| self.find(name))
            .filter(|&at| {
                self.module.symbols[at]
                    .rule
                    .as_ref()
                    .is_some_and(|r| r.kind == RuleKind::Assignment)
            })
            .collect()
    }

    /// Delphi §6.3: assignment rules must not depend on themselves, directly
    /// or through others. A rate rule may use its own subject: it is a
    /// differential equation, not a cycle. A depth-first search with three
    /// colours, written as a loop so a long chain of rules cannot exhaust the
    /// stack. One report per module.
    fn check_circular_rules(&mut self) {
        const DONE: u8 = 2;
        const OPEN: u8 = 1;
        let mut colour = vec![0u8; self.module.symbols.len()];
        for start in 0..self.module.symbols.len() {
            let is_assignment = self.module.symbols[start]
                .rule
                .as_ref()
                .is_some_and(|r| r.kind == RuleKind::Assignment);
            if !is_assignment || colour[start] == DONE {
                continue;
            }
            // (symbol, its dependencies, how many have been followed)
            let mut stack = vec![(start, self.rule_dependencies(start), 0usize)];
            colour[start] = OPEN;
            while let Some((symbol, dependencies, next)) = stack.last_mut() {
                if *next < dependencies.len() {
                    let target = dependencies[*next];
                    *next += 1;
                    match colour[target] {
                        DONE => {}
                        OPEN => {
                            // Name the whole path: far more use than "a is
                            // circular".
                            let mut path: Vec<&str> = stack
                                .iter()
                                .map(|(s, _, _)| self.module.symbols[*s].name.as_str())
                                .collect();
                            path.push(&self.module.symbols[target].name);
                            let message = messages::circular_rule(&path.join(" -> "));
                            let span = self.module.symbols[target].span;
                            self.diagnostics.add(Code::CircularRule, span, message);
                            return;
                        }
                        _ => {
                            colour[target] = OPEN;
                            let dependencies = self.rule_dependencies(target);
                            stack.push((target, dependencies, 0));
                        }
                    }
                } else {
                    colour[*symbol] = DONE;
                    stack.pop();
                }
            }
        }
    }

    /// Delphi §6.3: a function may not call itself (SBML functions cannot
    /// recurse), and a user function must be called with as many arguments
    /// as it takes. Built-ins such as `piecewise` are variadic.
    fn check_functions(&mut self) {
        let mut found = Vec::new();
        for function in &self.module.functions {
            let mut inside = Vec::new();
            calls(function.body.as_ref(), &mut inside);
            if inside.iter().any(|(name, _, _)| *name == function.name) {
                found.push((
                    Code::RecursiveFunction,
                    function.span,
                    messages::recursive_function(&function.name),
                ));
            }
        }

        let mut all = Vec::new();
        for function in &self.module.functions {
            calls(function.body.as_ref(), &mut all);
        }
        for symbol in &self.module.symbols {
            calls(symbol.initial_value.as_ref(), &mut all);
            calls(symbol.rule.as_ref().and_then(|r| r.expr.as_ref()), &mut all);
        }
        for reaction in &self.module.reactions {
            calls(reaction.kinetic_law.as_ref(), &mut all);
        }
        for event in &self.module.events {
            calls(Some(&event.trigger), &mut all);
            calls(event.delay.as_ref(), &mut all);
            for assignment in &event.assignments {
                calls(assignment.expr.as_ref(), &mut all);
            }
        }
        for (name, supplied, span) in all {
            if let Some(function) = self.module.function(&name) {
                if function.args.len() != supplied {
                    found.push((
                        Code::ArityMismatch,
                        span,
                        messages::arity_mismatch(&name, function.args.len(), supplied),
                    ));
                }
            }
        }
        for (code, span, message) in found {
            self.diagnostics.add(code, span, message);
        }
    }
}

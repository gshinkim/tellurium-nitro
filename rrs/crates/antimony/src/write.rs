//! The Antimony writer (Delphi `uAntWriter`).
//!
//! **Round-trip closure is a hard invariant** (Delphi §9.1.4 #1): anything
//! written here reads back to the same model. Two things follow:
//!
//! - every degenerate form the parser accepts is written correctly too: an
//!   empty participant list, a kinetic law written empty, an empty
//!   right-hand side, an event with no assignments;
//! - nothing is written that would read back differently: `--` (two minus
//!   signs read as the DNA operator), a bare `-` before a participant (which
//!   is rejected), a quote or line break inside a string.
//!
//! Statements are grouped by kind, so the first write of a file reorders
//! it; after that, writing again gives the same text. Comments are kept, and
//! no section headings are invented.

use crate::expr::{Expr, ExprKind};
use crate::model::{
    Constness, ConstraintKind, CvTerm, Document, LayoutValue, Metadata, Module, Participant,
    ReactionKind, RuleKind, SymbolKind,
};
use crate::number;

/// How the text is laid out. Every choice reads back to the same model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WriteOptions {
    /// Indent inside a `model … end` block and a function body.
    pub indent: String,
    /// Write `species` and `compartment` declarations even where use would
    /// imply them. On by default: a species used nowhere else would
    /// otherwise not come back as a species.
    pub explicit_declarations: bool,
    /// A blank line between groups of statements. Cosmetic.
    pub blank_line_between_groups: bool,
}

impl Default for WriteOptions {
    fn default() -> WriteOptions {
        WriteOptions {
            indent: "  ".into(),
            explicit_declarations: true,
            blank_line_between_groups: true,
        }
    }
}

/// Every module of a document as Antimony, with the default layout.
pub fn write_document(document: &Document) -> String {
    write_document_with(document, &WriteOptions::default())
}

/// One module as Antimony, with the default layout.
pub fn write_module(module: &Module) -> String {
    write_module_with(module, &WriteOptions::default())
}

pub fn write_document_with(document: &Document, options: &WriteOptions) -> String {
    document
        .modules
        .iter()
        .map(|m| write_module_with(m, options))
        .collect::<Vec<_>>()
        .join("\n")
}

pub fn write_module_with(module: &Module, options: &WriteOptions) -> String {
    let mut writer = Writer {
        options,
        out: String::new(),
        indent: String::new(),
    };
    writer.module(module);
    writer.out
}

/// A string as Antimony source: quoted, with `\`, `"` and line breaks
/// escaped, the exact inverse of what the lexer decodes. A real annotation
/// holds `Institute … "López-Neyra", CSIC.` over two lines.
pub fn quoted(text: &str) -> String {
    let escaped = text
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace("\r\n", "\\n")
        .replace('\n', "\\n")
        .replace('\r', "\\r");
    format!("\"{escaped}\"")
}

fn expr_text(expr: Option<&Expr>) -> String {
    expr.map(Expr::to_string).unwrap_or_default()
}

/// A participant: `[stoichiometry] [$]species`. Only a number or a name may
/// stand as a stoichiometry unbracketed, and never one starting with `-`,
/// which the parser rejects.
fn participant(p: &Participant) -> String {
    let mut out = String::new();
    if let Some(s) = &p.stoichiometry {
        let written = s.to_string();
        let plain = matches!(s.kind, ExprKind::Number { .. } | ExprKind::Symbol(_));
        if plain && !written.starts_with('-') {
            out.push_str(&written);
        } else {
            out.push_str(&format!("({written})"));
        }
        out.push(' ');
    }
    if p.is_boundary {
        out.push('$');
    }
    out.push_str(&p.species);
    out
}

fn side(participants: &[Participant]) -> String {
    participants
        .iter()
        .map(participant)
        .collect::<Vec<_>>()
        .join(" + ")
}

fn arrow(kind: ReactionKind) -> &'static str {
    // `->` is REVERSIBLE and `=>` IRREVERSIBLE (Delphi §9.1.1).
    match kind {
        ReactionKind::Reversible => "->",
        ReactionKind::Irreversible => "=>",
        ReactionKind::BiDirectional => "<->",
        ReactionKind::Inhibition => "-|",
        ReactionKind::Activation => "-o",
    }
}

struct Writer<'o> {
    options: &'o WriteOptions,
    out: String,
    indent: String,
}

impl Writer<'_> {
    fn line(&mut self, text: &str) {
        if !text.is_empty() {
            self.out.push_str(&self.indent);
            self.out.push_str(text);
        }
        self.out.push('\n');
    }

    fn blank(&mut self) {
        if self.options.blank_line_between_groups {
            self.out.push('\n');
        }
    }

    /// A blank line after a group, if anything was written.
    fn end_group(&mut self, wrote: bool) {
        if wrote {
            self.blank();
        }
    }

    /// The comments above a statement, each line as it was.
    fn comment(&mut self, comment: Option<&String>) {
        for line in comment.into_iter().flat_map(|c| c.split('\n')) {
            let line = line.trim();
            if !line.is_empty() {
                self.line(line);
            }
        }
    }

    /// `<subject> notes …`. A note with a line break or a quote is written in
    /// a ```` ``` ```` fence: the quoted form is one line. Its lines are written
    /// unindented and untrimmed, since the fence is raw and anything added
    /// would become part of the note.
    fn notes(&mut self, subject: &str, notes: &str) {
        if !notes.contains('\n') && !notes.contains('"') {
            self.line(&format!("{subject} notes {};", quoted(notes)));
            return;
        }
        self.line(&format!("{subject} notes ```"));
        for line in notes.split('\n') {
            self.out.push_str(line);
            self.out.push('\n');
        }
        self.out.push_str("```\n");
    }

    fn module(&mut self, module: &Module) {
        if let Some(name) = &module.name {
            let star = if module.is_main { "*" } else { "" };
            let line = format!("model {star}{name}({})", module.formal_args.join(", "));
            self.line(&line);
            self.indent = self.options.indent.clone();
        }
        if let Some(notes) = &module.notes {
            self.notes("model", notes);
        }
        if let Some(factor) = &module.conversion_factor {
            self.line(&format!("model conversionFactor {factor};"));
        }
        if let Some(sbo) = module.sbo_term {
            self.line(&format!("model.sboTerm = {sbo}"));
        }

        self.declarations(module);
        self.functions(module);
        self.reactions(module);
        self.interactions(module);
        self.events(module);
        self.algebraic_rules(module);
        self.objectives(module);
        self.constraints(module);
        self.rules_and_values(module);
        self.qualifiers(module);
        // Units before attachments: `time_unit is "time"` names a unit only
        // if the unit has been read already; otherwise it would read back as
        // a parameter.
        self.units(module);
        self.attachments(module);
        self.annotations(module);
        self.raw_sections(module);
        self.layout(module);

        if let Some(name) = &module.name {
            self.indent.clear();
            self.line("end");
            // The module's display name is written outside the block, where
            // the parser knows its subject is the module.
            if let Some(display) = &module.display_name {
                self.line("");
                self.line(&format!("{name} is {}", quoted(display)));
            }
        }
    }

    fn declarations(&mut self, module: &Module) {
        if !self.options.explicit_declarations {
            return;
        }
        let mut wrote = false;
        for s in &module.symbols {
            let keyword = match s.kind {
                SymbolKind::Species => "species",
                SymbolKind::Compartment => "compartment",
                _ => continue,
            };
            self.comment(s.comment.as_ref());
            let mut text = String::new();
            if s.is_substance_only {
                text.push_str("substanceOnly ");
            }
            text.push_str(keyword);
            text.push(' ');
            if s.is_boundary {
                text.push('$');
            }
            text.push_str(&s.name);
            if let Some(c) = &s.compartment {
                text.push_str(&format!(" in {c}"));
            }
            text.push(';');
            self.line(&text);
            wrote = true;
        }
        self.end_group(wrote);
    }

    fn functions(&mut self, module: &Module) {
        for f in &module.functions {
            self.comment(f.comment.as_ref());
            self.line(&format!("function {}({})", f.name, f.args.join(", ")));
            self.indent.push_str(&self.options.indent);
            self.line(&expr_text(f.body.as_ref()));
            let keep = self.indent.len() - self.options.indent.len();
            self.indent.truncate(keep);
            self.line("end");
        }
        self.end_group(!module.functions.is_empty());
    }

    fn reactions(&mut self, module: &Module) {
        for r in &module.reactions {
            self.comment(r.comment.as_ref());
            let label = if r.name_generated {
                String::new()
            } else {
                format!("{}: ", r.name)
            };
            let mut text = format!(
                "{label}{} {} {};",
                side(&r.reactants),
                arrow(r.kind),
                side(&r.products)
            );
            // A law written empty is `; ;`, not the same as no law.
            match &r.kinetic_law {
                Some(law) => text.push_str(&format!(" {law};")),
                None if r.has_empty_kinetic_law => text.push_str(" ;"),
                None => {}
            }
            self.line(&text);
        }
        self.end_group(!module.reactions.is_empty());
    }

    fn interactions(&mut self, module: &Module) {
        for i in &module.interactions {
            self.comment(i.comment.as_ref());
            let label = i
                .name
                .as_ref()
                .map(|n| format!("{n}: "))
                .unwrap_or_default();
            self.line(&format!(
                "{label}{} {} {};",
                i.species,
                i.kind.arrow(),
                i.reaction
            ));
        }
        self.end_group(!module.interactions.is_empty());
    }

    fn events(&mut self, module: &Module) {
        for e in &module.events {
            self.comment(e.comment.as_ref());
            let mut text = e
                .name
                .as_ref()
                .map(|n| format!("{n}: "))
                .unwrap_or_default();
            text.push_str("at ");
            // The delay comes before the trigger, joined by `after`.
            if let Some(delay) = &e.delay {
                text.push_str(&format!("{delay} after "));
            }
            text.push_str(&e.trigger.to_string());
            for (word, value) in [
                ("priority", &e.priority),
                ("persistent", &e.persistent),
                ("t0", &e.initial_value),
                ("fromTrigger", &e.from_trigger),
            ] {
                if let Some(value) = value {
                    text.push_str(&format!(", {word} = {value}"));
                }
            }
            text.push_str(": ");
            let assignments: Vec<String> = e
                .assignments
                .iter()
                .map(|a| format!("{} = {}", a.target, expr_text(a.expr.as_ref())))
                .collect();
            text.push_str(&assignments.join(", "));
            // `E0: at time > 5.5: ;` must read back too.
            text.push(';');
            self.line(&text);
        }
        self.end_group(!module.events.is_empty());
    }

    fn algebraic_rules(&mut self, module: &Module) {
        for r in &module.algebraic_rules {
            self.comment(r.comment.as_ref());
            let label = if r.name.is_empty() {
                String::new()
            } else {
                format!("{}: ", r.name)
            };
            // `_alg0: 0 = ;` (case 01244) is written empty too.
            self.line(&format!("{label}0 = {};", expr_text(r.expr.as_ref())));
        }
        self.end_group(!module.algebraic_rules.is_empty());
    }

    fn objectives(&mut self, module: &Module) {
        for o in &module.objectives {
            self.comment(o.comment.as_ref());
            let label = o
                .name
                .as_ref()
                .map(|n| format!("{n}: "))
                .unwrap_or_default();
            let word = if o.maximise { "maximize" } else { "minimize" };
            self.line(&format!("{label}{word} {};", o.expr));
        }
        self.end_group(!module.objectives.is_empty());
    }

    fn constraints(&mut self, module: &Module) {
        for c in &module.constraints {
            self.comment(c.comment.as_ref());
            let mut text = c
                .name
                .as_ref()
                .map(|n| format!("constraint {n}: "))
                .unwrap_or_default();
            match &c.kind {
                ConstraintKind::Expr(e) => text.push_str(&e.to_string()),
                // The range is written as the chained comparison it was.
                ConstraintKind::Range {
                    lower,
                    subject,
                    upper,
                    lower_inclusive,
                    upper_inclusive,
                } => {
                    let op = |inclusive: bool| if inclusive { "<=" } else { "<" };
                    text.push_str(&format!(
                        "{lower} {} {subject} {} {upper}",
                        op(*lower_inclusive),
                        op(*upper_inclusive)
                    ));
                }
            }
            self.line(&text);
        }
        self.end_group(!module.constraints.is_empty());
    }

    fn rules_and_values(&mut self, module: &Module) {
        let mut wrote = false;
        for s in &module.symbols {
            // Written by their own statements.
            if matches!(
                s.kind,
                SymbolKind::Reaction
                    | SymbolKind::Event
                    | SymbolKind::Function
                    | SymbolKind::UnitDef
            ) {
                continue;
            }
            if let Some(value) = &s.initial_value {
                // A species' or compartment's comment went with its
                // declaration.
                if !matches!(s.kind, SymbolKind::Species | SymbolKind::Compartment) {
                    self.comment(s.comment.as_ref());
                }
                self.line(&format!("{} = {value};", s.name));
                wrote = true;
            }
            // `p = 3;` then `p = ;` is the corpus's own shape (case 01234):
            // both are written, even with a value.
            if s.has_empty_value {
                self.line(&format!("{} = ;", s.name));
                wrote = true;
            }
            if let Some(rule) = &s.rule {
                let expr = expr_text(rule.expr.as_ref());
                let line = match rule.kind {
                    RuleKind::Assignment => format!("{} := {expr};", s.name),
                    RuleKind::Rate => format!("{}' = {expr};", s.name),
                };
                self.line(&line);
                wrote = true;
            }
        }
        self.end_group(wrote);
    }

    /// `const a, b;` and `var x;`, one statement each, as libantimony writes.
    /// A `$` written on something that is not a species (`$x = 4`) is kept
    /// as `const $x`; on a species the declaration carries it.
    fn qualifiers(&mut self, module: &Module) {
        let names = |constness: Constness| -> Vec<String> {
            module
                .symbols
                .iter()
                .filter(|s| s.constness == constness)
                .map(|s| {
                    if s.is_boundary && s.kind != SymbolKind::Species {
                        format!("${}", s.name)
                    } else {
                        s.name.clone()
                    }
                })
                .collect()
        };
        let (consts, vars) = (names(Constness::Const), names(Constness::Var));
        if !consts.is_empty() {
            self.line(&format!("const {};", consts.join(", ")));
        }
        if !vars.is_empty() {
            self.line(&format!("var {};", vars.join(", ")));
        }
        self.end_group(!consts.is_empty() || !vars.is_empty());
    }

    fn units(&mut self, module: &Module) {
        for u in &module.units {
            self.comment(u.comment.as_ref());
            // The definition is the source text, so this is exact.
            self.line(&format!("unit {} = {};", u.name, u.definition));
        }
        self.end_group(!module.units.is_empty());
    }

    fn attachments(&mut self, module: &Module) {
        let mut wrote = false;
        for s in &module.symbols {
            // A species' or compartment's place is in its declaration.
            if let Some(c) = &s.compartment {
                if !matches!(s.kind, SymbolKind::Species | SymbolKind::Compartment) {
                    self.line(&format!("{} in {c};", s.name));
                    wrote = true;
                }
            }
            if let Some(units) = &s.units {
                self.line(&format!("{} has {units};", s.name));
                wrote = true;
            }
            if let Some(factor) = &s.conversion_factor {
                self.line(&format!("{} conversionFactor {factor};", s.name));
                wrote = true;
            }
            if let Some(sbo) = s.sbo_term {
                self.line(&format!("{}.sboTerm = {sbo}", s.name));
                wrote = true;
            }
            if let Some(notes) = &s.notes {
                self.notes(&s.name, notes);
                wrote = true;
            }
        }
        for s in &module.symbols {
            if let Some(display) = &s.display_name {
                self.line(&format!("{} is {};", s.name, quoted(display)));
                wrote = true;
            }
        }
        for u in &module.units {
            if let Some(display) = &u.display_name {
                self.line(&format!("{} is {};", u.name, quoted(display)));
                wrote = true;
            }
        }
        self.end_group(wrote);
    }

    /// One subject's CV terms and Dublin Core. A module's subject is the word
    /// `model`, as the corpus and the reference write it.
    fn annotations_for(&mut self, subject: &str, terms: &[CvTerm], metadata: &Metadata) -> bool {
        let mut wrote = false;
        // Consecutive terms with the same qualifier and spelling were one
        // comma-separated statement, and are written as one.
        let mut i = 0;
        while i < terms.len() {
            let first = &terms[i];
            let keyword = if first.spelling.is_empty() {
                first.qualifier.keyword()
            } else {
                first.spelling.as_str()
            };
            let mut resources = vec![quoted(&first.resource)];
            let mut j = i + 1;
            while j < terms.len()
                && terms[j].qualifier == first.qualifier
                && terms[j].spelling == first.spelling
            {
                resources.push(quoted(&terms[j].resource));
                j += 1;
            }
            self.line(&format!("{subject} {keyword} {};", resources.join(", ")));
            wrote = true;
            i = j;
        }

        let mut field = |writer: &mut Self, qualifier: String, value: &str| {
            writer.line(&format!("{subject} {qualifier} {};", quoted(value)));
            wrote = true;
        };
        if let Some(created) = &metadata.created {
            field(self, "created".into(), created);
        }
        // One time is `modified`; several are numbered.
        for (n, time) in metadata.modified.iter().enumerate() {
            let word = if metadata.modified.len() == 1 {
                "modified".to_string()
            } else {
                format!("modified{}", n + 1)
            };
            field(self, word, time);
        }
        for creator in &metadata.creators {
            // `creator.email`, or `creator2.email` when the file numbered them.
            let prefix = if creator.index > 0 {
                format!("creator{}", creator.index)
            } else {
                "creator".to_string()
            };
            for (name, value) in [
                ("name", &creator.name),
                ("givenName", &creator.given_name),
                ("familyName", &creator.family_name),
                ("organization", &creator.organization),
                ("email", &creator.email),
            ] {
                if let Some(value) = value {
                    field(self, format!("{prefix}.{name}"), value);
                }
            }
        }
        wrote
    }

    fn annotations(&mut self, module: &Module) {
        let mut wrote = self.annotations_for("model", &module.cv_terms, &module.metadata);
        for s in &module.symbols {
            wrote |= self.annotations_for(&s.name, &s.cv_terms, &s.metadata);
        }
        // Units are not symbols; leaving them out would lose a unit's
        // annotation on the first save.
        for u in &module.units {
            wrote |= self.annotations_for(&u.name, &u.cv_terms, &u.metadata);
        }
        self.end_group(wrote);
    }

    /// Layout statements, in the order read: `S1.pos = {28, 35}`.
    fn layout(&mut self, module: &Module) {
        for statement in &module.layout {
            self.comment(statement.comment.as_ref());
            let value = match &statement.value {
                LayoutValue::Number(n) => number::format(*n),
                LayoutValue::Point(x, y) => {
                    format!("{{{}, {}}}", number::format(*x), number::format(*y))
                }
                LayoutValue::List(items) => format!("{{{}}}", items.join(", ")),
                LayoutValue::Word(word) => word.clone(),
                LayoutValue::Text(text) => quoted(text),
            };
            self.line(&format!("{} = {value}", statement.path.join(".")));
        }
        self.end_group(!module.layout.is_empty());
    }

    /// Unknown annotations, verbatim.
    fn raw_sections(&mut self, module: &Module) {
        for section in &module.unparsed_sections {
            self.line(&section.text);
        }
        self.end_group(!module.unparsed_sections.is_empty());
    }
}

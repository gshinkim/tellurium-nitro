//! The statement parser (Delphi `uAntParser`): recursive descent over the
//! lexer, building a [`Document`].
//!
//! Names are not resolved. Every reference is recorded by name, because
//! Antimony allows forward references (`S1 = k*2` before `k = 3`); the
//! semantic pass binds them.
//!
//! Parsing never stops at the first error. A statement that fails
//! resynchronises at the next `;` or line break (a reaction at the next line
//! break, a block at its `end`), so every problem in a file is reported at
//! once.

use std::collections::HashMap;

use crate::diagnostics::{Code, Diagnostics, Severity};
use crate::expr::{BinaryOp, Expr, ExprKind, ExprParser};
use crate::lexer::{Keyword, Lexer, Token, TokenKind};
use crate::messages;
use crate::model::{
    AlgebraicRule, Constness, Constraint, ConstraintKind, CvTerm, Document, Event, EventAssignment,
    Function, Interaction, InteractionKind, LayoutStatement, LayoutValue, Metadata, Module,
    Objective, Participant, Qualifier, RawKind, RawSection, Reaction, ReactionKind, Rule, RuleKind,
    Symbol, SymbolKind, UnitDef,
};
use crate::number;
use crate::span::Span;

/// Reads `source` into a document. Never fails: problems are in the
/// diagnostics, in the order they were found.
pub fn parse_document(source: &str) -> (Document, Diagnostics) {
    Parser::new(source).parse()
}

/// The keywords that open a `… end` block, for skipping one.
const BLOCK_OPENERS: [Keyword; 3] = [Keyword::Model, Keyword::Module, Keyword::Function];

/// The keywords that may open a type declaration.
const DECLARATION_KEYWORDS: [Keyword; 4] = [
    Keyword::Species,
    Keyword::Compartment,
    Keyword::Formula,
    Keyword::SubstanceOnly,
];

/// How far ahead the statement classifiers look for an arrow or an `=`.
const LOOKAHEAD: usize = 4096;

/// Where an annotation is filed.
#[derive(Clone, Copy)]
enum Subject {
    Module,
    Unit(usize),
    Symbol(usize),
}

/// The relational operators, which a constraint may chain.
fn is_relational(op: BinaryOp) -> bool {
    matches!(
        op,
        BinaryOp::Less | BinaryOp::LessEqual | BinaryOp::Greater | BinaryOp::GreaterEqual
    )
}

fn is_inclusive(op: BinaryOp) -> bool {
    matches!(op, BinaryOp::LessEqual | BinaryOp::GreaterEqual)
}

/// The reaction kind an arrow token stands for. `-|` and `-o` count, so a
/// reaction-shaped statement using one is still recognised; interactions
/// are checked first.
fn arrow_kind(kind: TokenKind) -> Option<ReactionKind> {
    Some(match kind {
        // `->` is REVERSIBLE and `=>` IRREVERSIBLE (Delphi §9.1.1).
        TokenKind::ArrowRev => ReactionKind::Reversible,
        TokenKind::ArrowIrrev => ReactionKind::Irreversible,
        TokenKind::ArrowBiDir => ReactionKind::BiDirectional,
        TokenKind::Inhibits => ReactionKind::Inhibition,
        TokenKind::Activates => ReactionKind::Activation,
        _ => return None,
    })
}

/// `creator1` is `creator` and 1, `created` is `created` and 0. The number
/// is a grouping key.
fn split_trailing_number(word: &str) -> (&str, u32) {
    let base = word.trim_end_matches(|c: char| c.is_ascii_digit());
    let index = word[base.len()..].parse().unwrap_or(0);
    (base, index)
}

fn non_empty(text: &str) -> Option<String> {
    (!text.is_empty()).then(|| text.to_string())
}

/// The state of one parse.
pub struct Parser<'a> {
    source: &'a str,
    lexer: Lexer<'a>,
    document: Document,
    /// The module statements land in, an index into `document.modules`.
    module: usize,
    /// Symbol name to position in the current module's `symbols`.
    index: HashMap<String, usize>,
    /// The comment above the statement being parsed, for whatever the
    /// statement produces.
    comment: Option<String>,
    /// The subject of the statement being parsed.
    subject_span: Span,
}

impl<'a> Parser<'a> {
    pub fn new(source: &'a str) -> Parser<'a> {
        Parser {
            source,
            lexer: Lexer::new(source),
            document: Document::default(),
            module: 0,
            index: HashMap::new(),
            comment: None,
            subject_span: Span::default(),
        }
    }

    /// Parses the whole source. A file need not declare a module: the
    /// document starts with one unnamed main module, which the first
    /// `model` header names.
    pub fn parse(mut self) -> (Document, Diagnostics) {
        self.skip_terminators();
        while self.kind(0) != TokenKind::Eof {
            self.parse_statement();
            self.skip_terminators();
        }
        (self.document, self.lexer.take_diagnostics())
    }

    // --- tokens ------------------------------------------------------------

    fn peek(&mut self, ahead: usize) -> &Token {
        self.lexer.peek(ahead)
    }

    fn kind(&mut self, ahead: usize) -> TokenKind {
        self.lexer.peek(ahead).kind
    }

    fn span(&mut self, ahead: usize) -> Span {
        self.lexer.peek(ahead).span
    }

    fn text(&mut self, ahead: usize) -> String {
        self.lexer.peek(ahead).text.clone()
    }

    fn is_keyword(&mut self, ahead: usize, keyword: Keyword) -> bool {
        self.lexer.peek(ahead).is_keyword(keyword)
    }

    fn is_any_keyword(&mut self, ahead: usize, keywords: &[Keyword]) -> bool {
        self.lexer.peek(ahead).is_any_keyword(keywords)
    }

    fn describe(&mut self, ahead: usize) -> String {
        self.lexer.peek(ahead).describe()
    }

    fn next(&mut self) -> Token {
        self.lexer.next_token()
    }

    fn try_consume(&mut self, kind: TokenKind) -> bool {
        self.lexer.try_consume(kind).is_some()
    }

    fn try_consume_keyword(&mut self, keyword: Keyword) -> bool {
        self.lexer.try_consume_keyword(keyword)
    }

    fn current_span(&self) -> Span {
        self.lexer.current().span
    }

    // --- diagnostics ---------------------------------------------------------

    fn error(&mut self, code: Code, span: Span, message: impl Into<String>) {
        self.lexer.diagnostics_mut().add(code, span, message);
    }

    /// `Unexpected <next token>, expected <what>`.
    fn unexpected(&mut self, expected: &str) {
        let (span, found) = (self.span(0), self.describe(0));
        self.error(
            Code::UnexpectedToken,
            span,
            messages::unexpected_token(&found, expected),
        );
    }

    fn expected_identifier(&mut self) {
        let (span, found) = (self.span(0), self.describe(0));
        self.error(
            Code::ExpectedIdentifier,
            span,
            messages::expected_identifier(&found),
        );
    }

    fn syntax_error(&mut self, span: Span, detail: &str) {
        self.error(Code::SyntaxError, span, messages::syntax_error(detail));
    }

    // --- recovery ------------------------------------------------------------

    fn skip_terminators(&mut self) {
        while matches!(self.kind(0), TokenKind::Semicolon | TokenKind::Newline) {
            self.next();
        }
    }

    fn at_statement_end(&mut self) -> bool {
        matches!(
            self.kind(0),
            TokenKind::Semicolon | TokenKind::Newline | TokenKind::Eof
        )
    }

    /// A statement ends with `;`, a line break or the end of the file. The
    /// corpus needs all three: `case00001 is "case00001"` ends almost every
    /// file with no semicolon.
    fn expect_statement_end(&mut self) {
        match self.kind(0) {
            TokenKind::Semicolon | TokenKind::Newline => {
                self.next();
            }
            TokenKind::Eof => {}
            _ => {
                self.unexpected("end of statement");
                self.resynchronise();
            }
        }
    }

    /// Statement-level recovery: on to the next terminator, so one bad
    /// statement costs one diagnostic.
    fn resynchronise(&mut self) {
        while !self.at_statement_end() {
            self.next();
        }
        if matches!(self.kind(0), TokenKind::Semicolon | TokenKind::Newline) {
            self.next();
        }
    }

    /// Recovery to the next line break. A reaction has `;` inside it
    /// (`J0: A -> B; law;`), so stopping at the first one would leave its law
    /// to be read as a statement of its own.
    fn resynchronise_line(&mut self) {
        while !matches!(self.kind(0), TokenKind::Newline | TokenKind::Eof) {
            self.next();
        }
        if self.kind(0) == TokenKind::Newline {
            self.next();
        }
    }

    /// Recovery to the `end` of a block, so its body is not read as
    /// statements. The opening keyword is consumed first, or it would count
    /// as a nested block.
    fn skip_to_block_end(&mut self) {
        self.next();
        let mut depth = 1;
        while self.kind(0) != TokenKind::Eof {
            if self.is_any_keyword(0, &BLOCK_OPENERS) {
                depth += 1;
            } else if self.is_keyword(0, Keyword::End) {
                depth -= 1;
                if depth == 0 {
                    self.next();
                    return;
                }
            }
            self.next();
        }
    }

    /// The source text from the start of `from` to the end of `to`, trimmed.
    fn raw_text_between(&self, from: Span, to: Span) -> String {
        from.union(to).slice_of(self.source).trim().to_string()
    }

    // --- the model -----------------------------------------------------------

    fn module(&self) -> &Module {
        &self.document.modules[self.module]
    }

    fn module_mut(&mut self) -> &mut Module {
        &mut self.document.modules[self.module]
    }

    fn find_symbol(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    fn symbol_mut(&mut self, at: usize) -> &mut Symbol {
        &mut self.document.modules[self.module].symbols[at]
    }

    /// Finds or makes a symbol in the current module. A new one is
    /// unclassified, which is what lets its kind be promoted later.
    fn declare_symbol(&mut self, name: &str, span: Span) -> usize {
        if let Some(at) = self.find_symbol(name) {
            return at;
        }
        let symbols = &mut self.document.modules[self.module].symbols;
        symbols.push(Symbol::new(name, span));
        let at = symbols.len() - 1;
        self.index.insert(name.to_string(), at);
        at
    }

    /// As [`Parser::declare_symbol`], and gives the symbol the statement's
    /// comment the first time it is seen: the comment belongs to the
    /// declaration, not to every later mention.
    fn declare(&mut self, name: &str, span: Span) -> usize {
        let at = self.declare_symbol(name, span);
        if self.symbol_mut(at).comment.is_none() {
            let comment = self.comment.clone();
            self.symbol_mut(at).comment = comment;
        }
        at
    }

    /// Delphi §6.2's promotion table: an unknown symbol may become anything,
    /// a parameter a species, compartment or stoichiometry, and nothing else
    /// moves. Anything else names both places.
    fn apply_declaration(&mut self, at: usize, kind: SymbolKind, span: Span) {
        let symbol = self.symbol_mut(at);
        if symbol.kind == kind {
            return;
        }
        if symbol.kind.can_promote_to(kind) {
            symbol.kind = kind;
            return;
        }
        let message = messages::type_conflict(
            &symbol.name,
            kind.name(),
            symbol.kind.name(),
            &symbol.decl_span.to_string(),
        );
        self.error(Code::TypeConflict, span, message);
    }

    /// Switches statements to a module, rebuilding the name index.
    fn enter_module(&mut self, at: usize) {
        self.module = at;
        self.index = self
            .module()
            .symbols
            .iter()
            .enumerate()
            .map(|(i, s)| (s.name.clone(), i))
            .collect();
    }

    // --- expressions -----------------------------------------------------------

    /// An expression, or `None` without a diagnostic if none starts here:
    /// an empty right-hand side is legal (`p = ;`).
    fn parse_expression(&mut self) -> Option<Expr> {
        let mut parser = ExprParser::new(&mut self.lexer);
        if !parser.at_expression_start() {
            return None;
        }
        parser.parse()
    }

    /// `A.B.C`, as one name.
    fn parse_symbol_ref(&mut self) -> (String, Span) {
        let first = self.next();
        let (mut name, mut span) = (first.text, first.span);
        while self.kind(0) == TokenKind::Dot && self.kind(1) == TokenKind::Identifier {
            self.next();
            let part = self.next();
            name.push('.');
            name.push_str(&part.text);
            span = span.union(part.span);
        }
        (name, span)
    }

    // --- statements --------------------------------------------------------------

    fn parse_statement(&mut self) {
        let token = self.peek(0).clone();
        self.comment = non_empty(&token.trivia);

        if self.looks_like_layout_statement() {
            self.parse_layout_statement();
            return;
        }

        match token.kind {
            TokenKind::DashDash | TokenKind::Ellipsis => {
                self.error(
                    Code::DnaNotSupported,
                    token.span,
                    messages::DNA_NOT_SUPPORTED,
                );
                self.resynchronise();
                return;
            }
            TokenKind::Identifier => {
                match token.keyword {
                    Some(Keyword::Model | Keyword::Module) => {
                        // `model` is a block header, but also the subject of
                        // a module attribute (`model notes "…"`), and of
                        // `model.sboTerm = 167`. What follows decides.
                        if self.kind(1) == TokenKind::Dot && self.is_keyword(2, Keyword::SboTerm) {
                            self.next();
                            self.next();
                            self.parse_module_sbo_term();
                        } else if self.is_any_keyword(
                            1,
                            &[Keyword::Notes, Keyword::ConversionFactor, Keyword::Is],
                        ) || self.looks_like_annotation_at(1)
                        {
                            self.parse_module_attribute();
                        } else {
                            self.parse_module();
                        }
                        return;
                    }
                    Some(Keyword::End) => {
                        // An `end` with no open module: consumed, not looped on.
                        self.next();
                        return;
                    }
                    Some(
                        Keyword::Species
                        | Keyword::Compartment
                        | Keyword::Formula
                        | Keyword::SubstanceOnly,
                    ) => {
                        self.parse_type_declaration();
                        return;
                    }
                    Some(Keyword::Const | Keyword::Var) => {
                        // `const species S1` is a declaration; `const c, k1;` a
                        // qualifier over a list (1432 and 654 corpus cases).
                        if self.is_any_keyword(1, &DECLARATION_KEYWORDS) {
                            self.parse_type_declaration();
                        } else {
                            self.parse_qualifier_declaration();
                        }
                        return;
                    }
                    Some(Keyword::Unit) => {
                        self.parse_unit_definition();
                        return;
                    }
                    Some(Keyword::Delete) => {
                        self.next();
                        self.parse_not_yet_supported("delete", "M10", false);
                        return;
                    }
                    Some(Keyword::Dna) => {
                        self.error(
                            Code::DnaNotSupported,
                            token.span,
                            messages::DNA_NOT_SUPPORTED,
                        );
                        self.resynchronise();
                        return;
                    }
                    Some(Keyword::Function) => {
                        self.parse_function();
                        return;
                    }
                    Some(Keyword::Import) => {
                        self.parse_not_yet_supported("import", "M10", false);
                        return;
                    }
                    Some(Keyword::Constraint) => {
                        self.parse_constraint();
                        return;
                    }
                    Some(Keyword::At) => {
                        self.parse_event(None, Span::default());
                        return;
                    }
                    Some(Keyword::Maximize | Keyword::Minimize) => {
                        self.parse_objective(None, Span::default());
                        return;
                    }
                    _ => {}
                }
                self.parse_symbol_statement();
                return;
            }
            TokenKind::Dollar
                if self.kind(1) == TokenKind::Identifier
                    && (self.kind(2) == TokenKind::Assign || self.is_keyword(2, Keyword::In)) =>
            {
                self.parse_dollar_declaration();
                return;
            }
            TokenKind::Dollar | TokenKind::Number | TokenKind::LParen | TokenKind::Minus => {
                // An unnamed reaction may begin with `$`, a coefficient or a
                // bracketed one; a bare flux bound with a number
                // (`-10 <= J2 <= 10`). The arrow decides.
                if self.looks_like_constraint() {
                    self.parse_bare_constraint();
                } else {
                    self.parse_reaction(None, Span::default());
                }
                return;
            }
            TokenKind::ArrowRev | TokenKind::ArrowIrrev | TokenKind::ArrowBiDir => {
                // `-> S1; k1;`: a reaction with no reactants.
                self.parse_reaction(None, Span::default());
                return;
            }
            _ => {}
        }
        let detail = format!("{} cannot start a statement", token.describe());
        self.syntax_error(token.span, &detail);
        self.resynchronise();
    }

    // --- layout statements ---------------------------------------------------------

    /// Whether a layout statement starts here (Antimony 3: `model.layout = on`,
    /// `species.color = azure`, `S1.pos = {28, 35}`, `J0.S1.b1 = {55, 31}`,
    /// `J0.--.position = {520, 483}`). `model.`, `species.`, `reaction.` and
    /// `compartment.` always start one, except `model.sboTerm`. Any other
    /// dotted name starts one when a later part is a layout word or `--`.
    /// In Antimony `A.x = 3` can also set a submodule's `x`; this library
    /// does not read submodules (M10), so such a statement is read as
    /// layout.
    fn looks_like_layout_statement(&mut self) -> bool {
        if self.kind(0) != TokenKind::Identifier || self.kind(1) != TokenKind::Dot {
            return false;
        }
        let first = self.text(0);
        if first == "model" {
            return !self.is_keyword(2, Keyword::SboTerm);
        }
        if matches!(first.as_str(), "species" | "reaction" | "compartment") {
            return true;
        }
        let mut at = 2;
        let mut layout = false;
        loop {
            match self.kind(at) {
                TokenKind::DashDash => layout = true,
                TokenKind::Identifier => layout |= is_layout_word(&self.text(at)),
                _ => return false,
            }
            match self.kind(at + 1) {
                TokenKind::Dot => at += 2,
                TokenKind::Assign => return layout,
                _ => return false,
            }
        }
    }

    /// `path = value`: the path is dotted ids (and `--`), the value a number,
    /// `{x, y}`, a braced list, words or a string.
    fn parse_layout_statement(&mut self) {
        let comment = self.comment.clone();
        let first = self.next();
        let mut span = first.span;
        let mut path = vec![first.text];
        while self.try_consume(TokenKind::Dot) {
            if !matches!(self.kind(0), TokenKind::Identifier | TokenKind::DashDash) {
                self.expected_identifier();
                self.resynchronise();
                return;
            }
            let part = self.next();
            span = span.union(part.span);
            path.push(part.text);
        }
        if !self.try_consume(TokenKind::Assign) {
            self.unexpected("\"=\"");
            self.resynchronise();
            return;
        }
        let Some(value) = self.parse_layout_value() else {
            self.unexpected("a layout value");
            self.resynchronise();
            return;
        };
        span = span.union(self.current_span());
        self.module_mut().layout.push(LayoutStatement {
            path,
            value,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    fn parse_layout_value(&mut self) -> Option<LayoutValue> {
        match self.kind(0) {
            TokenKind::Number | TokenKind::Minus => {
                self.parse_layout_number().map(LayoutValue::Number)
            }
            TokenKind::String => Some(LayoutValue::Text(self.next().text)),
            TokenKind::Identifier => {
                let mut words = vec![self.next().text];
                while self.kind(0) == TokenKind::Identifier {
                    words.push(self.next().text);
                }
                Some(LayoutValue::Word(words.join(" ")))
            }
            TokenKind::LBrace => {
                self.next();
                // Each item as written, and its value if it is a number.
                let mut items: Vec<(String, Option<f64>)> = Vec::new();
                loop {
                    self.skip_newlines();
                    match self.kind(0) {
                        TokenKind::RBrace => {
                            self.next();
                            break;
                        }
                        TokenKind::Number | TokenKind::Minus => {
                            let value = self.parse_layout_number()?;
                            items.push((number::format(value), Some(value)));
                        }
                        TokenKind::Identifier | TokenKind::String => {
                            items.push((self.next().text, None));
                        }
                        _ => return None,
                    }
                    self.skip_newlines();
                    if !self.try_consume(TokenKind::Comma) && self.kind(0) != TokenKind::RBrace {
                        return None;
                    }
                }
                Some(match items.as_slice() {
                    [(_, Some(x)), (_, Some(y))] => LayoutValue::Point(*x, *y),
                    _ => LayoutValue::List(items.into_iter().map(|(text, _)| text).collect()),
                })
            }
            _ => None,
        }
    }

    /// A number with an optional sign: `28`, `-3.5`.
    fn parse_layout_number(&mut self) -> Option<f64> {
        let negative = self.try_consume(TokenKind::Minus);
        if self.kind(0) != TokenKind::Number {
            return None;
        }
        let value = self.next().number;
        Some(if negative { -value } else { value })
    }

    /// Line breaks inside a braced list.
    fn skip_newlines(&mut self) {
        while self.kind(0) == TokenKind::Newline {
            self.next();
        }
    }

    /// A deferred construct gets one diagnostic naming its milestone, never
    /// a confusing syntax error; a block's whole body is skipped.
    fn parse_not_yet_supported(&mut self, construct: &str, milestone: &str, is_block: bool) {
        let span = self.span(0);
        self.error(
            Code::NotYetSupported,
            span,
            messages::not_yet_supported(construct, milestone),
        );
        if is_block {
            self.skip_to_block_end();
        } else {
            self.resynchronise();
        }
    }

    // --- model … end -------------------------------------------------------------

    fn parse_module(&mut self) {
        let start = self.span(0);
        self.next();
        // `model *case00001()`: the `*` marks the main module.
        let is_main = self.try_consume(TokenKind::Star);

        if self.kind(0) != TokenKind::Identifier {
            self.expected_identifier();
            self.resynchronise();
            return;
        }
        let name = self.next().text;

        // The first header names the module the document started with,
        // rather than adding a second, as long as that module is still
        // unused. Function definitions do not count as use: Antimony writes
        // them before the `model` block that calls them (BIOMD0000000064).
        let current = self.module();
        let unused = current.name.is_none()
            && current.reactions.is_empty()
            && current.symbols.len() == current.functions.len();
        if unused {
            let module = self.module_mut();
            module.name = Some(name.clone());
            module.span = start;
        } else {
            self.document
                .modules
                .push(Module::new(Some(name.clone()), start));
            let at = self.document.modules.len() - 1;
            self.enter_module(at);
        }
        let only = self.document.modules.len() == 1;
        self.module_mut().is_main = is_main || only;
        if is_main {
            self.document.main = self.module;
        }

        if self.try_consume(TokenKind::LParen) {
            while self.kind(0) == TokenKind::Identifier {
                let arg = self.next();
                if self.module().formal_args.contains(&arg.text) {
                    self.error(
                        Code::DuplicateFormalArg,
                        arg.span,
                        messages::duplicate_formal_arg(&arg.text),
                    );
                } else {
                    self.module_mut().formal_args.push(arg.text);
                }
                if !self.try_consume(TokenKind::Comma) {
                    break;
                }
            }
            if !self.try_consume(TokenKind::RParen) {
                self.unexpected("\")\"");
            }
        }

        self.skip_terminators();
        while self.kind(0) != TokenKind::Eof && !self.is_keyword(0, Keyword::End) {
            self.parse_statement();
            self.skip_terminators();
        }

        if self.is_keyword(0, Keyword::End) {
            self.next();
        } else {
            let span = self.span(0);
            self.error(
                Code::ExpectedEnd,
                span,
                messages::expected_end(&format!("model {name}"), &start.to_string()),
            );
        }
    }

    fn parse_module_attribute(&mut self) {
        self.subject_span = self.span(0);
        self.next();
        self.parse_module_attribute_for(self.module);
    }

    /// `notes`, `is`, `conversionFactor` or an annotation, on a module given
    /// as `model` or by name.
    fn parse_module_attribute_for(&mut self, target: usize) {
        if self.try_consume_keyword(Keyword::Notes) {
            if self.kind(0) == TokenKind::String {
                self.document.modules[target].notes = Some(self.next().text);
            } else {
                self.unexpected("a string");
            }
            self.expect_statement_end();
            return;
        }
        if self.try_consume_keyword(Keyword::Is) {
            if self.kind(0) == TokenKind::String {
                self.document.modules[target].display_name = Some(self.next().text);
            } else {
                self.unexpected("a display name");
            }
            self.expect_statement_end();
            return;
        }
        if self.try_consume_keyword(Keyword::ConversionFactor) {
            // Delphi §5.8: the model-level default.
            if self.kind(0) == TokenKind::Identifier {
                self.document.modules[target].conversion_factor = Some(self.parse_symbol_ref().0);
            } else {
                self.expected_identifier();
            }
            self.expect_statement_end();
            return;
        }
        let name = self.module().name.clone().unwrap_or_default();
        let subject = self.subject_span;
        if self.try_parse_annotation(&name, subject, true) {
            return;
        }
        let (span, found) = (self.span(0), self.describe(0));
        self.syntax_error(span, &format!("{found} cannot follow a module name"));
        self.resynchronise();
    }

    // --- annotations ---------------------------------------------------------------

    /// A qualifier, possibly dotted, then a string: `created "…"`,
    /// `creator1.givenName "Brett"`.
    fn looks_like_annotation_at(&mut self, at: usize) -> bool {
        self.kind(at) == TokenKind::Identifier
            && (self.kind(at + 1) == TokenKind::String
                || (self.kind(at + 1) == TokenKind::Dot
                    && self.kind(at + 2) == TokenKind::Identifier
                    && self.kind(at + 3) == TokenKind::String))
    }

    fn subject_annotations(&mut self, subject: Subject) -> (&mut Vec<CvTerm>, &mut Metadata) {
        let module = &mut self.document.modules[self.module];
        match subject {
            Subject::Module => (&mut module.cv_terms, &mut module.metadata),
            Subject::Unit(at) => {
                let unit = &mut module.units[at];
                (&mut unit.cv_terms, &mut unit.metadata)
            }
            Subject::Symbol(at) => {
                let symbol = &mut module.symbols[at];
                (&mut symbol.cv_terms, &mut symbol.metadata)
            }
        }
    }

    /// Files a recognised annotation on its subject. `false` means the
    /// qualifier is not one we know, and the caller keeps the raw text.
    fn file_annotation(
        &mut self,
        word: &str,
        field: &str,
        resources: &[String],
        subject: Subject,
        span: Span,
    ) -> bool {
        // A CV term, `S1 identity "…"`, never dotted. One qualifier may carry
        // several resources; each is its own term under that qualifier.
        if field.is_empty() {
            if let Some(qualifier) = Qualifier::from_word(word) {
                // A model qualifier means nothing on an element. Rejecting it
                // would stop the file loading; it is kept, with a warning.
                if qualifier.is_model_only() && !matches!(subject, Subject::Module) {
                    self.error(
                        Code::InformationLoss,
                        span,
                        messages::model_only_qualifier(word),
                    );
                }
                let (terms, _) = self.subject_annotations(subject);
                for resource in resources {
                    terms.push(CvTerm {
                        qualifier,
                        resource: resource.clone(),
                        spelling: word.to_string(),
                        span,
                    });
                }
                return true;
            }
        }

        // Dublin Core takes one value; if a list is given, `modified` keeps
        // them all and the others take the last.
        let last = resources.last().cloned().unwrap_or_default();
        let (base, index) = split_trailing_number(word);
        let (_, metadata) = self.subject_annotations(subject);
        match (base, field) {
            ("created", "") => metadata.created = Some(last),
            ("modified", "") => metadata.modified.extend(resources.iter().cloned()),
            ("creator", _) => {
                let creator = metadata.creator_mut(index);
                match field {
                    // `creator1 "…"` with no field is the whole name.
                    "" | "name" => creator.name = Some(last),
                    "givenName" => creator.given_name = Some(last),
                    "familyName" => creator.family_name = Some(last),
                    "organization" => creator.organization = Some(last),
                    "email" => creator.email = Some(last),
                    _ => return false,
                }
            }
            _ => return false,
        }
        true
    }

    /// `<subject> <qualifier> "<string>"`, the qualifier possibly dotted, and
    /// the string possibly a comma-separated list that may wrap across lines:
    ///
    /// ```text
    /// ATP identity "http://identifiers.org/chebi/CHEBI:15422",
    ///              "http://identifiers.org/kegg.compound/C00002"
    /// ```
    fn try_parse_annotation(
        &mut self,
        subject_name: &str,
        subject: Span,
        for_module: bool,
    ) -> bool {
        if self.kind(0) != TokenKind::Identifier {
            return false;
        }
        let word = self.text(0);
        let field = if self.kind(1) == TokenKind::String {
            String::new()
        } else if self.kind(1) == TokenKind::Dot
            && self.kind(2) == TokenKind::Identifier
            && self.kind(3) == TokenKind::String
        {
            self.text(2)
        } else {
            return false;
        };

        // The span starts at the SUBJECT, so a raw section written back says
        // what it annotates.
        let start = if subject.is_unknown() {
            self.span(0)
        } else {
            subject
        };
        while self.kind(0) != TokenKind::String {
            self.next();
        }
        let first = self.next();
        let mut resources = vec![first.text];
        let mut stop = first.span;

        // A comma holds the statement open across a line break, but only once
        // a string is known to follow, so a malformed line still reports at
        // the comma.
        while self.kind(0) == TokenKind::Comma {
            let mut ahead = 1;
            while self.kind(ahead) == TokenKind::Newline {
                ahead += 1;
            }
            if self.kind(ahead) != TokenKind::String {
                break;
            }
            while self.kind(0) != TokenKind::String {
                self.next();
            }
            let resource = self.next();
            resources.push(resource.text);
            stop = resource.span;
        }

        // `model …` is the module; so is a subject naming the module itself
        // (`MODEL_01186 identity "…"`). A unit is not a symbol, so it is
        // checked before declaring one: that would invent a parameter.
        let is_module = for_module
            || (!subject_name.is_empty() && self.module().name.as_deref() == Some(subject_name));
        let target = if is_module {
            Subject::Module
        } else if let Some(at) = self
            .module()
            .units
            .iter()
            .position(|u| u.name == subject_name)
        {
            Subject::Unit(at)
        } else {
            Subject::Symbol(self.declare(subject_name, subject))
        };

        // A vocabulary we do not know is kept as raw text, so a file that
        // loaded before never stops loading.
        let span = start.union(stop);
        if !self.file_annotation(&word, &field, &resources, target, span) {
            let text = self.raw_text_between(start, stop);
            self.module_mut().unparsed_sections.push(RawSection {
                kind: RawKind::Annotation,
                text,
                span,
            });
        }
        self.expect_statement_end();
        true
    }

    // --- declarations ----------------------------------------------------------------

    /// `[const|var|substanceOnly]* (species|compartment|formula) item, …` where
    /// an item is `[$]name [in compartment] [= value]`.
    fn parse_type_declaration(&mut self) {
        let mut constness = Constness::Default;
        let mut substance_only = false;
        while self.is_any_keyword(0, &[Keyword::Const, Keyword::Var, Keyword::SubstanceOnly]) {
            match self.next().keyword {
                Some(Keyword::Const) => constness = Constness::Const,
                Some(Keyword::Var) => constness = Constness::Var,
                _ => substance_only = true,
            }
        }

        let kind = if self.is_keyword(0, Keyword::Species) {
            SymbolKind::Species
        } else if self.is_keyword(0, Keyword::Compartment) {
            SymbolKind::Compartment
        } else {
            SymbolKind::Parameter
        };
        self.next();

        loop {
            // The `$` is per item: `species G in Cell, $T in Cell;`.
            let is_boundary = self.try_consume(TokenKind::Dollar);
            if self.kind(0) != TokenKind::Identifier {
                self.expected_identifier();
                self.resynchronise();
                return;
            }
            let token = self.next();
            let at = self.declare(&token.text, token.span);
            self.apply_declaration(at, kind, token.span);
            let symbol = self.symbol_mut(at);
            if is_boundary {
                symbol.is_boundary = true;
            }
            if substance_only {
                symbol.is_substance_only = true;
            }
            if constness != Constness::Default {
                symbol.constness = constness;
                // A const species is a boundary species.
                if constness == Constness::Const && kind == SymbolKind::Species {
                    symbol.is_boundary = true;
                }
            }

            if self.try_consume_keyword(Keyword::In) {
                if self.kind(0) == TokenKind::Identifier {
                    let compartment = self.parse_symbol_ref().0;
                    self.symbol_mut(at).compartment = Some(compartment);
                } else {
                    self.expected_identifier();
                }
            }
            if self.try_consume(TokenKind::Assign) {
                self.set_value(at);
            }
            if !self.try_consume(TokenKind::Comma) {
                break;
            }
        }
        self.expect_statement_end();
    }

    /// The right-hand side of `=`. An empty one is recorded and leaves any
    /// value already declared alone (Delphi §9.1.4 #3).
    fn set_value(&mut self, at: usize) {
        match self.parse_expression() {
            Some(value) => self.symbol_mut(at).initial_value = Some(value),
            None => self.symbol_mut(at).has_empty_value = true,
        }
    }

    /// `const c, k1;` and `var x, y;`, where an item may also say where it
    /// is and what it is worth, as libantimony allows: `const k1 = 0.5`,
    /// `var S1 in C = 3` (Delphi rejected these; SPEC.md R36).
    fn parse_qualifier_declaration(&mut self) {
        let constness = if self.is_keyword(0, Keyword::Const) {
            Constness::Const
        } else {
            Constness::Var
        };
        self.next();
        self.parse_qualified_items(constness);
    }

    /// `$x = 4` and `$x in C = 4`: libantimony reads the `$` as `const`,
    /// so this is `const $x = 4` (SPEC.md R36).
    fn parse_dollar_declaration(&mut self) {
        self.parse_qualified_items(Constness::Const);
    }

    /// The items of a `const` or `var` declaration: `[$]name [in c] [= value]`.
    fn parse_qualified_items(&mut self, constness: Constness) {
        loop {
            // A `$` on a species is redundant with `const`. On anything else
            // it is kept, so a later use as a species makes a boundary one
            // and the writer can write it back.
            let dollar = self.try_consume(TokenKind::Dollar);
            if self.kind(0) != TokenKind::Identifier {
                self.expected_identifier();
                self.resynchronise();
                return;
            }
            let token = self.next();
            let at = self.declare_symbol(&token.text, token.span);
            let symbol = self.symbol_mut(at);
            symbol.constness = constness;
            if dollar || (constness == Constness::Const && symbol.kind == SymbolKind::Species) {
                symbol.is_boundary = true;
            }
            if self.try_consume_keyword(Keyword::In) {
                if self.kind(0) == TokenKind::Identifier {
                    let compartment = self.parse_symbol_ref().0;
                    self.symbol_mut(at).compartment = Some(compartment);
                } else {
                    self.expected_identifier();
                }
            }
            if self.try_consume(TokenKind::Assign) {
                self.set_value(at);
            }
            if !self.try_consume(TokenKind::Comma) {
                break;
            }
        }
        self.expect_statement_end();
    }

    /// `unit name = definition`, the definition kept as text.
    fn parse_unit_definition(&mut self) {
        self.next();
        if self.kind(0) != TokenKind::Identifier {
            self.expected_identifier();
            self.resynchronise();
            return;
        }
        let name = self.next();
        if !self.try_consume(TokenKind::Assign) {
            self.unexpected("\"=\"");
            self.resynchronise();
            return;
        }
        let from = self.span(0);
        let mut to = from;
        while !self.at_statement_end() {
            to = self.span(0);
            self.next();
        }
        let definition = self.raw_text_between(from, to);
        let comment = self.comment.clone();
        self.module_mut().units.push(UnitDef {
            name: name.text,
            definition,
            display_name: None,
            cv_terms: Vec::new(),
            metadata: Metadata::default(),
            comment,
            span: name.span.union(to),
        });
        self.expect_statement_end();
    }

    // --- events ----------------------------------------------------------------------

    /// `[name:] at [delay after] trigger {, modifier = value} : [target = value, …]`.
    fn parse_event(&mut self, name: Option<String>, name_span: Span) {
        let start = if name_span.is_unknown() {
            self.span(0)
        } else {
            name_span
        };
        self.next();

        let Some(first) = self.parse_expression() else {
            let span = self.span(0);
            self.syntax_error(span, "an event needs a trigger");
            self.resynchronise_line();
            return;
        };
        // The delay may come first, joined by `after` (328 corpus cases).
        let (trigger, mut delay) = if self.try_consume_keyword(Keyword::After) {
            let Some(trigger) = self.parse_expression() else {
                let span = self.span(0);
                self.syntax_error(span, "\"after\" must be followed by a trigger");
                self.resynchronise_line();
                return;
            };
            (trigger, Some(first))
        } else {
            (first, None)
        };
        let (mut priority, mut persistent, mut t0, mut from_trigger) = (None, None, None, None);

        while self.try_consume(TokenKind::Comma) {
            let modifiers = [
                Keyword::Priority,
                Keyword::Persistent,
                Keyword::T0,
                Keyword::FromTrigger,
                Keyword::Delay,
            ];
            if !self.is_any_keyword(0, &modifiers) {
                let (span, found) = (self.span(0), self.describe(0));
                self.syntax_error(span, &format!("{found} is not an event modifier"));
                self.resynchronise_line();
                return;
            }
            let modifier = self.next().keyword;
            if !self.try_consume(TokenKind::Assign) {
                self.unexpected("\"=\"");
                self.resynchronise_line();
                return;
            }
            let Some(value) = self.parse_expression() else {
                let span = self.span(0);
                self.syntax_error(span, "an event modifier needs a value");
                self.resynchronise_line();
                return;
            };
            match modifier {
                Some(Keyword::Priority) => priority = Some(value),
                Some(Keyword::Persistent) => persistent = Some(value),
                Some(Keyword::T0) => t0 = Some(value),
                Some(Keyword::FromTrigger) => from_trigger = Some(value),
                _ => delay = Some(value),
            }
        }

        if !self.try_consume(TokenKind::Colon) {
            self.unexpected("\":\" before the event assignments");
            self.resynchronise_line();
            return;
        }

        // The list may be empty: `E0: at time > 5.5: ;` is a challenge case.
        let mut assignments = Vec::new();
        while self.kind(0) == TokenKind::Identifier {
            let (target, span) = self.parse_symbol_ref();
            if !self.try_consume(TokenKind::Assign) {
                self.unexpected("\"=\"");
                self.resynchronise_line();
                return;
            }
            let expr = self.parse_expression();
            assignments.push(EventAssignment { target, expr, span });
            if !self.try_consume(TokenKind::Comma) {
                break;
            }
        }

        let span = start.union(self.current_span());
        if let Some(name) = &name {
            let at = self.declare_symbol(name, name_span);
            self.apply_declaration(at, SymbolKind::Event, name_span);
        }
        let comment = self.comment.clone();
        self.module_mut().events.push(Event {
            name,
            trigger,
            delay,
            priority,
            persistent,
            initial_value: t0,
            from_trigger,
            assignments,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    // --- functions ---------------------------------------------------------------------

    /// `function name(args) body end`. The corpus has none: libantimony
    /// inlines SBML functions, so only fixtures cover this.
    fn parse_function(&mut self) {
        let start = self.span(0);
        self.next();
        if self.kind(0) != TokenKind::Identifier {
            self.expected_identifier();
            self.skip_to_block_end();
            return;
        }
        let name = self.next();
        let mut args: Vec<String> = Vec::new();
        if self.try_consume(TokenKind::LParen) {
            while self.kind(0) == TokenKind::Identifier {
                let arg = self.next();
                if args.contains(&arg.text) {
                    self.error(
                        Code::DuplicateFormalArg,
                        arg.span,
                        messages::duplicate_formal_arg(&arg.text),
                    );
                } else {
                    args.push(arg.text);
                }
                if !self.try_consume(TokenKind::Comma) {
                    break;
                }
            }
            if !self.try_consume(TokenKind::RParen) {
                self.unexpected("\")\"");
            }
        }

        self.skip_terminators();
        let body = self.parse_expression();
        if body.is_none() {
            let span = self.span(0);
            self.syntax_error(span, "a function needs a body");
        }
        self.skip_terminators();

        if self.is_keyword(0, Keyword::End) {
            self.next();
        } else {
            let span = self.span(0);
            self.error(
                Code::ExpectedEnd,
                span,
                messages::expected_end(&format!("function {}", name.text), &start.to_string()),
            );
        }

        let span = start.union(self.current_span());
        let at = self.declare_symbol(&name.text, name.span);
        self.apply_declaration(at, SymbolKind::Function, name.span);
        let comment = self.comment.clone();
        self.module_mut().functions.push(Function {
            name: name.text,
            args,
            body,
            comment,
            span,
        });
    }

    // --- rules, constraints and objectives -----------------------------------------------

    /// `name: 0 = expr`. The name is a label, not a symbol.
    fn parse_algebraic_rule(&mut self, name: String, name_span: Span) {
        self.next(); // the 0
        self.next(); // the =
        // `_alg0: 0 = ;` is real corpus text (case 01244): kept as empty.
        let expr = self.parse_expression();
        let span = name_span.union(self.current_span());
        let comment = self.comment.clone();
        self.module_mut().algebraic_rules.push(AlgebraicRule {
            name,
            expr,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    /// The body of a constraint: an expression, or a chained comparison,
    /// which is a range here and only here (Delphi decision 18). The
    /// expression parser reads `a <= b <= c` as `(a <= b) <= c`; this takes
    /// that apart.
    fn parse_constraint_body(&mut self) -> Option<ConstraintKind> {
        let Some(first) = self.parse_expression() else {
            let span = self.span(0);
            self.syntax_error(span, "a constraint needs an expression");
            return None;
        };
        if let ExprKind::Binary {
            op: upper_op,
            left,
            right: upper,
        } = &first.kind
        {
            if let ExprKind::Binary {
                op: lower_op,
                left: lower,
                right: subject,
            } = &left.kind
            {
                if is_relational(*upper_op) && is_relational(*lower_op) {
                    return Some(ConstraintKind::Range {
                        lower: (**lower).clone(),
                        subject: (**subject).clone(),
                        upper: (**upper).clone(),
                        lower_inclusive: is_inclusive(*lower_op),
                        upper_inclusive: is_inclusive(*upper_op),
                    });
                }
            }
        }
        Some(ConstraintKind::Expr(first))
    }

    /// `constraint [name:] body`.
    fn parse_constraint(&mut self) {
        let start = self.span(0);
        self.next();
        let mut name = None;
        if self.kind(0) == TokenKind::Identifier && self.kind(1) == TokenKind::Colon {
            name = Some(self.next().text);
            self.next();
        }
        let Some(kind) = self.parse_constraint_body() else {
            self.resynchronise_line();
            return;
        };
        self.push_constraint(name, kind, start);
    }

    fn push_constraint(&mut self, name: Option<String>, kind: ConstraintKind, start: Span) {
        let span = start.union(self.current_span());
        let comment = self.comment.clone();
        self.module_mut().constraints.push(Constraint {
            name,
            kind,
            message: None,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    /// A flux bound written without `constraint`: `0 <= J0`, `J1 <= 1000`.
    fn parse_bare_constraint(&mut self) {
        let start = self.span(0);
        let Some(kind) = self.parse_constraint_body() else {
            self.resynchronise_line();
            return;
        };
        self.push_constraint(None, kind, start);
    }

    /// `[name:] maximize|minimize expr`. The name is a label, not a symbol:
    /// libantimony reports no `OBJF` symbol for case 01612.
    fn parse_objective(&mut self, name: Option<String>, name_span: Span) {
        let start = if name_span.is_unknown() {
            self.span(0)
        } else {
            name_span
        };
        let maximise = self.is_keyword(0, Keyword::Maximize);
        self.next();
        let Some(expr) = self.parse_expression() else {
            let span = self.span(0);
            self.syntax_error(span, "an objective needs an expression");
            self.resynchronise();
            return;
        };
        let span = start.union(self.current_span());
        let comment = self.comment.clone();
        self.module_mut().objectives.push(Objective {
            name,
            maximise,
            expr,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    /// `= <number>` after `sboTerm`, or `None` having reported the problem.
    fn parse_sbo_value(&mut self) -> Option<u32> {
        self.next(); // `sboTerm`
        if !self.try_consume(TokenKind::Assign) {
            self.unexpected("\"=\"");
            self.resynchronise();
            return None;
        }
        if self.kind(0) != TokenKind::Number {
            self.unexpected("an SBO term number");
            self.resynchronise();
            return None;
        }
        Some(self.next().number.trunc() as u32)
    }

    /// `fb_0.sboTerm = 625` (98 corpus cases) or `A sboTerm = 236`. Without
    /// this a symbol named "fb_0.sboTerm" would appear.
    fn parse_sbo_term(&mut self, target: &str, span: Span) {
        let Some(value) = self.parse_sbo_value() else {
            return;
        };
        if let Some(module) = self.document.module_index(target) {
            self.document.modules[module].sbo_term = Some(value);
        } else {
            let at = self.declare_symbol(target, span);
            self.symbol_mut(at).sbo_term = Some(value);
        }
        self.expect_statement_end();
    }

    /// `model.sboTerm = 167`, with `model.` consumed. Set on the module: a
    /// symbol called "model" would be a phantom in every such file.
    fn parse_module_sbo_term(&mut self) {
        let Some(value) = self.parse_sbo_value() else {
            return;
        };
        self.module_mut().sbo_term = Some(value);
        self.expect_statement_end();
    }

    // --- reactions and interactions --------------------------------------------------------

    /// An identifier, an interaction arrow, an identifier. `-(` is two tokens
    /// (a minus before a bracket is far more common), so position decides.
    fn looks_like_interaction(&mut self) -> bool {
        if self.kind(0) != TokenKind::Identifier {
            return false;
        }
        if matches!(self.kind(1), TokenKind::Inhibits | TokenKind::Activates) {
            return self.kind(2) == TokenKind::Identifier;
        }
        self.kind(1) == TokenKind::Minus
            && self.kind(2) == TokenKind::LParen
            && self.kind(3) == TokenKind::Identifier
    }

    /// `[label:] species (-| | -o | -() reaction`. The species is declared
    /// by the mention; the reaction is not, since an interaction naming a
    /// reaction that does not exist is for the semantic pass to report.
    fn parse_interaction(&mut self, name: Option<String>, name_span: Span) {
        let start = self.span(0);
        let species = self.next().text;
        let kind = match self.kind(0) {
            TokenKind::Inhibits => InteractionKind::Inhibits,
            TokenKind::Activates => InteractionKind::Activates,
            _ => {
                self.next(); // the `-` of `-(`
                InteractionKind::Unknown
            }
        };
        self.next();

        if self.kind(0) != TokenKind::Identifier {
            self.expected_identifier();
            self.resynchronise();
            return;
        }
        let reaction = self.next().text;

        let at = self.declare(&species, start);
        self.apply_declaration(at, SymbolKind::Species, start);
        // A label is a symbol; an unlabelled interaction declares nothing.
        if let Some(name) = &name {
            let at = self.declare(name, name_span);
            self.apply_declaration(at, SymbolKind::Interaction, name_span);
        }
        let span = start.union(self.current_span());
        let comment = self.comment.clone();
        self.module_mut().interactions.push(Interaction {
            name,
            species,
            reaction,
            kind,
            comment,
            span,
        });
        self.expect_statement_end();
    }

    /// Whether an arrow comes before the end of the statement, outside
    /// brackets and before any `=` or `:=`.
    fn looks_like_reaction(&mut self) -> bool {
        let mut depth = 0i32;
        for i in 0..=LOOKAHEAD {
            let kind = self.kind(i);
            match kind {
                TokenKind::LParen => depth += 1,
                TokenKind::RParen => depth -= 1,
                TokenKind::Semicolon
                | TokenKind::Newline
                | TokenKind::Eof
                | TokenKind::Assign
                | TokenKind::ColonAssign => return false,
                _ => {}
            }
            if depth == 0 && arrow_kind(kind).is_some() {
                return true;
            }
        }
        false
    }

    /// A statement that is a comparison and nothing else is a flux bound.
    /// An arrow makes it a reaction and an `=` an assignment.
    fn looks_like_constraint(&mut self) -> bool {
        let mut depth = 0i32;
        for i in 0..=LOOKAHEAD {
            let kind = self.kind(i);
            match kind {
                TokenKind::LParen => depth += 1,
                TokenKind::RParen => depth -= 1,
                TokenKind::Semicolon | TokenKind::Newline | TokenKind::Eof => return false,
                _ => {}
            }
            if depth == 0 {
                if matches!(
                    kind,
                    TokenKind::Assign
                        | TokenKind::ColonAssign
                        | TokenKind::Colon
                        | TokenKind::Prime
                ) || arrow_kind(kind).is_some()
                {
                    return false;
                }
                if matches!(
                    kind,
                    TokenKind::Less
                        | TokenKind::LessEqual
                        | TokenKind::Greater
                        | TokenKind::GreaterEqual
                        | TokenKind::Equal
                ) {
                    return true;
                }
            }
        }
        false
    }

    /// Whether a stoichiometry comes next: a number, name or bracketed
    /// expression followed by a species name.
    fn stoichiometry_follows(&mut self) -> bool {
        let species_next =
            |kind: TokenKind| matches!(kind, TokenKind::Identifier | TokenKind::Dollar);
        match self.kind(0) {
            TokenKind::Number | TokenKind::Identifier => return species_next(self.kind(1)),
            TokenKind::LParen => {}
            _ => return false,
        }
        let mut depth = 0i32;
        for i in 0..=LOOKAHEAD {
            match self.kind(i) {
                TokenKind::LParen => depth += 1,
                TokenKind::RParen => {
                    depth -= 1;
                    if depth == 0 {
                        return species_next(self.kind(i + 1));
                    }
                }
                TokenKind::Semicolon | TokenKind::Newline | TokenKind::Eof => return false,
                _ => {}
            }
        }
        false
    }

    /// One side of a reaction, into `into`. Either side may be empty. `false`
    /// means a problem was reported.
    fn parse_participants(&mut self, into: &mut Vec<Participant>) -> bool {
        if self.at_statement_end() || arrow_kind(self.kind(0)).is_some() {
            return true;
        }
        loop {
            let start = self.span(0);

            // A negative stoichiometry is rejected by name (Delphi §5.7,
            // decision 14): it has no chemical meaning, and libantimony writes
            // it and then cannot read it back.
            if self.kind(0) == TokenKind::Minus {
                let mut species = String::new();
                for i in 1..=3 {
                    if self.kind(i) == TokenKind::Identifier {
                        species = self.text(i);
                        break;
                    }
                }
                self.error(
                    Code::NegativeStoichiometry,
                    start,
                    messages::negative_stoichiometry(&species),
                );
                return false;
            }

            // Units are off here, or `2 S2` would be 2 in units of S2.
            let mut stoichiometry = None;
            if self.stoichiometry_follows() {
                let mut parser = ExprParser::new(&mut self.lexer);
                parser.parse_unit_suffixes = false;
                let Some(value) = parser.parse() else {
                    return false;
                };
                // A name in the coefficient's place is a named
                // stoichiometry, not a parameter (Delphi decision 17).
                if let ExprKind::Symbol(name) = &value.kind {
                    let at = self.declare_symbol(name, value.span);
                    self.apply_declaration(at, SymbolKind::Stoichiometry, value.span);
                }
                stoichiometry = Some(value);
            }

            let is_boundary = self.try_consume(TokenKind::Dollar);
            if self.kind(0) != TokenKind::Identifier {
                self.expected_identifier();
                return false;
            }
            let token = self.next();
            into.push(Participant {
                species: token.text.clone(),
                stoichiometry,
                is_boundary,
                span: start.union(token.span),
            });
            let at = self.declare_symbol(&token.text, token.span);
            self.apply_declaration(at, SymbolKind::Species, token.span);
            if is_boundary {
                self.symbol_mut(at).is_boundary = true;
            }
            if !self.try_consume(TokenKind::Plus) {
                return true;
            }
        }
    }

    /// `_J` and the reaction's 1-based position among all reactions, or the
    /// next free number if an author already used that id.
    fn generate_reaction_id(&self) -> String {
        let module = self.module();
        let mut next = module.reactions.len() + 1;
        loop {
            let id = format!("_J{next}");
            let taken = self.find_symbol(&id).is_some()
                || module
                    .reactions
                    .iter()
                    .any(|r| r.name.eq_ignore_ascii_case(&id));
            if !taken {
                return id;
            }
            next += 1;
        }
    }

    /// `[name:] participants arrow participants [; [law]]`.
    fn parse_reaction(&mut self, name: Option<String>, name_span: Span) {
        let start = if name_span.is_unknown() {
            self.span(0)
        } else {
            name_span
        };
        let mut reaction = Reaction::new(name.clone().unwrap_or_default(), start);
        reaction.comment = self.comment.clone();

        if !self.parse_participants(&mut reaction.reactants) {
            self.resynchronise_line();
            return;
        }
        let Some(kind) = arrow_kind(self.kind(0)) else {
            self.unexpected("a reaction arrow");
            self.resynchronise_line();
            return;
        };
        self.next();
        reaction.kind = kind;
        if !self.parse_participants(&mut reaction.products) {
            self.resynchronise_line();
            return;
        }

        if reaction.reactants.is_empty() && reaction.products.is_empty() {
            // `J0: -> ; 2;` is text libantimony writes and cannot read back
            // (a challenge case). Accepted, so only a hint.
            self.lexer.diagnostics_mut().add_with(
                Code::EmptyReaction,
                Severity::Hint,
                reaction.span,
                messages::empty_reaction(name.as_deref().unwrap_or_default()),
            );
        }

        // A `;` may start the law or end the statement:
        //   reaction1: S1 => S2; k*S1;   a law
        //   R16: J => M; ;               a law written and empty
        //   R17: J => M;                 no law
        if self.kind(0) == TokenKind::Semicolon
            && (crate::expr::can_start_expression(self.kind(1))
                || self.kind(1) == TokenKind::Semicolon)
        {
            self.next();
            match self.parse_expression() {
                Some(law) => reaction.kinetic_law = Some(law),
                None => reaction.has_empty_kinetic_law = true,
            }
        }

        reaction.span = start.union(self.current_span());
        match &name {
            Some(name) => {
                let at = self.declare_symbol(name, name_span);
                self.apply_declaration(at, SymbolKind::Reaction, name_span);
            }
            None => {
                reaction.name = self.generate_reaction_id();
                reaction.name_generated = true;
            }
        }
        self.module_mut().reactions.push(reaction);
        self.expect_statement_end();
    }

    // --- statements that start with a name ---------------------------------------------------

    fn parse_symbol_statement(&mut self) {
        if self.kind(0) == TokenKind::Identifier && self.kind(1) == TokenKind::Colon {
            self.parse_labelled_statement();
            return;
        }

        // Interactions before reactions: `-|` and `-o` are reaction arrows
        // too, and `G1 -| J1` read as a reaction would make J1 a species.
        if self.looks_like_interaction() {
            self.parse_interaction(None, Span::default());
            return;
        }
        if self.looks_like_reaction() {
            self.parse_reaction(None, Span::default());
            return;
        }
        if self.looks_like_constraint() {
            self.parse_bare_constraint();
            return;
        }

        // `fb_0.sboTerm = 625`, caught before the dotted name is read whole.
        if self.kind(0) == TokenKind::Identifier
            && self.kind(1) == TokenKind::Dot
            && self.is_keyword(2, Keyword::SboTerm)
        {
            let token = self.next();
            self.next();
            self.parse_sbo_term(&token.text, token.span);
            return;
        }

        let (name, name_span) = self.parse_symbol_ref();

        if self.is_keyword(0, Keyword::SboTerm) {
            self.parse_sbo_term(&name, name_span);
            return;
        }

        // `case00001 is "case00001"`: the subject is a module, and declaring
        // a symbol for it would put a phantom parameter in every model.
        if let Some(module) = self.document.module_index(&name) {
            if self.is_any_keyword(0, &[Keyword::Is, Keyword::Notes, Keyword::ConversionFactor]) {
                self.parse_module_attribute_for(module);
                return;
            }
            // `m identity "…"`: an annotation on the model, caught before a
            // symbol named after it is declared.
            if self.looks_like_annotation_at(0)
                && self.try_parse_annotation(&name, name_span, false)
            {
                return;
            }
        }

        // `time_unit is "time"`: a unit is not a symbol either.
        if let Some(unit) = self.module().units.iter().position(|u| u.name == name) {
            if self.is_keyword(0, Keyword::Is) {
                self.next();
                if self.kind(0) == TokenKind::String {
                    self.module_mut().units[unit].display_name = Some(self.next().text);
                } else {
                    self.unexpected("a display name");
                }
                self.expect_statement_end();
                return;
            }
            // An annotation on a unit, checked here because declaring a
            // symbol comes next.
            if self.looks_like_annotation_at(0)
                && self.try_parse_annotation(&name, name_span, false)
            {
                return;
            }
        }

        let at = self.declare(&name, name_span);

        // `S1' = expr`, a rate rule.
        if self.try_consume(TokenKind::Prime) {
            if !self.try_consume(TokenKind::Assign) {
                self.unexpected("\"=\"");
                self.resynchronise();
                return;
            }
            let expr = self.parse_expression();
            self.symbol_mut(at).rule = Some(Rule {
                kind: RuleKind::Rate,
                expr,
            });
            self.expect_statement_end();
            return;
        }

        // `Vtot := expr`, an assignment rule.
        if self.try_consume(TokenKind::ColonAssign) {
            let expr = self.parse_expression();
            self.symbol_mut(at).rule = Some(Rule {
                kind: RuleKind::Assignment,
                expr,
            });
            self.expect_statement_end();
            return;
        }

        // `S1 in cell` and `S1 in cell = 5`.
        if self.try_consume_keyword(Keyword::In) {
            if self.kind(0) == TokenKind::Identifier {
                let compartment = self.parse_symbol_ref().0;
                self.symbol_mut(at).compartment = Some(compartment);
            } else {
                self.expected_identifier();
            }
            if self.try_consume(TokenKind::Assign) {
                self.set_value(at);
            }
            self.expect_statement_end();
            return;
        }

        // `S1 has substance_per_volume` (933 corpus cases).
        if self.try_consume_keyword(Keyword::Has) {
            if self.kind(0) == TokenKind::Identifier {
                let units = self.parse_symbol_ref().0;
                self.symbol_mut(at).units = Some(units);
            } else {
                self.expected_identifier();
            }
            self.expect_statement_end();
            return;
        }

        // `S1 is "Glucose"`, or `A.S1 is B.S2`, a module mapping.
        if self.try_consume_keyword(Keyword::Is) {
            match self.kind(0) {
                TokenKind::String => {
                    let display = self.next().text;
                    self.symbol_mut(at).display_name = Some(display);
                }
                TokenKind::Identifier => {
                    let (target, span) = self.parse_symbol_ref();
                    self.error(
                        Code::NotYetSupported,
                        span,
                        messages::not_yet_supported(
                            &format!("The \"is\" mapping onto {target}"),
                            "M10",
                        ),
                    );
                }
                _ => self.unexpected("a display name or a symbol"),
            }
            self.expect_statement_end();
            return;
        }

        // `S1 notes "…"`.
        if self.try_consume_keyword(Keyword::Notes) {
            if self.kind(0) == TokenKind::String {
                let notes = self.next().text;
                self.symbol_mut(at).notes = Some(notes);
            } else {
                self.unexpected("a string");
            }
            self.expect_statement_end();
            return;
        }

        // `S1 conversionFactor cf` (Delphi §5.8).
        if self.try_consume_keyword(Keyword::ConversionFactor) {
            if self.kind(0) == TokenKind::Identifier {
                let factor = self.parse_symbol_ref().0;
                self.symbol_mut(at).conversion_factor = Some(factor);
            } else {
                self.expected_identifier();
            }
            self.expect_statement_end();
            return;
        }

        // `S1 = 10`: a value or initial assignment; the semantic pass decides.
        if self.try_consume(TokenKind::Assign) {
            self.set_value(at);
            self.expect_statement_end();
            return;
        }

        // A name on its own declares it by mention.
        if self.at_statement_end() {
            self.expect_statement_end();
            return;
        }

        // `S1 identity "http://…"`.
        if self.try_parse_annotation(&name, name_span, false) {
            return;
        }

        let (span, found) = (self.span(0), self.describe(0));
        self.syntax_error(span, &format!("{found} cannot follow \"{name}\""));
        self.resynchronise();
    }

    /// A statement starting `name:`: an interaction, reaction, event,
    /// algebraic rule, submodule, objective, or DNA strand.
    fn parse_labelled_statement(&mut self) {
        let label = self.next();
        self.next(); // the `:`
        let (name, span) = (label.text, label.span);

        if self.looks_like_interaction() {
            self.parse_interaction(Some(name), span);
        } else if self.looks_like_reaction() {
            self.parse_reaction(Some(name), span);
        } else if self.is_keyword(0, Keyword::At) {
            self.parse_event(Some(name), span);
        } else if self.kind(0) == TokenKind::Number && self.kind(1) == TokenKind::Assign {
            // `_alg0: 0 = -1*k1 + S1 + S2;` (124 corpus cases).
            self.parse_algebraic_rule(name, span);
        } else if self.kind(0) == TokenKind::Identifier && self.kind(1) == TokenKind::LParen {
            self.parse_not_yet_supported("Submodule instantiation", "M10", false);
        } else if self.is_any_keyword(0, &[Keyword::Maximize, Keyword::Minimize]) {
            self.parse_objective(Some(name), span);
        } else if matches!(self.kind(0), TokenKind::DashDash | TokenKind::Ellipsis) {
            // `dna1: --P1--G1--`, named rather than a stray operator.
            let span = self.span(0);
            self.error(Code::DnaNotSupported, span, messages::DNA_NOT_SUPPORTED);
            self.resynchronise();
        } else {
            let (at, found) = (self.span(0), self.describe(0));
            self.syntax_error(at, &format!("{found} cannot follow \"{name}:\""));
            self.resynchronise();
        }
    }
}

/// A part of a dotted name that only a layout statement uses: a position, a
/// size, a line's ends and control points, a look, or `arcN`/`segN`. Read
/// without regard to case: the reference writes `linecolor`, libAntimony's
/// writer `lineColor`.
fn is_layout_word(part: &str) -> bool {
    let lower = part.to_ascii_lowercase();
    let numbered = |prefix: &str| {
        lower
            .strip_prefix(prefix)
            .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
    };
    numbered("arc")
        || numbered("seg")
        || matches!(
            lower.as_str(),
            "pos"
                | "position"
                | "x"
                | "y"
                | "size"
                | "width"
                | "height"
                | "species_end"
                | "rxn_end"
                | "b1"
                | "b2"
                | "color"
                | "fillcolor"
                | "linecolor"
                | "strokecolor"
                | "linewidth"
                | "strokewidth"
                | "fontcolor"
                | "font"
                | "fontsize"
                | "fontstyle"
                | "fontweight"
                | "shape"
        )
}

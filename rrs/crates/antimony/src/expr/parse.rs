//! The expression parser (Delphi `uAntExprParser`): precedence climbing over
//! the shared lexer, following Delphi §5.3's table exactly. Two of its levels
//! are unusual and deliberate:
//!
//! - `!` / `not` sits BETWEEN `&&` and `==`, so `!a == b` is `!(a == b)`;
//! - unary `+` / `-` bind TIGHTER than `^`, so `-2^2` is 4, not -4. The
//!   corpus has `^-` 42 times, which depends on it.
//!
//! The parser never resolves names, which is what lets `S1 = k*2` come
//! before `k = 3`.
//!
//! **Nesting is limited** (SPEC.md §5.3). The Delphi parser had no limit;
//! here a pathological input would overflow the stack and abort, so it is
//! reported as a syntax error instead. The Delphi parser had one function
//! per precedence level, a dozen stack frames for every bracket; this one
//! climbs the precedence table in a loop, which takes a few.

use super::{BinaryOp, Constant, Expr, ExprKind, UnaryOp, precedence};
use crate::diagnostics::{Code, Diagnostics};
use crate::lexer::{Keyword, Lexer, TokenKind};
use crate::messages;
use crate::span::Span;

// The limits are about four times what real models use, and small enough
// that every recursive walk of the deepest tree allowed fits a 1 MB stack
// (a Windows main thread, WebAssembly) even unoptimised. Measured in A2: the
// deepest tree in 689 BioModels is 192 (BIOMD0000000574) and the deepest
// brackets 15 (BIOMD0000000589); in the 1788-case corpus, 100 and 4. An
// unoptimised build overflowed 1 MB at about 1000-2000 levels of tree and
// 100-200 of brackets; an optimised one at 4000 and 400-800.

/// The deepest an expression tree may be. A sum of terms is as deep as it
/// has terms, so this is also the longest such sum.
pub const MAX_DEPTH: usize = 800;

/// How deeply brackets, calls, unary operators, `?:` and `^` may nest.
pub const MAX_NESTING: usize = 64;

/// Whether a token of this kind can begin an expression. A string cannot:
/// SBML has no string values, so a string where an expression could be is
/// some other statement.
pub fn can_start_expression(kind: TokenKind) -> bool {
    matches!(
        kind,
        TokenKind::Number
            | TokenKind::Identifier
            | TokenKind::LParen
            | TokenKind::Plus
            | TokenKind::Minus
            | TokenKind::Not
    )
}

/// A node and its depth, so the depth limit costs nothing to check.
struct Node {
    expr: Expr,
    depth: usize,
}

/// Reads expressions from a lexer it borrows. It stops at the end of the
/// expression, so the statement parser carries on from there. Problems go to
/// the lexer's diagnostics list.
pub struct ExprParser<'l, 'a> {
    lexer: &'l mut Lexer<'a>,
    /// Whether a number may carry a unit (`40 foo`). The statement parser
    /// turns this off in a reaction's participant list, where `2 S1` is a
    /// coefficient and a species, not 2 in units of S1.
    pub parse_unit_suffixes: bool,
    nesting: usize,
    /// Set once a limit has been reported, so a deep input gives one
    /// diagnostic, not one per level.
    too_deep: bool,
}

impl<'l, 'a> ExprParser<'l, 'a> {
    pub fn new(lexer: &'l mut Lexer<'a>) -> ExprParser<'l, 'a> {
        ExprParser {
            lexer,
            parse_unit_suffixes: true,
            nesting: 0,
            too_deep: false,
        }
    }

    /// One expression, or `None` having reported at least one problem.
    /// An empty right-hand side (`p = ;`) is legal Antimony, so the statement
    /// parser checks [`ExprParser::at_expression_start`] before asking.
    pub fn parse(&mut self) -> Option<Expr> {
        self.nesting = 0;
        self.too_deep = false;
        self.parse_conditional().map(|node| node.expr)
    }

    /// Whether the next token could begin an expression.
    pub fn at_expression_start(&mut self) -> bool {
        can_start_expression(self.peek_kind(0))
    }

    // --- helpers -----------------------------------------------------------

    fn peek_kind(&mut self, ahead: usize) -> TokenKind {
        self.lexer.peek(ahead).kind
    }

    fn peek_span(&mut self) -> Span {
        self.lexer.peek(0).span
    }

    fn peek_is_keyword(&mut self, keyword: Keyword) -> bool {
        self.lexer.peek(0).is_keyword(keyword)
    }

    fn diagnostics(&mut self) -> &mut Diagnostics {
        self.lexer.diagnostics_mut()
    }

    /// Reports a syntax error at the next token.
    fn fail(&mut self, detail: &str) -> Option<Node> {
        let span = self.peek_span();
        self.diagnostics()
            .add(Code::SyntaxError, span, messages::syntax_error(detail));
        None
    }

    /// Consumes a `kind`, or reports what was there instead.
    fn expect(&mut self, kind: TokenKind) -> bool {
        if self.lexer.try_consume(kind).is_some() {
            return true;
        }
        let found = self.lexer.peek(0);
        let (span, found) = (found.span, found.describe());
        let expected = format!("\"{}\"", kind.spelling());
        self.diagnostics().add(
            Code::UnexpectedToken,
            span,
            messages::unexpected_token(&found, &expected),
        );
        false
    }

    fn report_too_deep(&mut self) {
        if !self.too_deep {
            self.too_deep = true;
            self.fail(&format!(
                "the expression is nested too deeply (the limits are {MAX_NESTING} levels of brackets, calls and operators, and a tree {MAX_DEPTH} deep)"
            ));
        }
    }

    /// Goes one nesting level down, or reports that it is too deep. A failed
    /// parse abandons the expression, so only a successful level comes back
    /// up through [`ExprParser::leave`]; `parse` starts from zero each time.
    fn enter(&mut self) -> Option<()> {
        if self.too_deep || self.nesting >= MAX_NESTING {
            self.report_too_deep();
            return None;
        }
        self.nesting += 1;
        Some(())
    }

    fn leave(&mut self, node: Node) -> Option<Node> {
        self.nesting -= 1;
        Some(node)
    }

    /// A new node over `children`, unless that makes the tree too deep.
    fn make(&mut self, kind: ExprKind, span: Span, children: &[usize]) -> Option<Node> {
        let depth = 1 + children.iter().copied().max().unwrap_or(0);
        if depth > MAX_DEPTH {
            self.report_too_deep();
            return None;
        }
        Some(Node {
            expr: Expr::new(kind, span),
            depth,
        })
    }

    fn binary(&mut self, op: BinaryOp, left: Node, right: Node) -> Option<Node> {
        let span = left.expr.span.union(right.expr.span);
        let depths = [left.depth, right.depth];
        let kind = ExprKind::Binary {
            op,
            left: Box::new(left.expr),
            right: Box::new(right.expr),
        };
        self.make(kind, span, &depths)
    }

    fn unary(&mut self, op: UnaryOp, op_span: Span, operand: Node) -> Option<Node> {
        let span = op_span.union(operand.expr.span);
        let depths = [operand.depth];
        let kind = ExprKind::Unary {
            op,
            operand: Box::new(operand.expr),
        };
        self.make(kind, span, &depths)
    }

    /// The binary operator the next token spells, if any. The words `or`
    /// and `and` are the operators wherever an operator can stand.
    fn binary_op_ahead(&mut self) -> Option<BinaryOp> {
        let token = self.lexer.peek(0);
        Some(match token.kind {
            TokenKind::OrOr => BinaryOp::Or,
            TokenKind::AndAnd => BinaryOp::And,
            TokenKind::Identifier if token.keyword == Some(Keyword::Or) => BinaryOp::Or,
            TokenKind::Identifier if token.keyword == Some(Keyword::And) => BinaryOp::And,
            TokenKind::Equal => BinaryOp::Equal,
            TokenKind::NotEqual => BinaryOp::NotEqual,
            TokenKind::Less => BinaryOp::Less,
            TokenKind::LessEqual => BinaryOp::LessEqual,
            TokenKind::Greater => BinaryOp::Greater,
            TokenKind::GreaterEqual => BinaryOp::GreaterEqual,
            TokenKind::Plus => BinaryOp::Add,
            TokenKind::Minus => BinaryOp::Subtract,
            TokenKind::Star => BinaryOp::Multiply,
            TokenKind::Slash => BinaryOp::Divide,
            TokenKind::Percent => BinaryOp::Modulo,
            TokenKind::Caret => BinaryOp::Power,
            _ => return None,
        })
    }

    // --- the grammar --------------------------------------------------------

    /// `condition ? then : otherwise`, the loosest level, right-associative.
    fn parse_conditional(&mut self) -> Option<Node> {
        self.enter()?;
        let condition = self.parse_binary(precedence::OR)?;
        if self.peek_kind(0) != TokenKind::Question {
            return self.leave(condition);
        }
        self.lexer.next_token();
        let then = self.parse_conditional()?;
        if !self.expect(TokenKind::Colon) {
            return None;
        }
        let otherwise = self.parse_conditional()?;
        let span = condition.expr.span.union(otherwise.expr.span);
        let depths = [condition.depth, then.depth, otherwise.depth];
        let kind = ExprKind::Conditional {
            condition: Box::new(condition.expr),
            then: Box::new(then.expr),
            otherwise: Box::new(otherwise.expr),
        };
        let node = self.make(kind, span, &depths)?;
        self.leave(node)
    }

    /// Everything that binds at least as tightly as `min`: an operand, then
    /// any binary operators of precedence `min` or tighter. Left-associative
    /// operators parse their right side one level tighter; `^`, which is
    /// right-associative, at its own level.
    fn parse_binary(&mut self, min: u8) -> Option<Node> {
        let mut node = self.parse_operand(min)?;
        loop {
            // `a--b` is the DNA operator, not two minus signs. Saying so
            // names the real problem rather than a stray token afterwards.
            if min <= precedence::ADDITIVE && self.peek_kind(0) == TokenKind::DashDash {
                let span = self.peek_span();
                self.diagnostics()
                    .add(Code::DnaNotSupported, span, messages::DNA_NOT_SUPPORTED);
                return None;
            }
            let Some(op) = self.binary_op_ahead() else {
                return Some(node);
            };
            let level = op.precedence();
            if level < min {
                return Some(node);
            }
            self.lexer.next_token();
            let right = if op == BinaryOp::Power {
                self.enter()?;
                let exponent = self.parse_binary(precedence::POWER)?;
                self.leave(exponent)?
            } else {
                self.parse_binary(level + 1)?
            };
            node = self.binary(op, node, right)?;
        }
    }

    /// An operand at level `min`. `!` and `not` may start one only where an
    /// operand at their own level may stand, so `a && !b` is fine and
    /// `a == !b` is not, as in Delphi §5.3. `not` followed by `(` is the
    /// function `not(…)`, as are `and(…)` and `or(…)`.
    fn parse_operand(&mut self, min: u8) -> Option<Node> {
        let is_not = self.peek_kind(0) == TokenKind::Not
            || (self.peek_is_keyword(Keyword::Not) && self.peek_kind(1) != TokenKind::LParen);
        if !is_not || min > precedence::NOT {
            return self.parse_unary();
        }
        let op_span = self.peek_span();
        self.lexer.next_token();
        self.enter()?;
        let operand = self.parse_binary(precedence::NOT)?;
        let node = self.unary(UnaryOp::Not, op_span, operand)?;
        self.leave(node)
    }

    /// Unary `+` and `-`, tighter than `^`.
    fn parse_unary(&mut self) -> Option<Node> {
        let op = match self.peek_kind(0) {
            TokenKind::Plus => UnaryOp::Plus,
            TokenKind::Minus => UnaryOp::Minus,
            _ => return self.parse_primary(),
        };
        let op_span = self.peek_span();
        self.lexer.next_token();
        self.enter()?;
        let operand = self.parse_unary()?;
        let node = self.unary(op, op_span, operand)?;
        self.leave(node)
    }

    /// Numbers, names, calls and brackets.
    fn parse_primary(&mut self) -> Option<Node> {
        let token = self.lexer.peek(0).clone();
        match token.kind {
            TokenKind::Number => {
                self.lexer.next_token();
                let mut span = token.span;
                let units = self.parse_unit_suffix(&mut span);
                self.make(
                    ExprKind::Number {
                        value: token.number,
                        units,
                    },
                    span,
                    &[],
                )
            }
            TokenKind::LParen => {
                self.lexer.next_token();
                let mut inner = self.parse_conditional()?;
                if !self.expect(TokenKind::RParen) {
                    return None;
                }
                // No node for the brackets: writing puts back exactly the
                // ones the meaning needs. The span covers them.
                inner.expr.span = token.span.union(self.lexer.current().span);
                Some(inner)
            }
            TokenKind::Identifier => {
                // A name directly followed by `(` is a call, so `and(a, b)` is
                // the function while `a and b` is the operator.
                if self.peek_kind(1) == TokenKind::LParen {
                    self.lexer.next_token();
                    return self.parse_call(token.text, token.span);
                }
                if let Some(constant) = Constant::from_word(&token.text) {
                    self.lexer.next_token();
                    return self.make(ExprKind::Constant(constant), token.span, &[]);
                }
                let (name, span) = self.parse_dotted_name();
                self.make(ExprKind::Symbol(name), span, &[])
            }
            TokenKind::DashDash => {
                self.diagnostics().add(
                    Code::DnaNotSupported,
                    token.span,
                    messages::DNA_NOT_SUPPORTED,
                );
                None
            }
            _ => self.fail(&format!("{} cannot start an expression", token.describe())),
        }
    }

    /// The arguments of a call whose name has been read; the next token is
    /// its `(`.
    fn parse_call(&mut self, name: String, name_span: Span) -> Option<Node> {
        self.lexer.next_token();
        let mut args = Vec::new();
        let mut depths = Vec::new();
        if self.lexer.try_consume(TokenKind::RParen).is_none() {
            loop {
                let arg = self.parse_conditional()?;
                depths.push(arg.depth);
                args.push(arg.expr);
                if self.lexer.try_consume(TokenKind::Comma).is_none() {
                    break;
                }
            }
            if !self.expect(TokenKind::RParen) {
                return None;
            }
        }
        let span = name_span.union(self.lexer.current().span);
        self.make(ExprKind::Call { name, args }, span, &depths)
    }

    /// `A.B.C`: a dotted name is one symbol. Dotted access belongs to modules,
    /// which are not built yet, but the name has to survive reading.
    fn parse_dotted_name(&mut self) -> (String, Span) {
        let first = self.lexer.next_token();
        let (mut name, mut span) = (first.text, first.span);
        while self.peek_kind(0) == TokenKind::Dot && self.peek_kind(1) == TokenKind::Identifier {
            self.lexer.next_token();
            let part = self.lexer.next_token();
            name.push('.');
            name.push_str(&part.text);
            span = span.union(part.span);
        }
        (name, span)
    }

    /// The unit after a number, which is exactly ONE identifier. The
    /// language reference's own example, `x = 40 foo/3 seconds`, gives `40`
    /// the unit `foo` and `3` the unit `seconds`: the `/` divides. Corpus case
    /// 00928, `S1 = 0.00015 mole/C`, divides by the compartment `C`.
    ///
    /// Not a unit: an identifier followed by `(` (a call), or a keyword, which
    /// the surrounding statement is waiting for (`at 5.7 after …`).
    fn parse_unit_suffix(&mut self, span: &mut Span) -> Option<String> {
        if !self.parse_unit_suffixes {
            return None;
        }
        let next = self.lexer.peek(0);
        if next.kind != TokenKind::Identifier || next.keyword.is_some() {
            return None;
        }
        if self.peek_kind(1) == TokenKind::LParen {
            return None;
        }
        let unit = self.lexer.next_token();
        *span = span.union(unit.span);
        Some(unit.text)
    }
}

/// Reads `text` as one complete expression. Anything after it is an error,
/// since hiding it would hide the mistake. Returns `None` on any error.
pub fn parse_expression(text: &str) -> (Option<Expr>, Diagnostics) {
    let mut lexer = Lexer::new(text);
    let mut parser = ExprParser::new(&mut lexer);
    let mut expr = parser.parse();
    if expr.is_some() && !lexer.peek(0).is_eof() {
        let rest = lexer.peek(0).describe();
        ExprParser::new(&mut lexer).fail(&format!("unexpected {rest} after the expression"));
        expr = None;
    }
    (expr, lexer.take_diagnostics())
}

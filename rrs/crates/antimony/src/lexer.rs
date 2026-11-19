//! The Antimony tokenizer (Delphi `uAntLexer`, Delphi §5.1).
//!
//! The Delphi lexer was itself a rewrite that closed nine defects of an
//! earlier one (Delphi §12.2). They stay closed here:
//!
//! - a leading `-` is never part of a number; the sign is a unary operator
//!   and belongs to the expression parser, or `a-1` would be two operands;
//! - every token is complete, including at the end, and the lexer parks on
//!   the end-of-file token instead of running past it;
//! - [`Lexer::peek`] looks any distance ahead, which telling a reaction from
//!   an assignment needs;
//! - comments are kept as trivia on the next token, so the writer can put
//!   them back;
//! - `inf` and `nan` are numbers.
//!
//! `layout` is an ordinary word. The Delphi lexer read `layout … end` and
//! `layout { … }` as one token (Delphi N1, §5.6), but Antimony has no such
//! block: its layout is dotted statements (`model.layout = on`), which the
//! parser reads (docs/antimony-layout-design.md).
//!
//! **There are no reserved words.** The lexer never makes a keyword token: it
//! makes an identifier and notes the keyword it spells, and the parser decides
//! from the position whether that means anything (`conversionFactor`, `at`,
//! `is` and the rest are ordinary names elsewhere).

use std::collections::VecDeque;
use std::fmt;

use crate::diagnostics::{Code, Diagnostics};
use crate::messages;
use crate::number;
use crate::span::{LineIndex, Pos, Span};

/// What a token is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TokenKind {
    // Structural.
    Eof,
    /// A statement terminator. A `\` before the line break suppresses it.
    Newline,
    /// May spell a keyword; see [`Token::keyword`].
    Identifier,
    Number,
    String,

    // Punctuation.
    Semicolon,
    Comma,
    Colon,
    Dot,
    LParen,
    RParen,
    LBrace,
    RBrace,
    LBracket,
    RBracket,

    // Markers.
    /// `$`, a boundary species.
    Dollar,
    /// `'`, a rate rule.
    Prime,
    At,
    /// `...`, DNA syntax: out of scope, but one token so the parser can name it.
    Ellipsis,

    // Assignment.
    /// `=`
    Assign,
    /// `:=`
    ColonAssign,

    // Arrows. `->` is REVERSIBLE and `=>` IRREVERSIBLE (Delphi §9.1.1).
    /// `->`, reversible.
    ArrowRev,
    /// `=>`, irreversible.
    ArrowIrrev,
    /// `<->`
    ArrowBiDir,
    /// `-|`
    Inhibits,
    /// `-o`
    Activates,
    /// `--`, DNA strand syntax: rejected by the parser (Delphi N2).
    DashDash,

    // Arithmetic.
    Plus,
    Minus,
    Star,
    Slash,
    Percent,
    Caret,

    // Comparison.
    /// `==`
    Equal,
    /// `!=`
    NotEqual,
    Less,
    LessEqual,
    Greater,
    GreaterEqual,

    // Logical.
    /// `&&`
    AndAnd,
    /// `||`
    OrOr,
    /// `!`
    Not,
    /// `?`, which with `:` forms the `?:` operator.
    Question,

    /// A character the language does not use; always reported.
    Unknown,
}

impl TokenKind {
    /// The fixed spelling of a punctuation or operator token, or `""` for
    /// kinds whose text varies.
    pub fn spelling(self) -> &'static str {
        use TokenKind::*;
        match self {
            Eof | Newline | Identifier | Number | String | Unknown => "",
            Semicolon => ";",
            Comma => ",",
            Colon => ":",
            Dot => ".",
            LParen => "(",
            RParen => ")",
            LBrace => "{",
            RBrace => "}",
            LBracket => "[",
            RBracket => "]",
            Dollar => "$",
            Prime => "'",
            At => "@",
            Ellipsis => "...",
            Assign => "=",
            ColonAssign => ":=",
            ArrowRev => "->",
            ArrowIrrev => "=>",
            ArrowBiDir => "<->",
            Inhibits => "-|",
            Activates => "-o",
            DashDash => "--",
            Plus => "+",
            Minus => "-",
            Star => "*",
            Slash => "/",
            Percent => "%",
            Caret => "^",
            Equal => "==",
            NotEqual => "!=",
            Less => "<",
            LessEqual => "<=",
            Greater => ">",
            GreaterEqual => ">=",
            AndAnd => "&&",
            OrOr => "||",
            Not => "!",
            Question => "?",
        }
    }
}

/// A word with a meaning in some statement position. Never reserved: see the
/// module documentation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Keyword {
    // Blocks and modules.
    Model,
    Module,
    End,
    Function,
    Import,
    // Declarations.
    Species,
    Compartment,
    Formula,
    Reaction,
    SubstanceOnly,
    Const,
    Var,
    Unit,
    // Attachments.
    Has,
    Is,
    In,
    Notes,
    ConversionFactor,
    // Events. `t0` says whether a trigger already true at the start counts
    // as having just become true.
    At,
    After,
    Priority,
    Delay,
    Persistent,
    FromTrigger,
    T0,
    // Other statements. `dna` is out of scope (Delphi N2), but a keyword so
    // the parser can say so.
    Delete,
    Constraint,
    Layout,
    Dna,
    // The fbc package.
    Maximize,
    Minimize,
    SboTerm,
    // Operators spelled as words.
    And,
    Or,
    Not,
    Xor,
    // Number literals spelled as words.
    Inf,
    NaN,
}

impl Keyword {
    pub const ALL: [Keyword; 38] = [
        Keyword::Model,
        Keyword::Module,
        Keyword::End,
        Keyword::Function,
        Keyword::Import,
        Keyword::Species,
        Keyword::Compartment,
        Keyword::Formula,
        Keyword::Reaction,
        Keyword::SubstanceOnly,
        Keyword::Const,
        Keyword::Var,
        Keyword::Unit,
        Keyword::Has,
        Keyword::Is,
        Keyword::In,
        Keyword::Notes,
        Keyword::ConversionFactor,
        Keyword::At,
        Keyword::After,
        Keyword::Priority,
        Keyword::Delay,
        Keyword::Persistent,
        Keyword::FromTrigger,
        Keyword::T0,
        Keyword::Delete,
        Keyword::Constraint,
        Keyword::Layout,
        Keyword::Dna,
        Keyword::Maximize,
        Keyword::Minimize,
        Keyword::SboTerm,
        Keyword::And,
        Keyword::Or,
        Keyword::Not,
        Keyword::Xor,
        Keyword::Inf,
        Keyword::NaN,
    ];

    /// The keyword `word` spells, case-sensitively.
    pub fn from_word(word: &str) -> Option<Keyword> {
        use Keyword::*;
        Some(match word {
            "model" => Model,
            "module" => Module,
            "end" => End,
            "function" => Function,
            "import" => Import,
            "species" => Species,
            "compartment" => Compartment,
            "formula" => Formula,
            "reaction" => Reaction,
            "substanceOnly" => SubstanceOnly,
            "const" => Const,
            "var" => Var,
            "unit" => Unit,
            "has" => Has,
            "is" => Is,
            "in" => In,
            "notes" => Notes,
            "conversionFactor" => ConversionFactor,
            "at" => At,
            "after" => After,
            "priority" => Priority,
            "delay" => Delay,
            "persistent" => Persistent,
            "fromTrigger" => FromTrigger,
            "t0" => T0,
            "delete" => Delete,
            "constraint" => Constraint,
            "layout" => Layout,
            "dna" => Dna,
            "maximize" => Maximize,
            "minimize" => Minimize,
            "sboTerm" => SboTerm,
            "and" => And,
            "or" => Or,
            "not" => Not,
            "xor" => Xor,
            "inf" => Inf,
            "nan" => NaN,
            _ => return None,
        })
    }

    /// How the keyword is spelled.
    pub fn spelling(self) -> &'static str {
        use Keyword::*;
        match self {
            Model => "model",
            Module => "module",
            End => "end",
            Function => "function",
            Import => "import",
            Species => "species",
            Compartment => "compartment",
            Formula => "formula",
            Reaction => "reaction",
            SubstanceOnly => "substanceOnly",
            Const => "const",
            Var => "var",
            Unit => "unit",
            Has => "has",
            Is => "is",
            In => "in",
            Notes => "notes",
            ConversionFactor => "conversionFactor",
            At => "at",
            After => "after",
            Priority => "priority",
            Delay => "delay",
            Persistent => "persistent",
            FromTrigger => "fromTrigger",
            T0 => "t0",
            Delete => "delete",
            Constraint => "constraint",
            Layout => "layout",
            Dna => "dna",
            Maximize => "maximize",
            Minimize => "minimize",
            SboTerm => "sboTerm",
            And => "and",
            Or => "or",
            Not => "not",
            Xor => "xor",
            Inf => "inf",
            NaN => "nan",
        }
    }
}

/// One token.
#[derive(Clone, Debug, PartialEq)]
pub struct Token {
    pub kind: TokenKind,
    /// The identifier, the operator's spelling, a string's decoded contents,
    /// or a number as written.
    pub text: String,
    /// The keyword an identifier spells, if any. Set on the `inf` and `nan`
    /// number tokens too.
    pub keyword: Option<Keyword>,
    /// The value of a number token; 0 otherwise.
    pub number: f64,
    pub span: Span,
    /// The comments just before this token, one per line. A comment belongs
    /// to the statement that follows it, so a comment at the end of a line
    /// carries over to the next line's first token, and the comments at the
    /// end of a file are on the end-of-file token.
    pub trivia: String,
}

impl Token {
    fn new(kind: TokenKind, text: impl Into<String>, span: Span) -> Token {
        Token {
            kind,
            text: text.into(),
            keyword: None,
            number: 0.0,
            span,
            trivia: String::new(),
        }
    }

    /// For messages: `"->"`, `identifier "S1"`, `number 2.5`, `end of line`.
    pub fn describe(&self) -> String {
        match self.kind {
            TokenKind::Eof => "end of file".into(),
            TokenKind::Newline => "end of line".into(),
            TokenKind::Identifier => format!("identifier \"{}\"", self.text),
            TokenKind::Number => format!("number {}", number::format(self.number)),
            TokenKind::String => format!("string \"{}\"", self.text),
            TokenKind::Unknown => format!("character \"{}\"", self.text),
            _ => format!("\"{}\"", self.text),
        }
    }

    /// Whether this is an identifier spelling `keyword`.
    pub fn is_keyword(&self, keyword: Keyword) -> bool {
        self.kind == TokenKind::Identifier && self.keyword == Some(keyword)
    }

    /// Whether this is an identifier spelling one of `keywords`.
    pub fn is_any_keyword(&self, keywords: &[Keyword]) -> bool {
        self.kind == TokenKind::Identifier && self.keyword.is_some_and(|k| keywords.contains(&k))
    }

    pub fn is_eof(&self) -> bool {
        self.kind == TokenKind::Eof
    }

    /// `;`, a line break or the end of the file (Delphi §5.1).
    pub fn is_terminator(&self) -> bool {
        matches!(
            self.kind,
            TokenKind::Semicolon | TokenKind::Newline | TokenKind::Eof
        )
    }

    pub fn has_trivia(&self) -> bool {
        !self.trivia.is_empty()
    }
}

impl fmt::Display for Token {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.describe())
    }
}

/// `A`–`Z`, `a`–`z` or `_`.
pub fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}

/// `A`–`Z`, `a`–`z`, `0`–`9` or `_`.
pub fn is_ident_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// `0`–`9`.
pub fn is_digit(c: char) -> bool {
    c.is_ascii_digit()
}

/// `e` or `E`.
pub fn is_exponent_marker(c: char) -> bool {
    c == 'e' || c == 'E'
}

/// Reads tokens from a source, on demand.
pub struct Lexer<'a> {
    source: &'a str,
    bytes: &'a [u8],
    /// Byte offset of the next character to read.
    pos: usize,
    lines: LineIndex,
    diagnostics: Diagnostics,
    lookahead: VecDeque<Token>,
    current: Token,
    pending_trivia: String,
    /// Tokens made so far, not counting end of file. A line break before the
    /// first one is not a terminator.
    produced: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(source: &'a str) -> Lexer<'a> {
        let pos = if source.starts_with('\u{feff}') { 3 } else { 0 };
        Lexer {
            source,
            bytes: source.as_bytes(),
            pos,
            lines: LineIndex::new(source),
            diagnostics: Diagnostics::new(),
            lookahead: VecDeque::new(),
            current: Token::new(TokenKind::Eof, "", Span::default()),
            pending_trivia: String::new(),
            produced: 0,
        }
    }

    /// The next token, consumed. At the end of the source it returns the
    /// end-of-file token, every time.
    pub fn next_token(&mut self) -> Token {
        let token = match self.lookahead.pop_front() {
            Some(token) => token,
            None => self.scan_token(),
        };
        self.current = token.clone();
        token
    }

    /// A token ahead, not consumed: `peek(0)` is what [`Lexer::next_token`]
    /// returns next. Past the end it is the end-of-file token.
    pub fn peek(&mut self, ahead: usize) -> &Token {
        while self.lookahead.len() <= ahead {
            let token = self.scan_token();
            self.lookahead.push_back(token);
        }
        &self.lookahead[ahead]
    }

    /// Consumes the next token only if it is a `kind`.
    pub fn try_consume(&mut self, kind: TokenKind) -> Option<Token> {
        (self.peek(0).kind == kind).then(|| self.next_token())
    }

    /// Consumes the next token only if it spells `keyword`.
    pub fn try_consume_keyword(&mut self, keyword: Keyword) -> bool {
        let matched = self.peek(0).is_keyword(keyword);
        if matched {
            self.next_token();
        }
        matched
    }

    /// Every remaining token, end of file included.
    pub fn tokens(&mut self) -> Vec<Token> {
        let mut tokens = Vec::with_capacity(self.source.len() / 5 + 16);
        loop {
            let token = self.next_token();
            let end = token.is_eof();
            tokens.push(token);
            if end {
                return tokens;
            }
        }
    }

    /// The token most recently returned by [`Lexer::next_token`].
    pub fn current(&self) -> &Token {
        &self.current
    }

    pub fn source(&self) -> &'a str {
        self.source
    }

    /// The problems found so far.
    pub fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics
    }

    /// The list the parsers built on this lexer add their problems to, so
    /// that every problem in a source is in one list, in the order found.
    pub fn diagnostics_mut(&mut self) -> &mut Diagnostics {
        &mut self.diagnostics
    }

    /// Takes the problems found so far, leaving none.
    pub fn take_diagnostics(&mut self) -> Diagnostics {
        std::mem::take(&mut self.diagnostics)
    }

    /// The position of a byte offset in this source.
    pub fn pos_at(&mut self, offset: usize) -> Pos {
        self.lines.pos(self.source, offset)
    }

    // --- characters --------------------------------------------------------

    fn byte(&self, ahead: usize) -> u8 {
        self.bytes.get(self.pos + ahead).copied().unwrap_or(0)
    }

    fn at_end(&self) -> bool {
        self.pos >= self.bytes.len()
    }

    fn is_line_break(b: u8) -> bool {
        b == b'\n' || b == b'\r'
    }

    /// Steps over a line break: `\r\n` is one.
    fn skip_line_break(&mut self) {
        if self.byte(0) == b'\r' && self.byte(1) == b'\n' {
            self.pos += 2;
        } else {
            self.pos += 1;
        }
    }

    fn here(&mut self) -> Pos {
        self.lines.pos(self.source, self.pos)
    }

    fn span_from(&mut self, start: Pos) -> Span {
        Span::new(start, self.here())
    }

    fn make(&mut self, kind: TokenKind, text: impl Into<String>, start: Pos) -> Token {
        let span = self.span_from(start);
        let mut token = Token::new(kind, text, span);
        token.trivia = std::mem::take(&mut self.pending_trivia);
        token
    }

    // --- trivia ------------------------------------------------------------

    fn append_trivia(&mut self, text: &str) {
        if !self.pending_trivia.is_empty() {
            self.pending_trivia.push('\n');
        }
        self.pending_trivia.push_str(text);
    }

    fn read_line_comment(&mut self) {
        let start = self.pos;
        while !self.at_end() && !Self::is_line_break(self.byte(0)) {
            self.pos += 1;
        }
        self.append_trivia(&self.source[start..self.pos]);
    }

    fn read_block_comment(&mut self) {
        let start = self.pos;
        let open = self.here();
        self.pos += 2;
        while !self.at_end() {
            if self.byte(0) == b'*' && self.byte(1) == b'/' {
                self.pos += 2;
                self.append_trivia(&self.source[start..self.pos]);
                return;
            }
            self.pos += 1;
        }
        // Ran off the end. Block comments do not nest, so there is nothing
        // to recover to but the end; report it and keep what there is.
        let span = self.span_from(open);
        self.diagnostics.add(
            Code::UnterminatedComment,
            span,
            messages::UNTERMINATED_COMMENT,
        );
        self.append_trivia(&self.source[start..self.pos]);
    }

    /// Skips spaces, comments and line continuations. Returns where the
    /// first line break crossed began, if one was.
    fn skip_trivia(&mut self) -> Option<Pos> {
        let mut newline = None;
        while !self.at_end() {
            match self.byte(0) {
                b' ' | b'\t' | 0x0b | 0x0c => self.pos += 1,
                // A `\` just before a line break continues the statement.
                // Anywhere else it is not part of the language.
                b'\\' if Self::is_line_break(self.byte(1)) => {
                    self.pos += 1;
                    self.skip_line_break();
                }
                b'\n' | b'\r' => {
                    if newline.is_none() {
                        newline = Some(self.here());
                    }
                    self.skip_line_break();
                }
                // The language reference gives both `#` and `//`. The corpus
                // never uses `#`, because libantimony always writes `//`.
                b'#' => self.read_line_comment(),
                b'/' if self.byte(1) == b'/' => self.read_line_comment(),
                b'/' if self.byte(1) == b'*' => self.read_block_comment(),
                _ => break,
            }
        }
        newline
    }

    // --- scanners ----------------------------------------------------------

    fn scan_token(&mut self) -> Token {
        let newline = self.skip_trivia();

        // A run of blank lines is one terminator, and a line break before
        // the first token is none at all. The pending comments are NOT
        // attached to the newline: a comment belongs to the next statement,
        // so it carries across the line break.
        if let (Some(at), true) = (newline, self.produced > 0) {
            self.produced += 1;
            return Token::new(TokenKind::Newline, "", Span::from_pos(at, 1));
        }

        let start = self.here();
        if self.at_end() {
            // Parked: every later call returns end of file too. Comments at
            // the end of the file are on this token, so none are lost.
            return self.make(TokenKind::Eof, "", start);
        }

        self.produced += 1;
        let b = self.byte(0);
        if b.is_ascii_alphabetic() || b == b'_' {
            return self.scan_word(start);
        }
        if b.is_ascii_digit() || (b == b'.' && self.byte(1).is_ascii_digit()) {
            return self.scan_number(start);
        }
        if b == b'"' {
            return self.scan_string(start);
        }
        if b == b'`' && self.byte(1) == b'`' && self.byte(2) == b'`' {
            return self.scan_fenced_string(start);
        }
        self.scan_operator(start)
    }

    /// An identifier, or one of the four word spellings of a number.
    fn scan_word(&mut self, start: Pos) -> Token {
        let from = self.pos;
        while is_ident_char(self.byte(0) as char) {
            self.pos += 1;
        }
        let text = &self.source[from..self.pos];

        // Exactly four spellings are numbers, and no others. libantimony
        // writes `INF` and `NaN` (case 00950), and case 01810 declares
        // `Inf`, `INF_` and `inf_` as ordinary parameters, so this is a
        // closed list, not a case-insensitive match.
        let value = match text {
            "inf" | "INF" => Some((f64::INFINITY, Keyword::Inf)),
            "nan" | "NaN" => Some((f64::NAN, Keyword::NaN)),
            _ => None,
        };
        if let Some((value, keyword)) = value {
            let mut token = self.make(TokenKind::Number, text, start);
            token.number = value;
            token.keyword = Some(keyword);
            return token;
        }

        let mut token = self.make(TokenKind::Identifier, text, start);
        token.keyword = Keyword::from_word(text);
        token
    }

    /// Entered on a digit, or a `.` followed by a digit. Never on a `-`.
    fn scan_number(&mut self, start: Pos) -> Token {
        let from = self.pos;
        while self.byte(0).is_ascii_digit() {
            self.pos += 1;
        }
        if self.byte(0) == b'.' {
            self.pos += 1;
            while self.byte(0).is_ascii_digit() {
                self.pos += 1;
            }
        }
        if is_exponent_marker(self.byte(0) as char) {
            let mark = self.pos;
            self.pos += 1;
            if matches!(self.byte(0), b'+' | b'-') {
                self.pos += 1;
            }
            if self.byte(0).is_ascii_digit() {
                while self.byte(0).is_ascii_digit() {
                    self.pos += 1;
                }
            } else {
                // Not an exponent after all: `2e` leaves the `e` for the
                // next token.
                self.pos = mark;
            }
        }

        let text = &self.source[from..self.pos];
        let mut token = self.make(TokenKind::Number, text, start);
        match number::parse(text) {
            Some(value) => token.number = value,
            None => {
                self.diagnostics.add(
                    Code::InvalidNumber,
                    token.span,
                    messages::invalid_number(text),
                );
            }
        }
        token
    }

    /// A `"…"` string, with `\n`, `\r`, `\t`, `\\` and `\"` decoded. Any other
    /// escape keeps both characters, so nothing is lost. A string ends at a
    /// line break, which is reported as unterminated.
    fn scan_string(&mut self, start: Pos) -> Token {
        self.pos += 1;
        let mut text = String::new();
        let mut run = self.pos;
        let mut terminated = false;
        while !self.at_end() {
            match self.byte(0) {
                b'"' => {
                    text.push_str(&self.source[run..self.pos]);
                    self.pos += 1;
                    terminated = true;
                    break;
                }
                b'\n' | b'\r' => break,
                b'\\' if self.pos + 1 < self.bytes.len() => {
                    text.push_str(&self.source[run..self.pos]);
                    self.pos += 1;
                    let escaped = self.source[self.pos..].chars().next().unwrap_or('\\');
                    match escaped {
                        'n' => text.push('\n'),
                        'r' => text.push('\r'),
                        't' => text.push('\t'),
                        '\\' => text.push('\\'),
                        '"' => text.push('"'),
                        other => {
                            text.push('\\');
                            text.push(other);
                        }
                    }
                    self.pos += escaped.len_utf8();
                    run = self.pos;
                }
                _ => self.pos += 1,
            }
        }
        if !terminated {
            text.push_str(&self.source[run..self.pos]);
        }

        let token = self.make(TokenKind::String, text, start);
        if !terminated {
            self.diagnostics.add(
                Code::UnterminatedString,
                token.span,
                messages::UNTERMINATED_STRING,
            );
        }
        token
    }

    /// A ```` ```…``` ```` string, the language reference's form for a
    /// multi-line `model notes`. Read raw, with no escapes: real notes are
    /// XHTML, and a `\` in one is a backslash. A line break straight after
    /// the opening fence and one straight before the closing fence are not
    /// part of the text, and `\r\n` becomes `\n`.
    fn scan_fenced_string(&mut self, start: Pos) -> Token {
        self.pos += 3;
        if Self::is_line_break(self.byte(0)) {
            self.skip_line_break();
        }
        let from = self.pos;
        while !self.at_end()
            && !(self.byte(0) == b'`' && self.byte(1) == b'`' && self.byte(2) == b'`')
        {
            self.pos += 1;
        }
        let mut text = &self.source[from..self.pos];
        let terminated = !self.at_end();
        if terminated {
            self.pos += 3;
        }
        text = text.strip_suffix('\n').unwrap_or(text);
        text = text.strip_suffix('\r').unwrap_or(text);
        let text = text.replace("\r\n", "\n");

        let token = self.make(TokenKind::String, text, start);
        if !terminated {
            self.diagnostics.add(
                Code::UnterminatedString,
                token.span,
                messages::UNTERMINATED_STRING,
            );
        }
        token
    }

    fn scan_operator(&mut self, start: Pos) -> Token {
        use TokenKind::*;
        let (b1, b2) = (self.byte(1), self.byte(2));
        let kind = match self.byte(0) {
            b';' => Some(Semicolon),
            b',' => Some(Comma),
            // A leading-dot number was taken before we got here.
            b'.' if b1 == b'.' && b2 == b'.' => Some(Ellipsis),
            b'.' => Some(Dot),
            b'(' => Some(LParen),
            b')' => Some(RParen),
            b'{' => Some(LBrace),
            b'}' => Some(RBrace),
            b'[' => Some(LBracket),
            b']' => Some(RBracket),
            b'$' => Some(Dollar),
            b'\'' => Some(Prime),
            b'@' => Some(At),
            b'+' => Some(Plus),
            b'*' => Some(Star),
            b'/' => Some(Slash),
            b'%' => Some(Percent),
            b'^' => Some(Caret),
            b'?' => Some(Question),
            b':' if b1 == b'=' => Some(ColonAssign),
            b':' => Some(Colon),
            // `=>` is IRREVERSIBLE (Delphi §9.1.1).
            b'=' if b1 == b'>' => Some(ArrowIrrev),
            b'=' if b1 == b'=' => Some(Equal),
            b'=' => Some(Assign),
            b'<' if b1 == b'-' && b2 == b'>' => Some(ArrowBiDir),
            b'<' if b1 == b'=' => Some(LessEqual),
            b'<' => Some(Less),
            b'>' if b1 == b'=' => Some(GreaterEqual),
            b'>' => Some(Greater),
            b'!' if b1 == b'=' => Some(NotEqual),
            b'!' => Some(Not),
            b'&' if b1 == b'&' => Some(AndAnd),
            b'|' if b1 == b'|' => Some(OrOr),
            // `->` is REVERSIBLE (Delphi §9.1.1).
            b'-' if b1 == b'>' => Some(ArrowRev),
            b'-' if b1 == b'|' => Some(Inhibits),
            // DNA strand syntax: named here so the parser can reject it.
            b'-' if b1 == b'-' => Some(DashDash),
            // `-o` only when the `o` stands alone, or `x -old` would be
            // `x`, `-o`, `ld`.
            b'-' if b1 == b'o' && !is_ident_char(b2 as char) => Some(Activates),
            // `-(` is the unknown-interaction arrow, but deliberately not
            // one token: a minus before a bracket is far more common and
            // looks the same. The parser recognises it by position.
            b'-' => Some(Minus),
            _ => None,
        };
        if let Some(kind) = kind {
            self.pos += kind.spelling().len();
            return self.make(kind, kind.spelling(), start);
        }

        // Anything else is reported and consumed, so the lexer cannot stall
        // on it and the parser can carry on.
        let c = self.source[self.pos..].chars().next().unwrap_or('\0');
        self.pos += c.len_utf8();
        let token = self.make(Unknown, c, start);
        self.diagnostics.add(
            Code::UnexpectedChar,
            token.span,
            messages::unexpected_char(c),
        );
        token
    }
}

/// Every token in `source`, end of file included, and the problems found.
pub fn tokenize(source: &str) -> (Vec<Token>, Diagnostics) {
    let mut lexer = Lexer::new(source);
    let tokens = lexer.tokens();
    (tokens, lexer.take_diagnostics())
}

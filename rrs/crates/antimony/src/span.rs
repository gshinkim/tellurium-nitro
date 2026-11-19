//! Source positions (SPEC.md §4).
//!
//! Every token, expression node and model item carries a [`Span`], so a
//! diagnostic can point at the text and an editor can find what produced an
//! item. Positions are made by [`LineIndex`] from a byte offset, never counted
//! by hand as a scanner moves: the Delphi version's line numbers drifted
//! whenever a scanner crossed a line break without counting it, and working
//! from the offset makes that impossible.

use std::fmt;

/// A position in the source text.
///
/// `Pos::default()` (line 0) means "no position", for items built in code
/// rather than read from text.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Pos {
    /// Byte offset into the source, 0-based.
    pub offset: usize,
    /// Line, 1-based.
    pub line: u32,
    /// Column, 1-based, counted in characters (Unicode scalar values).
    pub column: u32,
}

impl Pos {
    pub fn new(line: u32, column: u32, offset: usize) -> Pos {
        Pos {
            offset,
            line,
            column,
        }
    }

    pub fn is_unknown(&self) -> bool {
        self.line == 0
    }
}

impl fmt::Display for Pos {
    /// `3:17`, or `<unknown>`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_unknown() {
            f.write_str("<unknown>")
        } else {
            write!(f, "{}:{}", self.line, self.column)
        }
    }
}

/// A stretch of source text: `start` up to, not including, `end`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Span {
    pub start: Pos,
    /// Exclusive: one past the last character.
    pub end: Pos,
}

impl Span {
    pub fn new(start: Pos, end: Pos) -> Span {
        Span { start, end }
    }

    /// The span of `len` bytes from `start`, all on one line. Only right for
    /// ASCII text, where a byte is a column.
    pub fn from_pos(start: Pos, len: usize) -> Span {
        let mut end = start;
        end.offset += len;
        end.column += len as u32;
        Span { start, end }
    }

    pub fn is_unknown(&self) -> bool {
        self.start.is_unknown()
    }

    /// Length in bytes; 0 for an unknown span.
    pub fn len(&self) -> usize {
        if self.is_unknown() || self.end.is_unknown() {
            0
        } else {
            self.end.offset.saturating_sub(self.start.offset)
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The smallest span covering both. An unknown span adds nothing.
    pub fn union(&self, other: Span) -> Span {
        if self.is_unknown() {
            return other;
        }
        if other.is_unknown() {
            return *self;
        }
        Span {
            start: if other.start.offset < self.start.offset {
                other.start
            } else {
                self.start
            },
            end: if other.end.offset > self.end.offset {
                other.end
            } else {
                self.end
            },
        }
    }

    /// The text this span covers in `source`, the text it was made from.
    /// Empty for an unknown span or one that does not fit `source`.
    pub fn slice_of<'a>(&self, source: &'a str) -> &'a str {
        if self.is_unknown() {
            return "";
        }
        source
            .get(self.start.offset..self.start.offset + self.len())
            .unwrap_or("")
    }
}

impl fmt::Display for Span {
    /// `3:4-7` on one line, `3:4..5:1` across lines, or `<unknown>`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_unknown() {
            f.write_str("<unknown>")
        } else if self.start.line == self.end.line {
            write!(
                f,
                "{}:{}-{}",
                self.start.line, self.start.column, self.end.column
            )
        } else {
            write!(f, "{}..{}", self.start, self.end)
        }
    }
}

/// Turns byte offsets into positions. Built once per source.
///
/// `\n`, `\r\n` and a lone `\r` each end a line. A leading byte-order mark is
/// not part of line 1.
#[derive(Clone, Debug)]
pub struct LineIndex {
    /// Byte offset of the first character of each line.
    starts: Vec<usize>,
    /// The last position computed, so that positions asked for in order
    /// (as the lexer does) cost only the text since the previous one.
    last: Pos,
}

impl LineIndex {
    pub fn new(source: &str) -> LineIndex {
        let bytes = source.as_bytes();
        let first = if source.starts_with('\u{feff}') { 3 } else { 0 };
        let mut starts = vec![first];
        for (i, &b) in bytes.iter().enumerate().skip(first) {
            match b {
                b'\n' => starts.push(i + 1),
                b'\r' if bytes.get(i + 1) != Some(&b'\n') => starts.push(i + 1),
                _ => {}
            }
        }
        LineIndex {
            starts,
            last: Pos::new(1, 1, first),
        }
    }

    /// The position of byte `offset`, which must be a character boundary in
    /// the source this index was built from.
    pub fn pos(&mut self, source: &str, offset: usize) -> Pos {
        let line = self.starts.partition_point(|&s| s <= offset).max(1);
        let line_start = self.starts[line - 1];
        let from = if line as u32 == self.last.line && self.last.offset <= offset {
            self.last
        } else {
            Pos::new(line as u32, 1, line_start)
        };
        let counted = source
            .get(from.offset..offset)
            .map_or(0, |text| text.chars().count());
        let pos = Pos::new(line as u32, from.column + counted as u32, offset);
        self.last = pos;
        pos
    }
}

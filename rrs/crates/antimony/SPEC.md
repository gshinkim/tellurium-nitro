# The `antimony` crate — Rust design

This crate is a Rust port of the author's Delphi Antimony implementation
("the Delphi version"). **The Delphi `SPECIFICATION.md` stays the normative
description of the language surface and its semantics**, and its decisions
record (§13, 25 decisions) stands unless this file says otherwise. This file
records only what is different in Rust: the shape of the code, and every
place the port departs from the Delphi design (§9). Where a later milestone
settles something left open here, it is written in here at that time.

The Antimony language reference,
<https://tellurium.readthedocs.io/en/latest/antimony.html>, outranks both on
what the *language* is (Delphi decision 19).

Section numbers written "Delphi §n" refer to the Delphi `SPECIFICATION.md`.

Written 2026-10-13 in milestone A0; §5.1, §5.2 and R14–R17 added in A1; §5.3 and R18–R21 in A2; §3.3's model as built, §5.4, the corpus index and differential tier in §8, and R22–R25 in A3; §5.5, §6 as built and R26–R27 in A4; the A5 note in §5.4; §5.6, R28–R29 and the end of §10 in A6; §5.7, R30–R32 and the A7 parts of §2, §6, §7 and §8 in A7; §6 as built, the BioModels tier in §8, and R33–R35 in A8; R36 in A9.

---

## 1. Name, place, licence

- **Package `antimony-rs`, library `antimony`.** The name `antimony` on
  crates.io belongs to an unrelated crate (a 2017 stream-processing topology
  builder). The library keeps the language's name, so users write
  `use antimony::…`; only the `Cargo.toml` line says `antimony-rs`. The folder
  is `crates/antimony`.
- **Licence: MIT OR Apache-2.0**, chosen by the author on 2026-10-13. This
  replaces Delphi decision 5 (BSD, to match libantimony) for this crate. The
  code is ported from the author's own Delphi source. If libantimony's C++
  source is ever consulted for a point of semantics, say so in the commit
  message and in this file; its BSD licence is compatible.
- **Movable.** Everything the crate needs is inside `crates/antimony/`:
  `Cargo.toml` uses no `workspace = true` keys, tests find their files from
  `CARGO_MANIFEST_DIR`, and nothing in the crate refers to the application
  that first uses it. The one exception is `PLAN.md`, the working plan, which
  is rewritten or dropped when the crate moves.

## 2. Dependencies and targets

- **The core has no dependencies.** Lexer, parser, semantics and the Antimony
  writer are plain `std`.
- **The `sbml` feature** (built in A7) adds `sbml-io` and nothing else, and
  the crate still builds for wasm32 with it (CI checks both). Until `sbml-io` is published the
  dependency is a path with a version (`{ path = "../sbml-io", version =
  "0.1", optional = true }`), so publishing this crate requires publishing
  `sbml-io` first.
- **Dev-dependencies** are allowed for tests and examples, since they never
  reach a user, but are kept few: `serde_json` (the corpus index) and
  `sbml-io` (reading the oracle's SBML).
- **Targets:** native (Windows, macOS, Linux) and `wasm32-unknown-unknown`.
  The core therefore never touches the file system, the clock or threads; it
  takes `&str` and returns values. Reading a file is the caller's job.
- `edition = "2024"`, `rust-version = "1.85"` (the first with that edition).
- `#![forbid(unsafe_code)]`.

## 3. The data model

The Delphi object model (Delphi §4) is a class hierarchy with owned object
lists and object references between items. In Rust it becomes **plain
structs with public fields**, `Vec` for lists, `Option` for things that may be
absent, and enums for alternatives. No trait objects, no `Rc`, no
interior mutability, and no getters or setters for plain data.

### 3.1 References between items are by name

Delphi binds names to objects (`TAntSymbol.Compartment: TAntSymbol`,
`TAntParticipant.Species`, the expression node's `Binding`). In Rust every
reference is the **identifier as a `String`** — `compartment: Option<String>`,
a participant's `species: String`, an expression's `Symbol(String)` — and is
looked up in its module when needed.

Why: Antimony allows forward references and type promotion, a caller may
edit the model freely between parse and write, and names are what both
writers emit. Pointers or indices would have to be kept in step with every
edit; names cannot go stale. Lookups are by name through methods on
`Module` (`module.symbol("S1")`); passes that do many lookups (the semantic
pass, the writers) build a temporary `HashMap<&str, usize>` for the duration
of the pass. The model itself stores only `Vec`s, so editing its public
fields can never leave an index out of date.

### 3.2 Absence is `Option`, not a sentinel

Delphi's sentinels become `Option`: an `SboTerm` of −1 becomes
`sbo_term: Option<u32>`; an empty `DisplayName`, `Notes`, `UnitsName` or
`CompartmentName` becomes `None`. Where the Delphi version uses a flag to
record that something was *written but empty* (`HasEmptyValue`,
`HasEmptyKineticLaw`, Delphi §9.1.4 #3), the flag stays a `bool` beside the
`Option`, because "written empty" and "not written" are different facts that
both have to survive a round trip.

### 3.3 One `Symbol`, with promotion

Delphi §4.1's unified symbol becomes one struct:

```rust
pub struct Symbol {
    pub name: String,
    pub kind: SymbolKind,          // Delphi TAntSymbolType
    pub constness: Constness,      // Default, Const, Var
    pub is_boundary: bool,
    pub is_substance_only: bool,
    pub compartment: Option<String>,
    pub initial_value: Option<Expr>,
    pub has_empty_value: bool,
    pub rule: Option<Rule>,        // assignment or rate rule, with its expression
    pub units: Option<String>,
    pub conversion_factor: Option<String>,
    pub display_name: Option<String>,
    pub notes: Option<String>,
    pub sbo_term: Option<u32>,
    pub cv_terms: Vec<CvTerm>,
    pub metadata: Metadata,
    pub comment: Option<String>,
    pub is_implicit: bool,         // never declared: §6.5's typo case
    pub source_path: String,       // seam 3 (Delphi §2.2.2)
    pub span: Span,
    pub decl_span: Span,
}
```

`SymbolKind` keeps all thirteen Delphi values (`Unknown`, `Species`,
`Parameter`, `Compartment`, `Reaction`, `Event`, `Function`, `Submodule`,
`UnitDef`, `Constraint`, `Interaction`, `Deleted`, `Stoichiometry`), and the
promotion table of Delphi §6.2 (with `parameter → stoichiometry`, decision
17) is `SymbolKind::can_promote_to`.

**As built in A3** (`src/model/items.rs`), the field list above is exact.
`Rule { kind: RuleKind, expr: Option<Expr> }`, where `RuleKind` is
`Assignment` or `Rate`: Delphi's `TAntSymbolRule` also had `srNone` (here
`rule: None`) and `srInitialAssign`, which nothing ever set (an initial
assignment is an `initial_value` that is not constant). `expr` is `None` for
a rule written empty, `p' = ;`.

The other model types follow Delphi §4.2 one for one, with these choices:

- `Participant { species, stoichiometry: Option<Expr>, is_boundary, span }`;
  `None` means 1.
- `Reaction` always has a `name`: an unlabelled one is named `_J1`, `_J2`, …
  by the parser, as in Delphi. `has_empty_kinetic_law` beside `kinetic_law`.
- `Interaction`, `Event`, `Constraint` and `Objective` have
  `name: Option<String>`; an unnamed one has no invented name.
  `AlgebraicRule.name` is a `String`, because only the labelled form exists.
- `Event.trigger` is an `Expr`, not an `Option`: the parser only adds an
  event once it has one.
- `Constraint.kind` is an enum: `ConstraintKind::Expr(Expr)`, or
  `ConstraintKind::Range { lower, subject, upper, lower_inclusive,
  upper_inclusive }` for a chained comparison (decision 18). Delphi had four
  nullable fields.
- `RawSection.kind` is `RawKind::Annotation`, not a string. (`RawKind::Layout`
  went with the `layout` block, R37.)
- `Module.layout` holds Antimony 3's layout statements in the order read, each
  a `LayoutStatement { path, value, comment, span }`: the dotted parts as
  written (`["J0", "S1", "arc2", "b1"]`, a null node's part `"--"`) and a
  `LayoutValue` (`Number`, `Point(x, y)` for a braced pair of numbers, `List`
  for any other braced list, `Word` for one or more words, `Text` for a
  string). The parser takes a statement as layout when it starts `model.`
  (but not `model.sboTerm`), `species.`, `reaction.` or `compartment.`, or is
  a dotted name with a layout word (`pos`, `size`, `b1`, `color`, `arc2`, …,
  any case) or `--` after its first part. The writer puts them after
  everything else, `path = value`, each with its comments.
- `Module.name` is `Option<String>`: `None` for the statements of a file
  with no `model` block.
- Every optional string is `Option<String>` (R4). So `S1 notes ""` keeps an
  empty note, where Delphi's empty string meant "no note" and dropped it.
- Lookups by name: `Module::symbol`, `symbol_mut`, `reaction`, `unit`,
  `unit_mut`, `function`, `symbols_of(kind)`; `Document::module`,
  `module_index`, `main_module`, `main_module_mut`.

### 3.4 Expressions are an owned `enum`

```rust
pub struct Expr {
    pub kind: ExprKind,
    pub span: Span,
}

pub enum ExprKind {
    Number { value: f64, units: Option<String> },   // 2.5, inf, nan; `40 foo`
    Constant(Constant),                             // pi, exponentiale, true, false
    Symbol(String),                                 // possibly dotted: A.S1
    Unary { op: UnaryOp, operand: Box<Expr> },
    Binary { op: BinaryOp, left: Box<Expr>, right: Box<Expr> },
    Conditional { condition: Box<Expr>, then: Box<Expr>, otherwise: Box<Expr> },
    Call { name: String, args: Vec<Expr> },
}
```

This mirrors `TAntExprKind`. Ownership makes Delphi §3.3's clone-on-assign
rule automatic: an `Expr` cannot be aliased, and copying one is `.clone()`.
The V8 defects Delphi §12.1 lists (no clone, leaks, unary operand in the
wrong field, an argument list on every leaf) cannot be written in this shape.

### 3.5 Documents and modules

`Document { modules: Vec<Module>, main: usize }` from the start, even though
it holds one module until modules are built (seam 1, Delphi §2.2.2); `main`
is the index of the `*` module, or else the first. `Document::default()` is
one unnamed main module, which is what parsing starts from.
`flatten(&Document) -> Module` exists from A4 and
returns a copy of the main module (seam 2); the SBML writer is written against
its result, never against the parse tree.

### 3.6 Determinism

Anything ordered is a `Vec`, kept in source order. `HashMap` is used only for
lookup and is never iterated to produce output, so the same input always
gives the same model, the same diagnostics in the same order, and the same
text.

## 4. Source positions

```rust
pub struct Pos {
    pub offset: usize,   // byte offset into the source, 0-based
    pub line: u32,       // 1-based
    pub column: u32,     // 1-based, counted in characters (Unicode scalar values)
}

pub struct Span {
    pub start: Pos,
    pub end: Pos,        // exclusive: one past the last character
}
```

Every token, expression node and model item carries a `Span` (Delphi G4). A
span's length in bytes is `end.offset - start.offset`, and
`&source[span.start.offset..span.end.offset]` is its text. `Span::default()`
(line 0) means "no position", for items built in code rather than parsed.

**Line and column are computed from the byte offset**, through a table of
line starts built once per source, rather than counted as the lexer moves.
Delphi recorded that only `ConsumeNewline` tracked lines and that any scanner
which crossed a newline by hand shifted every later diagnostic (its
non-negotiables). Deriving positions from the offset makes that defect
impossible to write. `\n`, `\r\n` and a lone `\r` each end a line. A leading
byte-order mark is skipped.

Departure: Delphi counts columns in UTF-16 code units; this crate counts
characters. They differ only for characters outside the Basic Multilingual
Plane, and the fixtures compare lines only.

## 5. Diagnostics

```rust
pub struct Diagnostic {
    pub severity: Severity,   // Hint, Warning, Error
    pub code: Code,
    pub message: String,
    pub span: Span,
}
```

- **`Code` keeps every Delphi `TAntDiagCode`**, without the `dc` prefix
  (`Code::TypeConflict`). `code.name()` gives the Delphi spelling,
  `"dcTypeConflict"`, so the `.expected` fixtures are shared unchanged and a
  host can match on the same names as before. Each code has the Delphi
  default severity.
- **Messages use the Delphi wording** (`uAntStrings`), so the two
  implementations can be compared on diagnostics as well as models.
- **Display:** `diagnostic.with_file("model.ant")` gives Delphi's one-line
  form, `model.ant(3:17): error: message [dcSyntaxError]`; `Display` gives the
  same without the file. The diagnostic does not store the file name: the
  parser takes a `&str`, and only the caller knows where it came from.
- **The list is `Diagnostics`**, a thin wrapper that reads as a
  `[Diagnostic]` slice (indexing, `iter`, `len`) and adds Delphi's
  `TAntDiagnostics` helpers: `add` (default severity), `add_with` (explicit
  severity, for strict mode), `count_of`, `has_errors`, `has_any`,
  `first_of`, `sort_by_span` (stable, by line then column), `summary`, and a
  `Display` of one line per diagnostic.
- Diagnostics are returned **in source order** (Delphi `SortBySpan`).
- The message texts are functions and constants in `messages`
  (`messages::negative_stoichiometry("A")`), one per Delphi resource string.

## 5.1 The lexer (built in A1)

- `Lexer::new(&str)` reads on demand: `next_token()`, `peek(n)` (any depth),
  `try_consume(kind) -> Option<Token>`, `try_consume_keyword`, `current()`,
  `tokens()` (all of them, end of file included), and `tokenize(&str)` for
  the whole source at once.
- A `Token` has `kind`, `text`, `keyword: Option<Keyword>`, `number: f64`,
  `span` and `trivia` (the comments before it, one per line, joined with
  `\n`). There are no reserved words: an identifier notes the keyword it
  spells and the parser decides by position.
- The lexer owns its `Diagnostics`; the parser takes them with
  `take_diagnostics()`. In Delphi the lexer and parser share one list.
- Behaviour is the Delphi lexer's in every detail, including its edge cases
  (a `\` before a line break inside a string escapes the break). The one
  exception is the `layout` block (R37). The corpus check is exact:
  the same 273,129 tokens over the 1788 files, and the same 44,001 in the
  5000-line benchmark.
- Scanning is by byte: every character the language gives a meaning to is
  ASCII, and UTF-8 never uses an ASCII byte inside a longer character, so
  slicing at those bytes is always at a character boundary. A non-ASCII
  character outside a string or comment is one `Unknown` token (Delphi made
  two for a character outside the BMP, one per UTF-16 unit).

## 5.2 Numbers (built in A1)

`number::parse` and `number::format` replace Delphi's `AntTryStrToFloat` and
`AntFloatToStr`. `format` writes the shortest text that reads back as the
identical value, in exponent form outside 1e-5 to 1e15 (`6.02214076e23`),
with `inf`, `-inf` and `nan` for the non-finite values: the same thresholds
as `sbml-io`. Delphi wrote `6.02214076E23`; both read back the same.

## 5.3 Expressions (built in A2)

- **The tree** is `Expr { kind: ExprKind, span }` as sketched in §3.4, with
  `UnaryOp`, `BinaryOp` and `Constant` enums. A number's unit is
  `Option<String>`. There is no `Binding` back-pointer (R3); a symbol is its
  name.
- **Methods** (Delphi in brackets): `to_string()` via `Display`, with the
  fewest brackets the meaning needs (`ToString`); `same_tree` (`EqualsTree`:
  ignores spans, `nan` equals `nan`; derived `==` compares spans too);
  `evaluate()` and `evaluate_with(&resolver)` (`TryEvaluate`, returning
  `Option<f64>`); `is_constant`; `symbols()` (`CollectSymbols`, first
  appearance order); `node_count`, `depth`, `children`. `Clone` is derived.
- **Evaluation** follows Delphi call by call: lazy `piecewise`,
  short-circuit `&&`/`||`, the n-ary MathML operators, a bare `log` as base
  10, `None` for a zero divisor or reciprocal. A domain error the Delphi code
  did not check for (`sqrt(-1)`, `ln(0)`) gives IEEE's NaN or infinity, as
  Rust's functions do.
- **The parser** is `ExprParser::new(&mut lexer)` with `parse()`,
  `at_expression_start()` and the `parse_unit_suffixes` flag (off inside a
  participant list), plus `parse_expression(&str)` for a whole string and
  `can_start_expression(kind)`. It adds its problems to the lexer's list.
  It accepts exactly the Delphi grammar and builds the same trees, but it is
  written as precedence climbing (one loop over the precedence table) rather
  than Delphi's one function per level: a bracket costs about 4 stack frames
  instead of 13.
- **Nesting is limited, where Delphi had no limit.** A recursive descent
  parser and every recursive walk of a tree (writing, evaluating, cloning,
  comparing, dropping) use stack in proportion to depth, and running out
  aborts the program, which §6 does not allow. So the parser rejects, with
  one `dcSyntaxError`, an expression nested more than `MAX_NESTING` = 64
  levels (brackets, calls, unary operators, `^`, `?:`) or a tree more than
  `MAX_DEPTH` = 800 deep (a sum of 801 terms). Measured in A2: the deepest
  tree in 689 BioModels is 192 and the deepest brackets 15; in the corpus,
  100 and 4. The limits are about four times that, and a test walks the
  deepest accepted trees on a 1 MB stack (a Windows main thread,
  WebAssembly) in an unoptimised build.
- **What is written is never nested deeper than what was read**, or a valid
  expression could be written in a form that cannot be read back. A unary
  operator applied to a unary operator is written with a space, `- -a`,
  where Delphi wrote a bracket, `-(-a)`; both keep `--` from being read as
  the DNA operator.
- Checked beyond the Delphi tests: all 96,019 expressions in the 689
  BioModels parse, read back as the same tree, and write back identically.

## 5.4 The statement parser (built in A3)

- `parser::parse_document(&str) -> (Document, Diagnostics)`. It never fails
  and never panics; problems are in the diagnostics in the order found. The
  facade of §6 (`antimony::parse`, with options and the semantic pass) wraps
  it from A4.
- **The whole Delphi parser is ported in A3**, not only the core statements
  the plan gave A3: events, functions, unit definitions, constraints, the
  fbc objective and `sboTerm`, interactions and annotations (Delphi M5 and
  M9) are all in the one unit, and splitting it would have meant writing
  half a dispatcher and returning to it. A5 is left with porting the M5
  tests and checking the 1629-case set.
- It follows Delphi statement by statement, including its recovery
  (resynchronise at the next terminator, a reaction at the next line, a
  block at its `end`), its lookahead limits (4096 tokens), its generated
  reaction ids (`_J` + position, skipping a taken id; reactions are compared
  ignoring ASCII case, as Delphi's `SameText` did), and its quirks. One
  example of a quirk kept: `unit x = ;` records the definition `;`.
- **A5** ported the Delphi M5 tests (`tests/m5.rs`: events, functions,
  algebraic rules, constraints, annotations, and the 1629-case M5 set). All
  38 passed against the A3 parser and the A4 semantic pass unchanged: no
  gap was found. Deliberately breaking the `t0` modifier and a range's
  strictness each failed a test.
- The parser keeps a `HashMap` from symbol name to index for the module it
  is filling, so declaring a symbol is not a linear search (SPEC §3.1: the
  model itself holds only `Vec`s).
- The `comp` package's syntax (`sub1: moddef()`, `A.x is B.y`,
  `delete …`, `import`) is reported as `dcNotYetSupported` naming Delphi's
  milestone M10, the Delphi wording; some of it (`sub1.p80 * conv`) is a
  syntax error, as in Delphi. Modules are not in this plan.
- Checked: the 836 core cases parse with no errors except the 5 rejected by
  design; the 10 challenge cases parse; all 1657 parseable non-`comp` cases
  parse (1663 less the 5 and case 01622, an fbc case with a negative
  stoichiometry that the index does not flag and Delphi also rejects); all
  689 BioModels parse with no errors; 10 parser-level invalid fixtures give
  exactly their expected codes and lines; and 1,800 mangled corpus files
  (truncated, characters deleted, junk inserted) give diagnostics, never a
  panic.

## 5.5 The semantic pass (built in A4)

- `semantics::analyse(&mut Document, &Options) -> Diagnostics` resolves and
  checks every module, in Delphi's order: resolve all modules, then
  validate all. `antimony::parse_with` runs it after the parser.
- **Resolution** follows `uAntSemantics` step by step and in the same
  order, so implicit parameters are created in the same order as Delphi
  creates them: symbol values and rules, then each reaction (participant by
  participant, then its law), events, function bodies (with the function's
  arguments in scope), algebraic rules, constraints and objectives, symbol
  values again, then `in` and `conversionFactor` attachments, and last,
  anything still unclassified becomes a parameter.
- **References are names** (R3), so Delphi's `Binding`, `Compartment` and
  `Species` pointers have no counterpart. `semantics::resolve(module, name)`
  answers "what does this name refer to", with Delphi's rules: a built-in
  function or csymbol refers to nothing, nor does a module's formal argument.
- **Validation** is Delphi §6.3's list as Delphi built it: const with a
  rule, a species in a non-compartment, a conversion factor that is not a
  parameter, §6.5's uninitialised and implicit symbols, circular assignment
  rules (one report per module, naming the path), recursive functions and
  arity. Delphi's checker recursed through the rule graph; this one is a loop
  with an explicit stack, so a chain of 20,000 rules cannot exhaust the
  stack (tested on a 1 MB thread). The order and wording of its report are
  the same.
- `flatten::flatten(&Document) -> Module` and `flatten_module(&Document,
  name)` are a copy, seam 2 of Delphi §2.2.2.
- Checked: the 12 invalid fixtures each give exactly their expected codes
  and lines and no other error; the core and whole non-`comp` corpus give no
  `dcImplicitParameter` and leave nothing unclassified; all 689 BioModels
  load with no errors. And **every diagnostic matches the Delphi
  converter's**, by file, line and code: all 818 over the 1788-case corpus
  (`comp` cases included) and all 26 over the BioModels.
- One Delphi test (`ParameterPromotesToCompartment`) never reached the
  promotion it names: `v = 2` leaves `v` unclassified until the end of the
  pass. A test with `formula v` does reach it; a deliberate breakage of that
  promotion showed the gap.

## 5.6 The Antimony writer (built in A6)

- `write_document(&Document) -> String` and `write_module(&Module)`, and
  `write_document_with` / `write_module_with` taking `WriteOptions` (indent,
  explicit declarations, blank lines between groups; Delphi's
  `TAntWriterOptions`). `write::quoted` is the string escaper.
- It follows `uAntWriter` group by group and line by line: declarations,
  functions, reactions, interactions, events, algebraic rules, objectives,
  constraints, values and rules, `const`/`var`, units, attachments,
  annotations, raw sections; the module's display name after `end`. Over
  the corpus its text is **identical to the Delphi writer's** for every file
  but in how numbers are spelled (`1e-6` where Delphi writes `1E-6`, §5.2):
  all 139 files that differ do so only there. Over the 689 BioModels the
  only other difference is 46 `X notes "";` lines, which Delphi dropped and
  this keeps (R24).
- **Round-trip closure holds**, checked with a port of Delphi's
  `uAntModelCompare` (`tests/common/compare.rs`) that also compares what
  Delphi's did not: a reaction's generated name, the module's SBO term,
  modifiers, formal arguments and the raw sections. Spans and symbol order
  are ignored (the writer groups statements, so the first write reorders a
  file). Every parseable non-`comp` corpus case (1657) reads back as the same
  model, both as parsed and after the semantic pass, and so do all 689
  BioModels; the second write of each equals the third.
- Comments come back on the item they were written above (tested). Delphi's
  writer drops the comment of a symbol it writes no statement for (a
  parameter with no value); so does this.
- Lines end with `\n` on every platform (R29).

## 5.7 SBML output (built in A7, `sbml` feature)

- `sbml::write_sbml(&Document) -> String` writes the main module, flattened,
  as SBML Level 3 Version 2 (Delphi decision 11: L3V2 only).
  `sbml::to_sbml(&Module) -> sbml_io::Document` gives the same as data, for
  a caller that wants to change it or validate it first; `sbml::to_math`
  converts one expression to MathML; `xml_escape` and `xml_escape_text` are
  the escapers; `DEFAULT_COMPARTMENT` is the id of the compartment made for
  species that name none. `sbml::information_loss(&Module) -> Diagnostics`
  names, as `InformationLoss` warnings, the annotations the SBML does not
  carry (below), so none is dropped silently.
- It follows `uAntSBMLWriter` element by element, but builds an
  `sbml_io::Document` instead of writing text (R10): `sbml-io` writes the
  XML, so the attribute order is `sbml-io`'s (`metaid` first) and lines end
  with `\n`. The mapping, in brief:
  - a literal value (a number, or a signed one) is the element's attribute;
    anything else, `3/2` included, is an initial assignment, so the author's
    expression is kept;
  - a species with no compartment goes in `default_compartment` (made with
    SBO:0000410, unless the model declares one); `substanceOnly` makes the
    value an amount;
  - rules are written assignment, then rate, then algebraic; a rule with no
    expression (`p' = ;`) is left out;
  - a constraint naming a reaction is a flux bound, which belongs to `fbc`,
    and is left out;
  - a named stoichiometry is the species reference's id, with the named
    symbol's value, constant unless a rule or event changes it;
  - modifiers are every species the kinetic law uses that is not a reactant
    or product (rule 21121), in order of name (R32), and an interaction is a
    modifier with its SBO term and label, replacing the bare one;
  - an unnamed event is given the id `_E1`, `_E2`, … by position, and its
    three required flags are always written, `true` unless the source says
    otherwise;
  - `time` and `avogadro` are csymbols, `delay` and `rateOf` csymbol calls;
    `sqrt` and `log10` are a root of degree 2 and a log of base 10; `root`
    and `log` with two arguments take the first as the degree or base (R31);
  - a unit definition is built from base units, `^` exponents and a leading
    multiplier; one it cannot represent is written empty, not wrong;
  - notes that are already XHTML go out as markup in a `<body>`, prose is
    escaped in a `<p>`; CV terms, creators and dates go in one RDF
    annotation, terms sharing a qualifier in one `rdf:Bag`, `bqmodel` or
    `bqbiol` by subject; a `meta_<id>` metaid is made only where there is
    something to anchor;
  - SBO terms go on the model and on every element a symbol becomes (an
    attribute, needing no metaid); a function definition takes the
    annotation, notes and SBO term of the symbol named after it. Delphi
    wrote no SBO terms.
- Not written, as in Delphi: the `fbc` package (objectives and flux bounds
  are read and kept, not emitted) and the `comp` package (flattening is
  still a copy of the main module). Annotations not written, each named by
  `information_loss`: an interaction's (it becomes a modifier, which has
  none of its own), and one in a vocabulary this library does not know
  (kept for Antimony only).
- A compartment or parameter changed by a rule or an event is
  `constant="false"` even when not declared `var` (R30).
- **Checked:** every parseable non-`comp` corpus case (1657) and all 689
  BioModels give SBML that passes `sbml_io::validate`, reads back with
  `sbml_io::read`, and has no errors in libSBML 5.20.5's full consistency
  check (python-libsbml; the Delphi output has none either). The
  differential tier finds our SBML **equal to Delphi's, element by
  element**, on all 1657 (§8). Over the 708 BioModels (A8) the only differences
  are R31 (3 models, `root(2, …)`) and R32 (43 models, modifier order).

## 6. Entry points and the error contract

Built: reading (A4), the Antimony writer (A6), SBML output (A7), and the
facade (A8). The crate-level documentation in `src/lib.rs` is the user's
guide (Delphi's `doc/api.md`, R35); this section is the design.

```rust
let loaded = antimony::parse(source);                      // never fails, never panics
let loaded = antimony::parse_with(source, &Options { strict: true, ..Default::default() });
let loaded = antimony::parse_bytes(&bytes, &options);      // a file's bytes; R33
loaded.document          // the model, as far as it could be read
loaded.diagnostics       // Diagnostics (§5): the parser's, then the semantic pass's
loaded.has_errors()
let doc = loaded.into_result()?;   // Result<Document, Diagnostics> for callers who want it
let flat = antimony::flatten::flatten(&doc);               // seam 2: a copy of the main module

antimony::write_document(&doc) -> String                  // A6; also write_module, WriteOptions
antimony::sbml::write_sbml(&doc) -> String                // A7, `sbml` feature; also to_sbml(&module)
antimony::sbml::write_sbml_module(&module) -> String       // a module that is not the main one

let m = doc.main_module();                                  // the queries (A8), on Module:
m.species() / floating_species() / boundary_species() / compartments() / parameters()
m.value(id) -> Option<f64>                                  // raw: None when nothing to compute
m.concentration(id) / amount(id) / compartment_size(id) / parameter_value(id) -> Option<f64>
m.stoichiometry(&participant) -> f64                        // 1 when silent
m.evaluate(&expr) -> Option<f64>                            // through this module's values
m.symbol(id) / reaction(id) / event(id) / function(id) / unit(id) -> Option<&T>
symbol.assignment_rule() / rate_rule() -> Option<&Expr>;  participant.stoichiometry_text()
```

- **A malformed file is not an error for the caller to unwrap.** It is a
  document plus error diagnostics (Delphi §7.1). `into_result()` exists for
  callers who want Delphi's stricter facade rule — no model after a failed
  load — as an ordinary `Result`.
- **No input may panic.** Every path through the lexer and parser recovers;
  a panic on any input is a bug. A fuzz-style test over mangled corpus files
  is added in A3.
- **`Options`** (a plain struct with `Default`, `semantics::Options`
  re-exported at the root): `strict` reports Delphi §6.5's two warnings
  (`dcUninitialisedSymbol`, `dcImplicitParameter`) as errors, as Delphi's
  `StrictMode` does; other warnings, such as `dcInformationLoss`, stay
  warnings. `initialise_missing` is Delphi §6.5.1's opt-in fill-in
  (`--init-zero`): off by default, still reported, zero except for a
  compartment, which gets 1.
- **Delphi's facade, re-shaped.** `TAntimony` is an object that is loaded
  and then queried; here `parse` returns the model and the queries are
  methods on `Module` (`src/model/values.rs`). Its three rules become:
  1. *A caller writes no tests* — kept. The value methods answer with each
     kind's default where the model is silent: a species 0, a compartment's
     size 1 (0 would divide every concentration by zero), a parameter 0
     (the semantic pass warned), a stoichiometry 1. `value` is the raw
     answer, and the fields keep their `Option`, for a caller who needs to
     know a value was absent.
  2. *A bad argument raises* — becomes `None`. An id that is not of the kind
     asked about (no such symbol, or a species asked for a parameter's
     value) answers `None`; a lookup by index is the `Vec`'s own `get`, or
     an index that panics like any other (R8).
  3. *A query before a load raises* — cannot arise: there is no model before
     `parse` returns one. A failed load is `Err` from `into_result`.
- **Values are worked out as Delphi's facade did**: a symbol's assignment
  rule if it has one (it holds at the start too), else its initial value,
  evaluated through the values of the symbols it names, ten deep at most
  (`a = b; b = a` is not a circular rule, so the semantic pass allows it).
  `time`, a function call or a symbol with no value makes it `None`. A
  species stores one number, an amount if `substanceOnly`, else a
  concentration; the other goes through its compartment's size (1 if it
  names none), and a size of 0 gives a concentration of 0, not infinity.
- **No index-based accessor family** (R9). `GetSpeciesCount` /
  `GetSpeciesId(i)` become the iterators; `GetReactantStoichiometry(r, i)`
  is `m.stoichiometry(&m.reactions[r].reactants[i])`; `GetKineticLaw` is
  the law's `to_string()`; `GetModelId` is `main_module().name`. Ids never
  carry the `$` of a boundary species. An unnamed event has no id to look
  up (Delphi's `GetEventIds` gave `''`); it is found in `events` by
  position.
- **Not ported:** `LoadFile` (R33) and `OnResolveImport`, which Delphi
  declared but never called; it belongs with `import`, when modules are
  built.
- **Tests:** `tests/facade.rs` ports `uAntFacadeTests` (27 of 30; the three
  about exceptions, reuse and missing files have nothing to test here and
  are listed at its top), plus six of its own.

## 7. Layout of the code

| Module | Ports | Notes |
|---|---|---|
| `span` | `TSourcePos`, `TSourceSpan` | §4; also `LineIndex` |
| `diagnostics` | `uAntDiagnostics` | §5 |
| `messages` | `uAntStrings` | §5 |
| `number` | `uAntTypes` (number conversion) | §5.2 |
| `builtins` | `uAntTypes` (built-in function lists) | |
| `lexer` | `uAntLexer`, `uAntTypes` (tokens, keywords, character classes) | public: an editor can use the tokens for highlighting |
| `expr` | `uAntExprNode`, `uAntExprParser` | public |
| `model` | `uAntTypes` (model enums), `uAntModel` | |
| `parser` | `uAntParser` | §5.4 |
| `semantics` | `uAntSemantics` | §5.5 |
| `flatten` | `uAntFlatten` | a copy, as in Delphi (§5.5) |
| `write` | `uAntWriter` | §5.6 |
| `sbml` | `uAntSBMLWriter`'s mapping | §5.7; `sbml` feature; built on `sbml-io` |
| `model::values` | `uAntFacade`'s query surface | §6; methods on `Module`, `Symbol`, `Participant` |
| crate root | `uAntFacade` (loading) | re-shaped (§6); the crate documentation is the guide (R35) |

Not ported: `uAntCompat` (a migration shim for Delphi callers; Delphi
decisions 23–25 do not apply), `Oracle/` (replaced by the differential tier,
§8), and the FPC concerns (Delphi decision 4).

## 8. Tests

- **Ported tests keep their Delphi names**, in snake case, in one integration
  test file per Delphi fixture (`tests/lexer.rs` for `uAntLexerTests`, …), so
  each Rust test can be traced to the test it ports. A Delphi test that does
  not apply in Rust is listed, with the reason, at the top of its file.
- **Before trusting a new test, make it fail once** (a deliberately broken
  expectation or input), then restore it.
- **`tests/invalid/`** holds the Delphi `Tests/corpus/invalid/` set: 12
  `.ant` files, each with a `.expected` file listing `<code> <line>` pairs
  (`#` lines are comments). These are the author's own files, copied in.
  The line may be one past the last line: that is end of file, where an
  unterminated block is reported.
- **External data** is found through environment variables. A test whose
  variable is unset prints a note and passes, so a fresh clone passes
  `cargo test`. Nothing from these locations is copied into the crate; the
  SBML Test Suite in particular is LGPL.

  | Variable | Points at | Used from |
  |---|---|---|
  | `ANTIMONY_CORPUS` | the folder holding `NNNNN/NNNNN-antimony.ant` (1788 cases) | A0 |
  | `ANTIMONY_CORPUS_INDEX` | the corpus `index.json` (milestone sets, `lossy`, `challenge`, `expect_reject`) | A3 |
  | `ANTIMONY_BIOMODELS` | the curated BioModels as Antimony, `<id>.ant`, made by `tools/biomodels_to_antimony.py` | A8 |
  | `ANTIMONY_ORACLE` | the Delphi `ant2sbml.exe` | A3 |
  | `ANTIMONY_SBML_OUT` | a folder to save each corpus case's SBML in, for checking with libSBML | A7 |

- **The corpus index** is read with `serde_json`, a dev-dependency that
  never reaches a user (`tests/common/mod.rs`). Its milestone sets are
  checked against Delphi's counts: 836 core, 1629 M5, 125 `comp`, 10
  challenge, 5 reject, 9 lossy.
- **The differential tier** (built in A3) is
  `cargo run --release -p antimony-rs --features sbml --example differential [case …]`,
  Windows only in practice. It runs `ANTIMONY_ORACLE` (`ant2sbml.exe`) on
  every non-`comp` case, caching its SBML under `target/antimony-oracle/`,
  reads that with `sbml-io` (also a dev-dependency) and compares it
  structurally with this crate's result for the same file, never as text.
  Until the semantic pass exists it compares only what the parser alone
  decides: each reaction's id, reversibility, participants, stoichiometry
  and whether it has a law; which symbols have assignment and rate rules;
  event ids and targets; function ids; the number of algebraic rules. Each
  later milestone widens it, and from A7 it compares SBML with SBML,
  classifying differences as Delphi decision 22 does (only a divergence
  fails). A3 result: 1657 cases agree, none diverge, 1 rejected by the
  oracle (01622), 7358 facts compared; with the arrows deliberately
  swapped, 399 cases diverge. **From A4** it covers every one of the 1788
  cases: each case's diagnostics are compared by line and code with the
  oracle's (from its stderr, also cached), and where the oracle wrote SBML
  outside `comp`, the model comparison adds which ids are species,
  compartments and parameters, with each species' compartment and boundary
  flag. The oracle's `default_compartment`, which its SBML writer adds for
  species without one, is the one difference it allows. A4 result: all 1788
  agree, 818 diagnostics and 17,943 facts over 1657 models compared; with
  the uninitialised-symbol warning deliberately suppressed, 188 cases
  diverge. **From A7** it compares SBML with SBML: this crate's
  `write_sbml` and the oracle's output, both read with `sbml-io`, element
  by element (the model's attributes and annotation, then each function,
  unit, compartment, species, parameter, initial assignment, rule,
  constraint, reaction and event), every attribute and every MathML node.
  The one allowance is whitespace in raw notes and annotations (Delphi
  indents them and ends their lines with CRLF; 50 cases differ only in
  that). Modifier order needs no allowance: the corpus agrees exactly. A7
  result: all 1788 agree, 818 diagnostics and 25,051 elements over 1657
  models compared; with `^` deliberately written as `*`, 81 cases diverge,
  with compartments written 2-dimensional 1338, with a creator's given
  name changed 11.
- **The SBML tests** (`tests/sbml.rs`, `sbml` feature) check that every
  document they write passes `sbml_io::validate` and reads back, and the
  corpus test does the same for all 1657 cases. libSBML is run by hand on
  the files `ANTIMONY_SBML_OUT` collects; it is not a test, since it needs
  Python.
- **The BioModels tier** (built in A8; Delphi `Tools/biomodels_tier.py`)
  is `tests/biomodels.rs`. The corpus is synthetic — its Antimony was written
  by libantimony from SBML made to test one feature at a time — so real
  published models are where prose notes, author lists, CV terms and odd
  vocabulary turn up. `tools/biomodels_to_antimony.py` converts the curated
  `<id>/*url.xml` with libantimony: 708 of the 717 convert, as for Delphi
  (the 9 others libantimony cannot read, or writes as invalid UTF-8). The
  test checks every one loads with no errors, round-trips (deep comparison;
  the second write equals the third), and, with `sbml`, gives SBML that
  passes `sbml_io::validate`, reads back, and is the same after a trip
  through this crate's Antimony writer. A8 result: all 708 pass; libSBML
  5.20.5 finds no error in their SBML (`tools/libsbml_check.py` on the files
  `ANTIMONY_SBML_OUT` collects); their diagnostics (28) match the Delphi
  oracle's by line and code; with the writer's CV terms deliberately left
  out, 562 fail. The Delphi tier's one "unsupported" model loads in both
  implementations now.
- **CI**: the workspace CI runs `cargo fmt --check`, clippy (native and
  `wasm32-unknown-unknown`, warnings as errors) and `cargo test` over every
  member, which covers this crate's core. Since A7 it also runs clippy
  (native and wasm32) and the tests for this crate with `--all-features`,
  and checks that `cargo tree` finds no dependency of the core. When the crate moves to its own
  repository it gets its own workflow doing the same.

## 9. Departures from the Delphi design

Every difference from the Delphi version is recorded here, with the reason.
Language behaviour is not changed by the port; these are about code shape,
licence and API.

| # | Delphi | Rust | Why |
|---|---|---|---|
| R1 | Licence BSD (decision 5) | MIT OR Apache-2.0 | The author's choice; the Rust ecosystem's usual pair |
| R2 | Class hierarchy, owned object lists | Plain structs, `Vec`, `Option`, `enum` (§3) | Rust idiom; ownership makes clone-on-assign automatic |
| R3 | Items bound to each other by object reference | References by name, looked up when needed (§3.1) | Survives editing and promotion; nothing to keep in step |
| R4 | Sentinels (−1, empty string) | `Option` (§3.2) | The type says what may be absent |
| R5 | Line counting in the lexer | Line and column computed from the byte offset (§4) | Makes the Delphi line-drift defect impossible |
| R6 | Columns in UTF-16 units | Columns in characters (§4) | Rust strings are UTF-8; differs only outside the BMP |
| R7 | Diagnostic stores a file name | File name supplied when displaying (§5) | The parser takes `&str` and never sees a file |
| R8 | Facade raises on bad ids and on queries after a failed load | `Option` for lookups; `Loaded::into_result()` for the strict form (§6) | Rust idiom; no exceptions |
| R9 | Index-based query API (`GetXxxCount`/`GetXxxId(i)`) | `Vec`s and iterator methods (§6) | Unnecessary with public `Vec`s |
| R10 | Native SBML XML writer (decision 1) | Built on `sbml-io`, following the same mapping | `sbml-io` already writes and validates SBML; one XML writer to maintain |
| R11 | `uAntCompat`, `Oracle/`, FPC (decisions 4, 23–25) | Not ported (§7) | Delphi-specific |
| R12 | Invariant `TFormatSettings` for every number conversion | `str::parse::<f64>` and Rust's own formatting | Both are locale-independent by definition, so the Delphi hazard cannot arise; never call C library number functions |
| R13 | `TStringList` must be made case-sensitive | `String` comparison is case-sensitive | The Delphi defect cannot arise; the BIOMD0000000557 `duspn`/`DUSPn` case is still tested |
| R14 | Comments in a token's trivia joined with the platform line break (CRLF on Windows) | Joined with `\n` | The same text on every platform |
| R15 | Lexer and parser share one diagnostics list | The lexer owns the list; parsers add to it through the lexer, and the owner takes it at the end (§5.1) | Ownership without sharing; still one list in the order found |
| R16 | A non-BMP character outside a string is two `ttUnknown` tokens | One `Unknown` token | Rust reads characters, not UTF-16 units (§5.1) |
| R17 | `AntVersionString` | `antimony::VERSION`, from `Cargo.toml` | One place for the version |
| R18 | No limit on expression nesting | 64 levels of nesting, a tree 800 deep (§5.3) | Running out of stack aborts; a diagnostic does not |
| R19 | `- -a` written as `-(-a)`, `!!a` as `!(!a)` | Written `- -a`, `! !a` (§5.3) | A bracket adds nesting, so the written form could exceed the limit |
| R20 | A bracketed expression's span stops before `)` | It includes the `)` | The span is the text of the expression |
| R21 | `--` after an operand is named as DNA syntax only after the first term of a sum (`a + b--c` was a generic syntax error) | Named wherever it follows an operand | A consequence of precedence climbing, and the better message |
| R22 | `TAntSymbolRule` with `srNone` and `srInitialAssign` | `rule: Option<Rule>`, `RuleKind` = `Assignment` or `Rate` (§3.3) | `srInitialAssign` was never set; `None` says "no rule" |
| R23 | A constraint holds `Expr` or `Lower`/`Subject`/`Upper`, all nullable | `ConstraintKind::Expr` or `ConstraintKind::Range` (§3.3) | Only the valid combinations can be written |
| R24 | An empty string meant "absent" for names, notes and display names | `Option<String>`, so `notes ""` keeps an empty note | R4; the empty string and no string are different facts |
| R25 | The plan split the parser between A3 and A5 | The whole parser is ported in A3 (§5.4) | It is one unit; A5 keeps the M5 tests and the 1629-case check |
| R26 | References bound by pointer (`Binding`) | `semantics::resolve(module, name)` (§5.5) | R3; the same rules for built-ins and formal arguments |
| R27 | Circular rules found by recursion | A loop with an explicit stack (§5.5) | A long chain of rules cannot exhaust the stack; same result and message |
| R28 | An unlabelled reaction is written with its generated label, `_J1: S1 -> S2`, and reading that back declares a symbol `_J1` the model did not have | `Reaction::name_generated` records that the name was made up, and the writer leaves the label out | Round-trip closure for unlabelled reactions; the author's file comes back as written. The author chose this on 2026-10-14 |
| R29 | The writer ends lines with the platform line break (CRLF on Windows) | `\n` | The same file on every platform, as R14 |
| R30 | A compartment's or parameter's SBML `constant` follows its declaration alone, so `x := 2*k` gives `constant="true"` with an assignment rule (libSBML 20903; an event assignment, 21213) | `constant="false"` when a rule or an event assignment changes it, as well as when declared `var` (§5.7) | Valid SBML. Neither the corpus nor BioModels has such a case, so the differential tier sees no difference |
| R31 | `root(n, x)` and `log(b, x)` written as `<root/>` or `<log/>` with two ordinary arguments | The first argument as the `<degree>` or `<logbase>` qualifier (§5.7) | That is MathML; `sbml-io` cannot read Delphi's form. BIOMD0000000056, 195 and 1044 use `root(2, …)` |
| R32 | Modifiers derived from a kinetic law sorted by a `TStringList`, in the locale's order (case ignored first: `Mpp` before `MpY`) | Sorted by character code (`MpY` before `Mpp`) | The same on every platform. Order means nothing in SBML; the corpus agrees exactly, and 43 BioModels differ in this, and in nothing else but R31 |
| R33 | `LoadFile` reads a file, reporting `dcFileNotFound` and `dcReadError` (with the decoder's message) | The crate opens no files; `parse_bytes` takes a file's bytes and reports invalid UTF-8 as `dcReadError` at the first bad byte | It runs in a browser, where there is no file system (§2); reading a file, and reporting a missing one, is the caller's part |
| R34 | The facade's queries raise on a bad index or id, with a message naming the method | `None`, or the `Vec`'s own indexing (§6) | R8: Rust's idiom; there is no method name to report |
| R35 | `doc/api.md`, a separate reference | The crate-level documentation (`src/lib.rs`) and the README, whose examples are run as doc tests | One place, shown by docs.rs, and checked by `cargo test` |
| R36 | `const`, `var` and `$` declare only: `const k1 = 0.5`, `var k = 3` and `$B = 4` are syntax errors | An item of `const` or `var` may have `in c` and `= value`, and a statement `$x = value` (or `$x in c …`) is `const $x = value`, as libantimony reads it. A `$` on something that is not a species is kept (`is_boundary`), so a later use as a species makes a boundary one, and the writer writes it back as `const $x`; it also writes `x in c;` for a symbol that is not a species or compartment, which it used to drop | libantimony accepts all three, and the language reference outranks the Delphi version (Delphi decision 19). Found when the crate's first user met `$B = 4` |
| R37 | `layout … end` and `layout { … }` are one token, kept verbatim (Delphi N1, §5.6) | `layout` is an ordinary word. Antimony 3's layout statements (`model.layout = on`, `S1.pos = {28, 35}`) are read into `Module.layout` and written back (§3) | Antimony has no `layout` block; its layout is dotted statements (the reference's "Layout and Render", libAntimony 3.0). The block rule swallowed `model.layout = on` and everything after it up to the next `end`, losing values and the model's own `end`. The author chose this on 2026-10-22 (docs/antimony-layout-design.md, D1) |

## 10. Left open, for the milestone that needs it

- Nothing at present.

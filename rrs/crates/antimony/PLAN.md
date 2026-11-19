# The `antimony` crate — porting plan

A pure-Rust reader and writer for the Antimony language, ported from the
user's complete Delphi implementation. It is a **reusable crate**: it lives in
the PathwayDesigner workspace for now and must be easy to move to its own
repository later. When it is done it replaces `crates/pd-antimony`, PathwayDesigner's
small Antimony subset.

Written 2026-10-14. This file is the starting point for the work; keep it up
to date as milestones finish, and record every departure from the Delphi
design in `SPEC.md` (created in A0).

---

## 1. The source to port

**`D:\Documents\Embarcadero\Studio\Projects\Antimony_Delphi_Version`** (a git
repository; Delphi 13 / BDS 37.0).

Read these first, in this order:

1. `CLAUDE.md` — current state, conventions, non-negotiables.
2. `SPECIFICATION.md` — the normative design (~1560 lines). In particular:
   §3 architecture, §4 object model, §5 language surface (and §5.2.1–5.2.2
   corrections), §6 semantics (two passes, type promotion, §6.5 uninitialised
   symbols), §8 SBML mapping, §9 testing, **§13 decisions record (25 settled
   decisions — do not re-litigate them)**.
3. `doc/grammar.ebnf` (the grammar), `doc/api.md` (the Delphi API and its
   "conventions that will bite you", §10), `doc/manual.md`.
4. The Antimony language reference,
   <https://tellurium.readthedocs.io/en/latest/antimony.html>, which outranks the
   Delphi spec on what the *language* is (decision 19).

What exists there:

| Delphi unit | Lines | Rust module | Port? |
|---|---|---|---|
| `uAntTypes` | 816 | `types` | yes |
| `uAntStrings` | 72 | (inline) | yes |
| `uAntDiagnostics` | 335 | `diagnostics` | yes |
| `uAntLexer` | 981 | `lexer` | yes |
| `uAntExprNode` | 1084 | `expr` | yes |
| `uAntExprParser` | 638 | `expr::parse` | yes |
| `uAntModel` | 1564 | `model` | yes |
| `uAntParser` | 2303 | `parser` | yes |
| `uAntSemantics` | 683 | `semantics` | yes |
| `uAntFlatten` | 183 | `flatten` | yes (identity for now, as in Delphi) |
| `uAntWriter` | 926 | `write` | yes |
| `uAntSBMLWriter` | 1568 | `sbml` (feature) | **re-done on top of `sbml-io`**, following its mapping |
| `uAntFacade` | 1333 | crate root API | re-shaped as idiomatic Rust (§3) |
| `uAntCompat` | 1566 | — | **no** (a Delphi migration shim) |
| `Oracle/` | — | — | **no**; the Rust differential tier replaces it (§5) |

Tests: `Tests/*.pas`, **491 DUnitX tests** (~9950 lines) in fixtures per unit:
lexer, expressions, types, parser, semantics, M5 (functions/events/units),
writer, SBML, facade, diagnostics, corpus. Port each fixture with its layer.

Status of the Delphi version: M0–M9 built and clean on every tier. Not built
(and not to be built in this port either, unless asked): M10 modules/`import`/
real flattening, `fbc` *output*, the unit algebra (unit definitions are kept as
raw text), interactions (§9.2 of its spec).

## 2. Goals and non-goals

**Goals**

- Parse Antimony into a typed model with source positions on everything, and
  diagnostics (with the Delphi codes) instead of panics.
- Write Antimony back: **round-trip closure** — anything we write, we read back
  to an identical model (Delphi §9.1.4 #1).
- Write SBML Level 3 Version 2 through `sbml-io` (optional `sbml` feature).
- Agree with the Delphi implementation on the whole corpus (§5).
- Replace `pd-antimony` in PathwayDesigner (milestone A9).

**Non-goals** — the same as the Delphi project's N1–N5: no layout/render syntax
(skip a `layout` block), no DNA strand syntax (a clear "not supported"), no SBML
→ Antimony, no simulation (only constant folding), no byte-identical formatting.

**Constraints**

- **Pure Rust, no dependencies** for the core (lexer → writer). The `sbml`
  feature adds `sbml-io` only. Must build for `wasm32-unknown-unknown`.
- **Movable.** Its `Cargo.toml` must not use `workspace = true` keys; its tests,
  fixtures and README live inside `crates/antimony/`; nothing in it may mention
  PathwayDesigner. Test data that can't be redistributed is found through
  environment variables (§5), never committed.
- **Deterministic.** `Vec` for anything ordered; maps only for lookup.

## 3. Rust design (settled in A0; `SPEC.md` is now the record)

- **Name:** folder `crates/antimony`. `antimony` is taken on crates.io (an
  unrelated 2017 crate), so the package is **`antimony-rs`** with library name
  `antimony` (`cargo test -p antimony-rs`; code says `use antimony::…`).
- **Licence:** **MIT OR Apache-2.0** — the user's choice, 2026-10-13.
- **Data:** plain structs with public fields, `Vec`, `Option`, enums — not a class
  hierarchy (the user disliked libSBML's). The unified symbol of Delphi §4.1
  becomes one `Symbol` struct with a `SymbolType` enum and promotion rules.
  Expressions are an owned `enum Expr` (clone instead of aliasing, as Delphi §3.3
  requires).
- **Spans:** byte offset + line + column + length on tokens, AST nodes and model
  items (G4).
- **Entry points** (sketch): `Document::parse(&str) -> Loaded { document,
  diagnostics }`; `write_antimony(&Module) -> String`; `to_sbml(&Module) ->
  sbml_io::Document` behind the `sbml` feature. Keep Delphi's error contract in
  spirit: a malformed file is not a panic or an `Err` for the caller to unwrap —
  it is a document plus error diagnostics; strict mode turns warnings into errors.
- **Diagnostics:** severity, code (keep the Delphi `dc…` names as a Rust enum),
  message, span. `InitialiseMissing` is an explicit, opt-in option (Delphi
  §6.5.1).

## 4. Milestones

One milestone at a time; the user checks each before it is committed. Each
milestone ports the matching Delphi tests first, then the code, and is done
when those tests pass, `cargo clippy` is clean and the crate builds for wasm.

| # | Milestone | Ports | Done when |
|---|---|---|---|
| A0 ✓ | Skeleton | — | Crate in the workspace (not yet used by pd-core), README, `SPEC.md` with the §3 decisions, licence files, CI runs its tests and a wasm build. **Done 2026-10-13** (the workspace CI covers it; no CI change was needed) |
| A1 ✓ | Types, diagnostics, lexer | `uAntTypes`, `uAntDiagnostics`, `uAntStrings`, `uAntLexer` | `uAntLexerTests`, `uAntTypesTests`, `uAntDiagnosticsTests` ported and passing. **Done 2026-10-13**: same 273,129 tokens as Delphi over the corpus |
| A2 ✓ | Expressions | `uAntExprNode`, `uAntExprParser` | `uAntExprTests` passing: precedence, functions, `ToString` formatting, evaluation, symbol collection. **Done 2026-10-13**; nesting limits added (SPEC §5.3) |
| A3 ✓ | Model and statement parser | `uAntModel`, `uAntParser` | `uAntParserTests` passing; every corpus case in the M3/M4 set (836, Delphi §9.1) parses with no errors. **Done 2026-10-13**: the whole parser was ported (SPEC R25); all 1657 parseable non-comp cases and 689 BioModels parse; differential tier built, 0 divergences |
| A4 ✓ | Semantics | `uAntSemantics`, `uAntFlatten` (identity) | `uAntSemanticsTests` passing; promotion table; §6.5 warnings; the `invalid/` set gives each expected code and line; zero `dcImplicitParameter` over the corpus. **Done 2026-10-14**: every diagnostic over the corpus and BioModels matches Delphi's; the differential tier now covers all 1788 cases |
| A5 ✓ | Functions, events, units, constraints, fbc model, annotations/CV terms | the parser already handles these (ported whole in A3, SPEC R25); A5 ports their tests and fills any gap they show | `uAntM5Tests` passing; the 1629-case set parses. **Done 2026-10-14**: all 38 tests pass with no change to the code; the 1629 cases load with no errors |
| A6 ✓ | Antimony writer | `uAntWriter` | `uAntWriterTests` passing; round-trip closure over all 1657 parseable non-`comp` corpus cases (deep model comparison, Delphi `uAntModelCompare`). **Done 2026-10-14**: closure holds over the corpus and all 689 BioModels; the text matches Delphi's writer except in number spelling; unlabelled reactions stay unlabelled (SPEC R28) |
| A7 ✓ | SBML output (`sbml` feature) | `uAntSBMLWriter`'s mapping, built on `sbml-io` | `uAntSBMLTests` passing; output validates with `sbml-io` and python-libsbml; the differential tier (§5) shows zero divergences over the 1629 SBML cases. **Done 2026-10-15**: all 68 tests pass; the 1657 corpus and 689 BioModels outputs are valid in `sbml-io` and libSBML; the differential tier compares SBML element by element and finds none of the 1657 different; three Delphi defects fixed (SPEC R30–R32) |
| A8 ✓ | Facade, docs, BioModels | `uAntFacade` (re-shaped), `doc/api.md` | README with examples; the BioModels tier (§5) runs clean. **Done 2026-10-15**: value queries on `Module` and `parse_bytes` (SPEC §6), the crate documentation as the guide and README examples, all run as doc tests; `uAntFacadeTests` ported; the BioModels tier passes on all 708 convertible models, with no libSBML errors and the same diagnostics as Delphi |
| A9 ✓ | Swap into PathwayDesigner | — | `pd-core::antimony` uses the new crate; PathwayDesigner's own Antimony tests pass unchanged; `pd-antimony` deleted; specs updated. **Done 2026-10-15**: the bridge reads the crate's model and writes S §11.3's text itself; every PathwayDesigner test passes, one input changed (a `model` with no `end`, which Antimony rejects; the author's choice); `$B = 4`, `const k = 1` and `var k = 1`, which the old reader took and Delphi rejected, are now read by the crate as libantimony reads them (SPEC R36) |

A9 keeps everything PathwayDesigner-specific in the bridge,
`crates/pd-core/src/antimony.rs` (the catalytic rewrite, null nodes, aliases,
warnings for what a diagram can't hold such as events and functions). Its
current input type is `pd_antimony::Description`; A9 makes it read the new
model instead. PathwayDesigner's behaviour must not change, apart from now
reading every valid Antimony file.

## 5. Test data and the differential tier

All of this lives outside the repo. Tests that need it read a path from an
environment variable and are **skipped with a note** when it is unset, so a
fresh clone still passes `cargo test`.

| Data | Where | Env var | Use |
|---|---|---|---|
| Corpus: 1788 `NNNNN-antimony.ant` files | `…\Antimony_Delphi_Version\SBML Test Suite\NNNNN\` | `ANTIMONY_CORPUS` | parse, round-trip, SBML tiers |
| Corpus index (milestone sets, `lossy`/`challenge` flags) | `…\Antimony_Delphi_Version\Tests\corpus\index.json` | `ANTIMONY_CORPUS_INDEX` | choose each milestone's acceptance set |
| Invalid inputs, 12 `.ant` + `.expected` | `…\Antimony_Delphi_Version\Tests\corpus\invalid\` | — | **copied into `crates/antimony/tests/invalid/`** in A0 (the user's own files; `M3.md` there is a session transcript and was not copied) |
| BioModels, 717 models | `…\Antimony_Delphi_Version\BioModels\<id>\`, converted to Antimony by `tools/biomodels_to_antimony.py` (708 convert) | `ANTIMONY_BIOMODELS` (the converted folder) | A8 tier, `tests/biomodels.rs` (Delphi `Tools/biomodels_tier.py`) |
| Reference implementation | `…\Antimony_Delphi_Version\CLI\Win64\Release\ant2sbml.exe` (already built) | `ANTIMONY_ORACLE` | differential tier |

The SBML Test Suite is LGPL and is **not** copied into the repo (the same rule
`sbml-io` follows).

**The differential tier** is the strongest check available: two
implementations that must agree. For each corpus case, run the Delphi
`ant2sbml.exe` (`--to-sbml`, and `--to-antimony` for the writer) and the Rust
crate on the same `.ant`, read both SBML files with `sbml-io`, and compare them
structurally (sorted elements, normalised MathML, numeric tolerance), never as
text. Classify differences as Delphi decision 22 does — equivalent, known gap,
ours better, divergence — and fail only on divergence. Delphi's
`Oracle/uAntSBMLCompare.pas` has the comparison rules to copy. Build it as a
`cargo run --example` (Windows only, since the oracle is a Windows exe) in A3,
not at the end: the Delphi spec's §9 says the early differential tier is what
found the defects nothing else saw.

If the Delphi code needs rebuilding: see the global CLAUDE.md (Delphi 13 is
BDS 37.0; source `rsvars.bat` before `msbuild`).

Things the corpus does not cover (Delphi §9.1.1) need hand-written fixtures:
`function … end` blocks, interaction arrows, `import`, `<->`, comments and
layout variety, and invalid input. Remember Delphi's finding that in
libantimony's output **`->` is reversible and `=>` irreversible**.

## 6. Working practices (from the PathwayDesigner sessions)

- One milestone per round: build it, run the tests, show the user, then commit
  and push only when the user says so.
- Record decisions and deviations in `crates/antimony/SPEC.md`; don't edit the
  Delphi spec.
- The user prefers plain, short explanations, one simple question at a time,
  and a recommendation rather than a menu.
- Before trusting a new test, make sure it can fail.
- On this machine, long Python edit scripts in bash heredocs break; write them
  to a scratchpad file. Remove files with `rm -f "${dir:?}"/*`.

## 7. Moving it out later

When the user wants its own repository: copy `crates/antimony/` (it is
self-contained), add a workspace-free CI, publish, and point PathwayDesigner's
`Cargo.toml` at the published version (or a git dependency). Nothing else in
PathwayDesigner should need to change, because only `pd-core::antimony` uses it.

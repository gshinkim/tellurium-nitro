# antimony

Read and write [Antimony](https://tellurium.readthedocs.io/en/latest/antimony.html),
the human-readable language for biochemical models, and convert it to
[SBML](https://sbml.org). Pure Rust with no dependencies, so it also works in
WebAssembly.

> **Status: complete, not yet published.** The package will be published as
> `antimony-rs` (the name `antimony` is taken on crates.io); the library is
> called `antimony`.

It is a port of a complete Antimony implementation written in Delphi, and is
checked against it case by case: over a corpus of 1788 models and 708 curated
BioModels it reports the same problems on the same lines and writes Antimony
that reads back as the same model, and over the corpus it writes the same
SBML element by element. Where it departs from the Delphi version, `SPEC.md`
says so and why.

## Using it

```toml
[dependencies]
antimony-rs = "0.1"                                    # Antimony only
antimony-rs = { version = "0.1", features = ["sbml"] } # and SBML output
```

Read a model, check it, and ask it questions:

```rust
let loaded = antimony::parse("
    compartment C = 4;
    species S1 in C, S2 in C;
    J1: S1 -> S2; k1*S1;
    S1 = 10; S2 = 0; k1 = 0.1;
");

// A malformed file is not a panic or an Err: you get the model and a list
// of diagnostics, each with a code, a message and a line and column.
for d in loaded.diagnostics.iter() {
    eprintln!("{d}"); // (3:17): warning: ... [dcUninitialisedSymbol]
}
if loaded.has_errors() {
    return;
}

let model = loaded.document.main_module();
assert_eq!(model.species().count(), 2);
assert_eq!(model.concentration("S1"), Some(10.0));
assert_eq!(model.amount("S1"), Some(40.0)); // through the size of C
let j1 = model.reaction("J1").unwrap();
assert!(j1.is_reversible()); // `->`; `=>` is irreversible
assert_eq!(j1.kinetic_law.as_ref().unwrap().to_string(), "k1 * S1");
```

Write it back as Antimony, or as SBML Level 3 Version 2:

```rust
let loaded = antimony::parse("J1: S1 -> S2; k1*S1; S1 = 10; S2 = 0; k1 = 0.1");
let text = antimony::write_document(&loaded.document);
let xml = antimony::sbml::write_sbml(&loaded.document); // the `sbml` feature
```

Where a model says nothing, the value methods answer with what is true of
that kind of quantity — a species 0, a compartment size 1, a stoichiometry 1 —
so a caller writes no test for a missing value; the fields themselves keep
`Option` for a caller who needs to know. The crate documentation is the full
guide: the error contract, the model, writing, and what is not built.

## What it does

- **Parses** Antimony into a typed model: species, parameters, compartments,
  reactions, interactions, rules, events, functions, units, constraints,
  objectives, annotations and notes. Every item carries its position in the
  source.
- **Reports problems as diagnostics**, never by panicking, so an editor can
  open a broken file and show every problem at once.
- **Writes Antimony back.** Anything it writes reads back as the same model.
- **Writes SBML Level 3 Version 2** (the `sbml` feature, using the `sbml-io`
  crate). The output passes libSBML's consistency checks over the whole
  corpus and the BioModels.

Not built: reading SBML, expanding submodules (`comp`) and `import`, writing
the `fbc` package to SBML, DNA strand syntax, unit algebra beyond base units,
and simulation. Antimony 3's layout statements (`model.layout = on`,
`S1.pos = {28, 35}`, `J0.S1.b1 = {55, 31}`) are read into
`Module::layout`, kept in order and written back; what they draw is the
caller's to decide.

## Design

`SPEC.md` records the design and every place it departs from the Delphi
version. In short: plain structs with public fields, `Vec` and `Option`, one
`Symbol` type whose kind is inferred as the file is read (as Antimony
requires), and expressions as an ordinary owned `enum`.

## Testing

```sh
cargo test
cargo test --features sbml
```

The tests that need the 1788-model test corpus, the BioModels, or the Delphi
reference implementation find them through environment variables and are
skipped with a note when those are not set (see `SPEC.md` §8):

| Variable | Points at |
|---|---|
| `ANTIMONY_CORPUS` | the folder holding `NNNNN/NNNNN-antimony.ant` |
| `ANTIMONY_CORPUS_INDEX` | the corpus `index.json` |
| `ANTIMONY_BIOMODELS` | the curated BioModels as Antimony, made by `tools/biomodels_to_antimony.py` |
| `ANTIMONY_ORACLE` | the reference converter, `ant2sbml.exe` |
| `ANTIMONY_SBML_OUT` | a folder to save the SBML in, for `tools/libsbml_check.py` |

With the corpus variables and the oracle set, the differential tier compares
this crate with the reference implementation over the whole corpus:

```sh
cargo run --release --features sbml --example differential
```

The two scripts in `tools/` need Python with the `antimony` and
`python-libsbml` packages (both come with Tellurium).

## Licence

MIT OR Apache-2.0, at your option. See `LICENSE-MIT` and `LICENSE-APACHE`.

# Tellurium Nitro

**Regular tellurium, with a new math engine that is much faster.** Your code and your results stay the same. You
just get them sooner: simulations run about 2 to 8 times faster, big models load about 100 times faster, and heavy
analyses (steady states, control coefficients, sensitivities) run hundreds to thousands of times faster.

- [The story](#the-story)
- [What is different from regular tellurium?](#what-is-different-from-regular-tellurium)
- [Install it (about 10 minutes)](#install-tellurium-nitro-step-by-step)
- [Use it](#use-it)
- [Use it without Python: the C library](#use-it-without-python-the-c-library-libroadrunners-c-api)
- [How much faster is it?](#how-much-faster-is-it)
- [How do we know the answers are right?](#how-do-we-know-the-answers-are-right)
- [For experts: how it works](#for-experts-how-it-works)
- [Every command, and what it does](#every-command-and-what-it-does)

## The story

[Tellurium](https://github.com/sys-bio/tellurium) is a Python tool that biologists use to simulate what happens
inside a cell: how the amounts of chemicals go up and down as reactions run. You describe a model ("S1 turns into S2
at this rate"), and tellurium tells you how it plays out over time, where it settles, and which reactions matter most.

The part of tellurium that does the actual math is called **roadrunner**. It's good, but it was built many years
ago, and it starts to hurt as models get bigger or questions get harder:

- **Loading a big model is slow.** A model with 1000 chemicals and 3400 reactions takes about 30 seconds to load
  before you can do anything with it.
- **Some questions take a long detour.** To find out which reaction controls a pathway, roadrunner nudges every
  number a little, re-solves the whole model, and repeats. On a 40-chemical model that is about 170,000 extra runs.
- **It uses one processor core**, even though most computers now have 8 to 20.

**Tellurium Nitro replaces that math engine with a new one**, written from scratch in a fast language called Rust.
It does the same jobs, but:

- it reads the model directly and turns it into fast machine code in a fraction of a second;
- it answers the hard questions with exact math instead of trial and error, so the answers are faster *and* more
  accurate;
- it uses every core of your computer when there is a lot of work to share.

Everything else is still tellurium. The only change to your code is the import line: `import tellurium_nitro as te`
instead of `import tellurium as te`. After that, `te.loada(...)`, `r.simulate(...)` and the rest work exactly as
before. libroadrunner and libsbml are not needed: the Rust engine does all the work.

## What is different from regular tellurium?

| | Regular tellurium | Tellurium Nitro |
|---|---|---|
| Package names | `tellurium`, `libroadrunner` | `tellurium_nitro`, `libroadrunner_nitro`, so both can be installed side by side |
| How you write code | `import tellurium as te` | `import tellurium_nitro as te`; every command after that is **the same** |
| Results | roadrunner's answers | **The same answers** (within the solver's normal error), and more exact for control coefficients |
| Math engine | roadrunner (C++) | New Rust engine (no libroadrunner, no libsbml) |
| Simulate a model over time | baseline | **about 2-8x faster** |
| Load a 1000-reaction model | about 30 s | **about 0.3 s** |
| Steady state, control coefficients, sensitivities | baseline | **600 to 130,000x faster** |
| Processor cores used | 1 | all of them, when it helps |
| Models roadrunner refuses (delays, algebraic rules, fast reactions) | error | **runs them** |
| SBML Test Suite (the standard correctness test) | — | **1923 of 1923 cases pass** |
| Use from C, C++, MATLAB, Julia, C#… | libroadrunner's C API (`libroadrunner_c_api`) | **the same C API** (`libroadrunner_nitro_c_api`, all 350 functions) |
| How you install it | `pip install tellurium` | from this GitHub page (see below) |

New things you can do that regular tellurium can't:

- **Run many simulations at once:** `r.simulateBatch(...)` runs one simulation per set of parameter values, on
  every core. Handy for parameter scans, fitting and "what if" questions.
- **Use it from plain roadrunner code:** if your scripts start with `import roadrunner`, change that one line to
  `import roadrunner_nitro as roadrunner` and the rest runs on the fast engine (`tellurium_nitro.nitro` still
  works for now, but is deprecated).
- **Flux balance analysis:** `te.fba(model)` solves flux balance problems (linear programs), which roadrunner
  doesn't do.

## Install Tellurium Nitro (step by step)

> **Its own names, so it can't clash with regular tellurium.** Nitro installs as the package `tellurium_nitro`
> (`import tellurium_nitro as te`), and its engine as `libroadrunner_nitro` (`import roadrunner_nitro`). Your
> existing `tellurium` and `libroadrunner` are never replaced or changed by name. Nitro needs neither libroadrunner
> nor libsbml; it shares a few packages with tellurium (antimony, libsedml, phrasedml), and pip may update their
> versions, so a fresh virtual environment (step 4) is still the safest way to keep your current setup as it is.

**1. Install Python 3.9 or newer, and git.** Check with `python3 --version` and `git --version`.

**2. Install Rust**, which compiles the fast engine. You only need to do this once.

- macOS / Linux: `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`, then open a new terminal
  (on macOS, `brew install rust` also works).
- Windows: download and run `rustup-init.exe` from https://rustup.rs. It will ask you to install the
  "Visual Studio C++ Build Tools"; say yes.

Check with `cargo --version`.

**3. Download Tellurium Nitro and its engine,** side by side in one folder:

```bash
git clone https://github.com/gshinkim/libroadrunner_nitro.git
git clone https://github.com/gshinkim/tellurium-nitro.git
cd tellurium-nitro
```

**4. Create and activate a fresh virtual environment:**

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
```

**5. Build and install the Rust engine** (`libroadrunner_nitro`). This takes 1-3 minutes the first time.

```bash
pip install --upgrade pip
pip install ../libroadrunner_nitro
```

**6. Install Tellurium Nitro and its dependencies:**

```bash
pip install .
```

**7. Check that it works:**

```bash
python -c "import tellurium_nitro as te; print('Nitro engine on:', te.turbo())"
```

It should print `Nitro engine on: True`. If the import fails, step 5 or 6 didn't finish; run it again and read
the error.

From now on, activate the environment (step 4, `source .venv/bin/activate`) whenever you want to use Nitro.

**8. (Optional) Build the C library**, if you want to use Nitro from C, C++, MATLAB, Julia or another
language instead of Python: see [Use it without Python](#use-it-without-python-the-c-library-libroadrunners-c-api).

Tested on macOS (Apple Silicon). Linux and Windows should work, but they're untested.

## Use it

### Your first simulation

Write it exactly as you would in regular tellurium; only the import line is different:

```python
import tellurium_nitro as te

r = te.loada('''
    S1 -> S2; k1*S1      # S1 turns into S2 at rate k1*S1
    k1 = 0.5; S1 = 10
''')
r.simulate(0, 10, 100)      # how S1 and S2 change from time 0 to 10 (100 points); runs on the fast engine
r.plot()
r.steadyState()             # where the system settles
print(r.getScaledFluxControlCoefficientMatrix())   # which reaction controls the flow
```

### If your code uses roadrunner directly

Change only the import line. Everything below it stays the same.

```python
import roadrunner_nitro as roadrunner          # was: import roadrunner
r = roadrunner.RoadRunner('my_model.xml')  # an SBML file, or SBML text
r.simulate(0, 10, 100)                      # runs on the fast engine
```

### Use it without Python: the C library (libroadrunner's C API)

Many programs use libroadrunner through its shared library and C API instead of Python: C and C++ programs,
MATLAB, Julia, R, C#, Delphi and others. Nitro's engine provides the same thing: a shared library,
`libroadrunner_nitro_c_api`, that exports **all 350 functions** of libroadrunner's C API (`rrc_api.h`,
`rrc_types.h`, `rrc_utilities.h`) with the same names, arguments and memory rules. It does not need Python,
libroadrunner or libsbml.

It is built from the engine's repository, [libroadrunner_nitro](https://github.com/gshinkim/libroadrunner_nitro),
whose README explains how to build it, use it from C, Python (ctypes), MATLAB, Julia and C#, switch an existing
libroadrunner program to it, and where it differs from libroadrunner's C API (only where libroadrunner's was wrong).

```bash
cd ../libroadrunner_nitro
cargo build --release -p roadrunner_nitro_c_api      # makes target/release/<library file>
```

**How we know it's right.** The whole SBML Test Suite runs through the C library alone, with no Python engine
involved. Reproduce from this repository (the library is found in `../libroadrunner_nitro/target/release`, or set
`RRC_LIB`):

```bash
python benchmarks/sts/run_capi.py <suite>/cases/semantic capi.jsonl --all
python benchmarks/sts/run_capi_stochastic.py <suite>/cases/stochastic capi_stoch.jsonl 10000
```

### Run many simulations at once

Give each setting you want to change a list of values. Nitro runs one simulation per position in the lists, all
at the same time:

```python
import tellurium_nitro as te
r = te.loada('S1 -> S2; k1*S1; k1 = 0.5; S1 = 10')
results = r.simulateBatch({'k1': [0.1, 0.5, 1.0, 2.0]}, 0, 10, 100)
print(results.shape)        # (4, 100, 3): 4 runs x 100 time points x 3 columns (time, S1, S2)
first_run = results[0]      # the simulation with k1 = 0.1
```

Each run starts from the model as it is now, and the model itself isn't changed. Species names mean
concentrations, the same as `r.S1`. Lists that change several settings are paired up position by position;
they don't make every possible combination.

### Switches

- `te.antimony(False)` reads Antimony with libantimony instead of Nitro's faster Rust reader (see below).
- `pip install ".[gpu]"` adds PyTorch for `te.gpu(True)` (see [GPU](#from-tellurium-turbo-batched-control-coefficients-and-gpu)).
- The engine is compiled for your own computer's processor (`libroadrunner_nitro/.cargo/config.toml`,
  `target-cpu=native`). To build a copy you can share with other machines, delete that file before step 5.

## How much faster is it?

All timings are on an Apple M5 Pro (18 cores, 48 GB), macOS, Python 3.12, libroadrunner 2.10.0.

### In short

| What you do | Regular tellurium | Tellurium Nitro | Faster by |
|---|---:|---:|---:|
| Load a 1000-species, 3400-reaction model | 30.9 s | 0.30 s | **103x** |
| Simulate that model (`simulate(0, 100, 1000)`) | 284 ms | 38 ms | **7.5x** |
| Simulate it again after `reset()` | 282 ms | 37 ms | **7.6x** |
| Simulate a 40-species model | 1.14 ms | 0.22 ms | **5.2x** |
| Steady state, 40-species model | 244 ms | 0.029 ms | **8,407x** |
| Control coefficients, 40-species model | 590 ms | 0.0173 ms | **34,228x** |
| Sensitivities, 86 parameters | 9.31 s | 14.8 ms | **631x** |
| Control coefficients for 200 parameter sets | 251 s | 1.87 ms | **134,024x** |

### Why 5x for some things and 5,000x for others?

**Simulating a model over time** ("start here, show me the next 100 seconds") is like driving from A to B.
Roadrunner was already a good driver on that trip: it uses a top-quality solver (CVODE) and turns each model into
fast machine code. Nitro drives the same road the same way, just with a lighter car: leaner code, fewer wasted steps,
less overhead. When the old engine was already efficient, 2-8x is what's left to gain.

**The big numbers come from questions where roadrunner takes a very long detour:**

- **"Where does the system settle?" (steady state).** Roadrunner can rebuild its whole copy of the model, sometimes
  twice, every time you ask. It's like rebuilding the kitchen every time you cook dinner. Nitro builds it once and
  reuses it.
- **"Which reaction is really in control?" (control coefficients).** Roadrunner answers by trial and error: nudge each
  number a tiny bit, re-solve the whole model, see what moved, and repeat for every number. That's about 170,000
  extra re-runs on the 40-species model, and the answers come out slightly fuzzy. Nitro works the answer out exactly
  with algebra, the way you'd take a derivative in calculus class, then solves it in one step.
- **"How sensitive is the result to each setting?" (sensitivities).** Roadrunner checks the 86 settings one at a time
  on a single processor core. Nitro uses exact math and checks them side by side on all 18 cores.

**Loading a big model** is faster because roadrunner compiles the whole model with a heavy compiler (LLVM) before
it can do anything. Nitro reads the model straight from the file and compiles only what it needs, only when it's
worth it.

**Accuracy.** Nitro's control coefficients are more accurate than roadrunner's (error 1e-11 vs 4e-8), because it
uses exact derivatives instead of nudging. Time courses are within roadrunner's own error range at the same
settings; on the 1000-species model the largest difference at any point is 2.9e-5.

<details>
<summary><b>Full benchmark tables</b></summary>

Three versions are compared:

| Version | What does the math |
|---|---|
| Regular tellurium (PyPI 2.2.13.1, libroadrunner 2.10.0) | roadrunner (C++) |
| [Tellurium Turbo](https://github.com/gshinkim/tellurium-turbo) | roadrunner (C++); adds batched, parallel `te.controlCoefficients` and optional GPU |
| **Tellurium Nitro** (this repo) | Rust engine `roadrunner_nitro` + JIT (no libroadrunner) |

Each task runs in a fresh process. Reproduce with `python benchmarks/run.py`, which writes `benchmarks/results.md`
and `benchmarks/results.json`.

- **Time** is the median wall time per call.
- **Cores** is CPU time ÷ wall time, so it shows how many cores were busy.
- **MB** is peak memory (RSS).
- **Error** is the maximum deviation from a high-precision reference, relative to the output's largest value. The
  reference is plain roadrunner at rtol 1e-12, atol 1e-14, with finite-difference step 1e-5. It shows who is
  closer to the true answer, not only whether the versions agree.

**Roadrunner's own benchmark set.** These are the `sosbench` models from the official
[roadrunner repo](https://github.com/sys-bio/roadrunner/tree/develop/data/sosbench), each simulated with its official
duration, steps and tolerances (`benchmarks/sosbench/`):

| Task | Regular tellurium | Tellurium Turbo | Tellurium Nitro | Nitro vs regular | Nitro vs Turbo | Error (regular / Turbo / Nitro) |
|---|---|---|---|---|---|---|
| sos:00009 | 1.41 ms · 1.0 cores · 248 MB | 1.47 ms · 1.0 cores · 249 MB | 0.752 ms · 1.0 cores · 195 MB | **1.9x** | **2.0x** | 2e-11 / 2e-11 / — |
| sos:00014 | 1.5 ms · 1.0 cores · 284 MB | 1.37 ms · 1.0 cores · 287 MB | 0.643 ms · 1.0 cores · 200 MB | **2.3x** | **2.1x** | 5e-06 / 5e-06 / — |
| sos:00022 | 5.32 ms · 1.0 cores · 248 MB | 5.31 ms · 1.0 cores · 248 MB | 2.83 ms · 1.0 cores · 196 MB | **1.9x** | **1.9x** | 2e-05 / 2e-05 / — |
| sos:00033 | 0.885 ms · 1.0 cores · 250 MB | 0.888 ms · 1.0 cores · 250 MB | 0.387 ms · 1.0 cores · 194 MB | **2.3x** | **2.3x** | 2e-09 / 2e-09 / — |
| sos:10000¹ | — | — | 118 ms · 1.0 cores · 215 MB | — | — | — / — / — |
| sos:10001 | 0.0219 ms · 1.0 cores · 247 MB | 0.0216 ms · 1.0 cores · 248 MB | 0.0104 ms · 1.0 cores · 193 MB | **2.1x** | **2.1x** | 9e-06 / 9e-06 / — |
| sos:repressilator¹ | — | — | 118 ms · 1.0 cores · 216 MB | — | — | — / — / — |

¹ Regular tellurium and Turbo fail on `10000` and `repressilator` with `CV_TOO_MUCH_WORK: The solver took mxstep (20000)
internal steps`. Nitro finishes them. A blow-up at atol 1e-4 makes it retry once with a tighter atol, and the
result matches the high-precision reference to 1e-8.

**The Tellurium Turbo benchmark (40-species model).** The model is from Turbo: a 40-species chain with saturating
side branches.

| Task | Regular tellurium | Tellurium Turbo | Tellurium Nitro | Nitro vs regular | Nitro vs Turbo | Error (regular / Turbo / Nitro) |
|---|---|---|---|---|---|---|
| ss | 244 ms · 1.0 cores · 265 MB | 243 ms · 1.0 cores · 266 MB | 0.029 ms · 1.0 cores · 195 MB | **8,407x** | **8,395x** | 0e+00 / 0e+00 / — |
| fcc | 590 ms · 1.0 cores · 266 MB | 581 ms · 1.0 cores · 265 MB | 0.0173 ms · 1.0 cores · 194 MB | **34,228x** | **33,686x** | 4e-08 / 4e-08 / — |
| sim | 1.14 ms · 1.0 cores · 261 MB | 1.14 ms · 1.0 cores · 260 MB | 0.22 ms · 1.0 cores · 195 MB | **5.2x** | **5.2x** | 3e-07 / 3e-07 / — |
| sens | 9.31 s · 1.0 cores · 274 MB | 9.17 s · 1.0 cores · 275 MB | 14.8 ms · 13.4 cores · 227 MB | **631x** | **621x** | — / 0e+00* / — |
| cc200 | 251 s · 1.0 cores · 267 MB | 8.2 s · 15.8 cores · 323 MB | 1.87 ms · 5.0 cores · 219 MB | **134,024x** | **4,389x** | 4e-08 / 4e-08 / — |
| load² | — | — | 75.4 ms · 1.0 cores · 340 MB | — | — | — / — / — |

Nitro's times above are from the latest engine (2026-11-24). Its error column shows — because this run
was made without the reference cache (`benchmarks/ref/`, computed once with plain roadrunner at tight
tolerances); the SBML Test Suite results above are the accuracy check for this version.

- `ss` = `r.resetAll(); r.steadyState()`
- `fcc` = `r.getScaledFluxControlCoefficientMatrix()`
- `sim` = `r.resetAll(); r.simulate(0, 100, 1001)`
- `sens` = `r.timeSeriesSensitivities(0, 50, 101)` (86 parameters)
- `load` = `te.loadSBMLModel` of a 1000-reaction model + first `simulate` + `steadyState`
- `cc200` = scaled flux and concentration control coefficients for 200 parameter sets: `te.controlCoefficients` in
  Turbo and Nitro, a `steadyState` + getters loop in regular tellurium (timed on 20 sets, ×10).

² `load` was run for Nitro only. Nitro now builds the model straight from SBML and doesn't compile it in
roadrunner. Before that change the same task took 13.4 s, and roadrunner's load alone takes about 10–11 s at
1000 reactions, so regular tellurium needs at least that.

\* There is no high-precision reference for `sens`, because roadrunner sensitivities at tight tolerance take
more than 15 minutes. That row shows agreement with regular tellurium instead.

**Big model (1000 species, 3400 reactions, 5952 parameters).** A synthetic kinetics model in Antimony
(`benchmarks/bigmodel/synthetic_1000_species.ant`), loaded with `te.loada` and run with `r.simulate(0, 100, 1000)`.
Regular tellurium is the `tellurium` package with libroadrunner. Reproduce with
`python benchmarks/bigmodel/run.py`.

| Task | Regular tellurium | Tellurium Nitro | Speedup |
|---|---:|---:|---:|
| load (`te.loada`) | 30.9 s | 0.30 s | **103x** |
| first `simulate` | 284 ms | 38 ms | **7.5x** |
| `reset()` + `simulate` again | 282 ms | 37 ms | **7.6x** |

The largest difference from regular tellurium at any point is 2.9e-5.

</details>

## How do we know the answers are right?

A faster engine is only useful if it gives the right answers. Nitro is checked three ways.

### 1. The SBML Test Suite: 1923 of 1923 pass

The [SBML Test Suite](https://github.com/sbmlteam/sbml-test-suite) is the standard set of test models that checks
whether a simulator gets SBML (the common file format for these models) right. Each case comes with the correct
answer. Every simulation result below comes from the Rust engine (Nitro has no other engine).

| Test cases | Cases | Pass | Rate |
|---|---:|---:|---:|
| Semantic, time courses (Rust engine) | 1,789 | 1,789 | **100%** |
| Semantic, flux balance (`te.fba`) | 34 | 34 | **100%** |
| Stochastic (10,000 runs each) | 100 | 100 | **100%** |
| **All semantic + stochastic cases** | **1,923** | **1,923** | **100%** |

- **Semantic:** includes algebraic rules, fast reactions, `delay()`, events, `comp` models and the `distrib` random
  draws. Roadrunner itself refuses or handles several of these differently.
- **Flux balance:** these are linear programs, which roadrunner doesn't do either, so `te.fba` solves them with SciPy.
- **Stochastic:** judged by the suite's own rule, "Z or Y failing 0-1 times per test". With the suite's standard
  deviation test applied strictly, 96 of 100 pass. The other 4 involve highly skewed distributions. Checked at up to
  a million runs, they match the expected values, but that test assumes bell-shaped data.

Full report, runners and per-case results: [`benchmarks/sts/`](benchmarks/sts/report.md).

### 2. Tellurium's own tests

The same upstream sys-bio/tellurium pytest suite (385 tests) was run against each version with
`python benchmarks/tests.py <upstream checkout>`:

| Engine | Wall time | Passed | Failed | Skipped |
|---|---|---|---|---|
| regular | 191.9 s | 363 | 14 | 8 |
| turbo | 191.5 s | 363 | 14 | 8 |
| nitro | 187.7 s | 363 | 14 | 8 |

All three give identical outcomes; no test differs from regular tellurium. The 14 failures also happen in the
official release, so they are upstream issues. Wall time is about the same for all three because the suite is
mostly SED-ML/OMEX parsing, file I/O and plotting on tiny models, not simulation.

**This repository's own tests** cover the tellurium layer (including SED-ML/OMEX). The engine's tests (Rust core,
C API, libroadrunner's own Python test suite and the SBML Test Suite) are in the engine's repository. Answers are
compared with libroadrunner's, recorded once in `tellurium_nitro/tests/reference/` and
`libroadrunner_nitro/tests/reference/`, so the tests never need libroadrunner installed.

```bash
python -m pytest tellurium_nitro/tests                 # this repository
cd ../libroadrunner_nitro && python -m pytest tests    # the engine (see its README for the rest)
```

### 3. Old tellurium bugs, now fixed

Nitro also fixes those 14 failures and a few other long-standing tellurium bugs. Its own test suite
(`pytest tellurium_nitro/tests`) passes in full: 417 passed, 0 failed, 8 skipped (the skips need internet access or
optional extras). What that means for you:

- **Reading data files in SED-ML now works.** Simulation descriptions that load experimental data from CSV, TSV or
  NuML files used to crash with an `AssertionError`, because newer libsedml versions return their own copy of a
  data type. That was 13 of the failing tests.
- **`te.executeSEDML("my_file.omex")` works.** Running an `.omex` archive file through `executeSEDML` used to crash
  with `NameError: CombineArchive`.
- **Turning plotting off works.** After `te.disablePlotting()`, making a tiled figure used to crash with a
  `NameError`.
- **No more silently missing features.** If the optional `phrasedml` package is missing, the inline OMEX tools
  now say so clearly instead of quietly disappearing.
- **No leftover temporary folders.** Running an inline OMEX no longer leaves a temporary folder behind each time.
- **Smaller fixes:** `distributed_sensitivity_analysis` no longer crashes on its bounds check, `run_case(...,
  py_code=True)` writes its Python script again, and two test files that needed the old `six` package now run.

You can also compare Nitro with roadrunner on your own machine: `python ../libroadrunner_nitro/check.py` prints both results and timings.

## For experts: how it works

A fork of [sys-bio/tellurium](https://github.com/sys-bio/tellurium) in which roadrunner's compute engine is
rewritten in Rust: steady states, time courses, events, rules, control analysis (MCA), sensitivities and
Gillespie. Each model's math is JIT-compiled to native code with Cranelift, and batches run on every CPU core.
The engine package `roadrunner_nitro` provides libroadrunner's Python API (`RoadRunner`, `NamedArray`, `Logger`,
`Config`, `SelectionRecord`, solver and option classes) on the Rust engine, so `roadrunner.RoadRunner(sbml)` code
runs without libroadrunner. Tellurium Nitro's `RoadRunner` is a subclass of it that adds tellurium's helpers
(plotting, export, `gillespie`, ...), as regular tellurium extends libroadrunner's.

Plain-language explainer of what changed and how it is built: [`docs/specJIT.pdf`](docs/specJIT.pdf).
Technical rewrite spec: [`docs/roadrunner-rs-spec.md`](docs/roadrunner-rs-spec.md).

- **`ExtendedRoadRunner`** (`tellurium_nitro/roadrunner/extended_roadrunner.py`): a subclass of
  `roadrunner_nitro.RoadRunner` (libroadrunner's Python API, in the engine's `python/roadrunner_nitro/`, over the
  Rust core in its `src/rr/`) that adds tellurium's helpers: plotting, export, `gillespie`, the short Jarnac names.
  `simulateBatch`, sensitivities and `te.controlCoefficients` use the batch model built by
  `roadrunner_nitro`'s `from_sbml` (`rrs_bridge.from_sbml` here).
- **Rust engine** ([libroadrunner_nitro](https://github.com/gshinkim/libroadrunner_nitro), Python module
  `roadrunner_nitro`):
  - Rate laws become a flat tape with exact symbolic derivatives (Jacobians, elasticities, parameter derivatives).
  - The tape is JIT-compiled to native code with Cranelift after 2000 evaluations, so short one-shot runs don't pay
    the compile cost. The JIT output is bit-identical to the interpreter.
  - **Steady state:** Newton with conservation laws.
  - **Time courses:** variable-order BDF with CVODE-style step and Newton control at roadrunner's default accuracy
    and tighter (relative tolerance 1e-5 or below), Rodas4 at looser tolerances, where it is faster on small or
    event-heavy models. Both use a sparse LU with fill-reducing ordering.
  - **Events:** delay, priority, persistence, `useValuesFromTriggerTime` and `initialValue=false`. Pending events
    carry across `simulate` calls.
  - **Rules:** rate, assignment and algebraic rules. Models with nonlinear algebraic constraints are solved as DAEs
    with [diffsol](https://github.com/martinjrobins/diffsol); everything else uses Nitro's own solvers.
  - **SBML features:** variable compartments, conversion factors, variable stoichiometry, function definitions,
    `avogadro`, `rateOf`, fast reactions, `delay()`, `comp` and `distrib` random draws.
  - **MCA:** closed form, `C = -L (N_R E L)^-1 N_R` and `F = I + E C`.
  - **Sensitivities:** forward, one parameter per core.
  - **Gillespie SSA.**
- **SBML reader** (`libroadrunner_nitro/src/sbml/`, on the pure-Rust `libroadrunner_nitro/crates/sbml-io`): reads SBML Level 1-3 straight into
  the engine, including hierarchical `comp` models (flattened in `sbml-io/src/comp.rs`). No libsbml, and no
  roadrunner compile: a 1000-species model loads in about 0.3 s.
- **Antimony:** read by a pure-Rust reader (`libroadrunner_nitro/crates/antimony`, from
  [PathwayDesigner-rs](https://github.com/sys-bio/PathwayDesigner-rs)), about 4x faster to load than libantimony.
  Anything it rejects goes to libantimony. `te.antimony(False)` always uses libantimony. Species are listed in
  declaration order, which can differ from libantimony's, so select result columns by name, not position.
- **RoadRunner core** (`libroadrunner_nitro/src/rr/`): libroadrunner's instance behaviour in Rust: selections, get/set values,
  reset, simulate, `oneStep`, steady state (with the nleq1/nleq2/newton settings, presimulation and
  approximation), Jacobians, eigenvalues (LAPACK's `zgeev`, ported, so the order matches libroadrunner's),
  elasticities and control coefficients, conserved moiety analysis (libstruct's QR, ported, so the species order
  and `_CSUM` totals match), Gillespie, integrator settings, model editing (`libroadrunner_nitro/src/sbml/edit.rs`), SBML level
  and version conversion, and save/load state. The Python `ExtendedRoadRunner` and the C library are both thin
  layers over this one core.
- **Integrators:** cvode (BDF / Rodas4, honours tolerances per species, step limits and variable step output),
  gillespie, and libroadrunner's euler, rk4 and rk45 (`libroadrunner_nitro/src/fixed.rs`).
- **C library** (`libroadrunner_nitro/capi/`): all 350 functions of libroadrunner's C API on the core; see
  [Use it without Python](#use-it-without-python-the-c-library-libroadrunners-c-api).
- **Where Nitro gives the correct answer and libroadrunner doesn't** (each case has a test):
  - Control coefficients and elasticities are exact, not finite differences.
  - euler, rk4 and rk45 apply events after the start time; libroadrunner ignores them.
  - rk45 stays accurate on very small or very large values (sosbench 00009, 00033).
  - Steady states that Newton can't reach from the current state are still found (e.g. sosbench 00009
    after `simulate(0, 10)`); models with no steady state still fail.
  - Roadrunner carries sensitivity state into the next `timeSeriesSensitivities` call, even after `reset()`, and
    ignores event changes inside sensitivities; Nitro starts each call fresh and applies them.
  - Roadrunner returns uninitialised values from `getDependentRatesOfChange` when moiety analysis is on.
  - C API bugs, listed in the C library section above.

## From Tellurium Turbo: batched control coefficients and GPU

### Quick start: using the GPU

Call `te.gpu(True)` **once** per session. After that, every `te.controlCoefficients` call uses the GPU
automatically; you don't pass anything else.

```python
import tellurium_nitro as te, numpy as np

if __name__ == '__main__':             # needed in .py scripts (uses all CPU cores); not needed in Jupyter
    te.gpu(True)                       # once: picks NVIDIA (CUDA) or Apple (MPS), else CPU

    r = te.loada('''
      J0: $X0 -> S1; k1*X0
      J1: S1 -> S2;  k2*S1
      J2: S2 -> $X1; k3*S2
      X0 = 10; k1 = 1; k2 = 2; k3 = 0.5
    ''')

    res = te.controlCoefficients(r, {'k1': np.linspace(0.5, 3, 500)})   # 500 parameter sets in one call
    print(res['fcc'][0])               # flux control coefficients for the first parameter set

    te.gpu(False)                      # optional: back to CPU
```

Only `te.controlCoefficients` uses the GPU. Regular roadrunner calls (`r.simulate()`, `r.steadyState()`,
`r.getScaledFluxControlCoefficientMatrix()`, ...) run exactly as before.

### Full API

```python
import tellurium_nitro as te, numpy as np

te.gpu("true")        # finds an NVIDIA (CUDA) or Apple Silicon (MPS) GPU, falls back to CPU
te.getGPU()           # -> "cuda", "mps" or "cpu"

r = te.loada(model)
res = te.controlCoefficients(r, {'k1': np.linspace(0.5, 3, 1000)})
res['fcc']            # (1000, reactions, reactions)  scaled flux control coefficients
res['ccc']            # (1000, species, reactions)    scaled concentration control coefficients
res['elasticities']   # (1000, reactions, species)    scaled elasticities
res['ufcc'], res['uccc']        # unscaled versions
res['species'], res['reactions'], res['ok']   # ids, and False where no steady state was found (NaN rows)
```

- `samples` can be a dict of equal-length arrays (zipped, not a grid), a list of dicts, or a pandas DataFrame.
  `None` evaluates the model's current values once.
- Steady states and elasticities run in parallel CPU worker processes (all cores when there are 64+ samples;
  override with `processes=`). Scripts need an `if __name__ == '__main__':` guard; notebooks don't.
- Control coefficients are computed as `C = -L (N_R E L)^-1 N_R` and `F = I + E C` in one batched solve,
  on the GPU when `te.gpu` is on. Results match roadrunner's own matrices to ~1e-15 on CPU.
  Apple MPS only supports float32, so results there agree to ~1e-5.
- GPU support needs PyTorch: `pip install -e ".[gpu]"`. Without it `te.gpu` warns and stays on CPU.

### Turbo benchmark report (roadrunner path, before the Rust engine)

Apple Silicon Mac, wall clock. Task: scaled flux + concentration control coefficients for every parameter set.

200 parameter sets, 40-species model:

| What you run (Tellurium Turbo) | Regular tellurium | Tellurium Turbo | Speedup |
|---|---|---|---|
| `te.controlCoefficients(r, samples, processes=1)` (1 CPU core, no GPU) | 321 s¹ | 105.5 s | **3.0×** |
| `te.controlCoefficients(r, samples)` (default: all CPU cores, no GPU) | 321 s¹ | 20.0 s | **16.1×** |
| `te.gpu(True)` then `te.controlCoefficients(r, samples)` (all CPU cores + MPS GPU) | 321 s¹ | 22.8 s | **14.1×** |

**So what do you get from just turning on `te.gpu(True)`?** Compared with regular tellurium, 14.1× faster.
But almost all of that comes from `te.controlCoefficients` itself, not the GPU: the same call *without*
`te.gpu` is 16.1× faster. On a model this small the GPU is slightly slower (22.8 s vs 20.0 s), because copying
data to the GPU costs more than the matrix math it speeds up.

Where the GPU does help, the matrix math alone on a large model (2000 parameter sets, 200 species):

| What you run | GPU off | `te.gpu(True)` (MPS) | Speedup |
|---|---|---|---|
| Matrix algebra only² | 3.86 s | 3.01 s | **1.3×** |

¹ Loop of `steadyState()` + `getScaledFluxControlCoefficientMatrix()` + `getScaledConcentrationControlCoefficientMatrix()`,
timed on 20 sets and extrapolated to 200.
² Regular tellurium has no batched equivalent, so this compares Tellurium Turbo with the GPU off (numpy) vs on. Steady states are not included.

**Takeaway:** ~95% of the time is roadrunner's steady-state solve, which only runs on CPU, so the big win is the
algorithm (steady state + elasticities once, instead of roadrunner's per-call MCA) and using every core.
The GPU only accelerates the matrix algebra: it slightly *slows down* small models (transfer overhead) and helps
for large models with many samples. CUDA was not benchmarked.

## Changed files (vs sys-bio/tellurium)

- `libroadrunner_nitro/` (new, since moved to its own repository,
  [libroadrunner_nitro](https://github.com/gshinkim/libroadrunner_nitro)): Rust engine `roadrunner_nitro`
  - `src/expr.rs` tape + derivatives, `src/jit.rs` Cranelift JIT, `src/model.rs` solvers/integrators/events/SSA/MCA,
    `src/sens.rs` sensitivities, `src/rand.rs` `distrib` random draws, `src/fixed.rs` euler/rk4/rk45
  - `src/sbml/` SBML reader, model editing, level/version conversion; `src/rr/` the RoadRunner core;
    `src/py.rs` Python bindings
  - `crates/sbml-io` (SBML reader/writer, `comp` flattening), `crates/antimony` (Antimony reader)
  - `capi/` (new): the C library `libroadrunner_nitro_c_api`, libroadrunner's headers, tests, `examples/simulate.c`
  - `tests/`: core, editing, integrator, conversion and limits tests against libroadrunner's recorded answers
- `tellurium_nitro/roadrunner/extended_roadrunner.py`: a subclass of `roadrunner_nitro.RoadRunner` (libroadrunner's
  Python API on the Rust core) adding tellurium's helpers
- `tellurium_nitro/analysis/rrs_bridge.py` (new): SBML -> Rust batch model for `simulateBatch`, sensitivities and
  `te.controlCoefficients`
- `tellurium_nitro/analysis/gpu.py` (new): `te.gpu()`, `te.getGPU()`, `te.controlCoefficients()` (Rust first)
- `tellurium_nitro/analysis/fba.py` (new): `te.fba()`, flux balance analysis with SciPy HiGHS
- `tellurium_nitro/__init__.py`: exports `te.turbo`, `te.gpu`, `te.getGPU`, `te.controlCoefficients`
- `tellurium_nitro/teconverters/`: inline OMEX works for archives that reference BioModels by URN
- `tellurium_nitro/tests/test_turbo.py` (+ `make_reference.py`, `reference/turbo.npz`), `tellurium_nitro/tests/test_gpu.py`,
  `tellurium_nitro/tests/test_bugfixes.py` (new): results vs libroadrunner's recorded answers, and bug regressions
- `tellurium_nitro/visualization/sbmldiagram.py`: removed the Python 2 `six` import
- `setup.py` / `requirements.txt`: `libroadrunner_nitro` instead of libroadrunner and python-libsbml; `phrasedml`;
  optional extra `gpu` (`torch>=2.0`)
- `benchmarks/` (new): `run.py` (speed/efficiency/accuracy), `tests.py` (pytest suite on all three), `sosbench/`,
  `sts/` (SBML Test Suite runners through Python, the core and the C API), `bigmodel/` (1000-species model)
- `docs/specJIT.pdf` (plain-language explainer), `docs/roadrunner-rs-spec.md` (technical spec)

> **Everything below is the original sys-bio/tellurium README.** Its `pip install tellurium` instructions install
> *regular* tellurium. To get Tellurium Nitro, follow [Install Tellurium Nitro](#install-tellurium-nitro-step-by-step) above.

---

<h1><img title="tellurium logo" src="./docs/images/tellurium_logo.png" height="50" />&nbsp;&nbsp;tellurium</h1>

<!-- These badges come from shield.io, zenodo or travis -->

<table style="width:100%">
  <tr>
    <!--<td><img alt="Github version", src="https://travis-ci.com/sys-bio/tellurium.svg?branch=master"></td> -->
    <td><img alt="Read the Docs" src="https://img.shields.io/readthedocs/tellurium"></td>
    <td><a href="https://badge.fury.io/gh/sys-bio%2Ftellurium"><img src="https://badge.fury.io/gh/sys-bio%2Ftellurium.svg" alt="GitHub version" height="18"></a></td>
  </tr>
</table> 

 <table style="width:100%">
  <tr>
    <td><img alt="Licence", src="https://img.shields.io/hexpm/l/tellurium"</td>
    <td><img alt="PyPI - Downloads", src="https://img.shields.io/pypi/dm/tellurium"></td>
    <td><img alt="Funding", src="https://img.shields.io/badge/Funding-NIGMS%20(GM123032)-blue"></td>
    <td><img alt="Funding", src="https://img.shields.io/badge/Funding-NIBIB%20(EB028887)-blue"></td>
    <td><a href="https://badge.fury.io/py/tellurium"><img src="https://badge.fury.io/py/tellurium.svg" alt="PyPI version" height="18"></a> </td>
   </tr>
</table> 

 <table style="width:100%">
  <tr>
    <td><a href="https://doi.org/10.5281/zenodo.2548944"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2548944.svg" alt="DOI"></a></td>
  </tr>
</table>

Copyright 2014-2026
Kiri Choi, J Kyle Medley, Matthias König, Kaylene Stocking, Caroline Cannistra, Michal Galdzicki, Ciaran Welsh, Lucian Smith, Adel Heydarabadipour, and Herbert Sauro

## What is it?

Tellurium is a comprehensive Python accessible SBML-compliant modeling platform for biochemical reaction systems. 

## Introduction

For installation instructions scroll down a few sections.

Tellurium is a python environment for reproducible dynamical modeling of biological networks. 
Tellurium provides the interface code to convert between standard formats and uses powerful 
libraries without requiring technical expertise, allowing teh user to focus on what’s important: 
building better models. Tellurium also provides first-class support for exchangeability via 
[COMBINE archives](http://co.mbine.org/documents/archive) and SBML, allowing you to share your models 
and simulations with many other tools.

Tellurium combines state-of-the-art scientific Python libraries, such 
as [NumPy](http://www.numpy.org/) and [SciPy](http://www.scipy.org/), 
and includes special-purpose systems biology Python tools. Out of the box, 
Tellurium includes the high performance simulation engine [libroadrunner](https://github.com/sys-bio/roadrunner), the model descrption language [antimony](http://antimony.sourceforge.net/), [phrasedml](http://phrasedml.sf.net/), [libsbml](http://sbml.org/Software/libSBML), and [libsedml](https://github.com/fbergmann/libSEDML).

The Tellurium (and libroadrunner project) project is funded from the NIH/NIGMS (GM081070) and NIH/NIBIB U24EB028887.

## Documentation 
* General: http://tellurium.readthedocs.org/en/latest/
* API: http://tellurium.readthedocs.io/en/latest/_apidoc/tellurium.html
* Webpage: http://tellurium.analogmachine.org/

## Usage Example

```{python}
import tellurium as te

r = te.loada('''
      S1 -> S2; k1*S1
      S1 = 10; S2 = 0
      k1 = 0.1
''')

result = r.simulate(0, 40, 500) 
r.plot()
```

<img src="./docs/images/tellurium-front-page-image.png" width=50% height=50%>

## Installation Instructions

Tellurium can be installed as a Python package, with a front-end (Notebook, IDE), or in a cloud environment (Google Colab):


A. [Python Pip](#Install-via-pip)

B. [IDE front-end](#Front-end-installation) based on [Spyder](https://www.spyder-ide.org/)  

C. [Notebook front-end](#front-end-2-tellurium-notebook)  

D. [Google Colab](#google-colab)  

We recommend first-time users choose one of the front-ends, while developers looking to integrate Tellurium use the pip package. All python desktop IDE graphical front-ends should work. We know for certain that telluirum will work from Spyder, pyCharm and Visual Studio Code.

The IDE front-end provides a MATLAB like experience with a code editor and Python console. 
The notebook front-end provides a notebook interface similar to [Jupyter](http://jupyter.org/), 
and features notebook cells for inline OMEX, a human-readable representation of COMBINE archives.

For any issues with installation and setup, please check [FAQ](https://github.com/sys-bio/tellurium/wiki/FAQ) or [contact us](#contact-us). 

-------

### Install via pip
[![PyPI version](https://badge.fury.io/py/tellurium.svg)](https://badge.fury.io/py/tellurium)

If you have a version of Python in an environment you're comfortable with, you can add Tellurium with pip:

```
pip install tellurium
```	

This will install on all platforms. For Colab see the end of this screen. 

For those using Anaconda, we currently do not have a separate Anaconda version of Tellurium, so you'll need to install it via pip, as above.

:exclamation:**For detailed instructions on how to setup Tellurium on Anaconda distributions, 
see [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).** 


For developers, the latest stable version from the repository can be installed via
```
pip install git+https://github.com/sys-bio/tellurium.git
```

-------

### Front-end installation

<img align="left" width="32px" id="windows" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/windows.png">
<h4>Windows</h4>
<br style="clear:both"/>

#### Approach 1: Spyder for Tellurium IDE
The Tellurium Spyder installers are tested with Windows 10 and 11, and come with Python 3.12.

:exclamation:**Attention to those upgrading to Tellurium Spyder IDE from an older version :** 
We strongly suggest you to completely remove the older version of Tellurium Spyder IDE prior to installing the latest version.

1. Download Tellurium Spyder for Windows:
    * [Python 3.12](https://sourceforge.net/projects/pytellurium/files/Tellurium-2.4/2.4.2/Tellurium-2.4.2-Python-3.12-win64-setup.exe/download)
2. Double-click the installer to start the installation
3. Follow the instructions

If you wish to use a different version of Python, you'll need to install Tellurium from the windows command line by using the command: `pip install tellurium`

NOTE: Installation requires administrative rights. It is recommended to accept the default settings.

#### Approach 2: Thonny IDE for Python
Tellurium can be used with [Thonny](https://en.wikipedia.org/wiki/Thonny), a free, open-source IDE for Python that is designed for beginners.

1. [Download the latest version of Thonny for Windows](https://thonny.org)
2. Double-click the installer `*.exe` to start the installation
3. Follow the instructions
4. Launch Thonny
5. Go to Tools > Manage packages...
6. Type `tellurium` in the search box and press the "Search on PyPI" button
7. Once the `tellurium` package is found in "Search results", click on it, and then press the "Install" button. Wait until the installation process is completed.

&nbsp;

<img align="left" width="32px" id="mac-osx" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/macos.png">
<h4>Mac OS X</h4>
<br style="clear:both"/>

For Mac users there a a number of options. By far, the easiest is to use the self-cotbained Thonny Python IDE.

#### Approach 1: Thonny IDE for Python
Tellurium can be used with [Thonny](https://en.wikipedia.org/wiki/Thonny), a free, open-source IDE for Python that is designed for beginners.

1. [Download the latest version of Thonny for Mac OS X](https://thonny.org)
2. Double-click the installer `*.pkg` to start the installation
3. Follow the instructions
4. Launch Thonny from Spotlight or directly from your Applications folder
5. Go to Tools > Manage packages...
6. Type `tellurium` in the search box and press the "Search on PyPI" button
7. Once the `tellurium` package is found in "Search results", click on it, and then press the "Install" button. Wait until the installation process is completed.

#### Approach 2: Spyder for Tellurium IDE

1. [Download the Anaconda distribution for Mac OS X.](https://www.anaconda.com/download). 

2. [Follow the installer instructions.](https://docs.anaconda.com/anaconda/install/linux). The installer will install Spyder as part of the installation process. If you customize the installation, you will have a chance to select which packages to install - ensure Spyder is selected.
3. When prompted, add Anaconda to your PATH (optional, but this will make the following steps easier).
4. [Open a Terminal](https://www.wikihow.com/Open-a-Terminal-Window-in-Mac) and run the command:
```
conda install msgpack-python
```
Next, install Tellurium itself:
```
pip install tellurium
```
If this fails, try using the default absolute path to Anaconda, which will probably be something like:
```
/Users/<your-user>/opt/anaconda3/bin/pip install tellurium
```
5. Launch the Anaconda Navigator via Launchpad. In Anaconda Navigator, run Spyder and then try `import tellurium` within the Spyder editor or console.

#### Approach 3: Spyder for Tellurium IDE (only for Mac OS X Legacy)

Currently, only a Legacy version of the Mac Spyder IDE is available. 

1. [Download Tellurium IDE for Mac OS X 10.10 or later](https://github.com/sys-bio/tellurium/releases/download/1.3.5-rc3/Tellurium-1.3.5-Spyder-2.3.8-OSX.dmg)
2. Double-click the .dmg file to open a new window
3. Double-click the Spyder icon


**For those who wish to use Tellurium with Spyder IDE on Mac OS X, we also suggest that you install through [Anaconda](https://www.anaconda.com/) by 
following the instructions on [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).**

&nbsp;

<img align="left" width="32px" id="redhat" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/redhat.png">
<h4>Linux (RedHat)</h3>
<br style="clear:both"/>

<img align="left" width="32px" id="debian" src="https://raw.githubusercontent.com/wiki/sys-bio/tellurium/img/debian.png">
<h4>Linux (Debian)</h4>
<br style="clear:both"/>

**For those who wish to use Tellurium with Spyder IDE on Linux, we also suggest that you install through [Anaconda](https://www.anaconda.com/) by 
following the instructions on [this page](https://github.com/sys-bio/tellurium/wiki/FAQ#i-would-like-to-use-tellurium-on-anaconda-what-should-i-do).**

:question:Looking for old releases? [You can find them here.](https://sourceforge.net/projects/pytellurium/files/notebook) To install an old release over a new release, you will need to [manually wipe the data directory.](https://github.com/sys-bio/tellurium/wiki/FAQ#how-do-i-uninstall-tellurium-notebook)

-------

### Google Colab
Tellurium can be used entirely in a browser using [Google Colab](https://colab.research.google.com/).  Ideally, it will work with the following steps (updated May 2023)
  
1. run a cell with ```!pip install -q tellurium```
2. (Very important) Restart the runtime (From the menu: 'Runtime / Restart session')
3. Test by typing ```import tellurium as te```
4. te.getVersionInfo()  # Check for versions if you want to 
   Try a model
5. r = te.loada ('S1 -> S2; k1*S1; k1 = 0.1; S1 = 10'); r.simulate(); r.plot()
  
The Python version behind Colab changes periodically, so what worked one day may stop working the next, but the following Colab notebook worked when used at ICSB 2022:  https://colab.research.google.com/drive/1wddLftHNhetbozZY29r2HRkzQLl1F_fs#scrollTo=l1bCgW46-adR and will hopefully be instructive.

-------
  
## Citing

If you use Tellurium in your research, we would appreciate following citations in any works you publish:

Medley et al. (2018). ["Tellurium notebooks—An environment for reproducible dynamical modeling in systems biology."](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006220) *PLoS Computational Biology*, 14(6), e1006220.

Choi et al. (2018). ["Tellurium: An extensible python-based modeling environment for systems and synthetic biology."](https://www.sciencedirect.com/science/article/pii/S0303264718301254) *Biosystems*, 171, 74-79.

## Contact Us

For general questions, requesting help, or reporting bugs, feel free to use the [GitHub issue tracker](https://github.com/sys-bio/tellurium/issues).  You can also post to the  [Tellurium-discuss](https://groups.google.com/d/forum/tellurium-discuss) mailing list.

## Legal

The source code for the Tellurium Python package is hosted at https://github.com/sys-bio/tellurium and is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0). Tellurium uses third-party dependencies which may be licensed under different terms. Consult the documentation for the respective third-party packages for more details.

TELLURIUM AND ALL SOFTWARE BUNDLED WITH TELLURIUM (HEREAFTER "THE SOFTWARE") IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THE SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## Every command, and what it does

Everything you can call in Tellurium Nitro. `te` is `import tellurium_nitro as te`. `r` is a model you loaded, for
example `r = te.loada('S1 -> S2; k1*S1; k1 = 0.5; S1 = 10')`.

⚡ = runs on Nitro's fast Rust engine. Everything else works exactly as in regular tellurium.
**New** = only in Nitro (or Tellurium Turbo), not in regular tellurium.

### `te.` — loading models

| Command | What it does |
|---|---|
| `te.loada(text)` | Load a model written in Antimony (the easy text format) and return `r`. Same as `te.loadAntimonyModel`. |
| `te.loadAntimonyModel(text)` | Same as `te.loada`. |
| `te.loadSBMLModel(sbml)` | Load a model from an SBML file path or SBML text and return `r`. Same as `te.loads`. |
| `te.loads(sbml)` | Same as `te.loadSBMLModel`. |
| `te.loadCellMLModel(cellml)` | Load a CellML model and return `r`. |
| `te.loadTestModel(name)` | Load one of roadrunner's built-in test models. |
| `te.getTestModel(name)` | Get a built-in test model as SBML text. |
| `te.listTestModels()` | List the built-in test models. |
| `te.model(name)` | Get back a model you loaded earlier, by its model name. |

### `te.` — converting between formats

| Command | What it does |
|---|---|
| `te.antimonyToSBML(text)` | Antimony text → SBML text. |
| `te.sbmlToAntimony(sbml)` | SBML text or file → Antimony text. |
| `te.antimonyToCellML(text)` | Antimony → CellML. |
| `te.sbmlToCellML(sbml)` | SBML → CellML. |
| `te.cellmlToAntimony(cellml)` | CellML → Antimony. |
| `te.cellmlToSBML(cellml)` | CellML → SBML. |
| `te.antimonyConverter()` | Converter object used by the notebook front end. |
| `te.getODEsFromModel(r)` | Print the model's equations (how each chemical changes) as text. |
| `te.getODEsFromSBMLFile(path)` | Same, from an SBML file. |
| `te.getODEsFromSBMLString(sbml)` | Same, from SBML text. |

### `te.` — Nitro switches (**new**)

| Command | What it does |
|---|---|
| `te.turbo()` | Kept for old scripts: always `True` (the Rust engine is the only engine; `te.turbo(False)` warns and changes nothing). |
| `te.antimony()` / `te.antimony(False)` | Read Antimony with Nitro's Rust reader (default, ~4x faster) or with the original libantimony. |
| `te.SBML()` | Kept for old scripts: always `True` (SBML is read by the Rust reader). |
| `te.gpu(True)` / `te.gpu(False)` | Use the GPU (NVIDIA or Apple) for the matrix math in `te.controlCoefficients`. |
| `te.getGPU()` | Which device `te.controlCoefficients` uses: `"cuda"`, `"mps"` or `"cpu"`. |

The `TELLURIUM_ENGINE` environment variable is no longer read (there is no roadrunner engine to switch to).

### `te.` — analysis

| Command | What it does |
|---|---|
| `te.controlCoefficients(r, samples)` ⚡ **new** | Control coefficients and elasticities for many parameter sets in one call, on all cores (and the GPU if on). Returns a dict with `fcc`, `ccc`, `elasticities`, `ufcc`, `uccc`, `species`, `reactions`, `ok`. |
| `te.fba(model)` **new** | Flux balance analysis: the best steady-state fluxes of an SBML `fbc` model (file path or text). Returns a dict of flux values and the objective. |
| `te.ParameterScan(r, ...)` | Run a time course for a range of values of one parameter and plot them. |
| `te.ParameterScan2D(r, ...)` | Same, for two parameters at once. |
| `te.SteadyStateScan(r, ...)` | Steady state for a range of values of one parameter. |
| `te.plotBifurcation(r, parameter, lower, upper)` | Draw a bifurcation diagram (how steady states change, split or vanish as a parameter changes). |
| `te.SensitivityAnalysis(...)` | Set up a sensitivity analysis to run with `te.distributed_sensitivity_analysis`. |
| `te.StochasticSimulationModel(...)` | Set up many random (Gillespie) runs for `te.distributed_stochastic_simulation`. |
| `te.ParameterEstimation(...)` | Fit parameters to data with a stochastic model. |
| `te.UncertaintySingleP(...)` | How much one parameter's uncertainty changes the output, with a confidence band plot. |
| `te.UncertaintyAllP(...)` | Same, for all parameters together. |
| `te.distributed_parameter_scanning(sc, ...)` | Parameter scan on a Spark cluster (`sc` = SparkContext). |
| `te.distributed_stochastic_simulation(sc, ...)` | Many stochastic runs on a Spark cluster. |
| `te.distributed_sensitivity_analysis(sc, ...)` | Sensitivity analysis on a Spark cluster. |
| `te.plot_distributed_stochastic(data)` / `te.plot_stochastic_result(result)` | Plot the results of those runs. |
| `te.sample_plot(result)` | Quick plot of a distributed result. |
| `te.getEigenvalues(matrix)` | Eigenvalues of a matrix. |
| `te.rank(A)` / `te.nullspace(A)` / `te.rref(A)` | Rank, null space and reduced row echelon form of a matrix. |

### `te.` — plotting

| Command | What it does |
|---|---|
| `te.plot(x, y)` | Make a line/scatter plot. |
| `te.show()` | Show the current figure. |
| `te.plotArray(result)` | Plot a result array, one line per column. |
| `te.plotWithLegend(r, result)` | Plot a result with a legend from `r`'s names. |
| `te.plotImage(img)` | Show an image. |
| `te.nextFigure()` | Start a new figure. |
| `te.tiledFigure()` / `te.newTiledFigure(...)` | Start a grid of small plots. |
| `te.newLowerTriFigure(...)` | Start a lower-triangle grid of plots. |
| `te.clearTiledFigure()` | Clear the grid. |
| `te.getPlottingEngine()` | Get the plotting engine (`matplotlib` or `plotly`). |
| `te.getDefaultPlottingEngine()` / `te.setDefaultPlottingEngine(name)` | Get or set the default plotting engine. |
| `te.disablePlotting()` | Turn all plotting off (useful on servers). |
| `te.setSavePlotsToPDF(True)` | Save plots to PDF as well. |
| `te.SBMLDiagram(sbml)` | Draw the reaction network as a diagram (needs graphviz). |
| `te.LatexExport(r, ...)` | Export a simulation plot as LaTeX. |

### `te.` — SED-ML, COMBINE archives and OMEX (sharing simulations)

| Command | What it does |
|---|---|
| `te.executeSEDML(file)` | Run the simulations described in a SED-ML file or `.omex` archive. |
| `te.sedmlToPython(file)` | Turn a SED-ML file into the equivalent Python script. |
| `te.executeCombineArchive(path)` | Run every SED-ML simulation in a COMBINE archive (`.omex`). |
| `te.convertCombineArchive(path)` | Read an archive and convert its contents to inline OMEX text. |
| `te.convertAndExecuteCombineArchive(path)` | Read an archive, convert it and run it. |
| `te.createCombineArchive(path, files, ...)` | Make a new archive from files. |
| `te.addFileToCombineArchive(...)` / `te.addFilesToCombineArchive(...)` | Add files to an archive (saved as a new archive). |
| `te.extractFileFromCombineArchive(path, entry)` | Get one file out of an archive as text. |
| `te.inlineOmex(...)` | Combine Antimony and phraSED-ML text into inline OMEX. |
| `te.executeInlineOmex(text)` | Run inline OMEX (Antimony model + phraSED-ML simulation). |
| `te.executeInlineOmexFromFile(path)` | Same, from a file. |
| `te.exportInlineOmex(text, path)` | Save inline OMEX as an archive file. |
| `te.getLastReport()` / `te.setLastReport(report)` | Get or set the last report a SED-ML run produced. |

### `te.` — files, versions and other helpers

| Command | What it does |
|---|---|
| `te.saveToFile(path, text)` / `te.readFromFile(path)` | Write or read a text file. |
| `te.listFiles(pattern)` | List files matching a pattern, such as `'*.xml'`. |
| `te.runTool(name)` | Run an external program. |
| `te.getVersionInfo()` / `te.printVersionInfo()` | Versions of tellurium and every library it uses. |
| `te.getTelluriumVersion()` / `te.__version__` | Tellurium's version number. |
| `te.noticesOff()` / `te.noticesOn()` | Hide or show tellurium's notices. |
| `te.getSupportedTestCases(begin, end)` | Which SBML test cases tellurium supports. |
| `te.getAppDir()` / `te.DumpJSONInfo()` | Internal helpers for the app folder and archive metadata. |

### `r.` — simulating

| Command | What it does |
|---|---|
| `r.simulate(start, end, points)` ⚡ | Run the model forward in time and return a table (time + each chemical). |
| `r.simulateBatch({'k1': [...]}, start, end, points)` ⚡ **new** | Many simulations at once on all cores, one per set of values. Returns an array of shape (runs, points, columns). |
| `r.gillespie(start, end, points)` ⚡ | Random (stochastic) simulation that counts individual molecules. |
| `r.oneStep(time, step)` ⚡ | Move the model forward by one step. |
| `r.plot()` | Plot the last simulation. |
| `r.show()` | Show the plot. |
| `r.lastRun` **new** | The last simulate/steady-state/MCA call: `{'method': ..., 'engine': 'rust', 'reason': None}`. |
| `r.timeCourseSelections = [...]` | Choose which columns `simulate` returns, for example `['time', 'S1', '[S2]']`. |
| `r.integrator.relative_tolerance = 1e-8` | Set accuracy (also `absolute_tolerance`, `maximum_num_steps`, `variable_step_size`, ...). |
| `r.setIntegrator('gillespie')` / `r.getIntegrator()` | Choose the time-stepping method (`cvode`, `gillespie`, `rk4`, `rk45`, `euler`). |
| `r.setIntegratorSetting(name, key, value)` | Change one setting of an integrator. |
| `r.getIntegratorByName(name)` / `r.makeIntegrator(name)` / `r.integratorExists(name)` | Look up or make integrators. |
| `r.getAvailableIntegrators()` / `r.getExistingIntegratorNames()` / `r.getRegisteredIntegratorNames()` | List integrators. |
| `r.setSeed(n)` / `r.getSeed()` / `r.resetSeed()` | Set, read or reset the random seed for stochastic runs. |
| `r.getSimulationData()` | The last simulation's results. |

### `r.` — steady state

| Command | What it does |
|---|---|
| `r.steadyState()` ⚡ | Find where the model settles (nothing changes any more); returns how close it got. |
| `r.steadyStateNamedArray()` ⚡ | Steady state values with names. |
| `r.getSteadyStateValues()` / `r.getSteadyStateValuesNamedArray()` ⚡ | Solve the steady state and return the selected values. |
| `r.mcaSteadyState()` ⚡ | Steady state prepared for control analysis. |
| `r.steadyStateSelections = [...]` / `r.getSteadyStateSelectionStrings()` | Choose or list what steady-state calls return. |
| `r.setSteadyStateSolver(name)` / `r.getSteadyStateSolver()` / `r.getSteadyStateSolverByName(name)` / `r.makeSteadyStateSolver(name)` / `r.steadyStateSolverExists(name)` | Choose or look up the steady-state solver. |
| `r.getExistingSteadyStateSolverNames()` / `r.getRegisteredSteadyStateSolverNames()` | List steady-state solvers. |
| `r.getSteadyStateThreshold()` / `r.setSteadyStateThreshold(x)` | How close to "no change" counts as steady. |

### `r.` — control analysis (MCA): which reaction controls what

| Command | What it does |
|---|---|
| `r.getScaledFluxControlCoefficientMatrix()` ⚡ | How much each reaction controls each flow, in percent terms. |
| `r.getUnscaledFluxControlCoefficientMatrix()` ⚡ | Same, in absolute units. |
| `r.getScaledConcentrationControlCoefficientMatrix()` ⚡ | How much each reaction controls each chemical level, in percent terms. |
| `r.getUnscaledConcentrationControlCoefficientMatrix()` ⚡ | Same, in absolute units. |
| `r.getScaledElasticityMatrix()` ⚡ | How sensitive each reaction is to each chemical, in percent terms. |
| `r.getUnscaledElasticityMatrix()` ⚡ | Same, in absolute units. |
| `r.getCC(variable, parameter)` ⚡ | One scaled control coefficient for a parameter. |
| `r.getuCC(variable, parameter)` ⚡ | One unscaled control coefficient. |
| `r.getEE(reaction, parameter)` / `r.getuEE(reaction, parameter)` | One scaled / unscaled elasticity. |
| `r.getScaledFloatingSpeciesElasticity(reaction, species)` | One scaled species elasticity. |
| `r.getUnscaledSpeciesElasticity(reaction, species)` | One unscaled species elasticity. |
| `r.getUnscaledParameterElasticity(reaction, parameter)` | One unscaled parameter elasticity. |
| `r.getDiffStepSize()` / `r.setDiffStepSize(x)` | Step size roadrunner uses to estimate slopes (Nitro uses exact slopes instead). |
| `r.getFluxThreshold()` / `r.setFluxThreshold(x)` | Flux below which `getCC` treats a flux as zero. |

### `r.` — Jacobians, eigenvalues and stability

| Command | What it does |
|---|---|
| `r.getFullJacobian()` / `r.fjac()` ⚡ | Table of how every rate changes when every chemical changes. |
| `r.getReducedJacobian()` ⚡ | Same, for the independent chemicals only. |
| `r.getFullEigenValues()` / `r.getFullEigenValuesNamedArray()` ⚡ | Eigenvalues of the full Jacobian (negative real parts = stable). |
| `r.getReducedEigenValues()` / `r.getReducedEigenValuesNamedArray()` ⚡ | Eigenvalues of the reduced Jacobian. |
| `r.getEigenValueIds()` | Names used to select eigenvalues. |
| `r.getFrequencyResponse(...)` | How the model responds to a signal wiggling at different frequencies. |

### `r.` — sensitivities

| Command | What it does |
|---|---|
| `r.timeSeriesSensitivities(start, end, points)` ⚡ | How much every chemical's time course shifts when each parameter is nudged. |
| `r.setSensitivitySolver(name)` / `r.getSensitivitySolver()` / `r.getSensitivitySolverByName(name)` / `r.makeSensitivitySolver(name)` / `r.sensitivitySolverExists(name)` | Choose or look up the sensitivity solver. |
| `r.getExistingSensitivitySolverNames()` / `r.getRegisteredSensitivitySolverNames()` | List sensitivity solvers. |

### `r.` — reading and changing values

| Command | What it does |
|---|---|
| `r.S1` / `r['[S1]']` / `r['S1']` | Read a value. `r.S1` and `r['[S1]']` give the concentration; `r['S1']` gives the amount. |
| `r.S1 = 5` / `r['k1'] = 2` | Change a value. |
| `r.getValue(name)` / `r.setValue(name, value)` | Read or change any value by name. |
| `r.setValues(names, values)` | Change several values at once. |
| `r.keys()` / `r.values()` / `r.items()` (and `iterkeys`, `itervalues`, `iteritems`) | Every name and value in the model, like a dictionary. |
| `r.getIds()` / `r.getSupportedIdTypes()` | Every name you can select or set. |
| `r.createSelection(name)` | Make a reusable selection object. |
| `r.getSelectedValues()` | Current values of the selected columns. |
| `r.resetSelectionLists()` | Put the selections back to the defaults. |
| `r.getGlobalParameterByName(name)` / `r.setGlobalParameterByName(name, value)` | Read or change a parameter. |
| `r.reset()` | Put chemicals back to their starting amounts (parameters stay as you changed them). |
| `r.resetAll()` | Put everything back to the current starting values, parameters too. |
| `r.resetParameter()` | Put parameters back to their starting values. |
| `r.resetToOrigin()` | Put the model back exactly as it was when first loaded. |
| `r.getCurrentTime()` | The model's current time. |

### `r.` — listing what is in the model

| Command | What it does |
|---|---|
| `r.getFloatingSpeciesIds()` / `r.fs()` | Chemicals that change. |
| `r.getBoundarySpeciesIds()` / `r.bs()` | Chemicals held fixed. |
| `r.getReactionIds()` / `r.rs()` | Reactions. |
| `r.getGlobalParameterIds()` / `r.ps()` | Parameters. |
| `r.getCompartmentIds()` / `r.vs()` | Compartments. |
| `r.getFloatingSpeciesConcentrations()` / `r.sv()` | Current concentrations. |
| `r.getFloatingSpeciesAmounts()` | Current amounts. |
| `r.getBoundarySpeciesConcentrations()` / `r.getBoundarySpeciesAmounts()` | Fixed chemicals' levels. |
| `r.getReactionRates()` / `r.rv()` | Current reaction speeds. |
| `r.getRatesOfChange()` / `r.dv()` | How fast each chemical is changing now. |
| `r.getIndependentRatesOfChange()` / `r.getDependentRatesOfChange()` | Same, split into independent and dependent chemicals. |
| `r.getGlobalParameterValues()` / `r.getCompartmentVolumes()` | Parameter values / compartment sizes. |
| `r.getConservedMoietyIds()` / `r.getConservedMoietyValues()` | Totals that never change (like total enzyme). |
| `r.getIndependentFloatingSpeciesIds()` / `r.getDependentFloatingSpeciesIds()` | Which chemicals are independent / determined by the totals. |
| `r.getFloatingSpeciesConcentrationIds()` / `r.getBoundarySpeciesConcentrationIds()` / `r.getFloatingSpeciesInitialConcentrationIds()` | Names for selecting concentrations or starting values. |
| `r.getAssignmentRuleIds()` / `r.getRateRuleIds()` / `r.getInitialAssignmentIds()` | Things set by rules or starting formulas. |
| `r.getNumFloatingSpecies()`, `r.getNumBoundarySpecies()`, `r.getNumReactions()`, `r.getNumGlobalParameters()`, `r.getNumCompartments()`, `r.getNumEvents()`, `r.getNumRateRules()`, `r.getNumConservedMoieties()`, `r.getNumIndFloatingSpecies()`, `r.getNumDepFloatingSpecies()`, `r.getNumberOfIndependentSpecies()`, `r.getNumberOfDependentSpecies()` | How many of each. |
| `...NamedArray` versions (e.g. `r.getFloatingSpeciesConcentrationsNamedArray()`, `r.getRatesOfChangeNamedArray()`, `r.getIndependentFloatingSpeciesAmountsNamedArray()`) | The same values with names attached. `...V` versions return plain lists. |
| `r.getKineticLaw(reaction)` | A reaction's rate formula. |
| `r.getHasOnlySubstanceUnits(species)` | Whether a chemical is measured as an amount only. |
| `r.getStoichiometryIds()` | Names of the stoichiometry values. |
| `r.getModelId()` / `r.getModelName()` / `r.setModelId(id)` / `r.setModelName(name)` | The model's id and name. |
| `r.getInfo()` | A summary of the loaded model and settings. |
| `r.isModelLoaded()` | Whether a model is loaded. |
| `r.getModel()` | The low-level model object. |

### `r.` — structure (stoichiometry and conservation)

| Command | What it does |
|---|---|
| `r.getFullStoichiometryMatrix()` / `r.sm()` | Which reactions make or use which chemicals. |
| `r.getReducedStoichiometryMatrix()` / `r.getNrMatrix()` | Same, for independent chemicals only. |
| `r.getExtendedStoichiometryMatrix()` | Same, including fixed chemicals. |
| `r.getReactantsStoichiometryMatrix()` / `r.getProductsStoichiometryMatrix()` | Only what reactions use / only what they make. |
| `r.getLinkMatrix()` / `r.getL0Matrix()` | How dependent chemicals follow from independent ones. |
| `r.getKMatrix()` | Null space of the reduced stoichiometry (steady-state flux patterns). |
| `r.getConservationMatrix()` | The totals that never change, as a matrix. |
| `r.conservedMoietyAnalysis = True` | Detect and use those totals (needed for some steady states). |

### `r.` — editing the model

| Command | What it does |
|---|---|
| `r.addSpecies(...)` / `r.addSpeciesConcentration(...)` / `r.addSpeciesAmount(...)` / `r.removeSpecies(id)` | Add or remove a chemical. |
| `r.addReaction(...)` / `r.removeReaction(id)` | Add or remove a reaction. |
| `r.addParameter(...)` / `r.removeParameter(id)` | Add or remove a parameter. |
| `r.addCompartment(...)` / `r.removeCompartment(id)` | Add or remove a compartment. |
| `r.addAssignmentRule(...)` / `r.addRateRule(...)` / `r.removeRules(id)` | Add or remove rules (formulas that set or change a value). |
| `r.addInitialAssignment(...)` / `r.removeInitialAssignment(id)` | Add or remove a starting-value formula. |
| `r.addEvent(...)` / `r.addTrigger(...)` / `r.addDelay(...)` / `r.addPriority(...)` / `r.addEventAssignment(...)` / `r.removeEvent(id)` / `r.removeEventAssignments(...)` | Build or remove events (things that happen when a condition becomes true). |
| `r.setPersistent(...)` / `r.setTriggerInitialValue(...)` | Event settings. |
| `r.setKineticLaw(reaction, formula)` | Change a reaction's rate formula. |
| `r.setReversible(reaction, True)` | Make a reaction reversible or not. |
| `r.setBoundary(species, True)` | Hold a chemical fixed or let it change. |
| `r.setConstant(id, True)` | Mark a value as constant. |
| `r.setHasOnlySubstanceUnits(species, True)` | Measure a chemical as an amount only. |
| `r.setInitAmount(...)` / `r.setInitConcentration(...)` / `r.setInitValue(...)` / `r.setInitStoichiometryValue(...)` | Change starting values. |
| `r.regenerateModel()` | Rebuild the model after edits. |
| `r.clearModel()` | Unload the model. |
| `r.load(sbml)` | Load a different model into `r`. |
| `r.validateCurrentSBML()` | Check the edited model is valid SBML. |

### `r.` — exporting and saving

| Command | What it does |
|---|---|
| `r.getSBML()` / `r.getCurrentSBML()` | SBML of the model as loaded / as it is now. |
| `r.getAntimony()` / `r.getCurrentAntimony()` | Antimony text, as loaded / as now. |
| `r.getCellML()` / `r.getCurrentCellML()` | CellML, as loaded / as now. |
| `r.getMatlab()` / `r.getCurrentMatlab()` | Matlab code, as loaded / as now. |
| `r.exportToSBML(path)` / `r.exportToAntimony(path)` / `r.exportToCellML(path)` / `r.exportToMatlab(path)` | Save the model in that format. |
| `r.getParamPromotedSBML(sbml)` | SBML with local parameters turned into global ones. |
| `r.saveState(path)` / `r.loadState(path)` | Save or restore the whole model and its state to a file. |
| `r.saveStateS()` / `r.loadStateS(data)` | Same, to and from memory. |
| `r.draw()` | Draw the reaction network (needs graphviz). |

### `r.` — engine information

| Command | What it does |
|---|---|
| `r.getCompiler()` | The roadrunner compiler object. |
| `r.getExtendedVersionInfo()` | roadrunner and library versions. |
| `r.getInstanceCount()` / `r.getInstanceID()` | How many roadrunner objects exist / this one's id. |
| `r.setOptions(...)` / `r.registerSolvers()` | Low-level roadrunner options. |

### Plain roadrunner code (**new**)

| Command | What it does |
|---|---|
| `import roadrunner_nitro as roadrunner` | Use instead of `import roadrunner`: libroadrunner's Python API (`RoadRunner`, `NamedArray`, `Logger`, `Config`, `SelectionRecord`, ...) on the fast engine, without libroadrunner. (`tellurium_nitro.nitro` is a deprecated alias.) |

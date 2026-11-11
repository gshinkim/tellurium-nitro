# Tellurium Nitro

A fork of [sys-bio/tellurium](https://github.com/sys-bio/tellurium) in which roadrunner's compute engine is
rewritten in Rust: steady states, time courses, events, rules, control analysis (MCA), sensitivities and
Gillespie. Each model's math is JIT-compiled to native code with Cranelift, and batches run on every CPU core.
**Your tellurium code does not change**: `te.loada(...)`, `r.simulate()`, `r.steadyState()`,
`r.getScaledFluxControlCoefficientMatrix()` and the rest return the same results, faster. Anything the Rust
engine can't do exactly the way roadrunner does goes to the original roadrunner without you doing anything.

| Version | What does the math |
|---|---|
| Regular tellurium (PyPI 2.2.13.1, libroadrunner 2.10.0) | roadrunner (C++) |
| [Tellurium Turbo](https://github.com/gshinkim/tellurium-turbo) | roadrunner (C++); adds batched, parallel `te.controlCoefficients` and optional GPU |
| **Tellurium Nitro** (this repo) | Rust engine `roadrunner_rs` + JIT; roadrunner only as a fallback |

Plain-language explainer of what changed and how it is built: [`docs/specJIT.pdf`](docs/specJIT.pdf).
Technical rewrite spec: [`docs/roadrunner-rs-spec.md`](docs/roadrunner-rs-spec.md).

## Benchmark: regular tellurium vs Tellurium Turbo vs Tellurium Nitro

Apple M5 Pro (18 cores, 48 GB), macOS. All three versions use the same Python 3.12 and libroadrunner 2.10.0, and each
task runs in a fresh process. Reproduce with `python benchmarks/run.py`, which writes `benchmarks/results.md` and
`benchmarks/results.json`.

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
| sos:00009 | 1.41 ms · 1.0 cores · 248 MB | 1.47 ms · 1.0 cores · 249 MB | 0.701 ms · 1.0 cores · 256 MB | **2.0x** | **2.1x** | 2e-11 / 2e-11 / 2e-11 |
| sos:00014 | 1.5 ms · 1.0 cores · 284 MB | 1.37 ms · 1.0 cores · 287 MB | 0.625 ms · 1.0 cores · 327 MB | **2.4x** | **2.2x** | 5e-06 / 5e-06 / 6e-06 |
| sos:00022 | 5.32 ms · 1.0 cores · 248 MB | 5.31 ms · 1.0 cores · 248 MB | 2.79 ms · 1.0 cores · 255 MB | **1.9x** | **1.9x** | 2e-05 / 2e-05 / 3e-05 |
| sos:00033 | 0.885 ms · 1.0 cores · 250 MB | 0.888 ms · 1.0 cores · 250 MB | 0.402 ms · 1.0 cores · 258 MB | **2.2x** | **2.2x** | 2e-09 / 2e-09 / 3e-10 |
| sos:10000¹ | — | — | 110 ms · 1.0 cores · 260 MB | — | — | — / — / 1e-08 |
| sos:10001 | 0.0219 ms · 1.0 cores · 247 MB | 0.0216 ms · 1.0 cores · 248 MB | 0.0176 ms · 1.0 cores · 254 MB | **1.2x** | **1.2x** | 9e-06 / 9e-06 / 3e-05 |
| sos:repressilator¹ | — | — | 110 ms · 1.0 cores · 260 MB | — | — | — / — / 1e-08 |

¹ Regular tellurium and Turbo fail on `10000` and `repressilator` with `CV_TOO_MUCH_WORK: The solver took mxstep (20000)
internal steps`. Nitro finishes them. A blow-up at atol 1e-4 makes it retry once with a tighter atol, and the
result matches the high-precision reference to 1e-8.

**The Tellurium Turbo benchmark (40-species model).** The model is from Turbo: a 40-species chain with saturating
side branches.

| Task | Regular tellurium | Tellurium Turbo | Tellurium Nitro | Nitro vs regular | Nitro vs Turbo | Error (regular / Turbo / Nitro) |
|---|---|---|---|---|---|---|
| ss | 244 ms · 1.0 cores · 265 MB | 243 ms · 1.0 cores · 266 MB | 0.055 ms · 1.0 cores · 268 MB | **4,433x** | **4,426x** | 0e+00 / 0e+00 / 2e-15 |
| fcc | 590 ms · 1.0 cores · 266 MB | 581 ms · 1.0 cores · 265 MB | 0.049 ms · 1.0 cores · 268 MB | **12,050x** | **11,859x** | 4e-08 / 4e-08 / 1e-11 |
| sim | 1.14 ms · 1.0 cores · 261 MB | 1.14 ms · 1.0 cores · 260 MB | 0.342 ms · 1.0 cores · 269 MB | **3.3x** | **3.3x** | 3e-07 / 3e-07 / 5e-08 |
| sens | 9.31 s · 1.0 cores · 274 MB | 9.17 s · 1.0 cores · 275 MB | 22.1 ms · 15.3 cores · 298 MB | **422x** | **415x** | — / 0e+00* / 5e-07* |
| cc200 | 251 s · 1.0 cores · 267 MB | 8.2 s · 15.8 cores · 323 MB | 90.9 ms · 1.1 cores · 293 MB | **2,755x** | **90x** | 4e-08 / 4e-08 / 3e-11 |

- `ss` = `r.resetAll(); r.steadyState()`
- `fcc` = `r.getScaledFluxControlCoefficientMatrix()`
- `sim` = `r.resetAll(); r.simulate(0, 100, 1001)`
- `sens` = `r.timeSeriesSensitivities(0, 50, 101)` (86 parameters)
- `cc200` = scaled flux and concentration control coefficients for 200 parameter sets: `te.controlCoefficients` in
  Turbo and Nitro, a `steadyState` + getters loop in regular tellurium (timed on 20 sets, ×10).

\* There is no high-precision reference for `sens`, because roadrunner sensitivities at tight tolerance take
more than 15 minutes. That row shows agreement with regular tellurium instead.

### Why is Nitro 2x faster in one table and thousands of times faster in the other?

The two tables measure different kinds of work.

**The first table (roadrunner's own benchmark) only plays the model forward in time**: "start here, show me what
happens over the next 100 seconds." Think of driving from A to B. Roadrunner was already a good driver on that trip:
it uses a top-quality solver (CVODE) and turns each model into fast machine code. Nitro drives the same road the
same way, just with a lighter car: leaner code, fewer wasted steps per move, less overhead. When both do the same
job and the old one was already efficient, the best you can do is around 2x.

**The big numbers in the second table come from questions where roadrunner takes a very long detour:**

- **"Where does the system settle?" (steady state).** Roadrunner can rebuild its whole copy of the model, sometimes
  twice, every time you ask. It's like rebuilding the kitchen every time you cook dinner. Nitro builds it once and
  reuses it.
- **"Which reaction is really in control?" (control coefficients).** Roadrunner answers by trial and error: nudge each
  number a tiny bit, re-solve the whole model, see what moved, and repeat for every number. That's about 170,000
  extra re-runs on this model, and the answers come out slightly fuzzy. Nitro works the answer out exactly with
  algebra, the way you'd take a derivative in calculus class, then solves it in one step.
- **"How sensitive is the result to each setting?" (sensitivities).** Roadrunner checks the 86 settings one at a time
  on a single processor core. Nitro uses exact math and checks them side by side on all 18 cores.

So the two tables aren't in conflict:

- **Same work, done better:** about 2-3x faster. That's the first table, plus the `sim` row in the second.
- **A smarter route that skips almost all the work:** 400 to 12,000x faster. That's steady state, control
  coefficients and sensitivities. Because Nitro computes these exactly instead of by trial and error, the answers
  are also more accurate.

**Accuracy in short.** Nitro's control coefficients are more accurate than roadrunner's (1e-11 vs 4e-8), because
it differentiates exactly instead of using finite differences. Time courses are within roadrunner's own error
range at the same tolerances.

**Tellurium's own test suite.** The same upstream sys-bio/tellurium pytest suite (385 tests) runs against each
install with `python benchmarks/tests.py <upstream checkout>`:

| Engine | Wall time | Passed | Failed | Skipped |
|---|---|---|---|---|
| regular | 191.9 s | 363 | 14 | 8 |
| turbo | 191.5 s | 363 | 14 | 8 |
| nitro | 187.7 s | 363 | 14 | 8 |

All three give identical outcomes; no test differs from regular tellurium. The 14 failures also happen in the
official release, so they are upstream issues. Wall time is about the same for all three because the suite is
mostly SED-ML/OMEX parsing, file I/O and plotting on tiny models, not simulation. The speed difference shows in the
compute benchmarks above.

## Install

```bash
pip install tellurium                       # dependencies (libroadrunner, antimony, ...)
git clone https://github.com/gshinkim/tellurium-nitro && cd tellurium-nitro
pip install -e .
brew install rust && pip install maturin    # or rustup (Linux/Windows) + pip install maturin
cd rrs && maturin develop --release         # builds the Rust engine into the active venv
python check.py                             # correctness vs roadrunner + timings
```

`rrs/.cargo/config.toml` builds for the local CPU (`target-cpu=native`); delete it to build a portable binary.
Without the Rust engine, Nitro behaves exactly like Tellurium Turbo. `te.turbo(False)` (or env
`TELLURIUM_ENGINE=roadrunner`) switches back to roadrunner globally. Only tested on macOS (Apple Silicon) so far.

## How it works

- **Rust engine** (`rrs/`, Python module `roadrunner_rs`):
  - Rate laws become a flat tape with exact symbolic derivatives (Jacobians, elasticities, parameter derivatives).
  - The tape is JIT-compiled to native code with Cranelift after 2000 evaluations, so short one-shot runs don't pay
    the compile cost. The JIT output is bit-identical to the interpreter.
  - **Steady state:** Newton with conservation laws.
  - **Time courses:** Rodas4 at default accuracy, and variable-order BDF with CVODE-style step and Newton control for
    tight tolerances. Both use a sparse LU with fill-reducing ordering.
  - **Events:** delay, priority, persistence, `useValuesFromTriggerTime` and `initialValue=false`. Pending events
    carry across `simulate` calls.
  - **Rules:** rate and assignment rules.
  - **SBML features:** variable compartments, conversion factors, variable stoichiometry, function definitions,
    `avogadro` and `rateOf`.
  - **MCA:** closed form, `C = -L (N_R E L)^-1 N_R` and `F = I + E C`.
  - **Sensitivities:** forward, one parameter per core.
  - **Gillespie SSA.**
- **Bridge** (`tellurium/analysis/rrs_bridge.py`): reads the model roadrunner loaded and builds the Rust model.
- **Facade** (`tellurium/roadrunner/extended_roadrunner.py`): for each call, it decides between Rust and roadrunner
  and keeps the roadrunner object in sync, so selections, getters and plotting keep working.
- **What Rust handles:**
  - `simulate`, including the gillespie integrator and `oneStep`
  - `steadyState` and its variants
  - every elasticity, control-coefficient, Jacobian and eigenvalue matrix, `getCC` and `getuCC`
  - `timeSeriesSensitivities` (k = 0, 1)
  - `te.controlCoefficients`
- **Still on roadrunner, and why:**
  - **Model loading, validation and export:** done once, and already fast.
  - **Scalar getters** (~1 µs in C++): the call into Rust would cost more than it saves.
  - **Structural matrices** (stoichiometry, link, conservation): computed once by libstruct, so the choice of
    independent species stays consistent.
  - **`getFrequencyResponse`:** no reference to verify a Rust version against.
  - **`k ≥ 2` sensitivities.**
  - **SBML `distrib` random functions.**
  - **Delay equations, algebraic rules and fast reactions:** roadrunner 2.10 refuses to load these models.
- **Known differences from roadrunner:**
  - Roadrunner carries sensitivity state into the next `timeSeriesSensitivities` call, even after `reset()`.
    Nitro starts each call fresh.
  - Roadrunner ignores event changes inside sensitivities.
  - Roadrunner returns uninitialised values from `getDependentRatesOfChange` when moiety analysis is on.
  - Nitro gives the correct values in all three cases. The tests document each one.

## From Tellurium Turbo: batched control coefficients and GPU

### Quick start: using the GPU

Call `te.gpu(True)` **once** per session. After that, every `te.controlCoefficients` call uses the GPU
automatically; you don't pass anything else.

```python
import tellurium as te, numpy as np

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
import tellurium as te, numpy as np

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

- `rrs/` (new): Rust engine `roadrunner_rs` (`src/expr.rs` tape + derivatives, `src/jit.rs` Cranelift JIT,
  `src/model.rs` solvers/integrators/events/SSA/MCA, `src/sens.rs` sensitivities, `src/py.rs` bindings, `check.py`)
- `tellurium/analysis/rrs_bridge.py` (new): roadrunner model -> Rust model
- `tellurium/roadrunner/extended_roadrunner.py`: facade routing roadrunner calls to Rust, `te.turbo()`
- `tellurium/analysis/gpu.py` (new): `te.gpu()`, `te.getGPU()`, `te.controlCoefficients()` (Rust first)
- `tellurium/__init__.py`: exports `te.turbo`, `te.gpu`, `te.getGPU`, `te.controlCoefficients`
- `tellurium/tests/test_turbo.py`, `tellurium/tests/test_gpu.py` (new): Rust / batched results vs plain roadrunner
- `tellurium/visualization/sbmldiagram.py`: removed the Python 2 `six` import
- `setup.py`: optional extra `gpu` (`torch>=2.0`)
- `benchmarks/` (new): `run.py` (speed/efficiency/accuracy), `tests.py` (pytest suite on all three), `sosbench/`
- `docs/specJIT.pdf` (plain-language explainer), `docs/roadrunner-rs-spec.md` (technical spec)

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

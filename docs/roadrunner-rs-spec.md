# roadrunner-rs: Rust rewrite of libRoadRunner — implementation spec

Audience: an AI coding agent implementing the engine. Read this whole file before writing code.
Source of truth for roadrunner behaviour: `sys-bio/roadrunner` @ `d16ff4aa` (2026-11-09). All `file:line`
references below are relative to that tree's root.

---

## 1. Goal, non-goals, baseline, targets

### Goal
A Rust SBML simulation engine, exposed to Python as `roadrunner_rs`, that:
1. Is numerically equivalent to libRoadRunner (within the tolerances in §6) for everything tellurium uses.
2. Is as fast as the hardware allows: no redundant work, exact derivatives, all CPU cores, SIMD where it pays, optional GPU for large batched linear algebra.
3. Can replace `roadrunner.RoadRunner` in tellurium behind a compatibility layer (§4.10), one API group at a time.

### Non-goals
- Bit-for-bit reproduction of roadrunner output. Note that roadrunner's own MCA uses finite differences and is only accurate to about 6–7 digits (§2.3).
- Matching roadrunner's RNG streams. roadrunner uses Mersenne Twister (`source/llvm/Random.cpp`). We match the distributions, not the seeds.
- C/C++ API, the C API wrappers (`wrappers/C`), plugins (`rrplugins/`), the Java/C# bindings, and the legacy C code generator (`source/c/`).
- SBML writing, validation, units checking, or the model-editing API (`addReaction`, `setKineticLaw`, …, `rrRoadRunner.cpp:6644-7905`). These stay with libsbml and antimony in Python. Editing triggers a full rebuild (§4.10).
- SBML features roadrunner does not support either: algebraic rules (`source/llvm/LLVMModelDataSymbols.cpp:186`), `fast` reactions (`LLVMModelDataSymbols.cpp:1472`), and the `delay()` csymbol (`source/llvm/ASTNodeCodeGen.cpp:623-640`). They are listed in §5 as later stretch goals.

### Measured baseline
Measured per parameter set on a 40-species model, Apple M5 Pro, libroadrunner from pip:

| call | time |
|---|---|
| `steadyState()` | 0.24 s |
| `getScaledFluxControlCoefficientMatrix()` | 0.57 s |
| `getScaledConcentrationControlCoefficientMatrix()` | 0.41 s |
| total per sample (current `te.controlCoefficients` roadrunner path) | ~1.2 s, ~95 % inside roadrunner |

### Targets
Single-thread time per parameter set, 40-species model. Assume R ≈ 50 reactions.

| operation | roadrunner | target | why it is reachable |
|---|---|---|---|
| RHS eval `f(x,p)` | ~1–3 µs (JIT, `-O0`, `rrConfig.cpp:141`) | ≤ 0.3 µs (Cranelift/AOT), ≤ 1 µs (tape) | about 200–500 flops |
| dense Jacobian (exact) | 2n FD RHS evals | ≤ 10 µs | forward AD with sparse coloring, or symbolic J tape |
| steady state (Newton, warm-ish start) | 0.24 s | ≤ 0.5 ms (≥ 400×) | ~10 iterations × (J + 40×40 LU ≈ 21 kflop) |
| steady state via presimulation fallback | n/a | ≤ 5 ms | BDF to t=1e4 with exact J, ~300 steps |
| unscaled elasticities E (R×S) | ~0.2–0.3 s | ≤ 10 µs (≥ 10⁴×) | same AD pass as J: `J = N·E` |
| CCC + FCC + scaled (given SS) | 0.41 + 0.57 s | ≤ 50 µs | one LU of the m×m reduced Jacobian plus two GEMMs |
| full MCA sample (SS + E + CCC + FCC) | ~1.2 s | ≤ 1 ms (≥ 1000×) | sum of the rows above |
| batch of B samples | B × 1.2 s / procs (process pool) | B × 1 ms / cores | rayon, no IPC, shared compiled model |

Scaling target: at least 0.85 × (physical performance cores) parallel efficiency for B ≥ 256.

The 1000× target is not a micro-optimisation claim. It holds because roadrunner does
O(R·S²) model resets with string lookups and recompiles the model several times for each call that a
reasonable engine does once (§2.3). Removing that redundant work is most of the win. Compiled
code and SIMD are secondary.

---

## 2. Current roadrunner architecture (map + where time goes)

### 2.1 Pipeline

| stage | files | notes |
|---|---|---|
| SBML load | `source/rrSBMLReader.cpp` (comp detection `:117`), `rrRoadRunner.cpp:1570` `load()` | libsbml parse; comp flattening via libsbml |
| conserved-moiety conversion | `source/conservation/ConservedMoietyConverter.cpp`, called from `source/llvm/ModelGeneratorContext.cpp:121-150` | rewrites the SBML document and adds `_CSUM` parameters; uses LibStructural (external dep `libroadrunner-deps`) |
| symbols / layout | `source/llvm/LLVMModelDataSymbols.cpp` (2241 lines), `LLVMModelSymbols.cpp`, `LLVMModelData.h` | flat arrays: compartments, floating/boundary species, global params, rate-rule values, reaction rates, CSR stoichiometry |
| codegen | `source/llvm/*CodeGen.cpp`: `EvalReactionRatesCodeGen`, `EvalRateRuleRatesCodeGen`, `EvalInitialConditionsCodeGen`, `EvalVolatileStoichCodeGen`, `EvalConversionFactorCodeGen`, `EventTriggerCodeGen`, `EventAssignCodeGen`, `GetValuesCodeGen`/`SetValuesCodeGen`, `GetPiecewiseTriggerCodeGen`; AST→IR in `ASTNodeCodeGen.cpp` (1583 lines); data access in `ModelDataIRBuilder.cpp` | each SBML entity becomes an LLVM function. Assignment rules are inlined at each use through symbol resolvers (`ModelDataSymbolResolver.cpp`, `SymbolForest.h`) |
| JIT | `source/llvm/LLJit.cpp` (`addModule` `:327-365`), `MCJit.cpp`, `JitFactory.cpp`; default LLJIT, optimisation level NONE (`rrConfig.cpp:140-141`) | compiled-model cache keyed by SBML MD5 (`LLVMModelGenerator.cpp:63-64, 441-537`), consulted only by `createModel` |
| executable model | `source/rrExecutableModel.h` (interface, 1053 lines), `source/llvm/LLVMExecutableModel.cpp` (3176 lines) | `getStateVectorRate` `:413` = volatile stoich + reaction rates + CSR `dgemv` (`rrSparse.cpp:171`) + conversion factors + rate rules |
| integrators | `source/CVODEIntegrator.cpp` (BDF/Adams `:899-905`, dense LU `:956-957`, **DQ Jacobian** `:976`, root finding for events `:928`), `RK4Integrator.cpp`, `RK45Integrator.cpp`, `EulerIntegrator.h`, `GillespieIntegrator.cpp` | factories: `IntegratorFactory.cpp`, `RegistrationFactory.cpp` |
| steady state | `NLEQ1Solver.cpp`, `NLEQ2Solver.cpp` → `rrNLEQ2Interface.cpp` (Fortran NLEQ2, **Jacobian = NULL ⇒ internal FD** `:214-216`), `KinsolSteadyStateSolver.cpp` + `NewtonIteration.cpp`/`BasicNewtonIteration.cpp`/`LinesearchNewtonIteration.cpp` (KINSOL, dense, DQ Jacobian: no `KINSetJacFn`), decorators `PresimulationProgramDecorator.cpp`, `PresimulationDecorator.cpp`, `ApproxSteadyStateDecorator.cpp` | default solver `nleq2` (`rrRoadRunner.cpp:604`) |
| sensitivities | `ForwardSensitivitySolver.cpp` (CVODES forward) | `timeSeriesSensitivities` `rrRoadRunner.cpp:2556` |
| events | `source/llvm/EventQueue.cpp`, `LLVMExecutableModel::applyEvents` `:2724`, `:2861` | priority, delay, persistence, `useValuesFromTriggerTime` |
| structural | LibStructural (`getLibStruct` `rrRoadRunner.cpp:774`): N, Nr, L0, L, K, conservation matrix | external C++ library; QR/LU-based |
| MCA | `rrRoadRunner.cpp:1916-2133` (`getEE`/`getuEE`), `:4407` `getuCC`, `:4546` `getCC`, `:4591-5091` matrices | see §2.3 |
| selections | `rrSelectionRecord.cpp` (517 lines), `RoadRunner::createSelection` `rrRoadRunner.cpp:5378`, `getValue(SelectionRecord)` `:1214` | string grammar: `S`, `[S]`, `S'`, `init(S)`, `init([S])`, `ec(J,S)`, `uec`, `cc`, `ucc`, `eigen(S)`, `eigenReal`, `eigenImag`, `stoich(S,J)`, `time`, `rateOf` |
| simulate | `rrRoadRunner.cpp:2153-2550` | fixed grid: one `integrator->integrate` per output point (`:2428-2450`), variable step at `:2266-2353`; `getSelectedValues` per row does a per-selection switch |
| Python | `wrappers/Python/roadrunner/roadrunner.i` (SWIG, 3277 lines) | `__getattr__` dynamic properties `:1279`, `simulate` `:1498`, `resetAll` `:1908`, `plot` `:1935`, `steadyStateSolver` `:2039`, `__getitem__` `:1175`; `NamedArray` in `named_array.py` |

### 2.2 Global and thread-safety constraints
These constraints are why `te.controlCoefficients` uses a process pool:
- NLEQ2 keeps a process-global `static ExecutableModel* callbackModel` (`rrNLEQ2Interface.cpp:25`) and throws if it is already set (`:201-203`). Only one NLEQ2 solve can run per process.
- `Config` is process-global. `getFullJacobian`/`getReducedJacobian` mutate `Config::ROADRUNNER_JACOBIAN_MODE` through a guard (`rrRoadRunner.cpp:3385-3386`, `3510`).
- The compiled-model cache is a global map behind a `Poco::Mutex` (`LLVMModelGenerator.cpp:63-64`).

The Rust design must be share-nothing per sample: an immutable `Arc<CompiledModel>` plus one owned `State` per thread.

### 2.3 Why `steadyState` and MCA are slow (concrete evidence)

**A. `steadyState()` recompiles the model twice per call when `conservedMoietyAnalysis` is off.**
Off is the default for `te.loada`, `tellurium/tellurium.py:539`.
- `auto_moiety_analysis` defaults to `true` (`source/SteadyStateSolver.cpp:44`).
- When it is true and moiety conversion is off, `steadyState` calls `setConservedMoietyAnalysis(true)` (`rrRoadRunner.cpp:1772`). If the model has 0 moieties it then calls `setConservedMoietyAnalysis(false)` (`:1775`). Otherwise it restores the old setting after the solve (`:1864`).
- Every `setConservedMoietyAnalysis` sets `RECOMPILE` and calls `regenerateModel(true)` (`rrRoadRunner.cpp:1890-1897`).
- `LLVMModelGenerator::regenerateModel` (`source/llvm/LLVMModelGenerator.cpp:278-294`) builds a fresh `ModelGeneratorContext` with no cache lookup. That means libsbml moiety conversion (`ModelGeneratorContext.cpp:121-150`), `doc->toSBML()` plus MD5, full LLVM codegen and JIT. It also clones the SBML document (`rrRoadRunner.cpp:1894`).
- Net effect: two full model builds per `steadyState()`. `getUnscaledConcentrationControlCoefficientMatrix` calls `steadyState()` once or twice (`:4884-4885`), so FCC and CCC pay the same cost again.
- *Verify first:* time `steadyState()` with `r.conservedMoietyAnalysis=True` set beforehand, or with `r.steadyStateSolver.auto_moiety_analysis=False`. If the 0.24 s falls by an order of magnitude, A dominates. This is also a zero-code speedup for the roadrunner fallback in `tellurium/analysis/gpu.py:_init_worker`.

**B. NLEQ2 uses a finite-difference Jacobian.**
- `NLEQ2(&n, &ModelFunction, NULL, …)` (`rrNLEQ2Interface.cpp:214-216`) passes `JAC = NULL`, so NLEQ2 builds J by numerical differencing: n extra RHS evaluations per Jacobian.
- Each RHS call goes through `setStateVector` plus `getStateVectorRate` (`:288-290`).
- KINSOL (`NewtonIteration.cpp:44-49`) and CVODE (`CVODEIntegrator.cpp:976`) also use difference-quotient Jacobians.
- Tolerance is `relative_tolerance=1e-12` with `maximum_iterations=100` (`NLEQSolver.cpp:35-38`). Presimulation retries at times `{0.1,1,10,100,1e3,1e4}` (`SteadyStateSolver.cpp:60`; loop at `PresimulationProgramDecorator.cpp:22-38`). Each retry constructs a new `CVODEIntegrator`.
- `steadyState()` returns `sqrt(Σ dxdt²)` (`rrNLEQ2Interface.cpp:509-519`).

**C. Each element of the elasticity matrix costs O(S) full model resets.**
`getUnscaledElasticityMatrix` (`rrRoadRunner.cpp:4766`) loops over all R×S pairs and calls `getUnscaledSpeciesElasticity(i,j)` (`:4794`). For each pair (`:4591-4763`):
- It reads all species values and all compartment volumes. Then `setAllFloatingSpeciesValues(setInitValuePtr, currentVals)` (`:4692`) sets S initial values one at a time.
- Each `setFloatingSpeciesInit{Concentrations,Amounts}` call triggers `reset(TIME|FLOATING)` (`LLVMExecutableModel.cpp:2952, 2971, 3012`).
- `reset` (`LLVMExecutableModel.cpp:778-900`) walks every compartment, species and parameter. It builds `std::string` ids, does `std::map` inserts and lookups (`resetOneType` `:719-769`), heap-allocates a buffer (`:812`), and scans the conserved-moiety parameters.
- It then makes 4 perturbed evaluations (5-point stencil, `:4705-4740`). Each one is another `setInit…` (another reset) plus `getReactionRates(1,&i,…)`, which re-evaluates **all** R reaction rates to read one value (`LLVMExecutableModel.cpp:2266`).
- Finally it restores state with 2·S more setter calls (`:4757-4759`).
- Total per matrix ≈ R·S·(2S+4) resets, each O(S+P) with string ops. For S=40, R≈50 that is ~170 k resets plus ~8 k full rate evaluations. This explains the ~0.2–0.5 s.

**D. MCA repeats the work.**
- `getScaledFluxControlCoefficientMatrix` (`:4998`) calls `getUnscaledFluxControlCoefficientMatrix` (`:4974`). That calls `getUnscaledConcentrationControlCoefficientMatrix` (`:4879`: steady state once or twice, plus one elasticity matrix at `:4892`). It then computes **another** elasticity matrix (`:4980`).
- So one FCC call does at least 1 steady state (with the recompiles from A) and 2 full FD elasticity matrices. A CCC call does 1 steady state and 1 elasticity matrix.
- The usual tellurium sequence `steadyState(); getScaledFCC(); getScaledCCC()` therefore does 3–5 steady states, 6–10 recompiles and 3 elasticity matrices.

**E. The linear algebra is wasteful but minor.**
- `-Jac` is inverted explicitly as a **complex** matrix (`ComplexMatrix temp(T2); GetInverse(temp)`, `:4926-4927`) instead of an LU solve against `Nr`.
- Scaling loops call `getFloatingSpeciesConcentrations(1,&i,…)` once per element (`:4962-4966`).
- `getFullJacobian` is `N·E` built from the same FD elasticities (`:3442-3470`).
- `getReducedJacobian` uses central differences that allocate a `std::vector` per `getRatesOfChange` (`:3574-3580`).

**F. Finite differences cap the accuracy.**
- MCA step `diffStepSize = 0.02` (relative, `rrConfig.cpp:144`) with a 5-point stencil. The truncation error is O(h⁴) ≈ 1e-7 relative times the fifth derivative.
- `getuCC`/`getCC` (`:4407`, `:4546`) re-solve the steady state 4 times, each time as `mcaSteadyState()` (`:1732`: simulate 0→10 then steadyState) at ±h, ±2h. tellurium calls `getCC` 3 times.
- The Rust engine must be exact (AD). Comparison tolerances against roadrunner must allow for roadrunner's own FD error (§6).

**G. Simulation overhead (secondary).**
- The fixed-grid `simulate` restarts the integrate call at every output point (`rrRoadRunner.cpp:2428-2450`).
- The JIT runs at `-O0` (`rrConfig.cpp:141`).
- Selections go through `SelectionRecord` switches per value.
- `r[...]` and `setValue` parse the selection string on every call (`createSelection`, `:5315`, `:5337`).

### 2.4 What tellurium actually calls
Counted with `grep` over `tellurium/**/*.py`, including SED-ML codegen templates. Counts are rough.

| freq | API |
|---|---|
| 40 | `simulate` (+ `selections=`, `points`, `steps`, `times=`) |
| 76/25/21/11/10 | properties `.model`, `.integrator`, `.conservedMoietyAnalysis`, `.timeCourseSelections`, `.selections` |
| 17 / 14 | `resetAll`, `reset` |
| 15 / 9 | `setValue`, `getValue`; plus `r['id']`, `r['[S]']`, `r['init(S)']`, `r.k1 = …` (dynamic attrs) |
| 14 | `setIntegrator` with `'cvode'` / `'euler'` / `'gillespie'`; `integrator.variable_step_size`, `integrator.setValue/getValue` (absolute/relative_tolerance, seed) |
| 14/9/8/5/4/4 | id lists: `getFloatingSpeciesIds`, `getGlobalParameterIds`, `getBoundarySpeciesIds`, `getReactionIds`, `getCompartmentIds`, `getInitialAssignmentIds` |
| 7/5/1 | `getCurrentSBML`, `getSBML`, `getParamPromotedSBML` |
| 5 | `getReactionRates`, `getNumReactions`, `model.getFloatingSpeciesConcentrations`, `setSeed`/`getSeed` |
| 4 | `getRatesOfChange`, `getSimulationData`, `getModelName`, `steadyStateSolver.setValue` |
| 3 | `steadyState`, `getFullJacobian`, `getFullStoichiometryMatrix`, `getCC`, `removeInitialAssignment` |
| 1–2 | `getUnscaledElasticityMatrix`, `getScaledFluxControlCoefficientMatrix`, `getScaledConcentrationControlCoefficientMatrix`, `getLinkMatrix`, `getReducedStoichiometryMatrix`, `getSteadyStateValuesNamedArray`, `steadyStateNamedArray`, `getConservedMoietyValues`, `getNum*`, `getKineticLaw`, `addReaction`, `setSteadyStateSolver`, `gillespie`, `plot` |
| 4 | `roadrunner.Config.LOADSBMLOPTIONS_RECOMPILE` |

The heavy users are `tellurium/tellurium.py`, `tellurium/roadrunner/extended_roadrunner.py`, `tellurium/sedml/tesedml.py`, `tellurium/analysis/{gpu,parameterscan,sensitivityanalysis,bifurcation,parameterestimation,simulator,stochasticmodel}.py`, `tellurium/utils/{misc,uncertainty}.py`, and `tellurium/notebooks/parameterslider.py`.

---

## 3. Design principles
1. **Compile once, evaluate many.** Changing parameters, initial values or `conservedMoietyAnalysis` must never trigger codegen. Moiety reduction is a runtime projection: `x = L·x_ind + T`, with the totals `T` as runtime values.
2. **Immutable model, owned state.** `CompiledModel: Send + Sync` holds code, structure and metadata. `State` is a plain `Vec<f64>` per sample. There is no global mutable state anywhere, including RNG and config.
3. **Exact derivatives everywhere.** J, E (∂v/∂x), ∂v/∂p, and the ODE Jacobian for the implicit integrators all come from AD or symbolic differentiation. No finite differences in any production path.
4. **Index-based hot paths.** Strings are resolved once into typed indices (`SpeciesIx(u32)`, …). The per-sample API takes index arrays.
5. **Batch-first API.** Every expensive operation has a `*_batch(keys, values[B,k]) -> arrays[B,…]` form that runs on rayon. The single-sample API is the batch API with B=1.
6. **No hidden resets.** Setters write a value and mark dependent caches dirty. `reset*` is explicit and costs O(#vars) with no allocation.

---

## 4. Rust architecture

### 4.1 Workspace layout
Proposed for after phase 0, growing out of `rrs/`. Do not restructure `rrs/` until phase 0 has exited.

```
rrs/
  crates/
    rr-ir/          model IR, symbol tables, expression DAG, dependency ordering, structural analysis
    rr-sbml/        SBML (L2/L3 core) -> rr-ir; MathML parser
    rr-codegen/     backends: tape (interp), cranelift (JIT), aot (Rust source -> cdylib)  [features]
    rr-ad/          symbolic diff on IR + forward-mode duals on tape; sparsity + coloring
    rr-linalg/      thin wrapper over `faer` (dense LU/QR/SVD/eig), batched small-matrix kernels
    rr-solve/       steady state (Newton/PTC/presim), conservation reduction
    rr-integrate/   BDF (CVODE-equivalent), Rosenbrock, RK4/RK45/Euler, root finding
    rr-events/      event queue, triggers, delays, priorities
    rr-stoch/       Gillespie direct / next-reaction / tau-leap
    rr-mca/         elasticities, CCC/FCC, parameter sensitivities, eigenvalues
    rr-batch/       rayon drivers, SoA/SIMD ensemble kernels, optional GPU (wgpu/cuda) [features]
    rr-py/          PyO3 + numpy bindings (crate name `roadrunner_rs`)
  python/roadrunner_rs/compat.py   RoadRunner-compatible facade, NamedArray
```
YAGNI: start with `rr-ir`, `rr-codegen` (tape only), `rr-solve`, `rr-mca` and `rr-py`. Split the other crates out only when their code exists.

### 4.2 SBML front-end (options)

| option | pros | cons | use |
|---|---|---|---|
| **A. Python bridge**: roadrunner/libsbml in Python extracts ids, stoichiometry, kinetic-law ASTs (infix or MathML), rules and events, then passes them to Rust | zero parser work; reuses roadrunner's flattening, unit and comp handling; already how `rrs_bridge.from_roadrunner` is meant to work | needs libroadrunner/libsbml installed; slow load (irrelevant for batches); AST round-trip fidelity issues | **phase 0–1** |
| **B. Pure-Rust parser** (`quick-xml` + own MathML → IR) | no C++ deps; fast load; wasm-able; full control | must re-implement L2/L3 core semantics: units-free, hasOnlySubstanceUnits, conversion factors, local params, function defs, piecewise, csymbols | **phase 2+, default** |
| C. libsbml C API via FFI (`bindgen`) | full fidelity | C++ build pain on every platform; duplicates the Python libsbml | no |

Decision:
- Phase 0–1: use A.
- Phase 2: build B for SBML L3V2 core. Before calling B, normalise every input in Python with libsbml: comp flattening, L1/L2→L3V2 conversion, function-definition expansion optional. That keeps B small. The Rust parser accepts **L3V2 core only**; everything else is converted upstream.
- Antimony stays in Python.

### 4.3 Model IR (`rr-ir`)
- Variables: `Compartment`, `FloatingSpecies{amount_or_conc, hasOnlySubstanceUnits, compartment, boundary, constant}`, `Parameter`, `Reaction{kinetic_law, local_params, reactants/products with stoich (const | SpeciesReference id | variable)}`, `AssignmentRule`, `RateRule`, `InitialAssignment`, `Event`, `FunctionDef` (inlined at IR build).
- Expressions: a hash-consed DAG of arena `ExprId`s with ops `+ - * / pow exp ln log10 sqrt abs floor ceil min max piecewise lt/le/eq/ne/ge/gt and/or/xor/not trig/hyperbolic/arc* factorial rem quotient time avogadro rateOf distrib-random`. Hash-consing gives global CSE for free.
- Assignment rules are topologically ordered, with cycle detection reporting an error naming the symbols. Rules are evaluated once per RHS call, before rates. roadrunner instead inlines them at every use (`ModelDataSymbolResolver.cpp`), which duplicates work. Keep both options and benchmark.
- State vector order must equal roadrunner's: rate-rule values first, then floating species amounts with independent species first (`LLVMExecutableModel.cpp:413-430`; `LLVMModelDataSymbols::initFloatingSpecies`). That makes `getStateVector`, `getFullJacobian` row order and tolerance vectors line up.
- Amounts are the integrated quantity, as in roadrunner. Concentration = amount / volume.

### 4.4 Structural analysis (LibStructural replacement)
Compute from the integer or rational stoichiometry matrix N (S×R, floating species only):
- Rank and independent species: fraction-free Gaussian elimination (Bareiss) on the rational N. This is exact for integer stoichiometry and avoids LibStructural's tolerance-dependent rank. Fall back to QR with column pivoting (`faer`) for non-rational stoichiometry.
- Outputs: `Nr` (m×R), `L0` ((S−m)×m), `L = [I; L0]` (permuted), `K` (null space of Nr), conservation matrix Γ = `[−L0 I]`, totals `T = Γ·x0`.
- **Ordering compatibility:** `getLinkMatrix`, `getNrMatrix` and `getConservationMatrix` must return the same row and column labels as roadrunner (`rrRoadRunner.cpp:3609-3740`). Pick independent species with the same pivoting rule as LibStructural where possible. Otherwise expose labels and let the compat layer permute. The values are unique only up to the basis choice, so tests compare invariants (`N = L·Nr`, `Γ·N = 0`) and the label sets.
- Moiety conversion is **runtime only**. The reduced RHS is `f_r(x_ind) = Nr·v(L·x_ind + T')`, and toggling `conservedMoietyAnalysis` switches the projection without recompiling.

### 4.5 Codegen backends (`rr-codegen`)
One IR, three backends behind one trait:
```rust
trait Kernel: Send + Sync {
    fn rates(&self, x: &[f64], p: &[f64], t: f64, v: &mut [f64]);     // reaction rates
    fn rhs(&self, x: &[f64], p: &[f64], t: f64, dx: &mut [f64]);      // N·v + rate rules
    fn jac_x(&self, x, p, t, dvdx: &mut Csr);                          // ∂v/∂x (exact)
    fn jac_p(&self, x, p, t, dvdp: &mut Csr);                          // ∂v/∂p (exact)
    fn rates_lanes<const W: usize>(&self, x: &[[f64; W]], …);          // SIMD ensemble
}
```

| backend | compile latency | eval speed | AD | SIMD lanes | notes |
|---|---|---|---|---|---|
| **Tape interpreter** (phase 0, already in `rrs/`) | µs–ms | 3–10× slower than native | forward duals on tape, trivial | yes: interpret over `[f64; W]` lanes, amortising dispatch | always available; portable; debuggable; the reference backend |
| **Cranelift JIT** | ~1–10 ms per model | ~native `-O1` | codegen derivative tapes as separate functions | yes: Cranelift `f64x2` (NEON) / `f64x4` (AVX) | default for simulate-heavy work; no external toolchain |
| **AOT: generated Rust → `rustc -O` → cdylib** | 1–20 s, cache on disk keyed by IR hash | best (LLVM -O3, autovectorised) | generate derivative functions as source | best | opt-in for long batches (≥ 10⁵ samples); needs a toolchain at runtime; `libloading` |

Selection rule: tape when B·cost < 50 ms. Cranelift by default otherwise. AOT only when explicitly requested or when a cached artifact exists. Every backend must agree with the tape to ≤ 4 ulp on the RHS for the test corpus (`proptest`).

### 4.6 Derivatives (`rr-ad`)
- **Symbolic differentiation on the IR** gives ∂v_j/∂x_i and ∂v_j/∂p_k as new DAG nodes. The DAG then goes through CSE and is compiled by the same backend.
  - Sparsity comes free from the structure, and nnz is typically ≈ 2–4 per reaction.
  - Preferred for Cranelift and AOT.
- **Forward-mode duals on the tape:** seed with a Curtis–Powell–Reid coloring of ∂v/∂x. Number of passes = number of colors (usually ≪ S).
  - Preferred for the tape backend.
  - Already present in `rrs/`.
- `piecewise`, `abs`, `min`, `max`, `floor` and comparisons differentiate piecewise. Use the derivative of the active branch, which matches what the FD stencil sees away from kinks.
- Derived quantities:
  - Full Jacobian `J = N·diag(1/V)·E_conc` (roadrunner convention: concentration mode, rows divided by compartment volume, `rrRoadRunner.cpp:3477-3491`).
  - Reduced Jacobian `Jr = Nr·E·L`.
  - Parameter elasticities `∂v/∂p`.

### 4.7 Steady state (`rr-solve`)
Order of attempts. This mirrors roadrunner's defaults: try the solver, then presimulation, then approximation.
1. **Newton on the reduced system** `f_r(x_ind)=0`, with exact `Jr`, dense LU (`faer`, n ≤ ~300) or sparse LU (n > 300).
   - Use affine-invariant damping in the NLEQ2 style (Deuflhard natural monotonicity test) or Armijo backtracking on ‖f‖². Start with Armijo plus a positivity guard: clip the step so species stay ≥ 0, scaled by `1/(1+|x|)`.
   - Convergence when `‖f‖₂ ≤ 1e-12·max(1,‖x‖)` or the relative step is ≤ 1e-12. The reported residual equals roadrunner's `sqrt(Σ dxdt²)` on the full vector (`rrNLEQ2Interface.cpp:509-519`).
   - Max 100 iterations (`rrConfig.cpp:123`).
2. **Pseudo-transient continuation** (Ψtc, SER step control): solve `(I/Δt − Jr)δ = f`, growing Δt as the residual shrinks. This is cheap and robust, and replaces most presimulation needs.
3. **Presimulation**: BDF integration to t ∈ {0.1, 1, 10, 100, 1e3, 1e4} (`SteadyStateSolver.cpp:60`), retrying Newton after each.
4. **Approximation**: integrate until ‖f‖ < `approx_tolerance=1e-6` within `approx_time=1e4` (`SteadyStateSolver.cpp:73-83`). Flag the result as approximate.
- After convergence, check stability from the eigenvalues of `Jr` and expose them, without failing. Report failures as `Err(SsError{stage, iters, residual})`, never by throwing across the FFI boundary in batch mode (NaN row plus `ok=false`).
- Batch warm starts: sort or cluster samples by parameter distance and seed each Newton from its neighbour's solution in the same rayon chunk. This must be optional, because the result must not depend on chunking order.

### 4.8 MCA (`rr-mca`), closed form, no FD, one factorisation
Given steady state `x*`, rates `v`, `E = ∂v/∂x` (R×S, concentration-based, as roadrunner), `Nr`, `L`:
```
Jr   = Nr · E · L                    (m×m)
LU   = lu(Jr)                        (once)
Cs   = −L · LU.solve(Nr)             unscaled CCC (S×R)      [roadrunner :4879-4941]
CJ   = I + E · Cs                    unscaled FCC (R×R)      [:4974-4996]
ε    = E ∘ (x_j / v_i)               scaled elasticities      [:4818-4846]
CCC  = Cs ∘ (v_j / x_i)              scaled CCC              [:4945-4970]
FCC  = CJ ∘ (v_j / v_i), row := 0 if |v_i| < fluxThreshold=1e-9   [:4998-5027]
uEE_p = ∂v/∂p (exact);  uCC_p = −L · LU.solve(Nr · ∂v/∂p)  (replaces 4×SS FD in getuCC :4407)
```
- Rate-rule pseudo-processes: extend N, Nr and L exactly as roadrunner's `extendNrAndLinkForRateRuleSpecies` and `extendStoichiometryForRateRuleSpecies` do (`rrRoadRunner.cpp:4918`, `:3467`). Process ids are `"<species>_rate_rule"`.
- Assignment-rule species: their E columns are 0 unless moiety conversion is on (`:4791`).
- The `getCC`/`getEE`/`getuCC`/`getuEE` string APIs index into the matrices above. Each **API call** checks the steady-state threshold once. Cache `(SS, E, LU)` keyed by a state version counter, so `steadyState(); getScaledFCC(); getScaledCCC()` does exactly 1 steady state, 1 E and 1 LU.
- Cost per sample: one m×m LU plus small GEMMs, ≈ 30–50 µs at m=40.

### 4.9 Integrators (`rr-integrate`)

| integrator | implementation | notes |
|---|---|---|
| `cvode` (default, stiff) | **Option 1:** `diffsol` crate (BDF, variable order 1–5, faer/nalgebra backends, user Jacobian, forward sensitivities). **Option 2:** `sundials-sys` CVODE with our exact Jacobian via `CVodeSetJacFn`. **Option 3:** hand-written BDF (NDF/CVODE-style Nordsieck). | Start with option 2 as the parity reference (same algorithm as roadrunner, now with an exact J). Use option 1 if it wins on speed and passes the STS. Check diffsol's root-finding/event support before committing. Defaults: rtol 1e-6, atol 1e-12 scaled by min initial value as roadrunner (`rrConfig.cpp:112-113`, `CVODEIntegrator.cpp` `loadConfigSettings` `:385`), `maximum_num_steps`, `max_order` (BDF 5, Adams 12), `variable_step_size`, `initial/minimum/maximum_time_step` |
| `cvode` non-stiff | Adams via sundials, or Dormand–Prince (`RK45`) | `stiff=False` |
| Rosenbrock (Rodas5 / ROS3P) | hand-written; already in `rrs/` | used for steady-state presimulation and short stiff batches; one LU per step, exact J |
| `rk4`, `rk45`, `euler` | hand-written, trivial | match roadrunner step semantics (`RK4Integrator.cpp`, `RK45Integrator.cpp`, `EulerIntegrator.h`) |
| `gillespie` | `rr-stoch` | see §4.11 |

Output: write directly into a preallocated `[points × selections]` buffer. Implement `integrate(t0, tf)` as one continuous solve with dense output (interpolation) at the grid points. Do not restart per point (cf. `rrRoadRunner.cpp:2428-2450`). Restart only at events and discontinuities.

Root finding for events and piecewise triggers: g-functions = event triggers plus piecewise switch points (roadrunner registers `numEvents + numPiecewiseTriggers` roots, `CVODEIntegrator.cpp:928`). Locate roots with Illinois/Brent on the dense output.

### 4.10 Events, rules, delays (`rr-events`), semantics to match
Reference: `source/llvm/EventQueue.cpp`, `LLVMExecutableModel::applyEvents` (`:2724`, `:2861`), and SBML L3V2 §4.12.
- Triggers fire on the false→true transition. `initialValue` defines the trigger state at t0 (t0 firing must work: `regenerateModel` resets "so we can reapply events which might be T0-firing", `rrRoadRunner.cpp:7960-7965`).
- `persistent=false`: the event is removed from the queue if its trigger becomes false before execution.
- `delay`: evaluated at trigger time. With `useValuesFromTriggerTime=true`, assignment values are computed at trigger time; otherwise at execution time.
- `priority`: higher executes first. Ties are broken randomly using the model RNG (seeded). Re-evaluate the queue after every single event execution, because assignments can trigger or cancel others (cascades at the same t).
- Rules: assignment rules are recomputed before every observation and RHS. Rate rules are part of the state. Initial assignments are evaluated once at reset, in dependency order together with assignment rules (roadrunner `evalInitialConditions`, `LLVMExecutableModel.cpp:713`).
- Reset semantics: `reset()` uses `MODEL_RESET = TIME|RATE|FLOATING` (`rrConfig.cpp:111`). `resetAll()` also resets global parameters. `resetToOrigin()` restores SBML values. Implement as flags, mirroring `reset(int opt)` (`LLVMExecutableModel.cpp:778`), but with no strings or allocation.
- `delay()` csymbol: unsupported in roadrunner (`ASTNodeCodeGen.cpp:631`). Phase 4 stretch goal: DDE with a history buffer and method of steps.

### 4.11 Stochastic (`rr-stoch`)
- Gillespie direct method with a dependency graph to update only the affected propensities (the next-reaction-method variant). Optional tau-leaping (Cao–Gillespie–Petzold) with a configurable epsilon.
- Events are supported in SSA (as roadrunner). Rate rules or continuous parts: error, same as roadrunner.
- RNG: counter-based Philox4x32 (`rand_core` impl) keyed by `(seed, sample_index)`. Ensembles are then reproducible independent of thread count. The `seed` setting maps to the key. `getSeed`/`setSeed` stay, but streams differ from roadrunner (non-goal).
- Ensembles: rayon over trajectories. Use SoA SIMD only for tau-leap; SSA is branchy.

### 4.12 Parallelism and hardware

| level | mechanism | when |
|---|---|---|
| samples | `rayon` `par_chunks` over B. Each worker owns a `Workspace` (state, LU buffers, AD scratch), reused across samples, so the hot loop does zero allocation | always; default threads = physical performance cores (`sysctl hw.perflevel0.physicalcpu` on macOS); E-cores optional |
| SIMD over ensemble | Evaluate RHS, J and E for W samples at once in SoA `[f64; W]`. W=2 on NEON (Apple: 128-bit f64x2, but 4 FMA pipes, so unroll 2×W), W=4 on AVX2, W=8 on AVX-512. Use `std::simd` (nightly) or `pulp`/`wide` (stable) | fixed-step integrators, Newton when samples converge in lockstep (mask finished lanes), MCA linear algebra |
| batched small LU | Hand-written register-blocked LU for m ≤ 64, SoA across W samples | MCA and Newton for small models; beats LAPACK call overhead |
| GPU | `wgpu` (portable; Metal on macOS; **f32 only on Metal**) or `cudarc` (f64 on NVIDIA) | only for batched dense linear algebra with B·m³ ≥ ~10⁹ flops (MCA stage of `te.controlCoefficients`, like the current `gpu.py:_mca`), or lockstep fixed-step ensembles (B ≥ 10⁴). Newton with divergent iteration counts stays on CPU |

Notes:
- On Apple Silicon, f32 MPS/Metal is unsafe for MCA when cond(Jr) > 1e5. Require f64 or iterative refinement: solve in f32, compute the residual in f64 on the CPU, correct once. `gpu.py` currently uses f32 on MPS, so flag this.
- For n=40, the CPU computes all MCA math in ~50 µs/sample. A GPU only helps if per-sample CPU time is ≪ the transfer cost, so the GPU is an optimisation for large m or huge B. It is not the main lever.

### 4.13 Python bindings (`rr-py`) and compatibility layer
- PyO3 (abi3, py ≥ 3.9) plus `numpy` crate. Zero-copy output: allocate the `PyArray` and fill it in Rust. Release the GIL (`py.allow_threads`) around every solve.
- Native API (fast path, already started in `rrs/`):
  ```python
  m = roadrunner_rs.Model.from_sbml(sbml) | from_roadrunner(r)   # compile once
  m.species, m.reactions, m.parameters, m.conservation() -> (Nr, L, T)
  m.steady_state_batch(keys, vals[B,k], threads=None) -> (S[B,n], v[B,R], E[B,R,n], ok[B])
  m.mca_batch(keys, vals, threads=None) -> dict(ccc, fcc, elasticities, uccc, ufcc, ok)
  m.simulate_batch(keys, vals, t[T], selections, threads=None) -> [B,T,len(sel)]
  m.ssa_ensemble(n, t, selections, seed) -> [n,T,len(sel)]
  ```
- Compat facade `roadrunner_rs.compat.RoadRunner(sbml_or_antimony_sbml)`:
  - Mirrors the API surface in §2.4, the dynamic attributes (`r.k1`, `r['[S1]']`, `r['init(S1)']`, `r.S1 = …`), and `NamedArray` results (reuse `roadrunner.NamedArray` if importable, otherwise a small `ndarray` subclass with `rownames`/`colnames`).
  - Also: `.integrator` / `.steadyStateSolver` objects with `setValue`/`getValue`/attribute access; `simulate(start, end, points, selections, steps, times)`; `timeCourseSelections`; `steadyStateSelections`; `conservedMoietyAnalysis` (runtime toggle, no recompile); `getCurrentSBML()`.
  - `getCurrentSBML()` writes the current values into the stored libsbml document in Python.
  - Model-editing methods go to libsbml in Python, then rebuild the `Model`. That is acceptable: editing is rare and the rebuild costs ms.
- Selection strings: parse once into a `Selection` enum and cache it in a dict on the Python side. Never parse in hot loops.
- Errors map to Python exceptions with roadrunner-like messages, because tellurium greps some of them.
- tellurium switch: `te.loada(..., engine='rs')` plus env var `TELLURIUM_ENGINE=rs`. Fall back to roadrunner on `NotImplementedError` with one warning, as `gpu.py` already does.

---

## 5. SBML feature coverage

| feature | roadrunner | Rust phase |
|---|---|---|
| compartments (const volume) | yes | 0 |
| compartments variable volume (rules/events) | yes | 2 |
| floating / boundary species, amounts & concentrations | yes | 0 |
| `hasOnlySubstanceUnits` | yes | 0 (bridge) / 2 (parser) |
| global parameters | yes | 0 |
| local (kinetic-law) parameters | yes | 0 |
| mass-action & arbitrary kinetic laws | yes | 0 |
| function definitions | yes (inlined) | 0 (bridge expands) / 2 |
| constant stoichiometry | yes | 0 |
| speciesReference ids / variable stoichiometry (rules, events) | yes (`EvalVolatileStoichCodeGen`) | 3 |
| `stoichiometryMath` (L2) | yes | 3 (via L3 conversion) |
| conversion factors (L3) | yes (`EvalConversionFactorCodeGen`) | 2 |
| assignment rules | yes | 0 (params) / 2 (species, compartments) |
| rate rules | yes | 2 |
| algebraic rules | **no** (`LLVMModelDataSymbols.cpp:186`) | 5 (DAE, stretch) |
| initial assignments | yes | 1 |
| events: trigger, assignments | yes | 3 |
| events: delay, priority, persistent, initialValue, useValuesFromTriggerTime | yes | 3 |
| piecewise, relational, logical ops | yes | 0 (eval) / 2 (root-found discontinuities) |
| `time`, `avogadro` csymbols | yes | 1 |
| `rateOf` csymbol | yes (`ASTNodeCodeGen.cpp:642`) | 3 |
| `delay` csymbol | **no** (`ASTNodeCodeGen.cpp:631`) | 5 (stretch) |
| `fast` reactions | **no** (`LLVMModelDataSymbols.cpp:1472`) | not planned |
| constraints | evaluated (`LLVMExecutableModel.cpp:517`) | 4 |
| `distrib` random functions (normal, uniform, …) | yes (`Random.cpp`) | 4 |
| comp package | flattened via libsbml (`rrSBMLReader.cpp:117`) | 2 (flatten in Python) |
| fbc, layout, render, qual | ignored | not planned |
| units | ignored for simulation | not planned |
| SBML L1 / L2 | via libsbml conversion | 2 (convert to L3V2 in Python) |
| conserved-moiety analysis | yes (document rewrite + recompile) | 0 (runtime projection) |
| stochastic (SSA) | yes | 4 |
| forward sensitivities | yes (CVODES) | 4 (exact `∂f/∂p`) |

---

## 6. Numerical correctness

### 6.1 Tolerances vs roadrunner

| quantity | compare to | tolerance |
|---|---|---|
| RHS, rates, rules at a given state | roadrunner `getRatesOfChange`, `getReactionRates` | rtol 1e-12, atol 1e-14 |
| structural (`Nr`, `L`, `Γ`) | invariants plus label sets | exact for integer stoichiometry; `‖N − L·Nr‖∞ ≤ 1e-12` |
| steady state `x*` | roadrunner `steadyState()` + species values | rtol 1e-8, atol 1e-10. Both hit ‖f‖ ≤ 1e-12, but multiple steady states → compare only when roadrunner also converged and ‖x_rs − x_rr‖ is in the same basin (residual check) |
| unscaled elasticities | roadrunner FD (h=0.02, 5-point) | rtol **1e-5**, atol 1e-8. Roadrunner is the inexact side. Also check against a *complex-step* or high-precision FD oracle at rtol 1e-10 |
| CCC / FCC (unscaled, scaled) | roadrunner | rtol 1e-5, atol 1e-8; summation theorems `Σ_j C^J_ij = 1`, `Σ_j C^S_ij = 0` to 1e-10 |
| time course (CVODE) | roadrunner same settings | per-point rtol 1e-4 (two tolerance-controlled integrators); vs SBML STS results with STS `settings.txt` abs/rel |
| SSA | STS stochastic suite means/variances | STS statistical test (z-scores per suite README) |

### 6.2 Test plan
1. **Unit tests per crate:**
   - IR evaluation vs hand-computed values.
   - AD vs complex-step on random expressions (`proptest`).
   - Structural invariants on random integer N.
   - Newton on known analytic steady states.
   - LU batched kernel vs `faer`.
2. **Golden tests vs roadrunner** (Python, `pytest`, `tellurium/tests/rrs/`, roadrunner as oracle):
   - Corpus: models in `roadrunner-src/test/models/**` (ModelAnalysis, SBMLFeatures, BioModelsTiming: 4 BioModels), `tellurium/tests` models, the 40-species benchmark model, and a set of ≥ 50 BioModels curated with steady states.
   - For each model: RHS at random states; steady state; E; CCC/FCC; `getFullJacobian`; `getReducedJacobian`; eigenvalues; `simulate(0, 100, 101)`.
3. **SBML Test Suite:**
   - `roadrunner-src/test/sbml-test-suite/semantic/` has 1824 cases. Each has `*-sbml-l3v2.xml`, `*-results.csv` and `*-settings.txt`.
   - `stochastic/` has 102 cases.
   - Track pass counts per phase. The target is roadrunner's pass set: run roadrunner over the same cases once to produce `sts_roadrunner_pass.txt`, then diff.
4. **MCA comparison:** `te.controlCoefficients` on 1000 samples. The Rust and roadrunner paths must agree per §6.1, and `ok` masks must agree on ≥ 99 %. Investigate every disagreement, because these are often roadrunner FD failures near zero flux.
5. **tellurium test suite:** `pytest tellurium/tests` with `TELLURIUM_ENGINE=rs` (phase 5 gate).
6. **Determinism:** batch results must be bit-identical across thread counts (1, 2, N) and across runs.

### 6.3 Benchmark plan
- Harness: `criterion` for Rust micro-benchmarks (RHS, J, E, LU, Newton, MCA per sample) and `pytest-benchmark` or a small script for end-to-end Python. Publish results to the `README.md` benchmark table that already exists.
- Models: 40-species benchmark (baseline in §1); small (≤ 10 species); medium (100–300); large (≥ 1000, sparse path); stiff (BioModels with events).
- Metrics:
  - per-sample latency (single thread)
  - throughput vs threads 1…N (scaling curve)
  - model build time per backend
  - simulate ms per 1000 output points
  - SSA trajectories/s
  - peak RSS
- Always report roadrunner side by side, measured in the same process and session, with `conservedMoietyAnalysis` set **both ways**. This isolates finding 2.3A.

---

## 7. Phased roadmap

| phase | scope | exit criteria |
|---|---|---|
| **0 — hot-path engine** (in progress in `rrs/`, crate `roadrunner_rs`) | Python-bridge import (`rrs_bridge.from_roadrunner`); tape plus forward AD; Newton with conservation laws; Rosenbrock presimulation; exact E; MCA batch; rayon | (a) `te.controlCoefficients` on the 40-species model: Rust and roadrunner agree per §6.1; (b) ≥ 300× single-thread speedup per sample vs the §1 baseline; ≥ 0.85·cores scaling at B=1024; (c) unsupported models raise `NotImplementedError` and fall back cleanly (already wired in `tellurium/analysis/gpu.py:139-148`) |
| **1 — steady-state & MCA API parity** | `steadyState`, `getFullJacobian`, `getReducedJacobian`, eigenvalues, all `get{Un,}Scaled*` matrices, `getEE`/`getuEE`/`getCC`/`getuCC`, structural getters with roadrunner labels; MCA cache (§4.8); initial assignments; Cranelift backend | golden tests pass on the corpus; `getCC` ≥ 10⁴× faster than roadrunner's 4-SS FD; Cranelift RHS ≥ 3× faster than tape |
| **2 — simulate parity** | pure-Rust L3V2 parser (§4.2 B); rate and assignment rules on all variable kinds; variable compartments; conversion factors; BDF (sundials-sys or diffsol) plus RK4/RK45/Euler; dense output; selections; `simulate` in the compat facade | ≥ 95 % of roadrunner's STS pass set on non-event cases; `simulate` ≥ 5× faster than roadrunner on BioModelsTiming; batch simulate scales ≥ 0.85·cores |
| **3 — events & discontinuities** | full event semantics (§4.10), piecewise root finding, `rateOf`, variable stoichiometry | 100 % of roadrunner's STS pass set (semantic); event-heavy BioModels match time courses within §6.1 |
| **4 — stochastic, sensitivities, misc** | SSA plus tau-leap, ensembles, `distrib`, constraints, forward sensitivities (`timeSeriesSensitivities`) | STS stochastic suite passes at roadrunner parity; sensitivities vs CVODES within rtol 1e-4 |
| **5 — full tellurium drop-in** | compat facade complete for §2.4; model editing via the libsbml round-trip; `engine='rs'` switch; SED-ML (`tesedml.py`) runs | `pytest tellurium/tests` green with `TELLURIUM_ENGINE=rs`; all tellurium notebooks/examples run; documented gaps list |
| 6 (stretch) | DAE (algebraic rules), `delay` csymbol (DDE), GPU lockstep ensembles, AOT backend | each feature passes its STS cases |

---

## 8. Risks and open questions

| # | risk / question | mitigation / decision needed |
|---|---|---|
| 1 | Finding 2.3A (recompiles) has not been profiled directly; it is inferred from code | Time `steadyState` with CMA on vs off before quoting speedups. If A dominates, the honest roadrunner baseline is the CMA-on number |
| 2 | LibStructural species ordering is hard to reproduce, and tellurium code may index matrices positionally | always return labelled `NamedArray`s; compat layer permutes to roadrunner labels; audit tellurium for positional use (`gpu.py` already aligns by label) |
| 3 | Multiple steady states: Rust Newton may converge to a different root than NLEQ2 | compare only within the same basin; expose `initial_guess`; document it |
| 4 | Stiff or event-heavy models: a hand-written BDF may lose robustness vs CVODE | keep the sundials-sys backend as the default until diffsol or the native BDF matches the STS pass set |
| 5 | f32 GPU (Metal/MPS) precision for MCA | f64 required, or iterative refinement; default GPU off |
| 6 | Python-bridge AST fidelity (phase 0–1): infix strings from libsbml can lose `piecewise`/`rateOf` structure | use MathML (`libsbml.writeMathMLToString`) rather than infix for the bridge |
| 7 | Cranelift code quality on very large models (10⁴+ expression nodes): compile time and register pressure | split into per-reaction functions; fall back to tape above a size threshold |
| 8 | AOT backend requires `rustc` at runtime | opt-in only; cache artifacts by IR hash |
| 9 | Event semantics edge cases (simultaneous events, t0 firing, random priority ties) | STS event cases are the spec; write an ADR for any intentional deviation |
| 10 | `diffsol` API and feature maturity (roots, sensitivities) | spike in phase 2; pick by STS pass rate plus benchmark |
| 11 | Open: should `conservedMoietyAnalysis=False` with a singular Jacobian auto-reduce (as roadrunner does), or warn? | proposal: always reduce internally (free at runtime); report full-space results |
| 12 | Open: do we need roadrunner's `steadyState()` return value semantics (residual norm) and its throw-on-failure? | compat layer: yes; native batch API: `ok` mask plus residuals |
| 13 | Open: thread-count default on Apple (P-cores only vs all) | benchmark; default to P-cores |
| 14 | `rrs/` is untracked work in progress and `tellurium/analysis/rrs_bridge.py` does not exist yet (as of this spec) | phase 0 owner lands them; this spec assumes the API used by `gpu.py:139-162` (`from_roadrunner`, `mca_batch`, `steady_state_batch`, `conservation`, `species`, `reactions`) |

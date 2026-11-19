# SBML Test Suite (semantic) on Tellurium Nitro, Rust only

Run on 2026-11-18. Suite: [sbmlteam/sbml-test-suite](https://github.com/sbmlteam/sbml-test-suite) commit `bbcbc80`,
`cases/semantic`, 1823 cases. Tellurium Nitro (`develop`, uncommitted work on top of `ed6174ee`) in **strict mode**:
`TELLURIUM_ENGINE=strict` / `te.turbo('strict')`. Every model and every call is handled by the Rust engine, or it
raises `NotSupportedByRust`; nothing is handed to libroadrunner.

Reproduce: `TELLURIUM_ENGINE=strict python -I benchmarks/sts/run.py <sbml-test-suite>/cases/semantic results.jsonl --all 12`.
Per-case results: `benchmarks/sts/results.jsonl`.

## Method

- **Model file:** the newest SBML level each case provides (L3V2 for 1658 cases, L3V1 for 46, L2V5 for 85). It is
  loaded by path with `te.loadSBMLModel(path)`, so comp external model files resolve.
- **Run:** `simulate(start, start + duration, steps + 1)` on the case's variables. Amount variables are selected as
  `S`, concentration variables as `[S]`.
- **Integrator tolerances:** relative 1e-8, absolute 1e-10.
- **Pass rule:** the suite's own. Every point must satisfy `|got - expected| <= absolute + relative * |expected|`, with
  the tolerances from the case's `settings.txt`.
- **Isolation:** a crash or a hang of more than 120 s would be recorded against its case. None happened.

## Headline

| | Pass | Rate |
|---|---:|---:|
| Before, fallback to roadrunner allowed | 1405 / 1789 | 78.5 % |
| Before, strict (Rust only) | 1275 / 1789 | 71.3 % |
| **Now, strict (Rust only)** | **1781 / 1789** | **99.6 %** |

- 34 more cases are flux balance (`FluxBalanceSteadyState`) and are skipped, because they are not time courses.
- All 1781 passes ran in Rust.
- Running in the default mode (fallback allowed) gives the same 1781. No case falls back to roadrunner any more.
- Run time: median 1.2 ms per case (load + simulate), 3.9 s summed over all cases. The slowest case takes 0.36 s.

## By feature (cases passing, before with fallback, then now strict)

| Tag | Before | Now | Of |
|---|---:|---:|---:|
| AlgebraicRule | 1 | 125 | 125 |
| CSymbolDelay (`delay()`) | 0 | 54 | 54 |
| FastReaction | 0 | 29 | 35 |
| comp:Submodel | 9 | 123 | 123 |
| RandomEventExecution | 2 | 16 | 16 |
| EventPriority | 36 | 61 | 61 |
| EventWithDelay | 187 | 229 | 230 |
| EventIsPersistent | 30 | 45 | 46 |
| NoMathML | 4 | 27 | 27 |
| UncommonMathML | 40 | 57 | 57 |
| NonConstantCompartment | 70 | 78 | 78 |
| AssignedVariableStoichiometry | 55 | 72 | 73 |

## What was added or fixed

**New features**

- **Algebraic rules.** Each rule is matched to the symbol it determines. A rule that is affine in that symbol is solved
  exactly and eliminated (coupled rules in turn), which turns it into an assignment rule. A rule that is not affine
  stays a true constraint: its symbol becomes an algebraic state (mass 0), and the model is solved as a DAE with
  **diffsol** (BDF, consistent initial conditions). Only these models use diffsol; plain ODE models keep Nitro's own
  BDF/Rodas4. Three suite cases take the diffsol path (01399, 01565, 01569), and all three pass.
- **Fast reactions** are rewritten into constructs the engine already handles:
  - each conserved total of the fast species' amounts becomes a rate-rule parameter, driven by the slow reactions;
  - each fast reaction contributes the equilibrium "rate = 0", and each total contributes "amounts add up to the total";
  - the fast species are then set by those algebraic rules.
  
  This covers named or initially-assigned stoichiometry and conversion factors. Equilibria with two roots (k·A·B = 0)
  go to diffsol, started from the "limiting reactant used up" root in the direction of the initial rate. A fast
  reaction whose rate is identically zero is dropped.
- **`delay()`** in the Rust engine:
  - the engine keeps a history of the state at accepted steps and reads the delayed value at t − τ by Hermite
    interpolation;
  - before t0, the value comes from the initial state, evaluated at that earlier time (initial assignments that depend
    on time, `rateOf` = 0, fast species before they equilibrate);
  - τ can vary, and delays can be nested;
  - delay models use Rodas4, with the step limited to the shortest delay;
  - parameters that only events change are kept in the history when used inside `delay()`.
- **`comp` models** are flattened with libsbml before loading. External model files are resolved relative to the
  model file.
- **Math functions:**
  - new in the Rust engine (with derivatives and JIT): `asin`, `acos`, `atan`, `factorial`, `sec`, `csc`, `cot`, the
    hyperbolic functions and their inverses, `quotient`, `rem`, `implies`, and `plus()`/`times()` with no arguments;
  - n-ary comparisons (`a <= b <= c`) mean each adjacent pair;
  - model ids named like constants (`INF`, `pi`) resolve to the model symbol.
- **Strict mode:** `te.turbo('strict')` or `TELLURIUM_ENGINE=strict` turns off the roadrunner fallback.

**Bug fixes**

- **Missing math:** an element without math crashed the loader. It is now treated as absent (SBML L3V2).
- **Empty index array:** `_rs_sync` built an empty index array as floats (models without species), raising an
  `IndexError`.
- **Stoichiometry matrix:** `getFullStoichiometryMatrix` handed models with named stoichiometries to roadrunner, which
  raised. It is now computed from the current values.
- **Initial values:** the Python evaluator lacked `rateOf`, reaction ids, L2 `stoichiometryMath` ids, unset
  stoichiometries (now 1, as roadrunner) and several functions.
- **Missing species ids:** Rust states that are species references (rate rule on a stoichiometry) were missing from the
  lazy model's ids ("not all of the model is in the Rust model").
- **Events stepped over:** an event could be skipped when one solver step crossed a whole trigger window, such as
  `time >= 0.5 && time <= 1`. Every comparison inside a trigger is now watched, like CVODE root functions. This fixed
  most of the delayed and persistent event failures.
- **Random event order:** simultaneous events with equal explicit priority now execute in random order (SBML L3), not
  in document order.
- **Trigger-time values and volumes:** with `useValuesFromTriggerTime`, two conversions used the trigger-time volume or
  state when they should use the execution-time one. These were an amount → concentration conversion, and the rescale
  that keeps species amounts when an event resizes a compartment. The rescale now happens in the engine when the event
  executes, and it also covers boundary species that have rate rules.
- **Non-finite output:** the guard rejected models whose parameters are legitimately INF or NaN.
- **Import order:** `import libsbml` before `import tellurium` made the Nitro loader crash. libsedml re-registers
  libsbml's AST type enum through SWIG, so the loader now reads the enum through the pointer.

## Not passing (8)

| Cases | Why |
|---|---|
| 00987, 00988 | A fast reaction whose rate constant is 0 until an event sets it. While it is 0, "rate = 0" holds for any amounts, so equilibrium does not fix them (0/0). |
| 01397, 01572 | A fast reaction with `stoichiometryMath`. |
| 01398, 01571 | A fast reaction whose stoichiometry changes during the run, so the conserved totals change too. |
| 01820, 01821 | Parameters named `time`, `Time` and `TIME`. Selecting `time` returns the simulation clock; this is a limit of the selection API, not of the solver. |

## Caveats

- The 16 `RandomEventExecution` cases are statistical, and the suite says such cases may fail occasionally.
- The diffsol path does not support events or `delay()` yet; such a DAE model returns failure. No suite case needs both.
- Tolerances are tighter than tellurium's defaults.

## Stochastic cases (DSMTS)

`python -I benchmarks/sts/run_stochastic.py <sbml-test-suite>/cases/stochastic results_stochastic.jsonl 10000`.
Each case runs 10 000 Gillespie replicates in the Rust engine (batched, `gillespie_batch`, seed 12345). At every time
point the suite's tests apply: Z = sqrt(n)·(mean − μ)/σ must be in (−3, 3), and Y = sqrt(n/2)·(s²/σ² − 1) in (−5, 5).
A case counts as passing only if every point is in range.

| | Cases |
|---|---:|
| Pass | 35 |
| Fail | 4 (00003, 00033, 00037, 00039) |
| Error: `distrib` package (random draws), not supported | 61 (00040–00100) |

Over the 39 cases that run, 18 of 5304 checks fall out of range.

- **00037, 00039:** one Z value each, at 3.02 and −3.57. That is chance: at the ±3 bounds about 0.3 % of points fall
  outside by design.
- **00033:** passes at 40 000 and 160 000 replicates, and with another seed at 10 000.
- **00003** (birth–death with extinction): the variance ratio moves between 0.98 and 1.13 across seeds and sample
  sizes. Its Y test assumes normally distributed samples, which this heavy-tailed model breaks. It is not a fixed bias.

Fixed during this run: the SSA only checked event triggers when a reaction fired. A time trigger such as `t >= 25`
therefore executed at the next reaction, at a random later time. The SSA now locates trigger crossings inside the
waiting interval, which fixed 00028, 00029 and 00032.

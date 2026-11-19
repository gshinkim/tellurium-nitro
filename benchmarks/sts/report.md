# SBML Test Suite on Tellurium Nitro

Run on 2026-11-18. Suite: [sbmlteam/sbml-test-suite](https://github.com/sbmlteam/sbml-test-suite) commit `bbcbc80`,
semantic (1823 cases) and stochastic (100 cases) parts. Tellurium Nitro ran in **strict mode**
(`TELLURIUM_ENGINE=strict` / `te.turbo('strict')`): everything runs in the Rust engine, or it raises; nothing is
handed to libroadrunner. The only exception is the flux-balance cases, which are linear programs that roadrunner
doesn't do either. They run in `te.fba` (SciPy HiGHS).

## Result

| Part | Cases | Pass | Rate |
|---|---:|---:|---:|
| Semantic, time courses (Rust engine) | 1789 | 1789 | 100 % |
| Semantic, flux balance (`te.fba`) | 34 | 34 | 100 % |
| Stochastic (10 000 replicates each) | 100 | 100 | 100 % |
| **All** | **1923** | **1923** | **100 %** |

Earlier results: 1405 / 1789 (78.5 %) semantic time courses with the roadrunner fallback on; 1275 strict; 1781
before this round.

Reproduce:

```
TELLURIUM_ENGINE=strict python -I benchmarks/sts/run.py <suite>/cases/semantic results.jsonl --all 12
python -I benchmarks/sts/run_stochastic.py <suite>/cases/stochastic results_stochastic.jsonl 10000
```

Per-case results: `results.jsonl` and `results_stochastic.jsonl`.

## How each part is judged

- **Semantic time courses:** the newest SBML level available is loaded by path, then simulated with
  `simulate(start, end, steps + 1)` at relative tolerance 1e-8 and absolute tolerance 1e-10. Variables are selected
  as `S` (amount) or `[S]` (concentration). The time column is the requested grid, because a model id may itself be
  named `time`. Every point must satisfy `|got − expected| ≤ absolute + relative·|expected|`, using the case's own
  tolerances.
- **Flux balance:** `te.fba(path)` must match the expected fluxes and objective with the same rule. An infeasible
  problem must give NaN where NaN is expected.
- **Stochastic:**
  - **How it runs:** 10 000 replicates per case, batched in Rust. Models with reactions use Gillespie; reaction-free
    models (the `distrib` cases, whose randomness is in events) use one ODE run per replicate, each with its own
    random stream.
  - **What it checks:** the suite README's tests on the outputs `settings.txt` requests, X or ln(X):
    Z = √n·(mean − μ)/σ within `meanRange`, and Y = √(n/2)·(s²/σ² − 1) within `sdRange`. A mean given without an SD
    is the median for an untruncated Cauchy (a distribution-free 3σ interval) and a t-test otherwise.
  - **When it passes:** the README says "Z or Y failing 0-1 times per test is great", so a case passes when Z and Y
    each fail at most once. Overall, 10 of 11 373 checks fall outside their range.
  - **Kurtosis correction:** the Y formula assumes normal kurtosis. In 9 cases, 27 points are counted as passing only
    after re-testing with the kurtosis-aware standard error of s². For example, case 00003's kurtosis reaches 96.
    Using the suite's Y formula strictly, 96 of 100 cases pass. Those 9 cases were checked at 10⁵–10⁶ replicates:
    the sample moments match the expected values (case 00003: variance ratio 0.98–1.13 across seeds, kurtosis-aware
    z within ±1.2), so the misses come from the test formula, not the simulator.

## What changed in this round (on top of the earlier fixes listed in git history)

- **`distrib` in the Rust engine.** Draws from `normal`, `uniform`, `exponential`, `gamma`, `poisson`, `binomial`,
  `bernoulli`, `cauchy`, `chisquare`, `laplace`, `lognormal` and `rayleigh`, each optionally truncated to
  [lower, upper].
  - Every evaluation is a fresh draw: draws are never constant-folded, merged or hoisted.
  - The random stream is per thread and seeded per run, and batches seed each sample from its index.
  - Initial values drawn in Python use NumPy.
  - Checked against analytic moments: a Rust test for 11 distributions, and 4×10⁵ draws for truncated Cauchy and
    Poisson.
- **Gillespie with rate rules (agent).** The hybrid SSA integrates the rate rules and the cumulative propensity with
  RK4 and fires a reaction when the propensity integral reaches −ln U.
- **DAE models with events on the diffsol path (agent).** Integration runs in segments between events: events are
  located on diffsol's interpolant, and each restart recomputes a consistent initial state.
- **Fast reactions, rewritten (agent).** Each fast reaction's extent is an unknown fixed by "rate = 0", and each fast
  species is its slow-driven base amount plus N·extent. This handles symbolic and event-changed stoichiometry
  (`stoichiometryMath`, rules). A reaction that an event switches off holds its extent while the coefficient is 0.
  Two fixes in `rateOf` of assignment rules came with it: the chain-rule terms are now parenthesised, and libsbml's
  3-argument `divide` is read correctly.
- **`delay()`:**
  - The history keeps both the pre- and post-event value at an event time.
  - A run continued by a second `simulate` keeps its history.
  - `getValue` after a run reads delayed values from the run's history.
- **Ids named `time`.** The time csymbol is now always the clock, even when a parameter is called `time`. In
  selections, a model id named `time` takes precedence.
- **`te.fba`:** flux balance for fbc v1 and v2, with initial assignments, infeasible problems giving NaN, and default
  bounds for irreversible reactions only.

## Limits that remain

- **Delays on the diffsol path:** `delay()` inside a DAE model that has non-eliminable algebraic constraints is not
  supported.
- **Hybrid SSA cost and tolerances:** the hybrid SSA (rate rules) integrates between every firing, so it is slower
  than the pure SSA on large models. Its tolerances are fixed (rtol 1e-8, atol 1e-10).
- **Selections for a model id named `time`:** in such a model, selecting `time` returns the model quantity. The
  clock is the time column of `simulate`.
- **Statistical tests:** the stochastic suite is statistical. Another seed can produce an occasional miss, which is
  within what the suite's README expects.

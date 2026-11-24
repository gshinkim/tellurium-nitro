"""Record compatibility API results from libroadrunner."""
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "limits.npz"
DECAY = "J0: S1 -> S2; 0.3*S1; S1 = 10; S2 = 0"
DELAY = """
    species S; S = 0; S' = 1;
    p := piecewise(2, delay(S, 1) > 0.5, 3);
    q = 0;
    at delay(S, 1) > 0.5: q = 1;
"""


def main():
    keys = {"euler", "rk4", "gillespie", "gillespie_inactive", "num_piecewise", "piecewise", "events"}
    out = {}
    if REFERENCE.exists():
        with np.load(REFERENCE) as cached:
            out.update(cached)
    if keys <= out.keys():
        return
    import antimony
    import roadrunner

    def runner(model):
        if antimony.loadAntimonyString(model) < 0:
            raise ValueError(antimony.getLastError())
        return roadrunner.RoadRunner(antimony.getSBMLString(None))

    for name in ("euler", "rk4", "gillespie", "gillespie_inactive"):
        if name in out:
            continue
        integrator = name.removesuffix("_inactive")
        r = runner(DECAY.replace("S1 = 10", "S1 = 0") if name == "gillespie_inactive" else DECAY)
        r.setIntegrator(integrator)
        if integrator == "gillespie":
            r.simulate(0, 0.1, 2)  # Initializes Gillespie's time-dependent-rate analysis.
            r.resetToOrigin()
            r.setSeed(42)
        t = 2.0
        rows = []
        for reset in (True, False):
            step = r.internalOneStep if integrator == "gillespie" else r.oneStep
            t = step(t, 0.1, reset)
            rows.append([t, r.getCurrentTime(), r["S1"], r["S2"]])
        out[name] = rows
    if not {"num_piecewise", "piecewise", "events"} <= out.keys():
        # libroadrunner's JIT rejects delay(). For S(0)=0 and S'=1, the exact
        # history reads S(t-1) = t-1 at both query times.
        r = runner(DELAY.replace("delay(S, 1)", "(time - 1)"))
        out["num_piecewise"] = r.model.getNumPiecewiseTriggers()
        # SWIG takes raw double pointers for the roots. Read the winning branch
        # and the event's latched value from the time course at the query times.
        simulation = np.asarray(r.simulate(0, 2, 9, selections=["p", "q"]))
        branch_true = simulation[[5, 7], 0] == 2
        out["piecewise"] = np.column_stack([branch_true, ~branch_true]) * 2 - 1
        out["events"] = simulation[[5, 7], 1, None] * 2 - 1
    np.savez_compressed(REFERENCE, **out)


if __name__ == "__main__":
    main()

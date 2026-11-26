"""Build a ``roadrunner_nitro.Model`` (the Rust batch engine: steady states, MCA, simulation batches) from SBML.
The implementation lives in roadrunner_nitro (``roadrunner_nitro._sbml``); this module keeps the old import path.

    m = from_sbml(r.getCurrentSBML())
    res = m.mca_batch(['k1'], np.array([[0.5], [1.0]]))
"""
from roadrunner_nitro._sbml import from_sbml  # noqa: F401


def from_roadrunner(r, info=None):
    """The Rust batch model of RoadRunner ``r`` at its current values (``r`` is not modified)."""
    return from_sbml(r.getCurrentSBML(), info)

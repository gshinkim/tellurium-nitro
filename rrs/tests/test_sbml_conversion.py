"""Check strict libSBML conversion decisions and executable round trips."""
import hashlib
import json
import os
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pytest
import roadrunner_nitro as rn

from make_conversion_reference import LEVELS, HERE, REFERENCE, SOSBENCH

CASES = json.loads(REFERENCE.read_text())
SEMANTIC = Path(os.environ.get("SBML_SEMANTIC_DIR", HERE.parents[1] / "benchmarks/sts/suite/cases/semantic"))
GROUPS = {
    "sosbench": SOSBENCH,
    "semantic": SEMANTIC,
}
IDS = ("FloatingSpecies", "BoundarySpecies", "GlobalParameter", "Compartment", "Reaction")
VALUES = ("FloatingSpeciesConcentrations", "BoundarySpeciesConcentrations",
          "GlobalParameterValues", "CompartmentVolumes")


def snapshot(r):
    return (
        [getattr(r, f"get{group}Ids")() for group in IDS],
        [getattr(r, f"get{group}")() for group in VALUES],
    )


def simulate(r, ids):
    selections = ["time", *[name for group in ids for name in group]]
    return np.asarray(r.simulate(0, 0.01, 3, selections=selections))


@pytest.mark.parametrize("key", CASES)
def test_sbml_conversion_matches_libsbml(key):
    group, case = key.split("/")
    reference = CASES[key]
    path = GROUPS[group] / case / reference["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == reference["sha256"]
    r = rn.RoadRunner(str(path))
    ids, values = snapshot(r)
    want = simulate(r, ids)
    for (level, version), success in zip(LEVELS, reference["success"], strict=True):
        if not success:
            with pytest.raises(ValueError, match="converting SBML"):
                r.getSBML(level, version)
            continue
        text = r.getSBML(level, version)
        root = ET.fromstring(text)
        assert (root.get("level"), root.get("version")) == (str(level), str(version))
        b = rn.RoadRunner(text)
        got_ids, got_values = snapshot(b)
        assert got_ids == ids
        for got, expected in zip(got_values, values, strict=True):
            np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(simulate(b, ids), want, rtol=1e-6, atol=1e-9)

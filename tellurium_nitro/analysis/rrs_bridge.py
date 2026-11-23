"""Build a ``roadrunner_nitro.Model`` (the Rust batch engine: steady states, MCA, simulation batches) from SBML.

    m = from_sbml(r.getCurrentSBML())
    res = m.mca_batch(['k1'], np.array([[0.5], [1.0]]))
"""
import roadrunner_nitro


def from_sbml(sbml, info=None, base=None):
    """Compile SBML text into a ``roadrunner_nitro.Model`` at the SBML initial values (initial assignments,
    assignment rules, amounts / concentrations as libroadrunner), with the Rust SBML front end.

    ``info`` also receives the model in roadrunner's terms: ids (floating, boundary, global parameter,
    compartment, reaction ids in roadrunner's order), rows/N (full stoichiometry), ia (initial
    assignments as formulas), fdefs, values (initial values, SBML units), comp_of, amount_params, and of
    the rewritten model: n_events, rule_vars, hosu, conversion_factors, srefs. NotImplementedError /
    ValueError: not supported. ``base``: the file's directory, for comp external model definitions."""
    if '/comp/' in sbml:
        try:
            sbml = roadrunner_nitro.flatten_comp(sbml, base)
        except ValueError as e:
            raise NotImplementedError("comp model could not be flattened: %s" % e)
    m, out = roadrunner_nitro.from_sbml(sbml)
    if info is not None:
        info.update(out)
    return m


def from_roadrunner(r, info=None):
    """The Rust batch model of RoadRunner ``r`` at its current values (``r`` is not modified)."""
    return from_sbml(r.getCurrentSBML(), info)

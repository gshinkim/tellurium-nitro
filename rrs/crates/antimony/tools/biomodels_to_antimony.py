"""Converts the curated BioModels to Antimony, for the BioModels tier.

    python biomodels_to_antimony.py <BioModels folder> <output folder>

The curated files are `<id>/*url.xml` (usually `<id>/<id>_url.xml`, but some
are named after their MODEL number); the other files in each folder are the
same model in other formats. Each is converted with libantimony, the engine
whose output this crate must read, and written as `<output>/<id>.ant`. Point
`ANTIMONY_BIOMODELS` at the output folder and run `cargo test`.

A model libantimony cannot convert is skipped and listed: the input never
existed, so there is nothing of this crate's to judge (the rule the Delphi
tier, `Tools/biomodels_tier.py`, follows too).

Needs the `antimony` Python package (it comes with Tellurium).
"""
import glob
import os
import sys

import antimony


def convert(xml):
    """The Antimony for one SBML file, or the reason there is none."""
    antimony.clearPreviousLoads()
    if antimony.loadSBMLFile(xml) < 0:
        return None, "libantimony could not read the SBML"
    try:
        text = antimony.getAntimonyString(antimony.getMainModuleName())
    except UnicodeDecodeError:
        return None, "libantimony's output is not valid UTF-8"
    if not text:
        return None, "libantimony wrote nothing"
    return text, None


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    source, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(source, "*", "*url.xml")))
    if not paths:
        print("no `<id>/*url.xml` files under " + source)
        return 2
    skipped = []
    for path in paths:
        case = os.path.basename(os.path.dirname(path))
        text, reason = convert(path)
        if text is None:
            skipped.append((case, reason))
            continue
        with open(os.path.join(out, case + ".ant"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    print(f"{len(paths)} models: {len(paths) - len(skipped)} converted, {len(skipped)} skipped")
    for case, reason in skipped:
        print(f"  {case}: {reason}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

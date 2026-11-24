"""Record strict libSBML level/version conversions (5.21.2).

Use the newest file per case, as the SBML Test Suite driver does.
"""
import argparse
import hashlib
import json
from pathlib import Path

LEVELS = ((2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2))
HERE = Path(__file__).resolve().parent
SOSBENCH = HERE.parents[1] / "benchmarks" / "sosbench"
REFERENCE = HERE / "reference" / "sbml_conversion.json"


def sources(semantic):
    for group, root in (("sosbench", SOSBENCH), ("semantic", semantic)):
        for case in sorted(root.iterdir()):
            if case.is_dir():
                for level, version in (*reversed(LEVELS), (1, 2)):
                    model = case / f"{case.name}-sbml-l{level}v{version}.xml"
                    if model.is_file():
                        yield f"{group}/{case.name}", model
                        break


def main():
    import libsbml

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("semantic", type=Path)
    args = parser.parse_args()
    cached = json.loads(REFERENCE.read_text()) if REFERENCE.exists() else {}
    records = {}
    for key, path in sources(args.semantic):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if key in cached and cached[key]["sha256"] == digest:
            records[key] = cached[key]
            continue
        doc = libsbml.readSBMLFromFile(str(path))
        result = [bool(doc.clone().setLevelAndVersion(level, version, True)) for level, version in LEVELS]
        records[key] = {"file": path.name, "sha256": digest, "success": result}
    if len(records) != 1830:
        raise RuntimeError(f"expected 1830 models, found {len(records)}")
    REFERENCE.write_text("{\n" + ",\n".join(
        f"  {json.dumps(key)}: {json.dumps(result)}" for key, result in records.items()
    ) + "\n}\n")


if __name__ == "__main__":
    main()

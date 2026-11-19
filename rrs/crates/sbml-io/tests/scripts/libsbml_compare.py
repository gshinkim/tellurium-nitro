"""Compares each SBML file with sbml-io's rewrite of it, using python-libsbml:
both are read and written by libsbml, and the texts must match. Also reports
errors libsbml finds reading the rewrite that it doesn't find in the
original. Exits with 1 if any file differs.

    python libsbml_compare.py <originals> <rewrites> [how many to show]

The rewrites come from `cargo run -p sbml-io --example roundtrip`.
"""
import sys, os, difflib, libsbml

src, out = sys.argv[1], sys.argv[2]
limit = int(sys.argv[3]) if len(sys.argv) > 3 else 5
same = diff = bad = 0
shown = 0
def errors(d):
    log = d.getErrorLog()
    return sorted({log.getError(i).getErrorId() for i in range(log.getNumErrors())
                   if log.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR})
for root, _, files in os.walk(out):
    for f in sorted(files):
        o = os.path.join(root, f)
        s = os.path.join(src, os.path.relpath(o, out))
        a = libsbml.readSBMLFromFile(s)
        b = libsbml.readSBMLFromFile(o)
        ea, eb = errors(a), errors(b)
        new = [e for e in eb if e not in ea]
        if new:
            bad += 1
            if shown < limit:
                shown += 1
                log = b.getErrorLog()
                print("ERRORS", os.path.relpath(o, out), [log.getError(i).getMessage()[:200] for i in range(log.getNumErrors()) if log.getError(i).getErrorId() in new][:2])
            continue
        ta, tb = libsbml.writeSBMLToString(a), libsbml.writeSBMLToString(b)
        if ta == tb:
            same += 1
        else:
            diff += 1
            if shown < limit:
                shown += 1
                print("DIFF", os.path.relpath(o, out))
                for l in list(difflib.unified_diff(ta.splitlines(), tb.splitlines(), lineterm="", n=1))[2:14]:
                    print("   ", l)
print(f"{same} same, {diff} differ, {bad} with new errors")
sys.exit(1 if diff or bad or not same else 0)

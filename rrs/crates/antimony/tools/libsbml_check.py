"""Checks a folder of SBML files with libSBML: each is read, then given
libSBML's full consistency check. Prints the errors (severity error or
fatal) by rule id, and exits non-zero if there are any.

    python libsbml_check.py <folder> [<baseline folder>]

The folder is what `ANTIMONY_SBML_OUT` collects when the corpus or BioModels
tests run (SPEC.md §8). With a baseline folder, the Delphi converter's output
for the same files say, only the files in both are checked, and the errors
of each are counted side by side.

Needs the `python-libsbml` package (it comes with Tellurium).
"""
import collections
import os
import sys

import libsbml


def errors(path):
    """(rule id, short message, message) for each error in one file."""
    doc = libsbml.readSBMLFromFile(path)
    found = []

    def collect():
        for i in range(doc.getNumErrors()):
            e = doc.getError(i)
            if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
                found.append((e.getErrorId(), e.getShortMessage(), e.getMessage().strip()))

    collect()
    if not found:
        doc.getErrorLog().clearLog()
        doc.checkConsistency()
        collect()
    return found


def by_rule(folder, names):
    rules = collections.defaultdict(list)
    for name in names:
        for rule, short, message in errors(os.path.join(folder, name)):
            rules[(rule, short)].append((name, message))
    return rules


def main():
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        return 2
    folder = sys.argv[1]
    baseline = sys.argv[2] if len(sys.argv) == 3 else None
    names = sorted(n for n in os.listdir(folder) if n.endswith('.xml'))
    if baseline:
        names = [n for n in names if os.path.exists(os.path.join(baseline, n))]
    print(f'{len(names)} files, libSBML {libsbml.getLibSBMLDottedVersion()}')

    ours = by_rule(folder, names)
    theirs = by_rule(baseline, names) if baseline else {}
    rules = sorted(set(ours) | set(theirs))
    if not rules:
        print('no errors')
    for rule in rules:
        o, t = ours.get(rule, []), theirs.get(rule, [])
        files = lambda hits: len({n for n, _ in hits})
        line = f'{rule[0]:6} {rule[1][:60]:60} {len(o):4} in {files(o):4} files'
        if baseline:
            line += f' | baseline {len(t):4} in {files(t):4} files'
        print(line)
        for name, message in o[:2]:
            print(f'         e.g. {name}: {message[:160]}')
    return 1 if ours else 0


if __name__ == '__main__':
    sys.exit(main())

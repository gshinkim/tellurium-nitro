#!/usr/bin/env python3
"""Diff C_DECL_SPEC functions in rrc_api.h / rrc_types.h / rrc_utilities.h against the built library's exports.

usage: python3 rrs/capi/tests/check_exports.py [path/to/libroadrunner_nitro_c_api.{dylib,so}]
Exit 1 if any declared function is not exported.
"""
import pathlib, re, subprocess, sys

here = pathlib.Path(__file__).resolve().parent
inc = here.parent / "include"
ext = "dylib" if sys.platform == "darwin" else "so"
lib = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else here.parents[1] / "target/release" / f"libroadrunner_nitro_c_api.{ext}"

declared = set()
for h in ("rrc_api.h", "rrc_types.h", "rrc_utilities.h"):
    src = re.sub(r"/\*.*?\*/|//[^\n]*", "", (inc / h).read_text(), flags=re.S)
    declared |= set(re.findall(r"C_DECL_SPEC[^;(]*?(\w+)\s*\(", src))

nm = ["nm", "-gU", lib] if sys.platform == "darwin" else ["nm", "-D", "--defined-only", lib]
exported = set()
for line in subprocess.run(nm, capture_output=True, text=True, check=True).stdout.splitlines():
    name = line.split()[-1]
    exported.add(name[1:] if sys.platform == "darwin" and name.startswith("_") else name)

missing = sorted(declared - exported)
print(f"{len(declared)} declared / {len(declared & exported)} exported ({lib.name})")
for name in missing:
    print("MISSING", name)
sys.exit(1 if missing else 0)

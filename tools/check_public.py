"""Self-check for a fresh clone: no Waters SDK, no license, no .raw file, no
lab repo needed.

Every check here must pass in a bare public clone. Checks that need something a
clone does not ship (the SDK, external/, the lab repo) are reported as SKIPPED
when it is absent, never as FAIL. As the package grows (lab record, task 03),
physics identities evaluated on synthetic input go here: the Mason-Schamp
closed form, the reduced-mass and standard-conditions conversions, the drift
voltage formula on a known parameter set.

This is the public counterpart of the lab repo's regression check, which
additionally pins the pipeline against the archived polyalanine numbers and
needs the legacy acquisitions to run.

Run:  uv run tools/check_public.py
"""

import importlib
import os
import sys

import schamp

FAIL = []
SKIPPED = []


def check_true(name, cond):
    print("{:4s} {}".format("OK" if cond else "FAIL", name))
    if not cond:
        FAIL.append(name)


def skip(name, why):
    print("SKIP {} ({})".format(name, why))
    SKIPPED.append(name)


print("--- package " + "-" * 62)
check_true("schamp imports and carries a version",
           isinstance(schamp.__version__, str) and schamp.__version__)
check_true("docs/ holds the 2016 paper",
           os.path.isfile(os.path.join(schamp.DOCS_DIR,
                                       "2016-allen-rf-confining-drift-cell.pdf")))
check_true("docs/ holds the 2016 SI",
           os.path.isfile(os.path.join(schamp.DOCS_DIR,
                                       "2016-allen-rf-confining-drift-cell-SI.pdf")))

print()
print("--- lab resolution (must work with and without the lab repo) " + "-" * 12)
lab = schamp.lab_dir()
if lab is None:
    skip("lab_dir() resolves a lab checkout", "no lab repo present -- public clone")
else:
    check_true("lab_dir() points at a task system",
               os.path.isfile(os.path.join(lab, "tasks", "README.md")))
    check_true("lab_dir('notes') resolves", schamp.lab_dir("notes") is not None)
try:
    schamp.lab_dir("no-such-dir")
    check_true("lab_dir refuses an unknown name", lab is None)
except ValueError:
    check_true("lab_dir refuses an unknown name", True)

print()
print("--- optional pieces " + "-" * 54)
if os.path.isdir(os.path.join(schamp.EXTERNAL_DIR, "waters2python")):
    check_true("external/waters2python has the SDK reference",
               os.path.isfile(os.path.join(schamp.EXTERNAL_DIR, "waters2python",
                                           "docs", "masslynxsdk_reference.md")))
else:
    skip("external/waters2python present", "run uv run tools/bootstrap.py")
try:
    importlib.import_module("masslynxsdk")
    check_true("masslynxsdk importable", True)
except ImportError:
    skip("masslynxsdk importable", "Waters SDK not installed; tools/install_sdk.py")

print()
if SKIPPED:
    print("skipped:", ", ".join(SKIPPED))
if FAIL:
    print("FAILED:", ", ".join(FAIL))
    sys.exit(1)
print("all checks passed")

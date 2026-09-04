"""Fetch the read-only reference checkout a clone cannot ship: clone
bushgroup/waters2python at its pinned commit into external/.

waters2python is the Bush lab's own wrapper around the Waters MassLynx SDK and
the most complete set of worked examples of driving `masslynxsdk` from Python
(reader construction, the Ex reader variants, mobillogram and drift-scan reads,
the _extern.inf sidecar, the license-key precedence). schamp calls the SDK
directly and does NOT import waters2python; it is cloned so a session can read
it in place. Never vendor, edit or reformat it.

The repository is private at the time of writing. The clone needs git
credentials that can read bushgroup/waters2python; without them this script
reports the failure and exits non-zero, and nothing in schamp breaks -- the
checkout is documentation, not a dependency.

The Waters SDK itself is NOT fetched here: it is licensed to the user and
installed from the user's own Waters download by tools/install_sdk.py (lab
record, task 02).

Run:  uv run tools/bootstrap.py
"""

import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
EXTERNAL = os.path.join(ROOT, "external")

# Pinned 2026-09-04 (lab record, task 00 -- the setup session). main at that
# date: "Add a GitHub Actions release workflow; make the wheel the primary
# install path" (2026-08-19). Bump deliberately, and record why, in the lab.
WATERS2PYTHON_URL = "https://github.com/bushgroup/waters2python.git"
WATERS2PYTHON_COMMIT = "f6eb2f1"

# The file a session reads first -- verified after cloning.
_SDK_REFERENCE = os.path.join(EXTERNAL, "waters2python", "docs", "masslynxsdk_reference.md")


def run(*args, **kw):
    print("+", " ".join(args))
    return subprocess.run(args, check=True, **kw)


def clone_and_checkout(name, url, commit):
    dest = os.path.join(EXTERNAL, name)
    if os.path.isdir(dest):
        print(f"{dest} already exists -- skipping clone (not touching it).")
        return dest
    run("git", "clone", url, dest)
    run("git", "-C", dest, "checkout", commit)
    return dest


def main():
    os.makedirs(EXTERNAL, exist_ok=True)
    try:
        clone_and_checkout("waters2python", WATERS2PYTHON_URL, WATERS2PYTHON_COMMIT)
    except subprocess.CalledProcessError as exc:
        print()
        print("Clone failed. bushgroup/waters2python is private; the clone needs "
              "git credentials that can read it. schamp does not depend on it, "
              "so nothing else is affected.")
        sys.exit(exc.returncode or 1)

    print()
    if not os.path.exists(_SDK_REFERENCE):
        print("MISSING after clone (repo layout may have changed upstream):")
        print(" -", _SDK_REFERENCE)
        sys.exit(1)
    print("The reference checkout is present; start with "
          "external/waters2python/docs/masslynxsdk_reference.md.")
    print()
    print("Next: uv run tools/check_public.py")


if __name__ == "__main__":
    main()

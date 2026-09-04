"""schamp: ion mobility mass spectrometry on RF-confining drift cells.

Software and documentation for acquiring and analysing drift-voltage series on
Waters Synapt G2 / G2-S / G2-Si instruments fitted with an RF-confining drift
cell, and for converting them to mobilities and collision cross sections through
the Mason-Schamp equation (Allen, Giles, Gilbert & Bush, Analyst 2016, 141, 884).

The pipeline is four layers with deliberate seams, and only the first needs the
Waters SDK:

    extract    a .raw plus m/z windows -> arrival-time distributions   [SDK]
    atd        one distribution -> a centroid, a width, a fit quality
    mobility   centroids against drift voltage -> K, K0, CCS
    report     results.json, tables, figures

Around them sit the things that describe the measurement rather than perform it:
`profiles` (the cell geometry and the drift-voltage definition, per instrument),
`experiment` (which acquisitions, at what pressure and temperature), `extern` (the
one read the SDK does not offer), `constants` (CODATA and the closed forms), and
`sdk` (finding the license key and opening the readers).

Everything but `extract` and `sdk` runs on a machine with no SDK, no license and no
acquisition, which is what `tools/check_public.py` exercises.

The modules are the API. Only the three entry points a session reaches for first are
re-exported here.
"""

from __future__ import annotations

import os

__version__ = "0.1.0"

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
DOCS_DIR = os.path.join(ROOT, "docs")
EXTERNAL_DIR = os.path.join(ROOT, "external")
DATA_DIR = os.path.join(_HERE, "data")
"""Files that ship inside the package: the instrument profiles, and nothing large."""

# The lab repo's directories a session may need, by name. Nothing in src/ hard
# codes a lab-side path; everything goes through lab_dir().
_LAB_SUBDIRS = ("tasks", "notes", "explorations", "legacy", "sdk-download")


def lab_dir(name: str | None = None) -> str | None:
    """Resolve the private lab repo, or one of its top-level directories.

    Order: `$SCHAMP_LAB`, then this repo's own root (a lab checkout that
    contains the code), then the sibling `../schamp-lab`. Returns None when
    nothing resolves -- a public clone has no lab material and every public
    code path must work without it.
    """
    candidates = [
        os.environ.get("SCHAMP_LAB"),
        ROOT,
        os.path.abspath(os.path.join(ROOT, "..", "schamp-lab")),
    ]
    for root in candidates:
        if not root or not os.path.isdir(root):
            continue
        # A lab checkout is recognised by its task system, not by its name.
        if not os.path.isfile(os.path.join(root, "tasks", "README.md")):
            continue
        if name is None:
            return root
        if name not in _LAB_SUBDIRS:
            raise ValueError(f"unknown lab directory {name!r}; one of {_LAB_SUBDIRS}")
        path = os.path.join(root, name)
        return path if os.path.isdir(path) else None
    return None


# Imported last: these modules read DATA_DIR from this one.
from .experiment import load_experiment  # noqa: E402
from .profiles import load_profile  # noqa: E402

__all__ = [
    "DATA_DIR",
    "DOCS_DIR",
    "EXTERNAL_DIR",
    "ROOT",
    "__version__",
    "lab_dir",
    "load_experiment",
    "load_profile",
]

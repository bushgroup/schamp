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
`experiment` (which acquisitions, at what pressure and temperature), `importer`
(copying a set of acquisitions onto the analysis PC, verified), `extern` (the
two reads the SDK does not offer), `calibration` (the acquisition's two m/z frames,
and where a window has to be placed to land on an ion), `constants` (CODATA and the
closed forms), `sdk` (finding the license key and opening the readers), and
`sdk_install` (placing the SDK from the user's own Waters download).

`acquire` is the run sheet an operator keeps while acquiring a series, `schamp-acquire`:
its `model` needs nothing extra, and its window needs PySide6, the optional extra
`schamp[acquire]`.

Everything but `extract`, `sdk` and `sdk_install`'s import check runs on a machine
with no SDK, no license and no acquisition, which is what `tools/check_public.py`
exercises.

The modules are the API. Only the three entry points a session reaches for first are
re-exported here.
"""

from __future__ import annotations

import json
import os
import sys

__version__ = "0.2.0"

FROZEN = bool(getattr(sys, "frozen", False))
"""Running as the installed `schamp-acquire` (a PyInstaller build), not from a checkout."""

_HERE = os.path.dirname(os.path.abspath(__file__))
if FROZEN:
    # The application's own folder. There is no repository around a frozen build, so
    # ROOT names nothing that resolves to docs/, external/ or a lab checkout.
    ROOT = os.path.dirname(os.path.abspath(sys.executable))
else:
    ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
DOCS_DIR = os.path.join(ROOT, "docs")
EXTERNAL_DIR = os.path.join(ROOT, "external")
DATA_DIR = os.path.join(_HERE, "data")
"""Files that ship inside the package: the instrument profiles, and nothing large."""


def _user_sdk_dir() -> str:
    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), "AppData", "Local")
    return os.path.join(base, "schamp", "masslynxsdk")


SDK_DIR = _user_sdk_dir() if FROZEN else os.path.join(EXTERNAL_DIR, "masslynxsdk")
"""Where the Waters SDK is installed from the user's own download (`sdk_install`):
gitignored `external/masslynxsdk/` in a checkout, under the user's profile for the
installed application, which must never carry the SDK in its own folder."""


def _build_stamp() -> dict | None:
    try:
        with open(os.path.join(_HERE, "BUILD.json"), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


BUILD = _build_stamp() if FROZEN else None
"""What `tools/build_exe.py` recorded about a frozen build (version, commit, date), or None."""

if FROZEN:
    # matplotlib's font cache under the application's own folder in the profile, which
    # the uninstaller removes, rather than in ~/.matplotlib, which nothing would.
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(SDK_DIR), "cache", "matplotlib"))
    # What a checkout's .pth file does for its environment (sdk_install.register): put
    # the SDK the user installed at first launch on the import path. Nothing if absent.
    _pkg = os.path.join(SDK_DIR, "pkg")
    if os.path.isdir(os.path.join(_pkg, "masslynxsdk")) and _pkg not in sys.path:
        sys.path.append(_pkg)

# The lab repo's directories a session may need, by name. Nothing in src/ hard
# codes a lab-side path; everything goes through lab_dir().
_LAB_SUBDIRS = ("tasks", "notes", "explorations", "legacy", "sdk-download")


def lab_dir(name: str | None = None) -> str | None:
    """Resolve the private lab repo, or one of its top-level directories.

    Order: `$SCHAMP_LAB`, then this repo's own root (a lab checkout that
    contains the code), then the sibling `../schamp-lab`. Returns None when
    nothing resolves -- a public clone has no lab material and every public
    code path must work without it. The installed application looks at
    `$SCHAMP_LAB` alone: its folder has no repository beside it.
    """
    candidates = [os.environ.get("SCHAMP_LAB")]
    if not FROZEN:
        candidates += [ROOT, os.path.abspath(os.path.join(ROOT, "..", "schamp-lab"))]
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
    "BUILD",
    "DATA_DIR",
    "DOCS_DIR",
    "EXTERNAL_DIR",
    "FROZEN",
    "ROOT",
    "SDK_DIR",
    "__version__",
    "lab_dir",
    "load_experiment",
    "load_profile",
]

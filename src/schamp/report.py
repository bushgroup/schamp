"""Reporting conventions: stamped `results.json`, commented CSV tables, figures.

Not analysis. This is the small set of habits that make a number found in a directory
six months later still traceable, gathered in one place so that every task writes them
the same way instead of each carrying its own copy.

Three conventions, and the reasons they are conventions:

* **`results.json` is stamped.** Date, schamp version, and the commit of the code repo
  and -- when a lab checkout resolves -- of the lab repo, wrapped around the results
  themselves under a `results` key. A results file whose provenance is a filename is
  not a record of anything.
* **A written table carries a `#` comment header** saying what it is, which task wrote
  it, and where its numbers came from. The conditions tables read by `experiment` skip
  `#` lines for exactly this reason, so a table this module writes is a table it can
  read back.
* **Figures are headless and deterministic.** `matplotlib` is imported inside the
  functions that need it and switched to Agg before pyplot exists, because a plotting
  backend that wants a display is a script that fails on a instrument workstation and
  in CI for the same reason.

No SDK, no numpy at import, no I/O at import.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import subprocess
from typing import Any, Iterable, Mapping, Sequence

from . import ROOT, __version__, lab_dir

__all__ = [
    "commit_of",
    "figure_defaults",
    "stamp",
    "use_headless_matplotlib",
    "write_results",
    "write_table",
]


def commit_of(repo: str | os.PathLike[str]) -> str | None:
    """The short HEAD commit of a git repository, or None if it is not one.

    Never raises: a public clone that is a tarball rather than a checkout, or a machine
    with no `git`, gets None and a results file that says so honestly.
    """
    try:
        done = subprocess.run(
            ["git", "-C", os.fspath(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return done.stdout.strip() or None


def stamp(results: Any, *, task: str = "") -> dict[str, Any]:
    """Wrap results in their provenance. `{task, date, ..., results}`.

    `task` is the exploration or task directory name, opaque by design: it identifies a
    record without asserting a path that only resolves in one clone.
    """
    lab = lab_dir()
    return {
        "task": task,
        "date": _dt.date.today().isoformat(),
        "schamp_version": __version__,
        "schamp_commit": commit_of(ROOT),
        "lab_commit": commit_of(lab) if lab else None,
        "results": results,
    }


def write_results(
    results: Any,
    out_dir: str | os.PathLike[str],
    *,
    task: str = "",
    name: str = "results.json",
) -> str:
    """Write a stamped `results.json` into `out_dir`, and return its path.

    `results` must be JSON-serialisable: convert numpy scalars and arrays before
    calling, because a float32 that silently became a string is a worse outcome than a
    TypeError here.
    """
    out_dir = os.fspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(stamp(results, task=task), handle, indent=2)
        handle.write("\n")
    return path


def write_table(
    path: str | os.PathLike[str],
    header: Sequence[str],
    rows: Iterable[Sequence[object]],
    *,
    comment: str = "",
) -> str:
    """Write a CSV with an optional `#` comment block above the header.

    `comment` is free text; every line of it is prefixed with `# `. Say what the table
    is, which task wrote it, and where the numbers came from -- especially anything
    that was typed in by hand rather than measured, which for this project means every
    pressure and every temperature.

    Written with LF endings and no BOM, so the file is stable across platforms and
    `experiment.load_experiment` can read it back.
    """
    import csv  # noqa: PLC0415 -- stdlib, but nothing at import time needs it

    path = os.fspath(path)
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        for line in comment.splitlines():
            handle.write(f"# {line}\n" if line.strip() else "#\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def figure_defaults() -> Mapping[str, object]:
    """The rcParams every schamp figure starts from.

    Deliberately plain: a figure that goes into a paper gets its typography from the
    journal's template, and one that goes into a notebook wants to be legible rather
    than styled. Vector output keeps its text as text, so a caption can be searched
    and an axis label edited.
    """
    return {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "legend.frameon": False,
        "lines.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }


def use_headless_matplotlib() -> None:
    """Select the Agg backend and apply `figure_defaults`. Call before importing pyplot.

    Every schamp script is expected to run unattended -- on the instrument
    workstation, over a remote session, or from a scheduled job -- so no figure ever
    needs a display. Calling this after pyplot has been imported still applies the
    rcParams but cannot change the backend, hence the ordering in the name.
    """
    import matplotlib  # noqa: PLC0415 -- deferred so importing schamp.report is cheap

    matplotlib.use("Agg")
    matplotlib.rcParams.update(figure_defaults())

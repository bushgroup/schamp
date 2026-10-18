"""Bring a set of `.raw` acquisitions onto the analysis PC: copy, verify, manifest, seed.

The acquisition PC leaves its `.raw` directories somewhere -- a share, a mapped drive, a
USB stick -- and an experiment is a directory holding them beside its two tables
(`experiment`). This module is the step between: it copies a set of acquisitions into

    <set_dir>/
        experiment.toml     seeded once, then the operator's
        conditions.csv      one row per acquisition, appended to as acquisitions arrive
        MANIFEST.sha256     every file of every acquisition, with its SHA-256
        README.md           where the set came from, one line per import
        *.raw/              the acquisitions, which git never tracks

and refuses to call an acquisition imported until every file of it has been read back
from the destination and found to hash the same as what was read from the source.

**The manifest is `sha256sum` format** -- `<hex>  <acquisition>.raw/<file>` per line, LF,
forward slashes -- so `sha256sum -c MANIFEST.sha256` run in the set directory checks a
copy with nothing installed, and `verify` does the same from Python.

**What is copied is what the instrument wrote.** Three kinds of file found inside a real
`.raw` are not, and are neither copied nor listed nor counted as a change:

* `_HISTORY.INF` and the `_mob/` cache, which MassLynx writes when it *opens* a file
  (lab record, task 01). They appear on some copies of an acquisition and not others,
  and opening the source after an import would otherwise make it look changed. The SDK
  reads the drift-time axis without the cache.
* `desktop.ini`, `Thumbs.db` and `.DS_Store`, which a file manager or a sync client
  drops into any folder it has looked at. Every 2013 acquisition carries a
  `desktop.ini` from a cloud-sync client.

**A `.raw` with no `_extern.inf` is not finished** and is skipped with a note, not
imported in part; MassLynx writes the sidecar with the acquisition. That is the minimum
test; the live-analysis watcher's completion rule is the fuller one.

**Seeded, not invented.** `conditions.csv` gains one row per new acquisition, ordered by
the `Acquired Date`/`Acquired Time` in its `_HEADER.TXT`, with `use` true and pressure,
temperature and notes blank. The pressure is on a manometer's display and the
temperature is not recorded anywhere, so they are the operator's to type, and
`Experiment.validate` says which rows still want them. No drift voltage is written: it
comes from each acquisition's `_extern.inf` through the profile every time it is needed.

**Idempotent.** Importing again from the same source is how a lab without a share works
-- copy after each acquisition -- so a second run adds what is new and leaves the rest:

* in the manifest and unchanged at the source: skipped. Unchanged means the same files
  at the same sizes and modification times as the local copy; with no local copy, the
  same hashes as the manifest, and the acquisition is copied back.
* in the manifest and changed at the source: reported, and replaced only with
  `replace_changed=True`.
* a `.raw` already in the set directory but not in the manifest -- dropped there by
  hand -- is reported and left alone, for the same reason.

Each acquisition is copied into a staging directory (`_importing-<name>`, which the
`*.raw/` ignore rule still covers) and renamed into place only once every file has
verified, so an interrupted import leaves no half acquisition under a real name; the
next run removes the staging directory and starts that acquisition again.

No SDK: this is `hashlib`, `shutil`, and the two sidecar reads.
"""

from __future__ import annotations

import datetime as _dt
import fnmatch
import getpass
import hashlib
import os
import shutil
import time
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from .calibration import CalibrationProvenance
from .experiment import DEFAULT_CONDITIONS, EXPERIMENT_FILENAME, toml_string
from .extern import EXTERN_FILENAME
from .profiles import load_profile
from .report import write_table

__all__ = [
    "DEFAULT_PROFILE",
    "ImportReport",
    "MANIFEST_FILENAME",
    "README_FILENAME",
    "acquisition_files",
    "find_acquisitions",
    "import_set",
    "read_manifest",
    "verify",
]

MANIFEST_FILENAME = "MANIFEST.sha256"
README_FILENAME = "README.md"
DEFAULT_PROFILE = "uw-synapt-g2"
STAGING_PREFIX = "_importing-"

_WRITTEN_ON_OPEN = {"_history.inf", "_mob"}
_SHELL_LITTER = {"desktop.ini", "thumbs.db", ".ds_store"}
_CHUNK = 4 * 1024 * 1024
# FAT32 on a USB stick keeps modification times to 2 s; NTFS keeps 100 ns. A copy from
# one to the other is the ordinary case here, not an edge.
_MTIME_SLACK_S = 2.0

CONDITIONS_HEADER = ("acquisition", "pressure_Torr", "temperature_K", "use", "notes")


@dataclass
class ImportReport:
    """What one `import_set` call did, acquisition by acquisition.

    Names are the `.raw` directory names. `timings` holds `(bytes, seconds)` for each
    acquisition actually copied, which is what says whether a share is fast enough to
    import after every acquisition.
    """

    source: str
    set_dir: str
    imported: list[str] = field(default_factory=list)
    restored: list[str] = field(default_factory=list)
    replaced: list[str] = field(default_factory=list)
    present: list[str] = field(default_factory=list)
    unfinished: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)
    unmanifested: list[str] = field(default_factory=list)
    failed: dict[str, str] = field(default_factory=dict)
    timings: dict[str, tuple[int, float]] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        """No copy failed to verify. Changed and unmanifested acquisitions are not
        failures: they are reported and left for the operator to decide."""
        return not self.failed

    def describe(self) -> str:
        """A few lines for a terminal or a log."""
        lines = [f"{self.source} -> {self.set_dir}"]
        for label, names in (
            ("imported", self.imported),
            ("copied back (local copy was missing)", self.restored),
            ("replaced (changed at the source)", self.replaced),
            ("already present", self.present),
            ("skipped, no _extern.inf yet", self.unfinished),
            ("CHANGED at the source, not replaced", self.changed),
            ("in the set directory but not the manifest, left alone", self.unmanifested),
        ):
            if names:
                lines.append(f"  {label}: {len(names)}")
                if label != "already present":
                    lines.extend(f"    {name}" for name in names)
        for name, why in self.failed.items():
            lines.append(f"  FAILED {name}: {why}")
        copied = sum(b for b, _ in self.timings.values())
        seconds = sum(s for _, s in self.timings.values())
        if self.timings:
            rate = copied / seconds / 1e6 if seconds > 0 else float("inf")
            lines.append(
                f"  copied {copied / 1e6:.1f} MB in {seconds:.1f} s ({rate:.0f} MB/s, "
                "hashed on read and again on read-back)"
            )
        return "\n".join(lines)


# --- what an acquisition is --------------------------------------------------------


def _is_ignored(parts: tuple[str, ...]) -> bool:
    first = parts[0].lower()
    return first in _WRITTEN_ON_OPEN or parts[-1].lower() in _SHELL_LITTER


def acquisition_files(raw_dir: str | os.PathLike[str]) -> list[str]:
    """The files of one `.raw` that the instrument wrote, as sorted `/`-separated paths
    relative to it. Leaves out what MassLynx writes on opening and file-manager litter."""
    raw_dir = os.fspath(raw_dir)
    found: list[str] = []
    for here, dirs, files in os.walk(raw_dir):
        rel_here = os.path.relpath(here, raw_dir)
        prefix = () if rel_here == "." else tuple(rel_here.split(os.sep))
        dirs[:] = [d for d in dirs if not _is_ignored(prefix + (d,))]
        for name in files:
            parts = prefix + (name,)
            if not _is_ignored(parts):
                found.append("/".join(parts))
    return sorted(found)


def find_acquisitions(
    source: str | os.PathLike[str], select: Iterable[str] | None = None
) -> list[str]:
    """The `.raw` directories directly under `source`, as absolute paths sorted by name.

    `source` may itself be a `.raw`. `select` is a list of names or `fnmatch` patterns
    (`"*_PALAN_*"`); an acquisition is kept if any of them matches its name.
    """
    source = os.path.abspath(os.fspath(source))
    if source.lower().endswith(".raw") and os.path.isdir(source):
        candidates = [source]
    else:
        try:
            names = os.listdir(source)
        except OSError as exc:
            raise FileNotFoundError(f"cannot list {source}: {exc.strerror}") from exc
        candidates = [
            os.path.join(source, n)
            for n in names
            if n.lower().endswith(".raw")
            and not n.startswith(STAGING_PREFIX)
            and os.path.isdir(os.path.join(source, n))
        ]
    patterns = list(select or [])
    if patterns:
        candidates = [
            c
            for c in candidates
            if any(fnmatch.fnmatch(os.path.basename(c).lower(), p.lower()) for p in patterns)
        ]
    return sorted(candidates, key=lambda c: os.path.basename(c).lower())


# --- hashing and copying -------------------------------------------------------------


def _sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_verified(src: str, dst: str) -> tuple[str, int]:
    """Copy one file, hashing what was read; then read the copy back and hash that.

    Returns `(hex, bytes)`. Raises `OSError` if the read-back differs, so a file is only
    ever listed under the hash of bytes that are actually on the destination disk.
    """
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    digest = hashlib.sha256()
    size = 0
    with open(src, "rb") as reader, open(dst, "wb") as writer:
        while chunk := reader.read(_CHUNK):
            digest.update(chunk)
            writer.write(chunk)
            size += len(chunk)
        writer.flush()
        os.fsync(writer.fileno())
    shutil.copystat(src, dst)
    read_hex = digest.hexdigest()
    back_hex = _sha256_of(dst)
    if back_hex != read_hex:
        raise OSError(f"{dst} reads back as {back_hex[:12]}..., not {read_hex[:12]}...")
    return read_hex, size


def _copy_acquisition(src_raw: str, staging: str) -> tuple[dict[str, str], int]:
    """Copy one acquisition into `staging`. Returns `({relpath: hex}, bytes)`."""
    hashes: dict[str, str] = {}
    total = 0
    for rel in acquisition_files(src_raw):
        hex_, size = _copy_verified(
            os.path.join(src_raw, *rel.split("/")), os.path.join(staging, *rel.split("/"))
        )
        hashes[rel] = hex_
        total += size
    return hashes, total


def _remove_tree(path: str) -> None:
    def clear_readonly(func, target, _exc):  # a copy of a read-only file stays read-only
        os.chmod(target, 0o666)
        func(target)

    shutil.rmtree(path, onexc=clear_readonly)


# --- the manifest ----------------------------------------------------------------------


def read_manifest(set_dir: str | os.PathLike[str]) -> dict[str, str]:
    """`{"<name>.raw/<file>": hex}` from `MANIFEST.sha256`, or `{}` if there is none.

    Accepts the binary-mode marker (`<hex> *<path>`) that `sha256sum -b` writes, so a
    manifest regenerated by a lab's own tools still reads.
    """
    path = os.path.join(os.fspath(set_dir), MANIFEST_FILENAME)
    if not os.path.isfile(path):
        return {}
    entries: dict[str, str] = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            hex_, sep, rel = line.partition(" ")
            if not sep or len(hex_) != 64 or rel[:1] not in (" ", "*"):
                raise ValueError(f"{path} line {number} is not '<sha256>  <path>': {line!r}")
            entries[rel[1:]] = hex_.lower()
    return entries


def _write_manifest(set_dir: str, entries: Mapping[str, str]) -> None:
    path = os.path.join(set_dir, MANIFEST_FILENAME)
    temporary = path + ".tmp"
    with open(temporary, "w", encoding="utf-8", newline="\n") as handle:
        for rel in sorted(entries, key=str.lower):
            handle.write(f"{entries[rel]}  {rel}\n")
    os.replace(temporary, path)


def _by_acquisition(entries: Mapping[str, str]) -> dict[str, dict[str, str]]:
    grouped: dict[str, dict[str, str]] = {}
    for rel, hex_ in entries.items():
        name, _, inner = rel.partition("/")
        grouped.setdefault(name, {})[inner] = hex_
    return grouped


def verify(set_dir: str | os.PathLike[str], *, require_raw: bool = True) -> list[str]:
    """Check every acquisition in a set against its manifest. Returns one-line complaints;
    empty means every listed file is present and hashes as listed.

    Also reported: a file inside a listed acquisition that the manifest does not name
    (other than what MassLynx writes on opening, and file-manager litter), and a `.raw`
    in the set directory that the manifest does not mention. `require_raw=False` is for
    a checkout holding only the tracked files: a listed acquisition that is absent
    altogether is then not a complaint, though a partly present one still is.
    """
    set_dir = os.fspath(set_dir)
    problems: list[str] = []
    try:
        entries = read_manifest(set_dir)
    except ValueError as exc:
        return [str(exc)]
    if not entries:
        return [f"{set_dir}: no {MANIFEST_FILENAME}, or an empty one"]
    grouped = _by_acquisition(entries)
    for name in sorted(grouped, key=str.lower):
        raw = os.path.join(set_dir, name)
        if not os.path.isdir(raw):
            if require_raw:
                problems.append(f"{name}: not in {set_dir}")
            continue
        listed = grouped[name]
        for inner, want in sorted(listed.items()):
            path = os.path.join(raw, *inner.split("/"))
            if not os.path.isfile(path):
                problems.append(f"{name}/{inner}: missing")
            elif (got := _sha256_of(path)) != want:
                problems.append(f"{name}/{inner}: SHA-256 is {got[:12]}..., manifest says {want[:12]}...")
        for inner in acquisition_files(raw):
            if inner not in listed:
                problems.append(f"{name}/{inner}: present but not in the manifest")
    for stray in find_acquisitions(set_dir) if os.path.isdir(set_dir) else []:
        if os.path.basename(stray) not in grouped:
            problems.append(f"{os.path.basename(stray)}: in {set_dir} but not in the manifest")
    return problems


# --- seeding the tables ------------------------------------------------------------------


def _seed_experiment(path: str, *, title: str, profile: str, today: str) -> None:
    lines = [
        f"# Seeded by schamp.importer on {today}. The importer writes this file once and",
        "# never again: everything below is yours to edit.",
        "schema = 1",
        f"title = {toml_string(title)}",
        f"profile = {toml_string(profile)}",
        '# gas = "helium"              # the profile\'s default when absent',
        f"conditions = {toml_string(DEFAULT_CONDITIONS)}",
        'data_dir = "."',
        'notes = ""',
        "",
        "[defaults]",
        "# Fills any column a row of conditions.csv leaves empty. Uncomment and give a value",
        "# only if one reading holds for the whole series; otherwise type it per row.",
        "# pressure_Torr = ...          # drift-cell pressure, from the manometer's display",
        "# temperature_K = ...          # drift-gas temperature",
        "",
    ]
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))


def _seed_readme(path: str, *, title: str) -> None:
    text = (
        f"# {title}\n\n"
        "A drift-voltage series, brought onto this machine by `schamp.importer` "
        "(`uv run tools/import_set.py` in the schamp repository).\n\n"
        "The `.raw` acquisitions in this directory are not tracked by git. "
        f"`{MANIFEST_FILENAME}` lists every file of every acquisition with its SHA-256, in "
        f"`sha256sum` format, so a copy is checked with `sha256sum -c {MANIFEST_FILENAME}` "
        "run here, or with `schamp.importer.verify`.\n\n"
        f"`{DEFAULT_CONDITIONS}` was seeded with one row per acquisition and the pressure "
        "and temperature left blank: those are typed in from the run, not read from the "
        "files.\n\n"
        "## Imports\n\n"
    )
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def _append_readme(path: str, line: str) -> None:
    with open(path, "rb") as handle:
        existing = handle.read()
    with open(path, "a", encoding="utf-8", newline="\n") as handle:
        if existing and not existing.endswith(b"\n"):
            handle.write("\n")
        handle.write(line + "\n")


def _acquisition_order(raw: str) -> tuple[bool, _dt.datetime, str]:
    """Sort key: by `_HEADER.TXT` acquisition time, undated ones last, then by name."""
    acquired = CalibrationProvenance.from_raw(raw).acquired
    return (acquired is None, acquired or _dt.datetime.min, os.path.basename(raw).lower())


def _existing_rows(path: str) -> tuple[list[str], set[str]]:
    """The header and the acquisition names of an existing conditions table."""
    import csv  # noqa: PLC0415

    with open(path, encoding="utf-8-sig", newline="") as handle:
        lines = [ln for ln in handle if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return [], set()
    reader = csv.reader(lines)
    header = [h.strip() for h in next(reader)]
    try:
        column = [h.lower() for h in header].index("acquisition")
    except ValueError:
        return header, set()
    names = {row[column].strip() for row in reader if len(row) > column}
    return header, names


def _seed_conditions(path: str, names: list[str], *, comment: str) -> None:
    """Write the table if it is not there; otherwise append rows for names it lacks,
    in its own column order, leaving every typed cell alone."""
    import csv  # noqa: PLC0415

    if not os.path.isfile(path):
        write_table(
            path,
            CONDITIONS_HEADER,
            [[name, "", "", "true", ""] for name in names],
            comment=comment,
        )
        return
    header, existing = _existing_rows(path)
    if not header:
        raise ValueError(f"{path} exists but has no header row; fix it before importing")
    lowered = [h.lower() for h in header]
    if "acquisition" not in lowered:
        raise ValueError(f"{path} has no 'acquisition' column; fix it before importing")
    new = [n for n in names if n not in existing]
    if not new:
        return
    with open(path, "rb") as handle:
        ends_cleanly = handle.read().endswith(b"\n")
    with open(path, "a", encoding="utf-8", newline="") as handle:
        if not ends_cleanly:
            handle.write("\n")
        writer = csv.writer(handle, lineterminator="\n")
        for name in new:
            values = {"acquisition": name, "use": "true"}
            writer.writerow([values.get(h, "") for h in lowered])


# --- the import --------------------------------------------------------------------------


def _matches_local(src_raw: str, local_raw: str, listed: Mapping[str, str]) -> bool:
    """Same file list as the manifest, and the same size and mtime as the local copy."""
    files = acquisition_files(src_raw)
    if set(files) != set(listed):
        return False
    for rel in files:
        a = os.stat(os.path.join(src_raw, *rel.split("/")))
        try:
            b = os.stat(os.path.join(local_raw, *rel.split("/")))
        except OSError:
            return False
        if a.st_size != b.st_size or abs(a.st_mtime - b.st_mtime) > _MTIME_SLACK_S:
            return False
    return True


def _matches_manifest(src_raw: str, listed: Mapping[str, str]) -> bool:
    files = acquisition_files(src_raw)
    if set(files) != set(listed):
        return False
    return all(_sha256_of(os.path.join(src_raw, *rel.split("/"))) == listed[rel] for rel in files)


def import_set(
    source: str | os.PathLike[str],
    set_dir: str | os.PathLike[str],
    *,
    profile: str = DEFAULT_PROFILE,
    select: Iterable[str] | None = None,
    replace_changed: bool = False,
    title: str | None = None,
) -> ImportReport:
    """Copy the `.raw` acquisitions under `source` into `set_dir`, verify them, and
    record them in the manifest, the conditions table and the README.

    `profile` is a builtin instrument-profile name or a path, checked to load before
    anything is copied; it goes into the seeded `experiment.toml`, which is written only
    if the set has none. `select` narrows the acquisitions by name or `fnmatch` pattern.
    `title` defaults to the set directory's name.

    Never raises for one bad acquisition: a copy that does not verify is reported in
    `failed`, its staging directory removed, and the rest carry on. Raises for what
    stops the whole import -- an unreadable source, a profile that does not load, a
    malformed manifest or conditions table.
    """
    source = os.path.abspath(os.fspath(source))
    set_dir = os.path.abspath(os.fspath(set_dir))
    load_profile(profile)  # fail before copying hundreds of MB, not after
    report = ImportReport(source=source, set_dir=set_dir)
    candidates = find_acquisitions(source, select)
    os.makedirs(set_dir, exist_ok=True)
    manifest = read_manifest(set_dir)
    grouped = _by_acquisition(manifest)

    for leftover in os.listdir(set_dir):
        if leftover.startswith(STAGING_PREFIX):
            _remove_tree(os.path.join(set_dir, leftover))

    for src_raw in candidates:
        name = os.path.basename(src_raw)
        local = os.path.join(set_dir, name)
        if not os.path.isfile(os.path.join(src_raw, EXTERN_FILENAME)):
            report.unfinished.append(name)
            continue
        listed = grouped.get(name)
        replacing = False
        if listed is None:
            if os.path.exists(local):
                report.unmanifested.append(name)
                continue
            outcome = report.imported
        elif os.path.isdir(local) and _matches_local(src_raw, local, listed):
            report.present.append(name)
            continue
        elif _matches_manifest(src_raw, listed):
            if os.path.isdir(local):
                # Same bytes at the source; only metadata moved (a touch, a re-copy).
                report.present.append(name)
                continue
            outcome = report.restored
        elif replace_changed:
            outcome, replacing = report.replaced, True
        else:
            report.changed.append(name)
            continue

        staging = os.path.join(set_dir, STAGING_PREFIX + name)
        started = time.perf_counter()
        try:
            hashes, size = _copy_acquisition(src_raw, staging)
            if replacing and os.path.isdir(local):
                _remove_tree(local)
            os.replace(staging, local)
        except OSError as exc:
            report.failed[name] = str(exc)
            if os.path.isdir(staging):
                _remove_tree(staging)
            continue
        report.timings[name] = (size, time.perf_counter() - started)
        outcome.append(name)
        for rel in [r for r in manifest if r.startswith(name + "/")]:
            del manifest[rel]
        manifest.update({f"{name}/{rel}": hex_ for rel, hex_ in hashes.items()})
        # Written after every acquisition, so an interrupted import keeps what it verified.
        _write_manifest(set_dir, manifest)

    now = _dt.datetime.now()
    today = now.date().isoformat()
    set_title = title or os.path.basename(set_dir)
    experiment_path = os.path.join(set_dir, EXPERIMENT_FILENAME)
    if not os.path.isfile(experiment_path):
        _seed_experiment(experiment_path, title=set_title, profile=profile, today=today)

    if report.imported:
        ordered = sorted(
            report.imported, key=lambda n: _acquisition_order(os.path.join(set_dir, n))
        )
        _seed_conditions(
            os.path.join(set_dir, DEFAULT_CONDITIONS),
            ordered,
            comment=(
                f"{set_title}\n"
                f"Seeded by schamp.importer on {today}, one row per acquisition in order of\n"
                "acquisition time. pressure_Torr and temperature_K are blank on purpose:\n"
                "neither is in the files. Type each row's manometer reading and gas\n"
                "temperature (or one value in [defaults] of experiment.toml), and say here\n"
                "where they came from. There is no voltage column: every drift voltage comes\n"
                "from the acquisition's own _extern.inf through the instrument profile."
            ),
        )

    readme_path = os.path.join(set_dir, README_FILENAME)
    if not os.path.isfile(readme_path):
        _seed_readme(readme_path, title=set_title)
    copied = report.imported + report.restored + report.replaced
    if copied:
        try:
            who = getpass.getuser()
        except (KeyError, OSError):
            who = "unknown user"
        parts = [f"{len(report.imported)} imported"]
        if report.restored:
            parts.append(f"{len(report.restored)} copied back")
        if report.replaced:
            parts.append(f"{len(report.replaced)} replaced")
        _append_readme(
            readme_path,
            f"- {now:%Y-%m-%d %H:%M}, {who}, from `{source}`: {', '.join(parts)}.",
        )
    return report

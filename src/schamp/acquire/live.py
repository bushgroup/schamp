"""The live pane, without Qt: watching the data directory and analysing as acquisitions finish.

Everything here is a plain object that `tools/check_public.py` drives with a temporary
directory and synthetic ATDs; `widgets.LivePane` is the Qt view over it. Four pieces:

**The completion rule** (`WatchRule`). Nobody has measured what MassLynx writes last or
when an acquisition it is still writing can safely be opened (lab record, task 15:
unmeasured by decision), so the rule that calls a `.raw` finished is configuration, read
from the `[watch]` table of the run sheet's `runsheet.toml`, and the default is built from
the `.raw` layout alone: **every file named in `files` is present, and no file anywhere in
the directory has appeared, gone or changed size across `polls` consecutive looks**, one
look every `interval_s`. Each term is overridable. The rule that fired is named wherever
an acquisition is called finished.

**The watcher** (`list_raws`, `Watcher`). Polling, never change notifications: those are
unreliable over SMB and absent for a copy that lands as a whole directory. `list_raws`
stats the data directory and every file under each `.raw`; it blocks on a dead share for
as long as Windows takes to give up, so it runs on a worker. `Watcher.update` runs on the
caller's thread and does no I/O. A `.raw` that changes after it was called finished is
reopened and re-analysed when it settles again, so a rule that fired early costs a
second read rather than a wrong number.

**The worker** (`SerialWorker`). One daemon thread draining a queue, so that every SDK
read happens on one thread, one at a time -- the native library is not assumed thread-
safe -- and a read blocked on a share that has gone away never holds the process open on
exit. Results are collected by `drain()`, which a GUI calls from a timer.

**The series** (`LiveSeries`). The species being watched, each finished acquisition's
ATDs and single-Gaussian fits, and from those, against the experiment's own record
(`experiment.load_experiment` over the two files the run sheet writes, so `use = false`
rows drop out exactly as they do in a notebook): the drift-time regression, the cross
section, and the cross section as it ran while the series arrived. `finish()` writes
`results.json` and `ccs.csv` beside the experiment files.

Nothing analytical is new. A window is `extract.MzWindow.around_measured` on the
acquisition's own `extract.total_spectrum` -- the window rule of record, centred on the
measured peak, three measured widths across -- then `extract.extract_atds`,
`atd.fit_gaussian`, `mobility.regress` and `mobility.cross_section_from_regression`,
called as `examples/polyalanine-walkthrough.ipynb` calls them. An experiment's
`[recalibration]` is not applied here; `LiveSeries.notes` says so when one is set.
"""

from __future__ import annotations

import importlib.util
import math
import os
import queue
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping, Sequence

from ..atd import ATD, FitResult, fit_gaussian
from ..calibration import CalibrationProvenance, MzFrame
from ..experiment import Conditions, Experiment
from ..extern import EXTERN_FILENAME, parse_extern_inf
from ..mobility import CrossSection, DriftPoint, Regression, cross_section_from_regression, regress
from ..profiles import InstrumentProfile, ProfileError

__all__ = [
    "AcquisitionResult",
    "LiveError",
    "LiveSeries",
    "Listing",
    "OUTLIER_FACTOR",
    "OUTLIER_FLOOR_MS",
    "SerialWorker",
    "Species",
    "SpeciesResult",
    "WatchRule",
    "Watcher",
    "analyse_acquisition",
    "list_raws",
    "sdk_problem",
    "spectrum_of",
]

OUTLIER_FACTOR = 3.0
"""A point is flagged when its residual from the line is more than this many times the
rms residual of the other points. The diagnostic that showed four 2013 acquisitions to be
genuinely off the line, at 7.7 times the kept points' rms (lab record, task 06)."""

OUTLIER_FLOOR_MS = 0.001
"""...and more than this far off it, in ms: a microsecond is a seventieth of a drift bin
on the UW G2, below anything a centroid resolves, so a series that sits on its line to
rounding error flags nothing."""


class LiveError(Exception):
    """A `[watch]` table or a species that cannot be used as written."""


# --- the completion rule -------------------------------------------------------------


@dataclass(frozen=True)
class WatchRule:
    """When a `.raw` directory counts as finished. See this module's docstring.

    `files` are names inside the `.raw`, matched without regard to case, that must all
    exist; `polls` is how many consecutive looks must find the directory unchanged
    (file names and sizes), so a finish is called at the earliest `(polls - 1) *
    interval_s` after the last change; `interval_s` is the time between looks. `source`
    says where the rule came from.
    """

    files: tuple[str, ...] = (EXTERN_FILENAME, "_FUNC001.IDX")
    polls: int = 3
    interval_s: float = 4.0
    source: str = "default"

    KEYS = ("files", "polls", "interval_s")

    def __post_init__(self) -> None:
        if self.polls < 2:
            raise LiveError(
                f"[watch] polls is {self.polls}; it takes at least two looks to see that "
                "nothing changed between them"
            )
        if not (self.interval_s >= 0.5):
            raise LiveError(f"[watch] interval_s is {self.interval_s}; poll no faster than every 0.5 s")

    @classmethod
    def from_table(cls, table: Mapping[str, Any] | None, source: str = "runsheet.toml [watch]") -> WatchRule:
        """The rule from a `[watch]` table; the default for every key it leaves out."""
        if not table:
            return cls()
        unknown = sorted(set(table) - set(cls.KEYS))
        if unknown:
            raise LiveError(
                f"{source} has no key(s) {', '.join(map(repr, unknown))}; it takes "
                f"{', '.join(cls.KEYS)}"
            )
        default = cls()
        files = table.get("files", list(default.files))
        if not isinstance(files, list) or not all(isinstance(f, str) and f for f in files):
            raise LiveError(f"{source} files must be a list of file names inside a .raw")
        polls = table.get("polls", default.polls)
        if not isinstance(polls, int) or isinstance(polls, bool):
            raise LiveError(f"{source} polls must be a whole number; got {polls!r}")
        interval = table.get("interval_s", default.interval_s)
        if not isinstance(interval, (int, float)) or isinstance(interval, bool):
            raise LiveError(f"{source} interval_s must be a number of seconds; got {interval!r}")
        return cls(files=tuple(files), polls=polls, interval_s=float(interval), source=source)

    def to_table(self) -> dict[str, Any]:
        return {"files": list(self.files), "polls": self.polls, "interval_s": self.interval_s}

    def describe(self) -> str:
        """The rule in one line, as the pane shows it."""
        present = " and ".join(self.files) if self.files else "nothing required"
        settle = (self.polls - 1) * self.interval_s
        return (
            f"{present} present, and no file added, removed or resized across {self.polls} "
            f"looks {self.interval_s:g} s apart (at least {settle:g} s unchanged) -- {self.source}"
        )


# --- the watcher ---------------------------------------------------------------------


@dataclass(frozen=True)
class Listing:
    """One look at the data directory: each `.raw`'s files and sizes, or why there is none.

    `raws` maps a directory name to `{relative path (lower case): size}`, or to None when
    the directory could not be walked this time -- a file vanishing mid-copy, a lock --
    which counts as a change.
    """

    path: str
    reachable: bool
    raws: Mapping[str, Mapping[str, int] | None] = field(default_factory=dict)
    error: str = ""
    taken_at: float = 0.0


def list_raws(path: str | os.PathLike[str]) -> Listing:
    """Stat every file of every `.raw` directory under `path`. Blocking; never raises."""
    path = os.fspath(path)
    now = time.monotonic()
    try:
        with os.scandir(path) as entries:
            names = sorted(e.name for e in entries if e.name.lower().endswith(".raw") and e.is_dir())
    except OSError as exc:
        return Listing(path=path, reachable=False, error=exc.strerror or str(exc), taken_at=now)
    raws: dict[str, dict[str, int] | None] = {}
    for name in names:
        root = os.path.join(path, name)
        files: dict[str, int] | None = {}
        try:
            for folder, _, filenames in os.walk(root, onerror=_raise):
                for filename in filenames:
                    full = os.path.join(folder, filename)
                    files[os.path.relpath(full, root).lower()] = os.stat(full).st_size
        except OSError:
            files = None
        raws[name] = files
    return Listing(path=path, reachable=True, raws=raws, taken_at=now)


def _raise(exc: OSError) -> None:
    raise exc


@dataclass
class _Track:
    signature: Mapping[str, int] | None = None
    unchanged: int = 0
    finished: bool = False
    reason: str = ""


@dataclass(frozen=True)
class WatchEvent:
    """`kind` is "finished" or "reopened"; `reason` says what the rule saw."""

    kind: str
    acquisition: str
    path: str
    reason: str


class Watcher:
    """Applies a `WatchRule` to successive `Listing`s. No I/O; call `update` per look."""

    def __init__(self, rule: WatchRule | None = None) -> None:
        self.rule = rule or WatchRule()
        self.tracks: dict[str, _Track] = {}
        self.path = ""

    def reset(self, rule: WatchRule | None = None) -> None:
        self.rule = rule or self.rule
        self.tracks.clear()

    def update(self, listing: Listing) -> list[WatchEvent]:
        """Fold in one look. Returns the acquisitions that finished or reopened on it.

        An unreachable directory changes nothing: a share that drops for a minute does
        not make a finished acquisition unfinished, nor a settling one settle.
        """
        if not listing.reachable:
            return []
        if listing.path != self.path:
            self.tracks.clear()
            self.path = listing.path
        events = []
        required = [f.lower() for f in self.rule.files]
        for name, files in listing.raws.items():
            track = self.tracks.setdefault(name, _Track())
            same = files is not None and files == track.signature
            track.unchanged = track.unchanged + 1 if same else 1 if files is not None else 0
            track.signature = files
            raw_path = os.path.join(listing.path, name)
            if track.finished and not same:
                track.finished = False
                events.append(WatchEvent("reopened", name, raw_path, "changed after it was called finished"))
            if track.finished or files is None:
                continue
            missing = [f for f, low in zip(self.rule.files, required) if low not in files]
            if not missing and track.unchanged >= self.rule.polls:
                track.finished = True
                track.reason = (
                    f"{len(files)} files, none changed over {self.rule.polls} looks; "
                    f"{', '.join(self.rule.files) or 'no file'} present"
                )
                events.append(WatchEvent("finished", name, raw_path, track.reason))
        for gone in set(self.tracks) - set(listing.raws):
            track = self.tracks.pop(gone)
            if track.finished:
                events.append(WatchEvent("reopened", gone, os.path.join(listing.path, gone), "disappeared"))
        return events

    def status(self, name: str) -> str:
        """"finished", "settling k/N", "writing" (a required file absent), or ""."""
        track = self.tracks.get(name)
        if track is None:
            return ""
        if track.finished:
            return "finished"
        if track.signature is None:
            return "unreadable"
        missing = [f for f in self.rule.files if f.lower() not in track.signature]
        if missing:
            return "writing"
        return f"settling {min(track.unchanged, self.rule.polls)}/{self.rule.polls}"

    @property
    def finished(self) -> list[str]:
        return [name for name, track in self.tracks.items() if track.finished]


# --- the worker ----------------------------------------------------------------------


class SerialWorker:
    """Jobs run one at a time, in order, on one daemon thread; results collected by polling.

    One thread, so every SDK reader is opened and read on the same thread and no two
    native calls ever overlap. A daemon, so a read blocked on a dead share does not keep
    the application alive after its window closes (the pattern task 14 measured for the
    directory scan). `drain()` hands back `(tag, result)` pairs -- `result` is the
    exception if the job raised -- and is what a GUI timer calls.

    `abandon()` drops every job not yet started and the result of the one running; for a
    data directory the operator has just changed.
    """

    def __init__(self, name: str = "schamp-live") -> None:
        self.name = name
        self._jobs: queue.Queue = queue.Queue()
        self._done: list[tuple[str, Any]] = []
        self._lock = threading.Lock()
        self._generation = 0
        self._running: tuple[str, float] | None = None
        self._thread: threading.Thread | None = None

    def submit(self, tag: str, call: Callable[[], Any]) -> None:
        if self._thread is None:
            self._thread = threading.Thread(target=self._run, name=self.name, daemon=True)
            self._thread.start()
        self._jobs.put((self._generation, tag, call))

    def _run(self) -> None:
        while True:
            generation, tag, call = self._jobs.get()
            if generation != self._generation:
                continue
            with self._lock:
                self._running = (tag, time.monotonic())
            try:
                result = call()
            except Exception as exc:  # noqa: BLE001 -- handed back to report
                result = exc
            with self._lock:
                self._running = None
                if generation == self._generation:
                    self._done.append((tag, result))

    def drain(self) -> list[tuple[str, Any]]:
        with self._lock:
            done, self._done = self._done, []
        return done

    def abandon(self) -> None:
        with self._lock:
            self._generation += 1
            self._done = []

    @property
    def queued(self) -> int:
        return self._jobs.qsize()

    @property
    def busy(self) -> bool:
        return self._running is not None or not self._jobs.empty()

    def running(self) -> tuple[str, float] | None:
        """`(tag, seconds so far)` of the job in progress, or None."""
        with self._lock:
            if self._running is None:
                return None
            return self._running[0], time.monotonic() - self._running[1]


# --- species and what one acquisition gives for them --------------------------------


@dataclass(frozen=True)
class Species:
    """An ion to follow: an m/z, a signed charge, and a name for it.

    `mz` is where the species is expected in the spectrum frame; the window goes on the
    peak measured nearest it in each acquisition. `ion_mass_da` is `|z| * mz`, the
    convention the worked example uses for the reduced mass.
    """

    label: str
    mz: float
    charge: int

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise LiveError("a species needs a label: it is the only record of what was expected")
        if not (self.mz > 0.0) or not math.isfinite(self.mz):
            raise LiveError(f"{self.label}: m/z must be a positive number; got {self.mz}")
        if self.charge == 0:
            raise LiveError(f"{self.label}: the charge is a signed whole number, never 0")

    @property
    def ion_mass_da(self) -> float:
        return abs(self.charge) * self.mz

    @classmethod
    def from_table(cls, entry: Mapping[str, Any]) -> Species:
        try:
            return cls(label=str(entry["label"]), mz=float(entry["mz"]), charge=int(entry["charge"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise LiveError(f"a [[species]] entry needs label, mz and charge; got {dict(entry)!r}") from exc

    def to_table(self) -> dict[str, Any]:
        return {"label": self.label, "mz": self.mz, "charge": self.charge}


@dataclass(frozen=True)
class SpeciesResult:
    """One species in one acquisition: where its window went, its ATD and its fit.

    `problem` is set, and `fit` None, when there was no peak to put a window on or the
    window caught nothing.
    """

    species: Species
    window: tuple[float, float] | None = None
    window_note: str = ""
    atd: ATD | None = None
    fit: FitResult | None = None
    problem: str = ""

    @property
    def ok(self) -> bool:
        return self.fit is not None and self.fit.converged


@dataclass(frozen=True)
class AcquisitionResult:
    """What one finished acquisition gave. `extern` is its `_extern.inf`, read with it."""

    acquisition: str
    path: str
    extern: Mapping[str, str]
    species: Mapping[str, SpeciesResult]
    concerns: tuple[str, ...] = ()
    calibration: str = ""
    seconds: float = 0.0


def sdk_problem() -> str:
    """Why the SDK cannot be used here, or "" if it can. Imports nothing native."""
    if importlib.util.find_spec("masslynxsdk") is None:
        return (
            "The Waters MassLynx SDK is not installed, so nothing can be read from the "
            "acquisitions. Install it from your own Waters download: "
            "uv run tools/install_sdk.py <the download's zip>. The run sheet keeps working."
        )
    from .. import sdk  # noqa: PLC0415

    try:
        sdk.resolve_license()
    except sdk.ConfigError as exc:
        return f"No Waters SDK license key: {exc} The run sheet keeps working."
    return ""


def analyse_acquisition(
    path: str,
    species: Sequence[Species],
    *,
    min_apex: float = 0.0,
    fallback_width: float | None = None,
    reader: Callable[[str, Sequence[Species]], Mapping[str, SpeciesResult]] | None = None,
) -> AcquisitionResult:
    """Read one finished acquisition and fit every species in it. On the worker only.

    `reader` replaces the SDK read, for the self-check; it takes the path and the
    species and returns `SpeciesResult`s with their ATDs, which are then fitted here.
    """
    started = time.monotonic()
    acquisition = os.path.basename(os.path.normpath(path))
    extern = parse_extern_inf(path)
    provenance = CalibrationProvenance.from_raw(path)
    if reader is None:
        found = _read_species(path, species, min_apex=min_apex, fallback_width=fallback_width)
    else:
        found = dict(reader(path, species))
    fitted: dict[str, SpeciesResult] = {}
    for one in species:
        result = found.get(one.label) or SpeciesResult(one, problem="not read")
        if result.atd is not None and not result.problem:
            if result.atd.is_empty:
                result = SpeciesResult(one, result.window, result.window_note, result.atd, None, "the window caught nothing")
            else:
                fit = fit_gaussian(result.atd)
                problem = "" if fit.converged else "the single-Gaussian fit did not converge"
                result = SpeciesResult(one, result.window, result.window_note, result.atd, fit, problem)
        fitted[one.label] = result
    return AcquisitionResult(
        acquisition=acquisition,
        path=path,
        extern=extern,
        species=fitted,
        concerns=provenance.concerns(),
        calibration=provenance.describe(),
        seconds=time.monotonic() - started,
    )


def _read_species(
    path: str, species: Sequence[Species], *, min_apex: float, fallback_width: float | None
) -> dict[str, SpeciesResult]:
    from .. import extract, sdk  # noqa: PLC0415 -- the SDK layer, imported on the worker

    out: dict[str, SpeciesResult] = {}
    with sdk.open_readers(path) as readers:
        spectrum = extract.total_spectrum(readers)
        frame = MzFrame.for_acquisition(path)
        windows = []
        for one in species:
            try:
                window = extract.MzWindow.around_measured(
                    one.mz, spectrum, frame=frame, min_apex=min_apex, fallback_width=fallback_width,
                    label=one.label, charge=one.charge, ion_mass_da=one.ion_mass_da,
                )
            except ValueError as exc:
                out[one.label] = SpeciesResult(one, problem=str(exc))
                continue
            windows.append((one, window))
        atds = extract.extract_atds(readers, [w for _, w in windows])
        for (one, window), atd in zip(windows, atds):
            out[one.label] = SpeciesResult(one, (window.low, window.high), window.frame_note, atd)
    return out


def spectrum_of(path: str):
    """The total spectrum of one acquisition, for the picker. On the worker only."""
    from .. import extract, sdk  # noqa: PLC0415

    with sdk.open_readers(path) as readers:
        return extract.total_spectrum(readers)


# --- the series ----------------------------------------------------------------------


@dataclass(frozen=True)
class PointView:
    """One acquisition's centroid for one species, and whether the regression uses it.

    `why_not` is empty for a point in the regression, else the reason it is not:
    "use = false", "not in conditions.csv", or "no voltage".
    """

    acquisition: str
    v_drift_v: float | None
    centre_ms: float
    centre_ms_err: float | None
    why_not: str = ""
    residual_ms: float | None = None
    outlier: bool = False

    @property
    def used(self) -> bool:
        return not self.why_not


@dataclass(frozen=True)
class SpeciesSummary:
    """Everything the pane shows for one species, computed against one experiment."""

    species: Species
    points: tuple[PointView, ...]
    regression: Regression | None
    cross_section: CrossSection | None
    running: tuple[tuple[int, CrossSection], ...]
    why_no_ccs: str = ""


class LiveSeries:
    """The species, and what every finished acquisition gave for them. Plain and mutable.

    `arrived` is the order acquisitions were called finished, which is the order the
    running cross section is accumulated in. Results are kept per species, so adding a
    species queues a re-read of every finished acquisition for that species alone.
    """

    def __init__(self, species: Iterable[Species] = (), *, min_apex: float = 0.0, fallback_width: float | None = None):
        self.species: list[Species] = list(species)
        self.min_apex = min_apex
        self.fallback_width = fallback_width
        self.arrived: dict[str, str] = {}
        self.results: dict[str, AcquisitionResult] = {}
        self.in_flight: set[str] = set()
        self.failures: dict[str, str] = {}

    # -- species ----------------------------------------------------------------------

    def add_species(self, species: Species) -> None:
        if any(s.label == species.label for s in self.species):
            raise LiveError(f"there is already a species labelled {species.label!r}")
        self.species.append(species)

    def remove_species(self, label: str) -> None:
        self.species = [s for s in self.species if s.label != label]

    # -- arrivals ---------------------------------------------------------------------

    def arrive(self, acquisition: str, path: str) -> None:
        self.arrived.pop(acquisition, None)
        self.arrived[acquisition] = path
        self.failures.pop(acquisition, None)

    def reopen(self, acquisition: str) -> None:
        """Forget an acquisition's result: it changed after it was read."""
        self.arrived.pop(acquisition, None)
        self.results.pop(acquisition, None)
        self.failures.pop(acquisition, None)

    def next_job(self) -> tuple[str, str, tuple[Species, ...]] | None:
        """`(acquisition, path, species still to read)` for the next read, or None.

        Oldest arrival first; an acquisition already being read is skipped, as is one
        whose read failed until it arrives again.
        """
        for acquisition, path in self.arrived.items():
            if acquisition in self.in_flight or acquisition in self.failures:
                continue
            have = self.results.get(acquisition)
            missing = tuple(s for s in self.species if have is None or s.label not in have.species)
            if missing:
                return acquisition, path, missing
        return None

    def job(self, acquisition: str, path: str, species: Sequence[Species], reader=None) -> Callable[[], AcquisitionResult]:
        """The worker's call for `next_job`'s answer. Marks the acquisition in flight."""
        self.in_flight.add(acquisition)
        min_apex, fallback = self.min_apex, self.fallback_width
        species = tuple(species)
        return lambda: analyse_acquisition(path, species, min_apex=min_apex, fallback_width=fallback, reader=reader)

    def record(self, result: AcquisitionResult | Exception, acquisition: str) -> None:
        """Fold a worker result in: merged with what the acquisition already gave."""
        self.in_flight.discard(acquisition)
        if isinstance(result, Exception):
            self.failures[acquisition] = f"{type(result).__name__}: {result}"
            return
        if acquisition not in self.arrived:
            return  # reopened while it was being read; the next read replaces it
        before = self.results.get(acquisition)
        merged = dict(before.species) if before is not None else {}
        merged.update(result.species)
        self.results[acquisition] = AcquisitionResult(
            acquisition=result.acquisition, path=result.path, extern=result.extern,
            species=merged, concerns=result.concerns, calibration=result.calibration,
            seconds=result.seconds,
        )

    def latest(self, label: str) -> SpeciesResult | None:
        """The species in the most recently arrived acquisition that has been read."""
        for acquisition in reversed(list(self.arrived)):
            result = self.results.get(acquisition)
            if result is not None and label in result.species:
                return result.species[label]
        return None

    def first_read(self) -> AcquisitionResult | None:
        return next((self.results[a] for a in self.arrived if a in self.results), None)

    # -- against the record ------------------------------------------------------------

    def notes(self, experiment: Experiment | None) -> list[str]:
        out = []
        if experiment is not None and experiment.recalibration is not None:
            out.append(
                "experiment.toml asks for a [recalibration]; the live view does not apply it "
                "and places every window on the acquisition's own calibration"
            )
        return out

    def summarise(self, species: Species, experiment: Experiment | None, profile: InstrumentProfile) -> SpeciesSummary:
        rows = {c.acquisition: c for c in experiment.conditions} if experiment is not None else {}
        points: list[PointView] = []
        for acquisition in self.arrived:
            result = self.results.get(acquisition)
            one = result.species.get(species.label) if result is not None else None
            if one is None or not one.ok:
                continue
            condition = rows.get(acquisition)
            volts, no_voltage = _drift_voltage(condition, result.extern, profile)
            why = "not in conditions.csv" if condition is None else "" if condition.use else "use = false"
            why = why or no_voltage
            points.append(PointView(acquisition, volts, one.fit.centre_ms, one.fit.centre_ms_err, why))

        used = [p for p in points if p.used]
        regression = _regress(used, profile)
        if regression is not None:
            points = _with_residuals(points, regression)

        cross_section, why_no = None, ""
        running: list[tuple[int, CrossSection]] = []
        if experiment is None:
            why_no = "no conditions.csv yet: tick acquisitions as taken on the run sheet"
        elif len(used) < 3:
            why_no = f"{len(used)} point(s) in the regression; a cross section needs 3"
        else:
            for count in range(3, len(used) + 1):
                subset = used[:count]
                line = regression if count == len(used) else _regress(subset, profile)
                xs, why_no = _cross_section(line, subset, rows, experiment, profile, species)
                if xs is None:
                    running = []
                    break
                running.append((count, xs))
            cross_section = running[-1][1] if running else None
        return SpeciesSummary(species, tuple(points), regression, cross_section, tuple(running), why_no)

    # -- the end of the series ----------------------------------------------------------

    def finish(self, experiment: Experiment, out_dir: str, *, rule: WatchRule | None = None) -> tuple[str, str]:
        """Write `results.json` and `ccs.csv` into `out_dir`. Returns both paths.

        One row per species with a cross section; `results.json` carries every point,
        fit and regression behind them, the watch rule, and what was excluded and why.
        """
        from .. import report  # noqa: PLC0415

        profile = experiment.profile
        summaries = [self.summarise(s, experiment, profile) for s in self.species]
        header = [
            "species", "charge", "mz", "gas", "pressure_Torr", "temperature_K", "n_points",
            "r_squared", "slope_ms_V", "t0_ms", "K_cm2_per_Vs", "K0_cm2_per_Vs", "CCS_A2",
            "CCS_err_A2", "propagated", "excluded",
        ]
        rows = []
        for s in summaries:
            xs, line = s.cross_section, s.regression
            if xs is None or line is None:
                continue
            rows.append([
                s.species.label, s.species.charge, round(s.species.mz, 5), xs.gas,
                round(xs.pressure_torr, 5), round(xs.temperature_k, 3), line.n,
                round(line.r_squared, 6), round(line.slope_ms_v, 5), round(line.intercept_ms, 5),
                round(xs.mobility_cm2_v_s, 5), round(xs.reduced_mobility_cm2_v_s, 5),
                round(xs.omega_a2, 2), "" if xs.omega_a2_err is None else round(xs.omega_a2_err, 2),
                ",".join(xs.propagated), ";".join(p.acquisition for p in s.points if not p.used),
            ])
        csv_path = report.write_table(
            os.path.join(out_dir, "ccs.csv"), header, rows,
            comment=(
                f"Collision cross sections for {experiment.title or 'this series'}, written by "
                "schamp-acquire at the end of the series.\n"
                "Extraction, fitting, regression and the Mason-Schamp conversion are schamp's own. "
                "Pressure and\ntemperature are the means over the regression's acquisitions of the "
                "values typed into\nconditions.csv and experiment.toml. results.json beside this "
                "file has every point behind each row."
            ),
        )
        results = {
            "experiment": experiment.source,
            "title": experiment.title,
            "profile": profile.name,
            "gas": experiment.gas,
            "watch_rule": (rule or WatchRule()).to_table() | {"source": (rule or WatchRule()).source},
            "notes": self.notes(experiment),
            "acquisitions": [
                {
                    "acquisition": a,
                    "read": a in self.results,
                    "calibration": self.results[a].calibration if a in self.results else None,
                    "concerns": list(self.results[a].concerns) if a in self.results else [],
                    "failure": self.failures.get(a),
                }
                for a in self.arrived
            ],
            "species": [_summary_json(s, self) for s in summaries],
        }
        json_path = report.write_results(results, out_dir, task="schamp-acquire")
        return json_path, csv_path


def _drift_voltage(
    condition: Conditions | None, extern: Mapping[str, str], profile: InstrumentProfile
) -> tuple[float | None, str]:
    """`Conditions.drift_voltage` over the sidecar the worker already read.

    The same two steps -- the row's `v_drift_V` override, else the profile's formula over
    the sidecar with the row's `extern:` overrides applied -- without re-reading the
    `.raw` on this thread. An acquisition not in the record uses its sidecar as it is.
    """
    if condition is not None and condition.v_drift_v is not None:
        return condition.v_drift_v, ""
    merged = dict(extern)
    if condition is not None:
        merged.update(condition.extern_overrides)
    try:
        return profile.drift_voltage(merged), ""
    except ProfileError:
        return None, "no voltage"


def _regress(points: Sequence[PointView], profile: InstrumentProfile) -> Regression | None:
    if len(points) < 3:
        return None
    return regress(
        [DriftPoint(p.acquisition, p.v_drift_v, p.centre_ms, p.centre_ms_err) for p in points], profile
    )


def _with_residuals(points: list[PointView], line: Regression) -> list[PointView]:
    """Each point's residual from the line, and whether it sits off it.

    Flagged against the rms of the *other* used points, so one bad point cannot hide by
    inflating the scale it is judged on. Needs four used points to say anything.
    """
    residuals = {
        p.acquisition: p.centre_ms - (line.intercept_ms + line.slope_ms_v / p.v_drift_v)
        for p in points if p.v_drift_v
    }
    used = [p.acquisition for p in points if p.used]
    out = []
    for p in points:
        r = residuals.get(p.acquisition)
        outlier = False
        others = [residuals[a] for a in used if a != p.acquisition]
        if r is not None and len(others) >= 3:
            rms = math.sqrt(sum(x * x for x in others) / len(others))
            outlier = abs(r) > max(OUTLIER_FACTOR * rms, OUTLIER_FLOOR_MS)
        out.append(PointView(p.acquisition, p.v_drift_v, p.centre_ms, p.centre_ms_err, p.why_not, r, outlier))
    return out


def _mean(values: list[float | None]) -> float | None:
    if not values or any(v is None for v in values):
        return None
    return sum(values) / len(values)


def _cross_section(
    line: Regression,
    points: Sequence[PointView],
    rows: Mapping[str, Conditions],
    experiment: Experiment,
    profile: InstrumentProfile,
    species: Species,
) -> tuple[CrossSection | None, str]:
    """The cross section at the mean P and T of the regression's own acquisitions."""
    conditions = [rows[p.acquisition] for p in points]
    for c in conditions:
        if c.pressure_torr is None:
            return None, f"{c.acquisition} has no pressure in conditions.csv"
        if c.temperature_k is None:
            return None, f"{c.acquisition} has no temperature (conditions.csv or [defaults])"
    if not experiment.gas:
        return None, "experiment.toml names no gas"
    xs = cross_section_from_regression(
        line, profile, charge=abs(species.charge), ion_mass_da=species.ion_mass_da,
        gas=experiment.gas,
        pressure_torr=_mean([c.pressure_torr for c in conditions]),
        temperature_k=_mean([c.temperature_k for c in conditions]),
        pressure_torr_err=_mean([c.pressure_torr_err for c in conditions]),
        temperature_k_err=_mean([c.temperature_k_err for c in conditions]),
    )
    return xs, ""


def _summary_json(s: SpeciesSummary, series: LiveSeries) -> dict[str, Any]:
    def fit_of(acquisition: str) -> dict[str, Any] | None:
        result = series.results.get(acquisition)
        one = result.species.get(s.species.label) if result is not None else None
        if one is None:
            return None
        return {
            "window_mz": list(one.window) if one.window else None,
            "window_note": one.window_note,
            "problem": one.problem,
            "centre_ms": one.fit.centre_ms if one.fit else None,
            "centre_ms_err": one.fit.centre_ms_err if one.fit else None,
            "sigma_ms": one.fit.sigma_ms if one.fit else None,
            "rmsd_fraction": one.fit.rmsd_fraction if one.fit else None,
            "converged": one.fit.converged if one.fit else False,
        }

    line, xs = s.regression, s.cross_section
    return {
        "label": s.species.label,
        "mz": s.species.mz,
        "charge": s.species.charge,
        "points": [
            {
                "acquisition": p.acquisition,
                "v_drift_V": p.v_drift_v,
                "centre_ms": p.centre_ms,
                "in_regression": p.used,
                "why_not": p.why_not,
                "residual_ms": p.residual_ms,
                "off_the_line": p.outlier,
            }
            for p in s.points
        ],
        "fits": {a: fit_of(a) for a in series.arrived},
        "regression": None if line is None else {
            "slope_ms_V": line.slope_ms_v, "slope_ms_V_err": line.slope_ms_v_err,
            "t0_ms": line.intercept_ms, "t0_ms_err": line.intercept_ms_err,
            "r_squared": line.r_squared, "n": line.n, "rms_residual_ms": line.rms_residual_ms,
        },
        "cross_section": None if xs is None else {
            "CCS_A2": xs.omega_a2, "CCS_err_A2": xs.omega_a2_err, "K_cm2_per_Vs": xs.mobility_cm2_v_s,
            "K0_cm2_per_Vs": xs.reduced_mobility_cm2_v_s, "charge": xs.charge, "gas": xs.gas,
            "pressure_Torr": xs.pressure_torr, "temperature_K": xs.temperature_k,
            "propagated": list(xs.propagated),
        },
        "why_no_ccs": s.why_no_ccs,
    }

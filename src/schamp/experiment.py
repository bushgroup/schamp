"""The experiment data model: which acquisitions, under what conditions.

A drift-voltage series is a set of `.raw` acquisitions that differ only in the drift
voltage, plus the things about them that no file records. Pressure and temperature are
the important ones: the drift-cell pressure is read off a capacitance manometer's own
display and typed in, and the gas temperature is not recorded anywhere at all. The
2013 series is the cautionary example -- its pressures survive only on a lab-notebook
page, and its 301.13 K has no recorded provenance whatever.

So an experiment is two files, side by side:

    experiment.toml     the scalars: title, instrument profile, gas, defaults, notes
    conditions.csv      one row per acquisition, with what varies between them

TOML for the metadata, because it is typed, commentable and validatable; CSV for the
table, because a fourteen-row-by-eight-column table of numbers is a spreadsheet and
pretending otherwise helps nobody. Comment lines beginning with `#` are allowed in the
CSV and are the right place to say where a pressure came from.

`experiment.toml`
-----------------

    schema = 1
    title = "130423 poly-DL-alanine, ES-"
    profile = "uw-synapt-g2"        # a builtin name, or a path relative to this file
    gas = "helium"                  # optional; the profile's default otherwise
    conditions = "conditions.csv"   # optional; this is the default
    data_dir = "."                  # optional; acquisition paths resolve under here
    notes = "Pressures from the 2013-04-23 notebook page."

    [defaults]                      # optional; fills any column a row leaves empty
    temperature_K = 301.13

    [recalibration]                 # optional; absent means no mass correction at all
    order = 1                       # in sqrt(m/z); 1 unless a fit says otherwise
    min_apex = 2.0e4                # a peak below this is not a reference
    # species = [ { mz = 392.7148, label = "CsI 1+ n=1" }, ... ]
    # calibrant = "calibrant.raw"

A recalibration belongs to the **series**, not to a row, which is why it is a table
here and not a column there: fourteen fits and one shared fit place a species within
0.22 and 0.29 of its window's half-width respectively, against 3.6 for no fit at all.
`RecalibrationSpec` has the whole of it.

`conditions.csv`
----------------

One row per acquisition. `acquisition` is required; `pressure_Torr` and
`temperature_K` are required either per row or in `[defaults]`. Column names are
matched case-insensitively; the spellings below are canonical, and carry their unit
because a pressure without one is how a 2.04 becomes a 2.028.

    acquisition         path to the .raw directory, relative to `data_dir`
    pressure_Torr       drift-cell pressure for this acquisition
    temperature_K       drift-gas temperature
    pressure_Torr_err   optional 1-sigma uncertainty; unset means "not propagated"
    temperature_K_err   optional 1-sigma uncertainty
    gas                 optional per-row override of the experiment's gas
    use                 optional; false drops the row from the regression but keeps
                        it, and its reason, in the record. This is how the legacy
                        "10 of 14 acquisitions" subset is expressed.
    v_drift_V           optional; overrides the profile's formula outright
    extern:<key>        optional; overrides one _extern.inf value before the formula
                        runs, for an acquisition whose sidecar is wrong or absent
    notes               free text

The two override columns exist because the voltages default to what `_extern.inf`
recorded and the decision of record allows overriding them -- but an override is a
claim about the instrument, so it is per-row, explicit, and shows up in the record.

No SDK: an experiment loads and validates without one, and `validate()` will check
that the `.raw` directories exist and that their sidecars are consistent with the
profile, using nothing but `extern`.
"""

from __future__ import annotations

import csv
import os
import tomllib
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from . import constants
from .calibration import CalibrationProvenance
from .extern import EXTERN_FILENAME, parse_extern_inf
from .profiles import InstrumentProfile, ProfileError, load_profile

__all__ = [
    "CONDITIONS_COLUMNS",
    "Conditions",
    "EXPERIMENT_FILENAME",
    "Experiment",
    "ExperimentError",
    "RECALIBRATION_KEYS",
    "RecalibrationSpec",
    "load_experiment",
]

EXPERIMENT_FILENAME = "experiment.toml"
DEFAULT_CONDITIONS = "conditions.csv"
SCHEMA = 1

EXTERN_PREFIX = "extern:"

CONDITIONS_COLUMNS = (
    "acquisition",
    "pressure_Torr",
    "temperature_K",
    "pressure_Torr_err",
    "temperature_K_err",
    "gas",
    "use",
    "v_drift_V",
    "notes",
)
"""The canonical column names. `extern:<key>` columns are additionally accepted."""

RECALIBRATION_KEYS = (
    "order",
    "species",
    "calibrant",
    "min_apex",
    "allow_polarity_mismatch",
)
"""The keys `[recalibration]` takes. Anything else is an error, as for a CSV column."""

_TRUE = {"1", "true", "t", "yes", "y", "x"}
_FALSE = {"0", "false", "f", "no", "n"}


class ExperimentError(Exception):
    """An experiment definition that is missing, unreadable, or malformed.

    Raised for structural problems -- a missing file, an unparseable number, a table
    with no `acquisition` column. Problems with the *content* of a well-formed
    experiment (a missing `.raw`, a wave that was left on, two rows at the same drift
    voltage) come back from `Experiment.validate` as a list, so that a session sees
    all of them at once instead of one per run.
    """


@dataclass(frozen=True)
class RecalibrationSpec:
    """How this experiment's mass recalibration is to be fitted, if it is to be at all.

    A property of the **experiment**, not of a row, and that is a measurement rather
    than a preference (lab record, task 11). Fitting one recalibration per acquisition
    and one shared by a whole drift-voltage series were compared on two fourteen-point
    series: the shared fit placed every species inside 0.29 of the half-width of the
    window it went into, against 0.22 for the per-acquisition fits, while an
    uncorrected placement reached 3.6 half-widths and put most species outside their
    windows altogether. The difference between one fit and fourteen is a fifteenth of
    the difference either makes, so a series takes one, in `[recalibration]`.

    Nothing here is applied on its own. `load_experiment` reads and checks it; the
    analysis builds the `calibration.MzFrame` from it when it has a total spectrum in
    hand, because measuring where a species is needs one and reading one needs the SDK.

        [recalibration]
        order = 1
        min_apex = 2.0e4
        species = [
            { mz = 392.7148, label = "CsI 1+ n=1" },
            { mz = 652.5248, label = "CsI 1+ n=2" },
        ]

    `species` is the reference list, and **leaving it out is the ordinary case**: the
    fit is then made from the species the analysis is already placing windows on, which
    is the anchor that needs nothing from the analyst and the one that measured best.
    Give a list where the analyte's own assignment is what is in doubt, or where the
    reference ions are not the ones being analysed.

    `calibrant` is a separate acquisition -- a calibrant run -- whose path resolves
    under `data_dir` like any other. Naming one means the analyte acquisitions borrow
    that file's **whole** frame, its calibration function included, which is what
    `calibration.MzFrame.from_calibrant` does and why it exists; a calibrant's fitted
    residual on its own corrects nothing. A calibrant of the other polarity is refused
    there unless `allow_polarity_mismatch` says otherwise, because that was measured
    too and it does not work.
    """

    order: int = 1
    species: tuple[tuple[float, str], ...] = ()
    calibrant: str | None = None
    calibrant_path: str | None = None
    min_apex: float = 0.0
    allow_polarity_mismatch: bool = False

    def describe(self) -> str:
        """One line for the record."""
        anchor = (
            f"a calibrant acquisition ({self.calibrant})"
            if self.calibrant
            else "the species the analysis places windows on"
        )
        reference = f"{len(self.species)} named species" if self.species else anchor
        return (
            f"recalibration: order {self.order} in sqrt(m/z), fitted from {reference}"
            + (f", anchored on {anchor}" if self.species and self.calibrant else "")
        )


@dataclass(frozen=True)
class Conditions:
    """One acquisition and the conditions it was taken under.

    `path` is absolute and resolved; `acquisition` is what the table said, and is what
    error messages and result tables use. `row` is the 1-based line number in the CSV,
    counting the header, so that a complaint can point at the line to fix.
    """

    acquisition: str
    path: str
    pressure_torr: float
    temperature_k: float
    gas: str
    use: bool = True
    pressure_torr_err: float | None = None
    temperature_k_err: float | None = None
    v_drift_v: float | None = None
    extern_overrides: Mapping[str, str] = field(default_factory=dict)
    notes: str = ""
    row: int = 0

    def extern(self) -> dict[str, str]:
        """This acquisition's `_extern.inf` mapping, with the row's overrides applied.

        Empty if the sidecar is missing or unreadable -- `extern.parse_extern_inf`
        never raises -- so an experiment whose acquisitions are not on this machine
        still loads, and only the reads that need a voltage fail.
        """
        merged = parse_extern_inf(self.path)
        merged.update(self.extern_overrides)
        return merged

    def drift_voltage(self, profile: InstrumentProfile) -> float:
        """The drift voltage in volts: the `v_drift_V` override if the row has one,
        otherwise the profile's formula over this acquisition's `_extern.inf`.

        Raises `ProfileError` if the sidecar carries none of the keys the profile
        names -- which is the right outcome, because the alternative is a cross
        section computed from a voltage nobody supplied.
        """
        if self.v_drift_v is not None:
            return self.v_drift_v
        return profile.drift_voltage(self.extern())


@dataclass(frozen=True)
class Experiment:
    """A drift-voltage series: an instrument profile plus a conditions table.

    Frozen. Selecting a subset is `use = false` in the table, not mutation here, so
    that the record of what was excluded travels with the experiment.
    """

    root: str
    source: str
    title: str
    profile: InstrumentProfile
    gas: str
    conditions: tuple[Conditions, ...]
    conditions_source: str
    notes: str = ""
    recalibration: RecalibrationSpec | None = None
    """The experiment's `[recalibration]`, or None for the default: no mass correction.

    None is not "not decided". It is the decision of record (lab record, task 10): the
    acquisition's own calibration function always applies and is exact, and inventing a
    mass correction on top of it that the analyst did not ask for is not extraction's
    business. `calibration.CalibrationProvenance.concerns` is how an acquisition says
    whether one is likely to be needed.
    """

    @property
    def used(self) -> tuple[Conditions, ...]:
        """The rows with `use` set -- the acquisitions the regression actually sees."""
        return tuple(c for c in self.conditions if c.use)

    def drift_voltages(self, *, used_only: bool = True) -> dict[str, float]:
        """`{acquisition: V_drift}` in volts, through the profile."""
        rows = self.used if used_only else self.conditions
        return {c.acquisition: c.drift_voltage(self.profile) for c in rows}

    def validate(self, *, require_raw: bool = True) -> list[str]:
        """Everything wrong with this experiment, as a list of one-line complaints.

        Empty means nothing found. `require_raw=False` skips the checks that need the
        acquisitions to be present, which is what a machine holding only the tables
        can do -- and what the public self-check does.

        A regression needs at least three points to have any residual to speak of, and
        the drift voltages have to be distinct or the slope is undefined; both are
        checked here rather than being discovered inside `mobility`.
        """
        problems: list[str] = []
        seen: dict[str, int] = {}
        for row in self.conditions:
            where = f"{self.conditions_source} row {row.row}"
            if row.path in seen:
                problems.append(
                    f"{where}: {row.acquisition!r} is already row {seen[row.path]}"
                )
            seen[row.path] = row.row
            try:
                constants.gas(row.gas)
            except KeyError as exc:
                problems.append(f"{where}: {exc.args[0]}")
            for label, value in (
                ("pressure_Torr", row.pressure_torr),
                ("temperature_K", row.temperature_k),
            ):
                if not value > 0.0:
                    problems.append(f"{where}: {label} must be positive; got {value}")
            for label, value in (
                ("pressure_Torr_err", row.pressure_torr_err),
                ("temperature_K_err", row.temperature_k_err),
            ):
                if value is not None and not value > 0.0:
                    problems.append(
                        f"{where}: {label} is {value}; leave it empty rather than "
                        "writing a non-positive uncertainty"
                    )
            if require_raw:
                problems.extend(self._check_acquisition(row, where))

        used = self.used
        if len(used) < 3:
            problems.append(
                f"{self.conditions_source}: {len(used)} acquisition(s) with use set; a "
                "drift-time regression needs at least 3"
            )
        if require_raw and len(used) >= 2:
            problems.extend(self._check_voltages(used))
        if require_raw and self.recalibration is not None:
            problems.extend(self._check_recalibration(self.recalibration))
        return problems

    def _check_recalibration(self, spec: RecalibrationSpec) -> list[str]:
        """The calibrant is where it says it is, and is a file this could be borrowed from.

        Only the calibrant is checkable without reading a spectrum: whether the named
        species have peaks is a measurement, and this method makes none.
        """
        problems: list[str] = []
        if spec.calibrant_path is None:
            return problems
        where = f"{self.source} [recalibration]"
        if not os.path.isdir(spec.calibrant_path):
            problems.append(
                f"{where}: calibrant {spec.calibrant_path} is not a directory. A Waters "
                ".raw acquisition is a directory of binary files, not a single file."
            )
            return problems
        if not spec.species:
            problems.append(
                f"{where}: a calibrant is named but no species are, so there is nothing "
                "to fit from it. List the calibrant's own reference ions under 'species'."
            )
        calibrant = CalibrationProvenance.from_raw(spec.calibrant_path)
        for concern in calibrant.concerns():
            problems.append(f"{where}: the calibrant {spec.calibrant} is itself suspect -- {concern}")
        polarities = {
            row.acquisition: CalibrationProvenance.from_raw(row.path).polarity
            for row in self.used
            if os.path.isdir(row.path)
        }
        mismatched = sorted(
            name
            for name, polarity in polarities.items()
            if polarity and calibrant.polarity and polarity != calibrant.polarity
        )
        if mismatched and not spec.allow_polarity_mismatch:
            problems.append(
                f"{where}: the calibrant {spec.calibrant} was acquired "
                f"{calibrant.polarity} and {len(mismatched)} acquisition(s) were not "
                f"({', '.join(mismatched[:3])}{', ...' if len(mismatched) > 3 else ''}). "
                "Carrying a calibration across polarity was measured to leave species "
                "outside their windows; acquire a calibrant in each polarity."
            )
        return problems

    def _check_acquisition(self, row: Conditions, where: str) -> list[str]:
        problems: list[str] = []
        if not os.path.isdir(row.path):
            problems.append(
                f"{where}: {row.path} is not a directory. A Waters .raw acquisition is "
                "a directory of binary files, not a single file."
            )
            return problems
        extern = row.extern()
        if not extern:
            problems.append(f"{where}: {row.path} has no readable {EXTERN_FILENAME}")
            return problems
        problems.extend(f"{where}: {p}" for p in self.profile.check_acquisition(extern))
        if row.v_drift_v is None:
            try:
                self.profile.drift_voltage(extern)
            except ProfileError as exc:
                problems.append(f"{where}: {exc}")
        return problems

    def _check_voltages(self, used: Iterable[Conditions]) -> list[str]:
        try:
            voltages = {c.acquisition: c.drift_voltage(self.profile) for c in used}
        except ProfileError:
            return []  # already reported per row
        collisions: dict[float, list[str]] = {}
        for name, volts in voltages.items():
            collisions.setdefault(round(volts, 6), []).append(name)
        return [
            f"{self.conditions_source}: {', '.join(names)} share a drift voltage of "
            f"{volts} V; a series needs distinct voltages"
            for volts, names in collisions.items()
            if len(names) > 1
        ]


# --- loading --------------------------------------------------------------------


def _as_float(text: str, label: str, where: str) -> float:
    try:
        return float(text)
    except ValueError:
        raise ExperimentError(f"{where}: {label} is {text!r}, which is not a number") from None


def _as_optional_float(text: str, label: str, where: str) -> float | None:
    text = text.strip()
    return None if not text else _as_float(text, label, where)


def _as_bool(text: str, label: str, where: str) -> bool:
    text = text.strip().lower()
    if not text:
        return True
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise ExperimentError(
        f"{where}: {label} is {text!r}; write true or false (empty means true)"
    )


def _read_conditions(
    path: str,
    *,
    data_dir: str,
    gas: str,
    defaults: Mapping[str, object],
) -> tuple[Conditions, ...]:
    """Parse a conditions CSV. `#` comment lines and blank lines are skipped."""
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            lines = [
                line
                for line in handle
                if line.strip() and not line.lstrip().startswith("#")
            ]
    except OSError as exc:
        raise ExperimentError(f"cannot read the conditions table {path}: {exc.strerror}") from exc
    if not lines:
        raise ExperimentError(f"{path} has no header row")

    reader = csv.DictReader(lines)
    if not reader.fieldnames:
        raise ExperimentError(f"{path} has no header row")

    canonical = {c.lower(): c for c in CONDITIONS_COLUMNS}
    mapping: dict[str, str] = {}
    extern_columns: dict[str, str] = {}
    unknown: list[str] = []
    for name in reader.fieldnames:
        stripped = (name or "").strip()
        if not stripped:
            continue
        if stripped.lower().startswith(EXTERN_PREFIX):
            extern_columns[stripped] = stripped[len(EXTERN_PREFIX) :].strip()
            continue
        target = canonical.get(stripped.lower())
        if target is None:
            unknown.append(stripped)
        else:
            mapping[target] = name
    if unknown:
        raise ExperimentError(
            f"{path}: unrecognised column(s) {', '.join(repr(u) for u in unknown)}. "
            f"Known columns are {', '.join(CONDITIONS_COLUMNS)}, plus "
            f"'{EXTERN_PREFIX}<_extern.inf key>' for a voltage override."
        )
    if "acquisition" not in mapping:
        raise ExperimentError(f"{path}: no 'acquisition' column")

    def cell(record: Mapping[str, str], column: str) -> str:
        source = mapping.get(column)
        if source is None:
            return ""
        return (record.get(source) or "").strip()

    def defaulted(record: Mapping[str, str], column: str, where: str) -> float:
        text = cell(record, column)
        if text:
            return _as_float(text, column, where)
        fallback = defaults.get(column)
        if fallback is None:
            raise ExperimentError(
                f"{where}: no {column}, and none in [defaults] of the experiment file"
            )
        return float(fallback)  # type: ignore[arg-type]

    rows: list[Conditions] = []
    for offset, record in enumerate(reader, start=2):
        where = f"{path} row {offset}"
        acquisition = cell(record, "acquisition")
        if not acquisition:
            raise ExperimentError(f"{where}: empty 'acquisition'")
        overrides = {
            key: (record.get(column) or "").strip()
            for column, key in extern_columns.items()
            if (record.get(column) or "").strip()
        }
        rows.append(
            Conditions(
                acquisition=acquisition,
                path=os.path.abspath(os.path.join(data_dir, acquisition)),
                pressure_torr=defaulted(record, "pressure_Torr", where),
                temperature_k=defaulted(record, "temperature_K", where),
                gas=cell(record, "gas") or gas,
                use=_as_bool(cell(record, "use"), "use", where),
                pressure_torr_err=_as_optional_float(
                    cell(record, "pressure_Torr_err"), "pressure_Torr_err", where
                ),
                temperature_k_err=_as_optional_float(
                    cell(record, "temperature_K_err"), "temperature_K_err", where
                ),
                v_drift_v=_as_optional_float(cell(record, "v_drift_V"), "v_drift_V", where),
                extern_overrides=overrides,
                notes=cell(record, "notes"),
                row=offset,
            )
        )
    if not rows:
        raise ExperimentError(f"{path} has a header but no acquisitions")
    return tuple(rows)


def _read_recalibration(block, source: str, data_dir: str) -> RecalibrationSpec | None:
    """`[recalibration]` from an `experiment.toml` table, or None if it has none.

    Strict about its keys for the same reason the conditions table is: a mistyped
    `min_apax` would otherwise be an experiment whose reference species are fitted from
    the baseline, with no complaint about it.
    """
    if block is None:
        return None
    if not isinstance(block, dict):
        raise ExperimentError(f"{source}: [recalibration] must be a table")
    unknown = [key for key in block if key not in RECALIBRATION_KEYS]
    if unknown:
        raise ExperimentError(
            f"{source}: [recalibration] has no key(s) {', '.join(repr(u) for u in unknown)}; "
            f"one of {', '.join(RECALIBRATION_KEYS)}"
        )

    order = block.get("order", 1)
    if not isinstance(order, int) or isinstance(order, bool) or order < 1:
        raise ExperimentError(
            f"{source}: [recalibration] order must be a whole number 1 or more; got {order!r}"
        )

    species: list[tuple[float, str]] = []
    raw_species = block.get("species", [])
    if not isinstance(raw_species, list):
        raise ExperimentError(
            f"{source}: [recalibration] species must be a list of tables, each with an "
            "'mz' and optionally a 'label'"
        )
    for index, entry in enumerate(raw_species, start=1):
        if not isinstance(entry, dict) or "mz" not in entry:
            raise ExperimentError(
                f"{source}: [recalibration] species {index} is {entry!r}; each needs an "
                "'mz' and may have a 'label'"
            )
        extra = set(entry) - {"mz", "label"}
        if extra:
            raise ExperimentError(
                f"{source}: [recalibration] species {index} has no key(s) "
                f"{', '.join(repr(e) for e in sorted(extra))}; only 'mz' and 'label'"
            )
        try:
            mz = float(entry["mz"])
        except (TypeError, ValueError) as exc:
            raise ExperimentError(
                f"{source}: [recalibration] species {index} has mz {entry['mz']!r}, "
                "which is not a number"
            ) from exc
        if not mz > 0.0:
            raise ExperimentError(
                f"{source}: [recalibration] species {index} has mz {mz}; an m/z is positive"
            )
        species.append((mz, str(entry.get("label") or "")))
    if species and len(species) < order + 2:
        raise ExperimentError(
            f"{source}: [recalibration] names {len(species)} species for an order-{order} "
            f"fit, which needs at least {order + 2} to leave a residual worth reading"
        )

    calibrant = block.get("calibrant")
    calibrant_path = None
    if calibrant is not None:
        calibrant = str(calibrant)
        calibrant_path = os.path.abspath(os.path.join(data_dir, calibrant))

    min_apex = block.get("min_apex", 0.0)
    try:
        min_apex = float(min_apex)
    except (TypeError, ValueError) as exc:
        raise ExperimentError(
            f"{source}: [recalibration] min_apex is {min_apex!r}, which is not a number"
        ) from exc
    if min_apex < 0.0:
        raise ExperimentError(
            f"{source}: [recalibration] min_apex is {min_apex}; an intensity floor is not negative"
        )

    mismatch = block.get("allow_polarity_mismatch", False)
    if not isinstance(mismatch, bool):
        raise ExperimentError(
            f"{source}: [recalibration] allow_polarity_mismatch is {mismatch!r}; it is true or false"
        )
    if mismatch and calibrant is None:
        raise ExperimentError(
            f"{source}: [recalibration] sets allow_polarity_mismatch with no calibrant to "
            "mismatch; the polarity check only applies to a borrowed calibrant frame"
        )

    return RecalibrationSpec(
        order=order,
        species=tuple(species),
        calibrant=calibrant,
        calibrant_path=calibrant_path,
        min_apex=min_apex,
        allow_polarity_mismatch=mismatch,
    )


def load_experiment(path: str | os.PathLike[str]) -> Experiment:
    """Load an experiment from its `experiment.toml`, or from the directory holding one.

    Raises `ExperimentError` for anything structurally wrong, naming the file and the
    field. Nothing is checked about the acquisitions themselves here -- call
    `Experiment.validate()` for that, and read the list it returns.
    """
    target = os.fspath(path)
    if os.path.isdir(target):
        target = os.path.join(target, EXPERIMENT_FILENAME)
    if not os.path.isfile(target):
        raise ExperimentError(f"no experiment definition at {target}")
    root = os.path.dirname(os.path.abspath(target))

    try:
        with open(target, "rb") as handle:
            table = tomllib.load(handle)
    except OSError as exc:
        raise ExperimentError(f"cannot read {target}: {exc.strerror}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ExperimentError(f"{target} is not valid TOML: {exc}") from exc

    schema = table.get("schema", SCHEMA)
    if schema != SCHEMA:
        raise ExperimentError(
            f"{target} declares experiment schema {schema!r}; this schamp understands {SCHEMA}"
        )

    profile_name = table.get("profile")
    if not profile_name:
        raise ExperimentError(
            f"{target}: no 'profile'. Name a builtin instrument profile, or give a path "
            "to your own; the drift voltage cannot be computed without one."
        )
    candidate = os.path.join(root, str(profile_name))
    try:
        profile = load_profile(candidate if os.path.isfile(candidate) else str(profile_name))
    except ProfileError as exc:
        raise ExperimentError(f"{target}: {exc}") from exc

    defaults = table.get("defaults", {})
    if not isinstance(defaults, dict):
        raise ExperimentError(f"{target}: [defaults] must be a table")
    canonical_defaults = {}
    for key, value in defaults.items():
        match = next((c for c in CONDITIONS_COLUMNS if c.lower() == str(key).lower()), None)
        if match is None:
            raise ExperimentError(
                f"{target}: [defaults] has no column {key!r}; one of "
                f"{', '.join(CONDITIONS_COLUMNS)}"
            )
        canonical_defaults[match] = value

    gas_name = str(table.get("gas") or profile.gas_default)
    try:
        constants.gas(gas_name)
    except KeyError as exc:
        raise ExperimentError(f"{target}: {exc.args[0]}") from exc

    conditions_path = os.path.join(root, str(table.get("conditions") or DEFAULT_CONDITIONS))
    data_dir = os.path.abspath(os.path.join(root, str(table.get("data_dir") or ".")))

    return Experiment(
        recalibration=_read_recalibration(table.get("recalibration"), target, data_dir),
        root=root,
        source=os.path.abspath(target),
        title=str(table.get("title") or os.path.basename(root)),
        profile=profile,
        gas=gas_name,
        conditions=_read_conditions(
            conditions_path, data_dir=data_dir, gas=gas_name, defaults=canonical_defaults
        ),
        conditions_source=os.path.abspath(conditions_path),
        notes=str(table.get("notes") or ""),
    )

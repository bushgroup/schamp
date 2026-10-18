"""The run sheet, without Qt: the planned series, the record it writes, the checks it runs.

An operator acquiring a drift-voltage series works down a list of settings and, for
each acquisition, notes the cell pressure off the manometer's display. This module is
that list and that record. It is everything `schamp-acquire` knows, and nothing in it
imports Qt, so all of it is checked by `tools/check_public.py` on a machine with no
display, no SDK and no acquisitions.

**The plan** is a series template (`load_template`): which `_extern.inf` keys the
series steps, the values of the others it holds fixed, and one row of settings per
acquisition in acquisition order. `uw-standard-14` ships with the package and is the
lab's standard order (lab record, task 14); a lab edits rows on the sheet, or loads a
template of its own by path. The nominal drift voltage of a row is the profile's own
`drift_voltage` over a mapping built from the planned settings -- the same function an
acquisition's recorded `_extern.inf` goes through, never a second formula.

**The record** is the experiment's two files and nothing more: `experiment.toml` holds
the session (title, profile, gas, data directory, notes, and the session temperature in
`[defaults]`) and `conditions.csv` one row per acquisition taken. An acquisition enters
`conditions.csv` when it is ticked as taken and leaves it if the tick is undone, so
**the tick is membership of that table** and is stored nowhere else. What the two files
cannot hold -- the planned settings, the rows not yet acquired and anything typed into
them, the operator, the file-name base, and the live pane's `[watch]` rule and
`[[species]]` (`live`) -- lives in a private `runsheet.toml` beside them, which no
analysis reads. Every save rewrites all three, each atomically.

**Sharing the table.** `schamp.importer` appends rows to the same `conditions.csv` as
acquisitions are copied in. A save first reads the table back and adopts any row it did
not write itself -- as a taken row, unplanned if its name matches nothing on the sheet --
so an import while the sheet is open is merged rather than overwritten. Cells the sheet
does not edit (`temperature_K`, `v_drift_V`, `extern:<key>`, uncertainties) and the
table's comment block are carried through untouched.

**The data directory** is where the `.raw` acquisitions appear: the experiment's
`data_dir`, which may be the set directory itself, a mapped drive or a UNC path, and may
be unreachable or empty for the whole session. `scan_data_dir` and `RunSheet.validate`
touch it and can block for as long as a dead share takes to time out, so a GUI calls
them off its own thread. Nothing else here touches it.
"""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from .. import DATA_DIR
from ..experiment import (
    CONDITIONS_COLUMNS,
    DEFAULT_CONDITIONS,
    EXPERIMENT_FILENAME,
    ExperimentError,
    RecalibrationSpec,
    _read_recalibration,
    load_experiment,
    read_conditions_table,
    toml_string,
    write_conditions_table,
    write_experiment,
)
from ..extern import EXTERN_FILENAME, as_float, parse_extern_inf
from ..profiles import InstrumentProfile, ProfileError, load_profile
from .live import WatchRule

__all__ = [
    "DEFAULT_TEMPLATE",
    "DataDirScan",
    "PlannedRow",
    "RUNSHEET_FILENAME",
    "RunSheet",
    "RunSheetError",
    "SeriesTemplate",
    "TOLERANCE_V",
    "VoltageCheck",
    "acquisition_name",
    "builtin_templates",
    "load_template",
    "scan_data_dir",
]

TEMPLATE_DIR = os.path.join(DATA_DIR, "series")
DEFAULT_TEMPLATE = "uw-standard-14"
RUNSHEET_FILENAME = "runsheet.toml"
SCHEMA = 1

TOLERANCE_V = 0.05
"""How far a recorded voltage may sit from the planned one before the sheet says so.

MassLynx records the setpoint as typed, to 0.1 V, so a knob that was set as planned
agrees exactly; anything over this is a knob that was not moved, or not moved far enough.
"""

EDITED_COLUMNS = ("pressure_Torr", "use", "notes")
"""The conditions columns the sheet edits per row. Every other column is carried through."""

DEFAULT_HEADER = ("acquisition", "pressure_Torr", "use", "notes")

CONDITIONS_COMMENT = """\
Written by schamp-acquire as the series was acquired: one row per acquisition taken,
in acquisition order. Pressures are typed from the drift-cell manometer's own display
at the time of each acquisition. The session temperature is in experiment.toml's
[defaults]. No voltage column: each drift voltage comes from the acquisition's own
_extern.inf through the instrument profile."""


class RunSheetError(Exception):
    """A template or run sheet that is missing, unreadable, or malformed."""


# --- the template ------------------------------------------------------------------


@dataclass(frozen=True)
class SeriesTemplate:
    """A planned drift-voltage series: what is stepped, what is held, and in what order.

    `vary` names the `_extern.inf` keys the series steps, and each of `rows` gives one
    value per key, in volts. `fixed` holds every other key the profile's drift voltage
    needs. `profile` is as the file wrote it; `profile_path` is what it resolved to, a
    path or a builtin name, ready for `load_profile`.
    """

    name: str
    title: str
    source: str
    profile: str
    profile_path: str
    vary: tuple[str, ...]
    fixed: Mapping[str, float]
    rows: tuple[tuple[float, ...], ...]


def builtin_templates() -> dict[str, str]:
    """`{name: path}` for every template shipped with the package."""
    if not os.path.isdir(TEMPLATE_DIR):
        return {}
    return {
        entry[: -len(".toml")]: os.path.join(TEMPLATE_DIR, entry)
        for entry in sorted(os.listdir(TEMPLATE_DIR))
        if entry.endswith(".toml")
    }


def _resolve_profile(name: str, relative_to: str) -> str:
    """A profile named in a file: a path relative to that file if one exists, else a
    builtin name -- the rule `load_experiment` applies."""
    candidate = os.path.join(relative_to, name)
    return os.path.abspath(candidate) if os.path.isfile(candidate) else name


def _check_keys(
    vary: Iterable[str], fixed: Iterable[str], profile: InstrumentProfile, where: str
) -> None:
    vary, fixed = list(vary), list(fixed)
    needed = set(profile.high_keys) | set(profile.low_keys)
    named = set(vary) | set(fixed)
    both = set(vary) & set(fixed)
    missing = sorted(needed - named)
    extra = sorted(named - needed)
    problems = []
    if both:
        problems.append(f"{', '.join(map(repr, sorted(both)))} both varied and fixed")
    if missing:
        problems.append(f"no value for {', '.join(map(repr, missing))}")
    if extra:
        problems.append(
            f"{', '.join(map(repr, extra))} not in the drift voltage of profile {profile.name!r}"
        )
    if problems:
        raise RunSheetError(
            f"{where}: {'; '.join(problems)}. `vary` and [fixed] together must name exactly "
            f"the profile's high_keys and low_keys: "
            f"{', '.join(map(repr, profile.high_keys + profile.low_keys))}."
        )


def load_template(name_or_path: str | os.PathLike[str]) -> SeriesTemplate:
    """Load a series template by builtin name or by path. Raises `RunSheetError`.

    Checked against the profile it names: `vary` and `[fixed]` together must be exactly
    the keys its drift voltage sums, so that every row has a nominal drift voltage.
    """
    text = os.fspath(name_or_path)
    path = builtin_templates().get(text, text)
    if not os.path.isfile(path):
        raise RunSheetError(
            f"no series template {text!r}: not a file, and not one of the builtin "
            f"templates ({', '.join(builtin_templates()) or 'none'})"
        )
    try:
        with open(path, "rb") as handle:
            table = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise RunSheetError(f"cannot read the series template {path}: {exc}") from exc

    if table.get("schema", SCHEMA) != SCHEMA:
        raise RunSheetError(f"{path}: template schema {table.get('schema')!r}; this schamp reads {SCHEMA}")
    profile_name = table.get("profile")
    if not isinstance(profile_name, str) or not profile_name:
        raise RunSheetError(f"{path}: no 'profile'; name the instrument profile the series is for")
    profile_path = _resolve_profile(profile_name, os.path.dirname(os.path.abspath(path)))
    try:
        profile = load_profile(profile_path)
    except ProfileError as exc:
        raise RunSheetError(f"{path}: {exc}") from exc

    vary = table.get("vary")
    if (
        not isinstance(vary, list)
        or not vary
        or not all(isinstance(k, str) and k for k in vary)
        or len(set(vary)) != len(vary)
    ):
        raise RunSheetError(f"{path}: 'vary' must be a list of distinct _extern.inf key names")
    fixed = table.get("fixed", {})
    if not isinstance(fixed, dict) or not all(
        isinstance(v, (int, float)) and not isinstance(v, bool) for v in fixed.values()
    ):
        raise RunSheetError(f"{path}: [fixed] must map _extern.inf keys to numbers of volts")
    _check_keys(vary, fixed, profile, path)

    rows = table.get("rows")
    if not isinstance(rows, list) or not rows:
        raise RunSheetError(f"{path}: 'rows' must be a non-empty list, one row per acquisition")
    parsed: list[tuple[float, ...]] = []
    for number, row in enumerate(rows, start=1):
        if (
            not isinstance(row, list)
            or len(row) != len(vary)
            or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in row)
        ):
            raise RunSheetError(
                f"{path}: row {number} is {row!r}; each row is {len(vary)} number(s) of volts, "
                f"one per key of 'vary' ({', '.join(vary)})"
            )
        parsed.append(tuple(float(v) for v in row))

    return SeriesTemplate(
        name=str(table.get("name") or os.path.splitext(os.path.basename(path))[0]),
        title=str(table.get("title") or ""),
        source=os.path.abspath(path),
        profile=profile_name,
        profile_path=profile_path,
        vary=tuple(vary),
        fixed={str(k): float(v) for k, v in fixed.items()},
        rows=tuple(parsed),
    )


def acquisition_name(base: str, number: int) -> str:
    """MassLynx's sample-list numbering: `<base>_NNN.raw`."""
    return f"{base}_{number:03d}.raw"


# --- the sheet ---------------------------------------------------------------------


@dataclass
class PlannedRow:
    """One line of the sheet: an acquisition, planned or not, taken or not.

    `settings` maps each varied key to its planned value in volts, and is empty for an
    acquisition that arrived without a plan -- an import, a row typed into the table by
    hand. `cells` holds the conditions-table cells by column name (`pressure_Torr`,
    `use`, `notes`, and any other column the table carries), as text, so a pressure is
    written exactly as it was typed.
    """

    acquisition: str
    settings: dict[str, float] = field(default_factory=dict)
    taken: bool = False
    cells: dict[str, str] = field(default_factory=dict)

    @property
    def planned(self) -> bool:
        return bool(self.settings)

    @property
    def use(self) -> bool:
        return self.cells.get("use", "").strip().lower() not in {"0", "false", "f", "no", "n"}


def _canonical(column: str) -> str:
    return next((c for c in CONDITIONS_COLUMNS if c.lower() == column.lower()), column)


@dataclass
class RunSheet:
    """A drift-voltage series being acquired: the session, the plan, and the record.

    Plain and mutable; a GUI edits the fields and calls `save()` after every change.
    `profile` names the instrument profile as `experiment.toml` will say it, a builtin
    name or a path relative to `set_dir`. `temperature_k` is the session temperature
    as typed, written to `[defaults]`, and `fixed` the held voltages for this session.
    """

    set_dir: str
    title: str = ""
    profile: str = "uw-synapt-g2"
    gas: str = ""
    data_dir: str = "."
    notes: str = ""
    temperature_k: str = ""
    operator: str = ""
    base_name: str = ""
    template: str = ""
    vary: tuple[str, ...] = ()
    fixed: dict[str, float] = field(default_factory=dict)
    rows: list[PlannedRow] = field(default_factory=list)
    header: list[str] = field(default_factory=lambda: list(DEFAULT_HEADER))
    comment: str = CONDITIONS_COMMENT
    other_defaults: dict[str, object] = field(default_factory=dict)
    recalibration: RecalibrationSpec | None = None
    watch: dict[str, object] = field(default_factory=dict)
    species: list[dict[str, object]] = field(default_factory=list)
    _written: set[str] = field(default_factory=set, repr=False)

    # -- construction ----------------------------------------------------------------

    @classmethod
    def from_template(
        cls,
        template: SeriesTemplate | str,
        set_dir: str | os.PathLike[str],
        *,
        base_name: str = "",
        title: str = "",
    ) -> RunSheet:
        """A new sheet for `set_dir` with the template's rows planned and named."""
        set_dir = os.path.abspath(os.fspath(set_dir))
        sheet = cls(
            set_dir=set_dir,
            title=title or os.path.basename(set_dir),
            base_name=base_name or os.path.basename(set_dir),
        )
        sheet.apply_template(template)
        return sheet

    def apply_template(self, template: SeriesTemplate | str) -> None:
        """Lay a template over the sheet: its keys, its fixed values, and its rows.

        Rows are matched by position. Where the sheet already has an unplanned row there
        -- an acquisition imported before any plan existed -- that row takes the
        template's settings and keeps its name; past the end of the sheet, new rows are
        added and named from `base_name`. A planned row already on the sheet keeps its
        own settings: re-applying a template does not overwrite what was planned.
        """
        if not isinstance(template, SeriesTemplate):
            template = load_template(template)
        profile = template.profile_path
        if os.path.isfile(profile):
            profile = os.path.relpath(profile, self.set_dir) if _same_drive(profile, self.set_dir) else profile
        self.profile = profile
        self.template = template.source if template.source not in builtin_templates().values() else template.name
        self.vary = template.vary
        self.fixed = dict(template.fixed)
        for index, values in enumerate(template.rows):
            settings = dict(zip(template.vary, values))
            if index < len(self.rows):
                if not self.rows[index].planned:
                    self.rows[index].settings = settings
            else:
                self.rows.append(
                    PlannedRow(
                        acquisition=acquisition_name(self.base_name or "acquisition", index + 1),
                        settings=settings,
                        cells={"use": "true"},
                    )
                )

    @classmethod
    def open(cls, set_dir: str | os.PathLike[str]) -> RunSheet:
        """Reopen a set: `experiment.toml`, `conditions.csv` and `runsheet.toml`.

        Any of the three may be missing -- an imported set has no `runsheet.toml`, and a
        sheet with nothing taken yet has no `conditions.csv` -- but not all three.
        `conditions.csv` decides what was taken; `runsheet.toml` supplies the plan.
        """
        set_dir = os.path.abspath(os.fspath(set_dir))
        experiment_path = os.path.join(set_dir, EXPERIMENT_FILENAME)
        runsheet_path = os.path.join(set_dir, RUNSHEET_FILENAME)
        if not os.path.isfile(experiment_path) and not os.path.isfile(runsheet_path):
            raise RunSheetError(f"{set_dir} holds neither {EXPERIMENT_FILENAME} nor {RUNSHEET_FILENAME}")
        sheet = cls(set_dir=set_dir, title=os.path.basename(set_dir))

        conditions = DEFAULT_CONDITIONS
        if os.path.isfile(experiment_path):
            table = _read_toml(experiment_path)
            sheet.title = str(table.get("title") or sheet.title)
            sheet.profile = str(table.get("profile") or sheet.profile)
            sheet.gas = str(table.get("gas") or "")
            sheet.data_dir = str(table.get("data_dir") or ".")
            sheet.notes = str(table.get("notes") or "")
            conditions = str(table.get("conditions") or DEFAULT_CONDITIONS)
            if conditions != DEFAULT_CONDITIONS:
                raise RunSheetError(
                    f"{experiment_path} names its table {conditions!r}; the run sheet keeps "
                    f"its record in {DEFAULT_CONDITIONS} and will not write elsewhere"
                )
            defaults = dict(table.get("defaults") or {})
            temperature = next((k for k in defaults if k.lower() == "temperature_k"), None)
            if temperature is not None:
                sheet.temperature_k = _number_text(defaults.pop(temperature))
            sheet.other_defaults = defaults
            sheet.recalibration = _read_recalibration(
                table.get("recalibration"), experiment_path, sheet.data_path
            )

        if os.path.isfile(runsheet_path):
            table = _read_toml(runsheet_path)
            if table.get("schema", SCHEMA) != SCHEMA:
                raise RunSheetError(f"{runsheet_path}: run-sheet schema {table.get('schema')!r}")
            sheet.template = str(table.get("template") or "")
            sheet.base_name = str(table.get("base_name") or "")
            sheet.operator = str(table.get("operator") or "")
            sheet.vary = tuple(str(k) for k in table.get("vary") or ())
            sheet.fixed = {str(k): float(v) for k, v in (table.get("fixed") or {}).items()}
            sheet.watch = dict(table.get("watch") or {})
            sheet.species = [dict(entry) for entry in table.get("species") or []]
            for entry in table.get("row") or []:
                values = entry.get("settings") or []
                if values and len(values) != len(sheet.vary):
                    raise RunSheetError(
                        f"{runsheet_path}: {entry.get('acquisition')!r} has {len(values)} "
                        f"setting(s) for {len(sheet.vary)} varied key(s)"
                    )
                sheet.rows.append(
                    PlannedRow(
                        acquisition=str(entry["acquisition"]),
                        settings=dict(zip(sheet.vary, (float(v) for v in values))),
                        cells={str(k): _number_text(v) for k, v in (entry.get("cells") or {}).items()},
                    )
                )
        if not sheet.base_name:
            sheet.base_name = os.path.basename(set_dir)

        sheet._adopt_table()
        return sheet

    # -- the record ------------------------------------------------------------------

    @property
    def conditions_path(self) -> str:
        return os.path.join(self.set_dir, DEFAULT_CONDITIONS)

    def _adopt_table(self) -> None:
        """Bring the sheet into line with `conditions.csv` as it is on disk.

        A row in the table is taken, and its cells are the table's. A row the sheet has
        never seen is adopted -- into the sheet row of the same name if there is one,
        else appended unplanned. A row the sheet itself wrote and has since untaken is
        not re-adopted: that is the sheet's own edit, not news.
        """
        if not os.path.isfile(self.conditions_path):
            return
        table = read_conditions_table(self.conditions_path)
        self.header = list(table.header)
        self.comment = table.comment
        by_name = {row.acquisition: row for row in self.rows}
        on_disk = set()
        for record in table.rows:
            name = table.get(record, "acquisition")
            if not name:
                continue
            on_disk.add(name)
            cells = {_canonical(k): v for k, v in record.items() if _canonical(k) != "acquisition"}
            row = by_name.get(name)
            if row is None:
                if name in self._written:
                    continue
                row = PlannedRow(acquisition=name)
                self.rows.append(row)
                by_name[name] = row
            if row.taken or name in self._written:
                continue
            row.taken = True
            for key, value in cells.items():
                if value.strip() or not row.cells.get(key, "").strip():
                    row.cells[key] = value
        self._written |= on_disk

    def save(self) -> None:
        """Write `conditions.csv`, `experiment.toml` and `runsheet.toml`, each atomically.

        Adopts any row another writer added to `conditions.csv` first. Raises `OSError`
        if a file cannot be replaced -- on Windows, typically because a spreadsheet has
        it open -- leaving that file as it was; the sheet stays as it is in memory, and
        the next save tries again.
        """
        os.makedirs(self.set_dir, exist_ok=True)
        self._adopt_table()
        taken = [row for row in self.rows if row.taken]
        if taken or os.path.isfile(self.conditions_path):
            header = list(self.header)
            lowered = {h.lower() for h in header}
            for row in taken:
                for key in list(EDITED_COLUMNS) + list(row.cells):
                    if key.lower() not in lowered:
                        header.append(key)
                        lowered.add(key.lower())
            records = []
            for row in taken:
                record = {"acquisition": row.acquisition}
                record.update({_canonical(k): v for k, v in row.cells.items()})
                records.append({h: record.get(_canonical(h), "") for h in header})
            write_conditions_table(self.conditions_path, header, records, comment=self.comment)
            self.header = header
        self._written = {row.acquisition for row in taken}

        defaults: dict[str, object] = dict(self.other_defaults)
        if self.temperature_k.strip():
            try:
                defaults["temperature_K"] = float(self.temperature_k)
            except ValueError:
                defaults["temperature_K"] = self.temperature_k.strip()
        write_experiment(
            self.set_dir,
            title=self.title,
            profile=self.profile,
            gas=self.gas or None,
            data_dir=self.data_dir or ".",
            notes=self.notes,
            defaults=defaults,
            recalibration=self.recalibration,
            comment=(
                "Written by schamp-acquire; edits made here while the run sheet is open are\n"
                f"overwritten. Its plan of the series is in {RUNSHEET_FILENAME}."
            ),
        )
        self._write_runsheet()

    def _write_runsheet(self) -> None:
        lines = [
            "# The schamp-acquire run sheet: the planned series and the rows not yet taken.",
            f"# The record of what was taken is {EXPERIMENT_FILENAME} + {DEFAULT_CONDITIONS};",
            "# no analysis reads this file.",
            f"schema = {SCHEMA}",
            f"template = {toml_string(self.template)}",
            f"base_name = {toml_string(self.base_name)}",
            f"operator = {toml_string(self.operator)}",
            f"vary = [{', '.join(toml_string(k) for k in self.vary)}]",
            "",
            "[fixed]",
        ]
        lines += [f"{toml_string(k)} = {float(v)!r}" for k, v in self.fixed.items()]
        lines += [
            "",
            "# The live pane calls a .raw finished by the [watch] rule (schamp.acquire.live.WatchRule).",
            "# Without a [watch] table the default applies; write one to change any of its terms:",
        ]
        lines += [f"#   {k} = {_toml_value(v)}" for k, v in WatchRule().to_table().items()]
        if self.watch:
            lines.append("[watch]")
            lines += [f"{toml_string(k)} = {_toml_value(v)}" for k, v in self.watch.items()]
        for entry in self.species:
            lines += ["", "[[species]]"]
            lines += [f"{toml_string(k)} = {_toml_value(v)}" for k, v in entry.items()]
        for row in self.rows:
            lines += ["", "[[row]]", f"acquisition = {toml_string(row.acquisition)}"]
            if row.planned:
                values = ", ".join(repr(float(row.settings[k])) for k in self.vary)
                lines.append(f"settings = [{values}]")
            if not row.taken:
                cells = ", ".join(
                    f"{toml_string(k)} = {toml_string(v)}" for k, v in row.cells.items() if v.strip()
                )
                if cells:
                    lines.append(f"cells = {{ {cells} }}")
        text = "\n".join(lines) + "\n"
        path = os.path.join(self.set_dir, RUNSHEET_FILENAME)
        temp = f"{path}.{os.getpid()}.tmp"
        try:
            with open(temp, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
            os.replace(temp, path)
        finally:
            if os.path.exists(temp):
                os.remove(temp)

    # -- editing ---------------------------------------------------------------------

    def propose_names(self) -> None:
        """Name every row not yet taken `<base_name>_NNN.raw`, NNN its position on the
        sheet -- MassLynx's sample-list numbering. Taken rows keep their names."""
        taken = {row.acquisition for row in self.rows if row.taken}
        for index, row in enumerate(self.rows, start=1):
            if not row.taken:
                name = acquisition_name(self.base_name or "acquisition", index)
                row.acquisition = name if name not in taken else row.acquisition

    def insert_row(self, index: int, settings: Mapping[str, float] | None = None) -> PlannedRow:
        """A new planned row at `index`, its settings copied from `settings` or zero."""
        row = PlannedRow(
            acquisition=acquisition_name(self.base_name or "acquisition", index + 1),
            settings={k: float((settings or {}).get(k, 0.0)) for k in self.vary},
            cells={"use": "true"},
        )
        self.rows.insert(index, row)
        return row

    # -- the physics -----------------------------------------------------------------

    def load_profile(self) -> InstrumentProfile:
        """The instrument profile, resolved the way `load_experiment` resolves it."""
        return load_profile(_resolve_profile(self.profile, self.set_dir))

    def planned_extern(self, row: PlannedRow) -> dict[str, str]:
        """The `_extern.inf` mapping the row's plan amounts to: fixed and varied keys."""
        mapping = {k: repr(float(v)) for k, v in self.fixed.items()}
        mapping.update({k: repr(float(v)) for k, v in row.settings.items()})
        return mapping

    def nominal_voltage(self, row: PlannedRow, profile: InstrumentProfile | None = None) -> float | None:
        """The planned drift voltage, through the profile's own formula; None if unplanned."""
        if not row.planned:
            return None
        return (profile or self.load_profile()).drift_voltage(self.planned_extern(row))

    def check_recorded(
        self,
        row: PlannedRow,
        extern: Mapping[str, str],
        profile: InstrumentProfile | None = None,
    ) -> VoltageCheck:
        """Compare what an acquisition recorded with what the row planned."""
        profile = profile or self.load_profile()
        try:
            recorded = profile.drift_voltage(extern)
            error = ""
        except ProfileError as exc:
            recorded, error = None, str(exc)
        if not row.planned:
            return VoltageCheck(None, recorded, (), error)
        planned = self.planned_extern(row)
        knobs = []
        for key, text in planned.items():
            want = float(text)
            got = as_float(extern, key)
            if got is None or abs(got - want) > TOLERANCE_V:
                knobs.append((key, want, got))
        return VoltageCheck(profile.drift_voltage(planned), recorded, tuple(knobs), error)

    # -- checks (these touch the data directory: call them off the GUI thread) ------

    @property
    def data_path(self) -> str:
        """`data_dir` resolved against the set directory."""
        return os.path.abspath(os.path.join(self.set_dir, self.data_dir or "."))

    def validate(self, *, require_raw: bool = False) -> list[str]:
        """`Experiment.validate` over the record as saved, plus the sheet's own checks.

        Reads the files `save()` wrote, so a GUI calls it after saving. With
        `require_raw`, it stats every acquisition under `data_dir`, which can block on a
        share that has gone away.
        """
        problems: list[str] = []
        names = [row.acquisition for row in self.rows]
        for name in sorted({n for n in names if names.count(n) > 1}):
            problems.append(f"{name} appears {names.count(name)} times on the run sheet")
        if not any(row.taken for row in self.rows):
            problems.append(
                f"no acquisition is ticked as taken yet, so {DEFAULT_CONDITIONS} has nothing to check"
            )
            return problems
        try:
            experiment = load_experiment(self.set_dir)
        except ExperimentError as exc:
            return problems + [str(exc)]
        return problems + experiment.validate(require_raw=require_raw)


def _toml_value(value: object) -> str:
    """A scalar or a list of scalars as TOML, for the run sheet's own tables."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    return toml_string(str(value))


def _same_drive(a: str, b: str) -> bool:
    return os.path.splitdrive(os.path.abspath(a))[0].lower() == os.path.splitdrive(os.path.abspath(b))[0].lower()


def _read_toml(path: str) -> dict:
    try:
        with open(path, "rb") as handle:
            return tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise RunSheetError(f"cannot read {path}: {exc}") from exc


def _number_text(value: object) -> str:
    """A TOML value as the text a cell would hold: 301.13, not 301.13000000000001."""
    if isinstance(value, bool):
        return "true" if value else "false"
    return repr(value) if isinstance(value, float) else str(value)


# --- planned against recorded ------------------------------------------------------


@dataclass(frozen=True)
class VoltageCheck:
    """One acquisition's drift voltage, planned and recorded, and the knobs that differ.

    `knobs` lists `(key, planned, recorded)` for every key of the plan whose recorded
    value is missing or more than `TOLERANCE_V` from it -- the knob that did not move.
    `error` is the profile's complaint when the sidecar lacks a key it needs.
    """

    planned_v: float | None
    recorded_v: float | None
    knobs: tuple[tuple[str, float, float | None], ...] = ()
    error: str = ""

    @property
    def delta_v(self) -> float | None:
        if self.planned_v is None or self.recorded_v is None:
            return None
        return self.recorded_v - self.planned_v

    @property
    def agrees(self) -> bool:
        delta = self.delta_v
        return delta is not None and abs(delta) <= TOLERANCE_V and not self.knobs


@dataclass(frozen=True)
class DataDirScan:
    """What one look at the data directory found.

    `externs` maps each `.raw` directory name to its parsed `_extern.inf`, empty for an
    acquisition MassLynx has not finished writing (no sidecar yet). `reachable` is False,
    with `error` saying why, when the directory could not be listed at all.
    """

    path: str
    reachable: bool
    externs: Mapping[str, Mapping[str, str]] = field(default_factory=dict)
    error: str = ""


def scan_data_dir(path: str | os.PathLike[str]) -> DataDirScan:
    """List the `.raw` directories under `path` and read each one's `_extern.inf`.

    Blocking, and on a dead share it blocks until Windows gives up: call it from a
    worker thread. Never raises.
    """
    path = os.fspath(path)
    try:
        with os.scandir(path) as entries:
            raws = sorted(
                entry.name
                for entry in entries
                if entry.name.lower().endswith(".raw") and entry.is_dir()
            )
    except OSError as exc:
        return DataDirScan(path=path, reachable=False, error=exc.strerror or str(exc))
    externs = {}
    for name in raws:
        raw = os.path.join(path, name)
        externs[name] = parse_extern_inf(raw) if os.path.isfile(os.path.join(raw, EXTERN_FILENAME)) else {}
    return DataDirScan(path=path, reachable=True, externs=externs)

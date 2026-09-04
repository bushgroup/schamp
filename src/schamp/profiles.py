"""Instrument profiles: the cell geometry and the drift-voltage definition, as data.

A drift voltage is not a property of an experiment, it is a property of an instrument.
Which `_extern.inf` keys add up to the applied high potential, how many resistor steps
the ladder has and how many of them fall outside the drift region, and how long the
drift region is -- all of that differs between the copies of this cell now running on
G2, G2-S and G2-Si instruments, and none of it belongs in analysis code. So it lives
in a small TOML file per instrument, and `mobility` asks the profile.

`src/schamp/data/profiles/uw-synapt-g2.toml` is shipped as the reference example and
is annotated to be read as the format specification. A lab adds its own instrument by
copying that file, editing it, and loading it by path.

The one rule the format enforces on itself: **the drift length and the divider ratio
are independent fields, never derived from each other.** 167 resistor steps at a
1.5 mm electrode pitch make the quoted 25.05 cm drift region and 168 make the quoted
25.2 cm mechanical length; Waters gives all of those numbers, they do not close, and
nothing in the physics needs them to. Reconstructing one from another would silently
choose which of the vendor's numbers to disbelieve.

No SDK and no data files needed: profiles load and evaluate on any machine, which is
what lets the public self-check exercise the drift-voltage formula.
"""

from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field
from typing import Mapping

from . import DATA_DIR
from .extern import as_float

__all__ = [
    "PROFILE_DIR",
    "InstrumentProfile",
    "ProfileError",
    "builtin_profiles",
    "load_profile",
]

PROFILE_DIR = os.path.join(DATA_DIR, "profiles")

SCHEMA = 1


class ProfileError(Exception):
    """A profile that is missing, unreadable, or does not describe an instrument."""


@dataclass(frozen=True)
class InstrumentProfile:
    """One instrument's cell geometry, drift-voltage definition and key names.

    Constructed by `load_profile`; the fields mirror the TOML one for one. Frozen,
    because a profile is a description of hardware and nothing downstream has any
    business editing it mid-analysis.
    """

    name: str
    title: str
    instrument: str
    source: str
    """Where this profile was loaded from -- a path, for the record and for results.json."""

    drift_length_cm: float
    resistor_steps: int
    untapped_steps: int
    high_keys: tuple[str, ...]
    low_keys: tuple[str, ...]
    gas_default: str

    drift_length_cm_err: float | None = None
    mechanical_length_cm: float | None = None
    electrode_count: int | None = None
    electrode_pitch_mm: float | None = None
    electrode_bore_mm: float | None = None
    excluded_keys: tuple[str, ...] = ()
    keys: Mapping[str, str] = field(default_factory=dict)
    reference: str = ""
    require_wave_off: bool = False
    expect_pushes_per_bin: int | None = None

    @property
    def divider_ratio(self) -> float:
        """`untapped_steps / resistor_steps`: the fraction of the applied drop lost.

        On the UW G2 this is 2/170. Neither plate bounding the drift region is tapped
        at the end of the ladder -- the first RF electrode sits one resistor below
        VappH and the exit plate one resistor above VappL -- so one step is lost at
        each end of a 170-step ladder.
        """
        return self.untapped_steps / self.resistor_steps

    def applied_voltage(self, extern: Mapping[str, str]) -> tuple[float, float]:
        """`(VappH, VappL)` in volts, summed over `high_keys` and `low_keys`.

        Signs are whatever `_extern.inf` recorded: `Helium Exit` reads -40.0 on the
        UW G2, so having it in `high_keys` correctly subtracts 40 V. Raises
        `ProfileError` naming every key the sidecar did not carry a number for --
        never falls back to zero, which is the failure mode that put every legacy
        centroid in the wrong units.
        """
        missing: list[str] = []
        totals: list[float] = []
        for keys in (self.high_keys, self.low_keys):
            total = 0.0
            for key in keys:
                value = as_float(extern, key)
                if value is None:
                    missing.append(key)
                else:
                    total += value
            totals.append(total)
        if missing:
            raise ProfileError(
                f"profile {self.name!r} needs _extern.inf keys that this acquisition "
                f"does not carry as numbers: {', '.join(repr(k) for k in missing)}. "
                "Key naming differs between instruments; check [drift_voltage] in "
                f"{self.source}."
            )
        return totals[0], totals[1]

    def drift_voltage(self, extern: Mapping[str, str]) -> float:
        """The drift voltage in volts, from one acquisition's `_extern.inf` mapping.

            V_drift = (VappH - VappL) * (1 - untapped_steps / resistor_steps)

        This is the expression the paper's ESI Fig. S2B prints and the one Waters' own
        guide to the cell gives. It is the decision of record, and the 2013 workbooks'
        variant -- which took the correction from VappH rather than from the drop, and
        applied it at one end rather than two -- is a defect rather than an
        alternative: it runs 0.58 % high, non-proportionally, so it does not cancel out
        of a slope.
        """
        high, low = self.applied_voltage(extern)
        return (high - low) * (1.0 - self.divider_ratio)

    def check_acquisition(self, extern: Mapping[str, str]) -> list[str]:
        """Sanity checks against one acquisition's `_extern.inf`; a list of complaints.

        Empty means nothing to complain about, including when the sidecar was
        unreadable -- an absent parameter cannot fail a check, only a present and wrong
        one can. Returns strings rather than raising so that a caller can report every
        problem in a series at once.
        """
        problems: list[str] = []
        if self.require_wave_off:
            key = self.keys.get("wave_height", "")
            height = as_float(extern, key) if key else None
            if height is not None and height != 0.0:
                problems.append(
                    f"{key} is {height} V, not 0: the traveling wave was on, so this is "
                    "not a linear-field drift measurement"
                )
        if self.expect_pushes_per_bin is not None:
            key = self.keys.get("pushes_per_bin", "")
            pushes = as_float(extern, key) if key else None
            if pushes is not None and int(pushes) != self.expect_pushes_per_bin:
                problems.append(
                    f"{key} is {int(pushes)}, not {self.expect_pushes_per_bin}: one drift "
                    "bin is no longer one pusher period"
                )
        return problems


def builtin_profiles() -> dict[str, str]:
    """`{name: path}` for the profiles that ship with schamp, by file stem."""
    try:
        names = sorted(os.listdir(PROFILE_DIR))
    except OSError:
        return {}
    return {
        os.path.splitext(n)[0]: os.path.join(PROFILE_DIR, n)
        for n in names
        if n.endswith(".toml")
    }


def _resolve(name_or_path: str | os.PathLike[str]) -> str:
    text = os.fspath(name_or_path)
    if os.path.isfile(text):
        return text
    builtin = builtin_profiles()
    if text in builtin:
        return builtin[text]
    raise ProfileError(
        f"no instrument profile {text!r}: not a file, and not one of the profiles "
        f"schamp ships ({', '.join(builtin) or 'none'}). Copy "
        f"{os.path.join(PROFILE_DIR, 'uw-synapt-g2.toml')} and load yours by path."
    )


def _require(table: Mapping[str, object], section: str, key: str, path: str) -> object:
    try:
        return table[key]
    except KeyError:
        raise ProfileError(f"{path}: [{section}] is missing {key!r}") from None


def load_profile(name_or_path: str | os.PathLike[str]) -> InstrumentProfile:
    """Load an instrument profile by builtin name (`"uw-synapt-g2"`) or by path.

    Raises `ProfileError` for anything wrong with the file, naming the file and the
    field. A profile is read once at the start of an analysis and is then the single
    source of the geometry, so failing loudly here is much cheaper than a plausible
    default surviving into a cross section.
    """
    path = _resolve(name_or_path)
    try:
        with open(path, "rb") as handle:
            table = tomllib.load(handle)
    except OSError as exc:
        raise ProfileError(f"cannot read the profile {path}: {exc.strerror}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ProfileError(f"{path} is not valid TOML: {exc}") from exc

    schema = table.get("schema", SCHEMA)
    if schema != SCHEMA:
        raise ProfileError(
            f"{path} declares profile schema {schema!r}; this schamp understands {SCHEMA}"
        )

    cell = table.get("cell", {})
    voltage = table.get("drift_voltage", {})
    gas_table = table.get("gas", {})
    checks = table.get("checks", {})
    if not isinstance(cell, dict) or not isinstance(voltage, dict):
        raise ProfileError(f"{path}: [cell] and [drift_voltage] must be tables")

    high = _require(voltage, "drift_voltage", "high_keys", path)
    low = voltage.get("low_keys", [])
    if not isinstance(high, list) or not high:
        raise ProfileError(f"{path}: [drift_voltage] high_keys must be a non-empty list")
    if not isinstance(low, list):
        raise ProfileError(f"{path}: [drift_voltage] low_keys must be a list")

    steps = int(_require(cell, "cell", "resistor_steps", path))  # type: ignore[arg-type]
    untapped = int(_require(cell, "cell", "untapped_steps", path))  # type: ignore[arg-type]
    if steps <= 0 or untapped < 0 or untapped >= steps:
        raise ProfileError(
            f"{path}: [cell] needs 0 <= untapped_steps < resistor_steps; got "
            f"{untapped} of {steps}"
        )
    length = float(_require(cell, "cell", "drift_length_cm", path))  # type: ignore[arg-type]
    if length <= 0.0:
        raise ProfileError(f"{path}: [cell] drift_length_cm must be positive; got {length}")

    stem = os.path.splitext(os.path.basename(path))[0]
    keys = table.get("keys", {})
    if not isinstance(keys, dict):
        raise ProfileError(f"{path}: [keys] must be a table of _extern.inf key names")

    return InstrumentProfile(
        name=str(table.get("name", stem)),
        title=str(table.get("title", stem)),
        instrument=str(table.get("instrument", "")),
        source=os.path.abspath(path),
        drift_length_cm=length,
        drift_length_cm_err=_optional_float(cell.get("drift_length_cm_err")),
        mechanical_length_cm=_optional_float(cell.get("mechanical_length_cm")),
        electrode_count=_optional_int(cell.get("electrode_count")),
        electrode_pitch_mm=_optional_float(cell.get("electrode_pitch_mm")),
        electrode_bore_mm=_optional_float(cell.get("electrode_bore_mm")),
        resistor_steps=steps,
        untapped_steps=untapped,
        high_keys=tuple(str(k) for k in high),
        low_keys=tuple(str(k) for k in low),
        excluded_keys=tuple(str(k) for k in voltage.get("excluded_keys", [])),
        gas_default=str(gas_table.get("default", "helium")),
        keys={str(k): str(v) for k, v in keys.items()},
        reference=str(table.get("reference", "")),
        require_wave_off=bool(checks.get("require_wave_off", False)),
        expect_pushes_per_bin=_optional_int(checks.get("expect_pushes_per_bin")),
    )


def _optional_float(value: object) -> float | None:
    return None if value is None else float(value)  # type: ignore[arg-type]


def _optional_int(value: object) -> int | None:
    return None if value is None else int(value)  # type: ignore[arg-type]

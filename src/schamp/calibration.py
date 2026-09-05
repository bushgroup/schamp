"""The two m/z frames of one acquisition, and the conversion between them.

An analysis knows a species' m/z from its formula. To read that species out it has to
say where the peak *is*, and on a Waters `.raw` that is two questions, not one, because
the SDK's two families of reader do not use the same m/z axis:

* the **spectrum frame** -- what `ReadScan`, `ReadDriftScan` and
  `MassLynxScanProcessor` return as the m/z axis;
* the **chromatogram frame** -- the axis `ReadMobillogram` accepts its `lowMass` and
  `highMass` bounds in. It returns no axis at all, so this one is only ever measured by
  what comes back for a given bound.

The difference between them is not an error in either reader. It is the acquisition's
own `$$ Cal Function N` from `_HEADER.TXT`: a polynomial in sqrt(m/z) whose square is
the corrected mass, which the spectrum readers apply and the chromatogram reader does
not. That is measured and not assumed (lab record, task 10):

* the order of sqrt(m/z) against time-of-flight channel equals the order of that
  polynomial on every acquisition tested, at the single-precision floor of the returned
  axis -- an identity function gives an axis linear in channel, a cubic gives a cubic;
* applying it to a peak's position in the chromatogram frame lands on the same peak's
  position in the spectrum frame, across the mass range, on acquisitions carrying both
  a linear and a cubic function, to the precision a chromatogram-frame centroid has;
* the decisive case, because it holds everything else fixed: one acquisition with two
  copies carrying **different** calibration functions has one chromatogram frame,
  identical species for species, and two spectrum frames differing by exactly the two
  polynomials.

So the conversion is a property of the file, not a fit, and `CalFunction` is it. On an
acquisition whose function is the identity the two frames coincide and every call here
is a no-op -- which is worth saying, because that case looks like "there is no problem
here" until an acquisition that carries a real function turns up.

**What this does not fix is the acquisition's own mass accuracy.** Both frames can be
displaced together, and on badly calibrated data they are: a calibration transferred
from another day and another polarity has been seen to put a known species 350 ppm out
at the bottom of the range and 155 ppm out at the top, while a calibrant acquired that
morning read within 4 ppm on the same instrument. `Recalibration` is the optional second
step for that, fitted from species whose m/z is known -- optional by decision, because
correcting an instrument's calibration is not something an extraction should do behind
the analyst's back, and never silent when it is on.

**Before fitting anything, ask the file why it is wrong.** Every acquisition records
which calibration it is carrying and which calibrant acquisition produced it, in
`$$ Cal MS1 Dynamic Params`, and records when each was taken. `CalibrationProvenance`
reads that, and on the acquisitions this was measured against it accounts for the whole
of the mass error without a single fit (lab record, task 11): the well behaved file was
calibrated from itself that morning and reads within 4 ppm, while the ones 120 to 350 ppm
out carry a calibration three months old -- one of them carried across polarity, which
MassLynx itself marks `<X+>` in that line. A correction not needed is worth more than a
correction fitted, so `concerns()` is the first thing to look at and `Recalibration` is
what is left when the answer is that the data is what it is.

`MzFrame` is the pair, and `MzFrame.reader_mz(mz)` is the whole point of the module:
the number to hand `ReadMobillogram` for a species of computed m/z. `describe()` on any
of these returns the one-line statement of what was applied, for the record, and
`MzFrame.note_for(mz)` adds what was true of that particular species -- notably whether
a fitted recalibration reached it at all.
"""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

from .extern import cal_function_coefficients, parse_extern_inf, parse_header_txt

__all__ = [
    "CalFunction",
    "CalibrationError",
    "CalibrationProvenance",
    "MzFrame",
    "Recalibration",
    "STALE_AFTER_DAYS",
]

# Newton on a function that is within a part in a thousand of the identity over any
# range these instruments acquire. Two passes converge to double precision; the extra
# ones cost nothing and cover a pathological calibration.
_INVERT_PASSES = 6
_INVERT_TOLERANCE_RELATIVE = 1.0e-12

STALE_AFTER_DAYS = 1.0
"""How old a calibration may be before `CalibrationProvenance.concerns` mentions it.

A convention, and the only number in this module that is not a measurement. What was
measured (lab record, task 11) is the two ends: a calibration fitted from the acquisition
itself, the same morning, leaves 4 ppm, and one carried for 92 days leaves 120 to 350 ppm.
Nothing in between was measured, because nothing in between was acquired. One day is the
interval an accurate-mass method is normally recalibrated on, so it is what `concerns`
defaults to; pass your own where your instrument's drift is known.
"""

_DYNAMIC_PARAMS = "Cal MS1 Dynamic Params"
# MassLynx prefixes both fields of that line with this when the calibration it applied
# was fitted in the other polarity. It is the one thing on the line that is a defect
# rather than a description, and it survives into the file unchanged.
_CROSS_POLARITY = "<X+>"
# The stamps are written in whatever format the workstation was set to. Acquisition
# stamps have been seen as dd-Mmm-yyyy and calibration stamps as mm/dd/yy on the same
# file, so both are tried against both, and a stamp matching none is kept as text
# rather than guessed at -- a date guessed the wrong way round is a silent wrong answer.
_DATE_FORMATS = ("%d-%b-%Y", "%m/%d/%y", "%d/%m/%Y", "%Y-%m-%d")
_TIME_FORMATS = ("%H:%M:%S", "%H:%M")
_NUMBER = re.compile(r"-?\d+(?:\.\d+)?")


class CalibrationError(Exception):
    """A calibration that cannot be applied, inverted, or fitted."""


def _as_array(mz) -> np.ndarray:
    values = np.asarray(mz, dtype=float)
    if not np.all(np.isfinite(values)):
        raise CalibrationError("m/z values must be finite")
    if np.any(values <= 0.0):
        raise CalibrationError("m/z values must be positive; sqrt(m/z) is the variable here")
    return values


def _restore_shape(result: np.ndarray, original) -> float | np.ndarray:
    return float(result) if np.isscalar(original) or np.ndim(original) == 0 else result


def _invert(forward, target: np.ndarray) -> np.ndarray:
    """The `mz` for which `forward(mz)` is `target`, by Newton with a numeric slope.

    Numeric rather than analytic because the forward map is a polynomial composed with
    a square and the derivative is not the interesting part of this module. Raises
    rather than returning a half-converged answer: a window placed at the wrong m/z is
    a silent wrong result, which is the failure mode this package exists to remove.
    """
    guess = target.copy()
    for _ in range(_INVERT_PASSES):
        step = np.maximum(1.0e-6, 1.0e-7 * guess)
        value = forward(guess)
        slope = (forward(guess + step) - value) / step
        if np.any(slope == 0.0):
            raise CalibrationError("the calibration function is not locally invertible")
        guess = guess - (value - target) / slope
    residual = np.abs(forward(guess) - target)
    if np.any(residual > _INVERT_TOLERANCE_RELATIVE * np.abs(target)):
        raise CalibrationError(
            "the calibration function did not invert; its coefficients are not a "
            "near-identity mass correction"
        )
    return guess


@dataclass(frozen=True)
class CalFunction:
    """One acquisition's `$$ Cal Function N`: the map from chromatogram to spectrum frame.

    `coefficients` are lowest order first, as `_HEADER.TXT` writes them, and are a
    polynomial in **sqrt(m/z)** whose square is the corrected mass. `function` is
    numbered from 1, as that file numbers it; the SDK numbers its functions from 0 and
    the caller converts.

    Build it with `from_raw`, which returns the identity for an acquisition that carries
    no calibration of its own and raises for one whose line cannot be read -- there is no
    third state, because a window has to be placed either way.
    """

    coefficients: tuple[float, ...]
    function: int = 1
    source: str = ""

    def __post_init__(self) -> None:
        if len(self.coefficients) < 2:
            raise CalibrationError(
                f"a calibration function needs at least a constant and a slope; got "
                f"{self.coefficients!r}"
            )
        if not all(np.isfinite(c) for c in self.coefficients):
            raise CalibrationError(f"non-finite calibration coefficients: {self.coefficients!r}")

    @classmethod
    def identity(cls, function: int = 1, source: str = "") -> CalFunction:
        """The function that leaves every mass where it is."""
        return cls(coefficients=(0.0, 1.0), function=function, source=source)

    @classmethod
    def from_raw(
        cls, raw_dir, function: int = 1, *, missing_is_identity: bool = True
    ) -> CalFunction:
        """The calibration function of a `.raw` acquisition directory.

        A file with no such line is the identity by default, which is what MassLynx
        means by writing `0,1` explicitly and what an acquisition with no calibration
        behaves as. Pass `missing_is_identity=False` where the difference matters and
        an absent line should raise.
        """
        coefficients = cal_function_coefficients(raw_dir, function)
        if coefficients is None:
            if missing_is_identity:
                return cls.identity(function=function, source=str(raw_dir))
            raise CalibrationError(
                f"{raw_dir} carries no readable '$$ Cal Function {function}' in its "
                "_HEADER.TXT"
            )
        return cls(coefficients=coefficients, function=function, source=str(raw_dir))

    @property
    def order(self) -> int:
        """The polynomial's order: 1 for a linear function, 3 for a cubic."""
        return len(self.coefficients) - 1

    @property
    def is_identity(self) -> bool:
        """Whether the two frames coincide, so that every conversion here is a no-op."""
        return self.coefficients == (0.0, 1.0)

    def to_spectrum_frame(self, mz):
        """Chromatogram frame in, spectrum frame out: the polynomial, squared."""
        values = _as_array(mz)
        root = np.sqrt(values)
        return _restore_shape(np.polyval(self.coefficients[::-1], root) ** 2, mz)

    def to_chromatogram_frame(self, mz):
        """Spectrum frame in, chromatogram frame out -- the frame a window is placed in."""
        values = _as_array(mz)
        if self.is_identity:
            return _restore_shape(values, mz)
        return _restore_shape(
            _invert(lambda m: np.asarray(self.to_spectrum_frame(m), dtype=float), values), mz
        )

    def describe(self) -> str:
        """One line for the record: what this is and whether it does anything."""
        if self.is_identity:
            return f"cal function {self.function}: identity (the two m/z frames coincide)"
        terms = ", ".join(f"{c:.12g}" for c in self.coefficients)
        return (
            f"cal function {self.function}: order {self.order} in sqrt(m/z), squared; "
            f"coefficients ({terms})"
        )


def _stamp(header: dict[str, str], date_key: str, time_key: str) -> _dt.datetime | None:
    """A `datetime` from one of `_HEADER.TXT`'s date/time pairs, or None."""
    text = f"{(header.get(date_key) or '').strip()} {(header.get(time_key) or '').strip()}"
    for date_format in _DATE_FORMATS:
        for time_format in _TIME_FORMATS:
            try:
                return _dt.datetime.strptime(text, f"{date_format} {time_format}")
            except ValueError:
                continue
    return None


@dataclass(frozen=True)
class CalibrationProvenance:
    """Where an acquisition's mass calibration came from, and how old it was.

    Read from `_HEADER.TXT` and `_extern.inf` alone -- no SDK, no license, no scan read
    -- so this is available before anything has been extracted, which is when it is
    useful. `calibration` and `calibrant` are the last two fields of
    `$$ Cal MS1 Dynamic Params`: what MassLynx called the calibration and which
    acquisition it was fitted from. `cross_polarity` is MassLynx's own `<X+>` marker on
    those fields, meaning the calibration was carried over from the other polarity.

    `concerns()` is the point of the class. It returns a tuple of one-line complaints
    rather than raising, the same convention `experiment.Experiment.validate` follows,
    because a stale calibration is not a reason to refuse to read a file -- it is a
    reason to expect a `Recalibration` to be needed, and to say so in the record.

    A field MassLynx did not write, or wrote in a format that does not parse, is None
    rather than a guess.
    """

    calibration: str | None = None
    calibrant: str | None = None
    acquired: _dt.datetime | None = None
    calibrated: _dt.datetime | None = None
    polarity: str | None = None
    mass_range: tuple[float, float] | None = None
    cross_polarity: bool = False
    source: str = ""

    @classmethod
    def from_raw(cls, raw_dir) -> CalibrationProvenance:
        """The provenance of a `.raw` acquisition directory. Never raises.

        A file with none of these lines gives an object whose fields are all None, and
        whose `concerns()` says the file records nothing about its calibration -- which
        is itself worth knowing and is not the same as a calibration that is fine.
        """
        header = parse_header_txt(raw_dir)
        raw = header.get(_DYNAMIC_PARAMS)
        calibration = calibrant = None
        mass_range = None
        if raw is not None:
            fields = [field.strip() for field in raw.split(",")]
            calibration = (fields[-2] if len(fields) >= 2 else "") or None
            calibrant = (fields[-1] if fields else "") or None
            if len(fields) >= 2 and _NUMBER.fullmatch(fields[0] or ""):
                mass_range = (float(fields[0]), float(fields[1]))
        return cls(
            calibration=calibration,
            calibrant=calibrant,
            acquired=_stamp(header, "Acquired Date", "Acquired Time"),
            calibrated=_stamp(header, "Cal Date", "Cal Time"),
            polarity=(parse_extern_inf(raw_dir).get("Polarity") or "").strip() or None,
            mass_range=mass_range,
            cross_polarity=bool(
                (calibration or "").startswith(_CROSS_POLARITY)
                or (calibrant or "").startswith(_CROSS_POLARITY)
            ),
            source=str(raw_dir),
        )

    @property
    def age_days(self) -> float | None:
        """How old the calibration was when this file was acquired, in days."""
        if self.acquired is None or self.calibrated is None:
            return None
        return (self.acquired - self.calibrated).total_seconds() / 86400.0

    def concerns(self, *, stale_after_days: float = STALE_AFTER_DAYS) -> tuple[str, ...]:
        """One line per reason to expect this acquisition's masses to be off.

        Empty is the good answer, and it means the file's own record gives no reason to
        doubt its calibration -- not that the calibration was checked, which only a
        species of known m/z can do.
        """
        found: list[str] = []
        if self.calibration is None and self.calibrant is None:
            found.append(
                "the file records no mass calibration at all: no "
                f"'$$ {_DYNAMIC_PARAMS}' line in its _HEADER.TXT"
            )
        if self.cross_polarity:
            found.append(
                f"the calibration ({self.calibration}) was carried over from the other "
                "polarity, which is MassLynx's own '<X+>' marker; expect the mass error "
                "to be large and to change sign across the range"
            )
        age = self.age_days
        if age is None and (self.calibration or self.calibrant):
            found.append("the calibration's date or the acquisition's did not parse, so its age is unknown")
        elif age is not None and age > stale_after_days:
            found.append(
                f"the calibration is {age:.1f} days old at acquisition (fitted from "
                f"{self.calibrant}), past the {stale_after_days:g} day this was asked with"
            )
        elif age is not None and age < 0.0:
            found.append(
                f"the calibration is dated {-age:.1f} days after the acquisition, which "
                "is not a thing that can happen; one of the two stamps is wrong"
            )
        return tuple(found)

    def describe(self) -> str:
        """One line for the record."""
        age = self.age_days
        parts = [f"calibration {self.calibration!r}" if self.calibration else "no calibration recorded"]
        if self.calibrant:
            parts.append(f"from {self.calibrant!r}")
        if age is not None:
            parts.append(f"{age:.1f} days old at acquisition")
        if self.polarity:
            parts.append(f"acquired {self.polarity}")
        if self.cross_polarity:
            parts.append("carried across polarity")
        return ", ".join(parts)


@dataclass(frozen=True)
class Recalibration:
    """A fitted correction for the acquisition's own mass error. Optional, and never silent.

    `coefficients` are a polynomial in sqrt(m/z), lowest order first, mapping a species'
    **true** sqrt(m/z) to the sqrt(m/z) it is **found at** in the spectrum frame -- that
    direction, because placing a window is the forward use and the residual then lands
    in the frame the window is placed in. `to_true` inverts it.

    `fit` takes the two m/z series and the order. Order 1 is the default and is not a
    placeholder: a time-of-flight instrument's own mass error, where it has one, is
    dominated by the time-zero and gain terms that a linear function in sqrt(m/z)
    absorbs, and on the acquisitions this was measured against orders 2 to 5 lowered the
    residual by less than the centroid precision. Raise the order when a fit says to,
    which is what `rms_th` and `max_th` are for.
    """

    coefficients: tuple[float, ...]
    order: int
    rms_th: float
    max_th: float
    species: int
    mz_range: tuple[float, float]
    labels: tuple[str, ...] = field(default=())
    offered: int = 0
    """How many species were put forward, where that differs from how many were fitted.

    `from_spectrum` sets it: a species with no measurable peak is skipped rather than
    fitted, and a fit from four of forty species is a different claim from a fit from
    forty. Zero means it was not recorded, which is what `fit` leaves it as.
    """

    @classmethod
    def fit(
        cls,
        true_mz,
        measured_mz,
        *,
        order: int = 1,
        labels: tuple[str, ...] = (),
        offered: int = 0,
    ) -> Recalibration:
        """Fit sqrt(m_measured) against sqrt(m_true) at `order`.

        `true_mz` is computed from formulae; `measured_mz` is where each species was
        found on the spectrum readers' axis. Both are in the same order and the same
        units. Needs at least two more species than the order, so that the residual
        means something rather than being zero by construction.
        """
        true = _as_array(true_mz)
        measured = _as_array(measured_mz)
        if true.shape != measured.shape or true.ndim != 1:
            raise CalibrationError(
                f"true_mz and measured_mz must be one-dimensional and the same length; "
                f"got {true.shape} and {measured.shape}"
            )
        if order < 1:
            raise CalibrationError(f"a recalibration needs order 1 or more; got {order}")
        if true.size < order + 2:
            raise CalibrationError(
                f"an order-{order} recalibration needs at least {order + 2} species to "
                f"leave a residual worth reading; got {true.size}"
            )
        if labels and len(labels) != true.size:
            raise CalibrationError("labels, where given, must match the species one for one")

        coefficients = np.polyfit(np.sqrt(true), np.sqrt(measured), order)
        residual_root = np.sqrt(measured) - np.polyval(coefficients, np.sqrt(true))
        # Back to Th, which is the unit a window is chosen in: d(m) = 2*sqrt(m)*d(sqrt m).
        residual_th = 2.0 * np.sqrt(measured) * residual_root
        return cls(
            coefficients=tuple(float(c) for c in coefficients[::-1]),
            order=int(order),
            rms_th=float(np.sqrt(np.mean(residual_th**2))),
            max_th=float(np.abs(residual_th).max()),
            species=int(true.size),
            mz_range=(float(true.min()), float(true.max())),
            labels=tuple(labels),
            offered=int(offered),
        )

    @classmethod
    def from_spectrum(
        cls,
        spectrum,
        species: Iterable[float | Sequence],
        *,
        order: int = 1,
        min_apex: float = 0.0,
    ) -> Recalibration:
        """Fit from where a total spectrum actually puts a set of species of known m/z.

        `species` is an iterable of m/z, or of `(mz, label)` pairs. Each is looked for in
        `spectrum` with `MassSpectrum.peak_near`, and the ones that have a peak clearing
        `min_apex` are what the fit is made from; the rest are skipped and counted in
        `offered`. That skipping is not a convenience -- a reference list for a real
        sample will name species that are not in every acquisition, and a fit is worth
        more than a refusal.

        **The reference species are normally the ones the analysis is already looking
        for** (decision of record, lab record task 11), which is what makes this need no
        list from the analyst: an oligomer series, a protein's charge-state envelope,
        whatever the windows are being placed on. That looks circular and is not, for
        one reason: `peak_near` only looks within its own search half-width, so a species
        whose peak is not nearly where it was expected contributes nothing rather than
        contributing a wrong anchor, and `rms_th` and `max_th` over the whole set are
        what say whether the assignment held. Read them. A residual far larger than the
        peaks are wide is a misassignment, not a mass error.

        A calibrant acquisition is the other anchor, and it is not used through this
        method alone -- see `MzFrame.from_calibrant`, which explains why.
        """
        pairs = [(float(s[0]), str(s[1])) if isinstance(s, Sequence) else (float(s), "") for s in species]
        if not pairs:
            raise CalibrationError("a recalibration needs species to fit from; none were given")
        true, measured, labels = [], [], []
        for mz, label in pairs:
            peak = spectrum.peak_near(mz, min_apex=min_apex)
            if peak is None:
                continue
            true.append(mz)
            measured.append(peak.mz)
            labels.append(label)
        if len(true) < order + 2:
            raise CalibrationError(
                f"only {len(true)} of {len(pairs)} species have a peak clearing "
                f"{min_apex:g} counts in "
                f"{getattr(spectrum, 'acquisition', '') or 'this spectrum'}; an order-"
                f"{order} recalibration needs at least {order + 2}"
            )
        return cls.fit(
            true, measured, order=order, labels=tuple(labels), offered=len(pairs)
        )

    def to_measured(self, mz):
        """True m/z in, the m/z it is found at in the spectrum frame out."""
        values = _as_array(mz)
        root = np.sqrt(values)
        return _restore_shape(np.polyval(self.coefficients[::-1], root) ** 2, mz)

    def to_true(self, mz):
        """The inverse: a measured spectrum-frame m/z back to what the species is."""
        values = _as_array(mz)
        return _restore_shape(
            _invert(lambda m: np.asarray(self.to_measured(m), dtype=float), values), mz
        )

    def covers(self, mz) -> bool:
        """Whether `mz` is inside the range this was fitted over.

        **It reports; it does not refuse** (decision of record, lab record task 11).
        Extrapolating one of these fits was measured rather than assumed, on a
        calibrant spanning 393-1172 Th applied to species reaching past both ends: the
        extrapolated species were placed *better* than the interpolated ones, 0.16
        against 0.34 of the half-width of the window they went into. An order-1 function
        in sqrt(m/z) is a time-zero and gain correction, and those do not stop being true
        outside the range they were fitted over.

        What does go wrong is an anchor that does not belong to the acquisition -- a
        calibrant of the wrong polarity put 8 % of species outside their window at every
        m/z, covered and extrapolated alike -- and no range check catches that. So
        refusing here would cost the cases that work and would not buy the case that
        fails. `MzFrame.note_for` puts the answer in the window's record instead, which
        is where a reader can weigh it.
        """
        low, high = self.mz_range
        values = np.asarray(mz, dtype=float)
        return bool(np.all((values >= low) & (values <= high)))

    def describe(self) -> str:
        """One line for the record."""
        offered = (
            f" of {self.offered} offered" if self.offered and self.offered != self.species else ""
        )
        return (
            f"recalibration: order {self.order} in sqrt(m/z), {self.species} species"
            f"{offered} over {self.mz_range[0]:.1f}-{self.mz_range[1]:.1f} Th, residual "
            f"{self.rms_th:.4f} Th rms and {self.max_th:.4f} Th max"
        )


@dataclass(frozen=True)
class MzFrame:
    """How a computed m/z becomes a bound `ReadMobillogram` can be given.

    Two steps, and the first is not optional:

    1. the acquisition's own mass error, if a `Recalibration` was fitted and passed --
       computed m/z to where the species is actually found;
    2. `CalFunction`, always -- spectrum frame to chromatogram frame.

    With no recalibration this trusts the acquisition's calibration, which is the right
    default: the correction is a conversion the file carries, and inventing a mass
    correction the analyst did not ask for is not extraction's business (decision of
    record, lab record task 10).
    """

    cal: CalFunction
    recalibration: Recalibration | None = None

    @classmethod
    def for_acquisition(cls, raw_dir, *, recalibration: Recalibration | None = None) -> MzFrame:
        """The frame of a `.raw` directory, from its own `_HEADER.TXT`."""
        return cls(cal=CalFunction.from_raw(raw_dir), recalibration=recalibration)

    @classmethod
    def from_calibrant(
        cls,
        calibrant_raw_dir,
        recalibration: Recalibration,
        *,
        for_acquisition=None,
        allow_polarity_mismatch: bool = False,
    ) -> MzFrame:
        """The calibrant acquisition's **whole** frame, to place windows in another file.

        This is how a calibrant run is used, and the reason it is a constructor of its
        own rather than a note in a docstring is that the obvious way is wrong (lab
        record, task 11, measured on a calibrant acquired twenty minutes after the
        series it was tested against).

        A `Recalibration` fitted from a calibrant measures what is left *after* that
        file's own calibration function, and a calibrant is normally calibrated from
        itself, so its residual is a few ppm of nothing. Handing that residual to an
        analyte acquisition -- `MzFrame(cal=analyte_cal, recalibration=calibrant_recal)`
        -- therefore corrects almost nothing, because the analyte's error is in its
        calibration and not in its residual. Measured: it left every species of a cation
        series 3.6 half-widths outside its window, which is where they started.

        Taking the calibrant's `CalFunction` **as well** is the transfer that works,
        because together the two are the whole map from the instrument's own axis to
        true m/z -- which is what calibrating the analyte from the calibrant means.
        Measured on the same series: 0.34 half-widths, every species inside its window,
        against 0.22 for a fit from the analyte itself.

        `for_acquisition`, where given, is the `.raw` the frame is going to be used on,
        and it is checked: **a calibrant of the other polarity is refused**, because that
        transfer was measured too and it fails -- 8 % of species outside their window,
        at every m/z, whether or not the calibrant's range reached them. It is the same
        defect MassLynx marks `<X+>`. Pass `allow_polarity_mismatch=True` to do it
        anyway, on an instrument where you have shown it holds.
        """
        cal = CalFunction.from_raw(calibrant_raw_dir)
        if for_acquisition is not None:
            calibrant = CalibrationProvenance.from_raw(calibrant_raw_dir)
            analyte = CalibrationProvenance.from_raw(for_acquisition)
            mismatched = (
                calibrant.polarity
                and analyte.polarity
                and calibrant.polarity != analyte.polarity
            )
            if mismatched and not allow_polarity_mismatch:
                raise CalibrationError(
                    f"the calibrant {calibrant_raw_dir} was acquired {calibrant.polarity} "
                    f"and {for_acquisition} was acquired {analyte.polarity}; carrying a "
                    "calibration across polarity was measured to leave species outside "
                    "their windows at every m/z. Acquire a calibrant in each polarity, "
                    "or pass allow_polarity_mismatch=True if you have shown otherwise"
                )
        return cls(cal=cal, recalibration=recalibration)

    @classmethod
    def identity(cls) -> MzFrame:
        """The frame that changes nothing -- for tests, and for data with no calibration."""
        return cls(cal=CalFunction.identity())

    def reader_mz(self, mz):
        """The m/z to hand `ReadMobillogram` for a species of computed m/z `mz`."""
        placed = self.recalibration.to_measured(mz) if self.recalibration else mz
        return self.cal.to_chromatogram_frame(placed)

    def computed_mz(self, mz):
        """The inverse: a position in the chromatogram frame back to a species' m/z."""
        spectrum = self.cal.to_spectrum_frame(mz)
        return self.recalibration.to_true(spectrum) if self.recalibration else spectrum

    @property
    def is_identity(self) -> bool:
        """Whether this frame moves nothing at all."""
        return self.cal.is_identity and self.recalibration is None

    def describe(self) -> str:
        """One line for the record, naming both steps."""
        parts = [self.cal.describe()]
        parts.append(
            self.recalibration.describe() if self.recalibration else "recalibration: none"
        )
        return "; ".join(parts)

    def note_for(self, mz) -> str:
        """`describe()`, plus what was true of this particular species.

        The only thing that varies species by species is whether a fitted recalibration
        reached it, and `Recalibration.covers` reports rather than refuses, so this is
        where that report goes. `extract.MzWindow` uses it for `frame_note` and `atd.ATD`
        carries it onward, which is the whole of the "never silent" rule: a window
        extrapolated past its fit says so in the record it travels with, and a reader
        who thinks that matters can find every one of them.
        """
        note = self.describe()
        if self.recalibration is not None and not self.recalibration.covers(mz):
            low, high = self.recalibration.mz_range
            note += (
                f"; {float(mz):.4f} Th is outside the {low:.1f}-{high:.1f} Th the "
                "recalibration was fitted over, so it is extrapolated"
            )
        return note

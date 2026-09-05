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

`MzFrame` is the pair, and `MzFrame.reader_mz(mz)` is the whole point of the module:
the number to hand `ReadMobillogram` for a species of computed m/z. `describe()` on any
of these returns the one-line statement of what was applied, for the record.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .extern import cal_function_coefficients

__all__ = [
    "CalFunction",
    "CalibrationError",
    "MzFrame",
    "Recalibration",
]

# Newton on a function that is within a part in a thousand of the identity over any
# range these instruments acquire. Two passes converge to double precision; the extra
# ones cost nothing and cover a pathological calibration.
_INVERT_PASSES = 6
_INVERT_TOLERANCE_RELATIVE = 1.0e-12


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

    @classmethod
    def fit(
        cls,
        true_mz,
        measured_mz,
        *,
        order: int = 1,
        labels: tuple[str, ...] = (),
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
        """Whether `mz` is inside the range this was fitted over. Extrapolation is not."""
        low, high = self.mz_range
        return bool(np.all((np.asarray(mz, dtype=float) >= low) & (np.asarray(mz) <= high)))

    def describe(self) -> str:
        """One line for the record."""
        return (
            f"recalibration: order {self.order} in sqrt(m/z), {self.species} species over "
            f"{self.mz_range[0]:.1f}-{self.mz_range[1]:.1f} Th, residual "
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

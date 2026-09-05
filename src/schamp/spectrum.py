"""A mass spectrum, and the one measurement an analysis needs from it: where a peak is
and how wide.

`extract.total_spectrum` produces these -- the acquisition summed over every retention
scan and every drift bin, which is one SDK call -- and everything here runs without the
SDK, so a spectrum read once can be measured, saved and re-measured anywhere.

**A `MassSpectrum` is always in the spectrum frame**, the m/z axis `ReadScan`,
`ReadDriftScan` and the scan processor return. It is not the frame `ReadMobillogram`
takes its bounds in; `calibration` converts, and `extract.MzWindow.around_measured` is
what does the converting in practice. Nothing here converts anything, because a spectrum
that had been quietly moved into another frame would be impossible to check afterwards.

Two things the peak measurement is shaped by, and both are measured properties of these
instruments rather than choices of style (lab record, task 10):

* **A centroid is taken over the points at or above half the apex, not over the whole
  peak.** The next isotopologue is 1.003 Th away at 1+ and 0.5 Th at 2+, and a centroid
  taken down to the baseline is pulled towards it -- which then looks like a mass error
  and is really an envelope shape.
* **Resolving power is not a constant.** On the 2013 acquisitions it doubles from about
  12,000 at 230 Th to 26,000 near 1000 Th and falls back to 13,500 by 1650 Th, so a peak
  is 2.7 time-of-flight channels wide over most of the range and 6.9 at the top. That is
  why a width in Th, a width in channels and a width in resolving-power units all fail
  at one end or the other, and why `SpectrumPeak.fwhm` exists: the honest width for a
  window is a multiple of the width the peak actually has.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["MassSpectrum", "SpectrumPeak"]

# How far either side of an expected m/z `peak_near` looks, relative to m/z. About
# 0.3 Th at 800: wider than any frame displacement seen on a working instrument and
# narrower than half the 2+ isotope spacing, so the tallest thing in the search window
# is the species being asked for or a neighbour that a caller needs to know about.
SEARCH_HALFWIDTH_RELATIVE = 3.75e-4


@dataclass(frozen=True)
class SpectrumPeak:
    """One peak of a mass spectrum, measured rather than fitted.

    `mz` is the intensity-weighted centroid of the points at or above half the apex, in
    the spectrum frame. `fwhm` is the full width at half maximum with both crossings
    interpolated between the points either side; where the profile never reaches half
    height on one side -- an unresolved shoulder, which the high-mass end of a
    polyalanine envelope is full of -- the last point above half is used instead and
    `resolved` is False, so the width is a lower bound and says so.

    No Gaussian is fitted. `atd.fit_gaussian` fits the arrival-time distribution, where
    the shape is the measurement; here the shape is not, and a fit would only add a
    model to a number that three points already give.
    """

    mz: float
    apex_mz: float
    apex: float
    fwhm: float
    area_above_half: float
    points_above_half: int
    resolved: bool = True

    @property
    def resolving_power(self) -> float:
        """m/z over the full width at half maximum. `inf` if the width came out zero."""
        return self.mz / self.fwhm if self.fwhm > 0.0 else float("inf")


@dataclass(frozen=True)
class MassSpectrum:
    """One profile mass spectrum in the spectrum frame, and where its peaks are.

    `mz` and `intensity` are parallel and ascending in m/z. `acquisition` is the `.raw`
    it came from and `summed_over` says what was added together to make it -- "all
    retention scans, all drift bins" for a total spectrum -- because a spectrum with no
    record of what it covers cannot be compared with another one.

    A spectrum may legitimately be empty or nearly so; `is_empty` says which, and nothing
    here raises for it.

    The axis is the acquisition's own: time-of-flight channels, so the spacing grows as
    the square root of m/z, and channels holding nothing are absent rather than zero. Do
    not resample it. A window narrower than one channel is not narrower than one channel,
    and a profile put on a uniform grid finer than the channel spacing produces peak
    widths that are an artefact of the grid.
    """

    mz: np.ndarray
    intensity: np.ndarray
    acquisition: str = ""
    summed_over: str = ""

    def __post_init__(self) -> None:
        if len(self.mz) != len(self.intensity):
            raise ValueError(
                f"mz and intensity must be the same length; got {len(self.mz)} and "
                f"{len(self.intensity)}"
            )
        if np.any(np.diff(self.mz) <= 0.0):
            raise ValueError("a spectrum's m/z axis must ascend; this one does not")

    @property
    def total(self) -> float:
        """Summed intensity over the whole spectrum."""
        return float(np.sum(self.intensity))

    @property
    def is_empty(self) -> bool:
        """True when there is nothing here to measure.

        A real and ordinary outcome, not a malformed object: one drift bin of one
        retention scan of a mobility acquisition can hold a handful of points or none at
        all, and a narrow combine is a legitimate thing to ask for. Callers check this
        the way they check `atd.ATD.is_empty`; `peak_near` returns None rather than
        raising, for the same reason.
        """
        return len(self.mz) == 0 or not np.any(self.intensity > 0.0)

    @property
    def mz_range(self) -> tuple[float, float]:
        """First and last m/z on the axis. Raises on a spectrum with no points."""
        if len(self.mz) == 0:
            raise ValueError("an empty spectrum has no m/z range; check is_empty first")
        return float(self.mz[0]), float(self.mz[-1])

    def mz_step_at(self, mz: float) -> float:
        """The acquisition's own m/z spacing near `mz`: one time-of-flight channel.

        The median spacing over the points within a search half-width, so a gap where a
        channel held nothing does not become the answer. This is the floor on any window
        width, and the quantity a bound is snapped to -- see `extract`.
        """
        halfwidth = abs(float(mz)) * SEARCH_HALFWIDTH_RELATIVE
        inside = (self.mz > mz - halfwidth) & (self.mz < mz + halfwidth)
        if np.count_nonzero(inside) < 2:
            raise ValueError(f"the spectrum has fewer than two points near {mz:g} Th")
        return float(np.median(np.diff(self.mz[inside])))

    def peak_near(
        self,
        mz: float,
        *,
        halfwidth: float | None = None,
        min_apex: float = 0.0,
    ) -> SpectrumPeak | None:
        """The tallest peak within `halfwidth` of `mz`, measured. None if there is none.

        `halfwidth` defaults to `SEARCH_HALFWIDTH_RELATIVE` of `mz`. `min_apex` is the
        apex intensity below which a peak is not worth measuring -- a summed spectrum has
        a baseline, and a centroid of the baseline is a number with no meaning. It
        returns None rather than raising, because "this species is not in this
        acquisition" is an ordinary answer, and the caller decides what that means.

        The walk out from the apex stops at half height **or** at a minimum, whichever
        comes first, so an unresolved shoulder ends the peak instead of joining it.
        """
        if halfwidth is None:
            halfwidth = abs(float(mz)) * SEARCH_HALFWIDTH_RELATIVE
        if halfwidth <= 0.0:
            raise ValueError(f"a search half-width must be positive; got {halfwidth}")

        inside = (self.mz > mz - halfwidth) & (self.mz < mz + halfwidth)
        if not inside.any():
            return None
        apex_index = int(np.flatnonzero(inside)[np.argmax(self.intensity[inside])])
        apex = float(self.intensity[apex_index])
        if apex <= 0.0 or apex < min_apex:
            return None

        half = 0.5 * apex
        last = len(self.mz) - 1
        low = apex_index
        while (
            low > 0
            and self.intensity[low - 1] >= half
            and self.intensity[low - 1] < self.intensity[low]
        ):
            low -= 1
        high = apex_index
        while (
            high < last
            and self.intensity[high + 1] >= half
            and self.intensity[high + 1] < self.intensity[high]
        ):
            high += 1

        counts = self.intensity[low : high + 1]
        centroid = float((self.mz[low : high + 1] * counts).sum() / counts.sum())
        left, left_resolved = self._half_crossing(low, low - 1, half)
        right, right_resolved = self._half_crossing(high, high + 1, half)
        return SpectrumPeak(
            mz=centroid,
            apex_mz=float(self.mz[apex_index]),
            apex=apex,
            fwhm=right - left,
            area_above_half=float(counts.sum()),
            points_above_half=int(high - low + 1),
            resolved=bool(left_resolved and right_resolved),
        )

    def _half_crossing(self, inner: int, outer: int, half: float) -> tuple[float, bool]:
        """Where the profile crosses `half` between two points, and whether it does.

        Only interpolates when the outer point is genuinely below half height and below
        the inner one. Where the walk stopped at a valley that never got there, the
        linear solution points the wrong way and returns a width smaller than the peak or
        a negative one; there is no crossing to find, so the last point above half is the
        answer and the caller is told the width is a lower bound.
        """
        if outer < 0 or outer > len(self.mz) - 1:
            return float(self.mz[inner]), False
        if self.intensity[outer] >= half or self.intensity[outer] >= self.intensity[inner]:
            return float(self.mz[inner]), False
        span = (half - self.intensity[inner]) / (self.intensity[outer] - self.intensity[inner])
        return float(self.mz[inner] + span * (self.mz[outer] - self.mz[inner])), True

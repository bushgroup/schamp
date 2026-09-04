"""Arrival-time distributions, and the single peak fitted to one.

An ATD is the intensity of one m/z window against drift time, summed over every
retention-time scan of an acquisition. It is what the extraction layer produces and
what the mobility layer consumes one centroid of, so it is the interface between the
only SDK-dependent module and everything else -- which is why the container lives
here, with no SDK import anywhere in sight.

The container carries **both axes**, drift bin and drift time in milliseconds, because
the 2013 analysis lost the distinction: it looked for a pusher-interval key that this
instrument does not write, silently defaulted the period to 1.0, and produced every
centroid in bins while labelling them milliseconds. Both axes present and named makes
that class of error impossible to make quietly, and lets a fit be compared against the
legacy numbers in the units the legacy actually used.

The fit itself -- a single Gaussian, its moments and its quality metrics -- is task 05
and is a stub here. Its result fields are fixed now, because task 06 and the figures
consume them: it reproduces the eight summary columns the legacy `FIT.py` wrote, so
that the modern pipeline can be checked against the surviving 2013 outputs column for
column before it is trusted to replace them.

No SDK.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ATD", "FitResult", "Moments", "fit_gaussian", "moments"]


@dataclass(frozen=True)
class ATD:
    """One m/z window's arrival-time distribution from one acquisition.

    `drift_bin`, `drift_time_ms` and `intensity` are parallel arrays of equal length,
    one entry per drift bin, in acquisition order. `drift_time_ms` comes from the SDK's
    own drift-time axis, never from a computed bin width.

    `mz_low` and `mz_high` are the window as the extraction asked for it, and `label`
    names what the window is meant to contain (`"n=22, 2-"`), because a bare m/z range
    is not a species. `acquisition` is the `.raw` this came from.
    """

    drift_bin: np.ndarray
    drift_time_ms: np.ndarray
    intensity: np.ndarray
    mz_low: float
    mz_high: float
    acquisition: str = ""
    label: str = ""

    def __post_init__(self) -> None:
        lengths = {len(self.drift_bin), len(self.drift_time_ms), len(self.intensity)}
        if len(lengths) != 1:
            raise ValueError(
                "drift_bin, drift_time_ms and intensity must be the same length; got "
                f"{len(self.drift_bin)}, {len(self.drift_time_ms)}, {len(self.intensity)}"
            )
        if not lengths.pop():
            raise ValueError("an ATD with no bins is not an ATD")
        if self.mz_high <= self.mz_low:
            raise ValueError(
                f"the m/z window must have positive width; got {self.mz_low} to {self.mz_high}"
            )

    @property
    def mz_centre(self) -> float:
        """The window's centre in Th. What `FitResult.mz` reports, and only that.

        It is the centre of the *window*, not a measured mass: a window is placed
        around an expected m/z, and how well the ion sits in it is a question for the
        spectrum, not the ATD.
        """
        return 0.5 * (self.mz_low + self.mz_high)

    @property
    def total(self) -> float:
        """Summed intensity over the window. The legacy `Area` column is this, not a fit."""
        return float(np.sum(self.intensity))

    @property
    def peak_bin(self) -> int:
        """The drift bin of the highest point -- the starting guess, not the answer.

        A centroid comes from `fit_gaussian`; this is for sanity checks and for the
        bin-for-bin comparisons against the legacy extraction.
        """
        return int(self.drift_bin[int(np.argmax(self.intensity))])

    @property
    def is_empty(self) -> bool:
        """True when nothing landed in the window.

        An inverted window -- a `low,high` pair written the wrong way round, which the
        2013 peaks file managed three times -- returns all zeros rather than failing,
        so an empty ATD is a real and silent outcome that a caller must check for.
        """
        return not np.any(self.intensity > 0.0)


@dataclass(frozen=True)
class Moments:
    """Intensity-weighted moments of an ATD: the initial guess, and a fit-free fallback.

    Times in milliseconds. `centre` is the intensity-weighted mean drift time, `sigma`
    the weighted standard deviation, `height` the maximum intensity, `area` the summed
    intensity. Computed without any fitting, so they exist even for a distribution too
    ugly to fit.
    """

    centre_ms: float
    sigma_ms: float
    height: float
    area: float


@dataclass(frozen=True)
class FitResult:
    """A single Gaussian fitted to one ATD, plus how badly it fits.

    The first eight fields are the legacy `FIT.py` summary columns, so that the modern
    fit can be compared against the 2013 outputs directly. Two differences are
    deliberate and both are in the names: times are in **milliseconds** here, with
    `centre_bins` carried alongside for the comparison, and `area` is the summed
    intensity of the data rather than the integral of the fit -- which is what the
    legacy column was, despite its name.

    `resolution` is `centre / FWHM` with `FWHM = sigma * sqrt(8 ln 2)`, the legacy
    definition. It is reported, not analysed: peak-width and resolving-power analysis
    is out of scope by decision of record.
    """

    mz: float
    centre_ms: float
    sigma_ms: float
    height: float
    area: float
    resolution: float
    rmsd_fraction: float
    area_error_fraction: float

    centre_bins: float
    centre_ms_err: float | None = None
    sigma_ms_err: float | None = None
    converged: bool = False


def moments(atd: ATD) -> Moments:
    """Intensity-weighted moments of an ATD, in milliseconds. Task 05 owns the fit;
    these are needed as its starting guess and are cheap enough to have now.

    Negative intensities are clipped to zero before weighting: the Waters ADC baseline
    can go slightly negative, and a negative weight makes a nonsense of a moment.
    """
    weights = np.clip(np.asarray(atd.intensity, dtype=float), 0.0, None)
    total = float(np.sum(weights))
    if total <= 0.0:
        raise ValueError(
            f"no positive intensity in the window {atd.mz_low}-{atd.mz_high} Th of "
            f"{atd.acquisition or 'this acquisition'}; nothing to take moments of"
        )
    times = np.asarray(atd.drift_time_ms, dtype=float)
    centre = float(np.sum(weights * times) / total)
    variance = float(np.sum(weights * (times - centre) ** 2) / total)
    return Moments(
        centre_ms=centre,
        sigma_ms=float(np.sqrt(max(variance, 0.0))),
        height=float(np.max(weights)),
        area=total,
    )


def fit_gaussian(atd: ATD) -> FitResult:
    """Fit a single Gaussian to an ATD and report it with its quality metrics.

    Not implemented: task 05. The signature and `FitResult` are fixed now so that
    nothing downstream has to move when it lands.

    The intended method is the legacy one, deliberately, because reproducing the 2013
    summary columns is how the modern fit gets validated: moments for the initial
    guess, then least squares on `height * exp(-(t - centre)^2 / (2 sigma^2))` with no
    baseline term. What must be added is convergence reporting and parameter
    uncertainties, neither of which the legacy fit produced.
    """
    raise NotImplementedError(
        "single-Gaussian ATD fitting is task 05; schamp.atd.moments works today"
    )

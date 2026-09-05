"""Mobility: centroids against drift voltage, to K, K0 and a cross section.

Pure physics, no SDK, no file reads. It takes a fitted arrival time per acquisition,
the drift voltage the instrument profile computes for each, and the pressure and
temperature from the conditions table, and it returns a mobility and a cross section
with their uncertainties. Everything it needs about the hardware arrives in an
`InstrumentProfile`; the drift length and the divider ratio are never written here.

The chain is eqns (2) to (4) of the paper:

    v = K E = K V / L                                             (2)
    t_D = (L^2 / K)(1 / V) + t_0                                  (3)
    Omega = (3 z e / 16 N) (2 pi / (mu k_B T))^1/2 (1 / K)         (4)

so a straight line through arrival time against reciprocal drift voltage has slope
`L^2 / K` and intercept `t_0`, the transport time from the cell exit to the TOF
analyser. `t_0` is a fitted nuisance parameter and is reported, because a value far
from the expected few hundred microseconds is the clearest sign that something is
wrong with a series.

The whole chain lives here: `regress` fits the line, `mobility_from_slope` inverts
its slope, and `cross_section_from_regression` carries the result through eqn (4) to a
`CrossSection` -- the object that exists to enforce one rule: **a cross section is
never quoted without its gas, its charge state, and the pressure and temperature it
was reduced with.** `CrossSection` cannot be constructed without all four.

What the error bar covers
-------------------------

Eqn (4) is a product of powers, so every uncertainty enters as a relative one and the
whole propagation is four logarithmic derivatives:

    Omega  ~  slope^1 . L^-2 . P^-1 . T^+1/2

The slope term is `L^2 / K`, so a fractional error on the slope is a fractional error
on Omega, one for one. `L` enters squared. `P` and `T` enter through the number
density `N = P / (k_B T)` and, for `T`, again through the square root of eqn (4) --
`T^-1 . T^-1/2` from N and the root gives `T^+1/2` once the `1/K` is written out.
`cross_section_from_regression` combines whichever of the four the caller supplies in
quadrature and records their names in `CrossSection.propagated`, because an error bar
whose contents are not written down is not a measurement of anything.

The 2013 workbooks propagated the slope alone, and the conditions table makes the
other three optional and reserved rather than required (lab record, task 03): the
2013 series has a temperature with no recorded provenance at all, so demanding a
`temperature_K_err` would block the very experiment this package was written to
reproduce. Supply what you measured; the error bar then says what it covers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from . import constants
from .profiles import InstrumentProfile

__all__ = [
    "CrossSection",
    "DriftPoint",
    "Regression",
    "cross_section",
    "cross_section_from_regression",
    "mobility_from_slope",
    "regress",
]


@dataclass(frozen=True)
class DriftPoint:
    """One acquisition's contribution to a drift-voltage series.

    `drift_time_ms` is the fitted centroid of the arrival-time distribution and
    `v_drift_v` the drift voltage from the instrument profile. `acquisition` identifies
    the `.raw` so that an outlying point can be traced back to a file.
    """

    acquisition: str
    v_drift_v: float
    drift_time_ms: float
    drift_time_ms_err: float | None = None

    def __post_init__(self) -> None:
        if not (self.v_drift_v > 0.0):
            raise ValueError(
                f"{self.acquisition or 'this acquisition'} has a drift voltage of "
                f"{self.v_drift_v} V; a linear-field series needs a positive drift voltage"
            )

    @property
    def inverse_v(self) -> float:
        """`1 / V_drift`: the abscissa of eqn (3), and of the published panel B."""
        return 1.0 / self.v_drift_v


@dataclass(frozen=True)
class Regression:
    """A straight-line fit of arrival time against reciprocal drift voltage.

    `slope_ms_v` is in ms V, the units of `t_D` against `1/V`, and `intercept_ms` is
    `t_0`. The errors are the standard errors of the fit, on Excel's `LINEST`
    convention -- the residual variance carries `n - 2` degrees of freedom -- so a
    series regressed here and the same series regressed in a spreadsheet agree to
    double precision, which is what makes the 2013 workbooks checkable at all.

    The fit error is the *only* uncertainty the legacy analysis had.
    `cross_section_from_regression` propagates P, T and L alongside it when the
    conditions table and the profile supply them, and `CrossSection.propagated` names
    whichever of the four actually went in.

    `excluded` records the acquisitions that were dropped before the fit, so the
    subset a number came from travels with the number. `rms_residual_ms` is the root
    mean square of the residuals in milliseconds: a series whose points lie on eqn (3)
    leaves residuals at the level of the centroid uncertainty, and one that does not
    is the first thing to look at.
    """

    slope_ms_v: float
    slope_ms_v_err: float
    intercept_ms: float
    intercept_ms_err: float
    r_squared: float
    n: int
    excluded: tuple[str, ...] = ()
    rms_residual_ms: float = 0.0
    weighted: bool = False


@dataclass(frozen=True)
class CrossSection:
    """A collision cross section, with everything needed to interpret it.

    Every field is required. A cross section in helium at 2.03 Torr and 301 K is not
    the same measurement as one in nitrogen at 1.5 Torr, and a number quoted without
    its gas, its charge state and its P and T cannot be compared with anything -- which
    is exactly how the 2013 compilation ended up with two columns from one set of
    acquisitions carrying two different pressures.

    `omega_a2` is in square angstroms, `mobility_cm2_v_s` is K at the stated P and T,
    and `reduced_mobility_cm2_v_s` is K0 at 273.15 K and 760 Torr.
    """

    omega_a2: float
    omega_a2_err: float | None
    mobility_cm2_v_s: float
    reduced_mobility_cm2_v_s: float
    charge: int
    ion_mass_da: float
    gas: str
    pressure_torr: float
    temperature_k: float
    propagated: tuple[str, ...] = ()
    """Which quantities the error bar includes: `("slope",)`, or more when the
    conditions table supplies uncertainties. Empty means there is no error bar."""

    def __str__(self) -> str:
        error = "" if self.omega_a2_err is None else f" +/- {self.omega_a2_err:.1f}"
        return (
            f"{self.omega_a2:.1f}{error} A^2, z={self.charge}, {self.gas}, "
            f"{self.pressure_torr:.4g} Torr, {self.temperature_k:.5g} K"
        )


def mobility_from_slope(slope_ms_v: float, drift_length_cm: float) -> float:
    """Eqn (3) rearranged: `K = L^2 / slope`, in cm^2 V^-1 s^-1.

    `slope_ms_v` is the slope of arrival time in milliseconds against reciprocal drift
    voltage, so it is in ms V and is divided by 1000 here. `drift_length_cm` is the
    *electrical drift region*, from the instrument profile -- 25.05 cm on the UW G2,
    not the 25.2 cm mechanical length.

    The K this returns is at the experiment's own pressure and temperature. Reduce it
    with `constants.reduce_mobility` to compare it with anything.
    """
    if not (slope_ms_v > 0.0):
        raise ValueError(
            f"the slope must be positive; got {slope_ms_v} ms V. Arrival time rises with "
            "1/V, so a non-positive slope means the axes are swapped or the fit failed."
        )
    if not (drift_length_cm > 0.0):
        raise ValueError(f"the drift length must be positive; got {drift_length_cm} cm")
    return drift_length_cm**2 / (slope_ms_v / 1000.0)


def cross_section(
    mobility_cm2_v_s: float,
    *,
    charge: int,
    ion_mass_da: float,
    gas: str,
    pressure_torr: float,
    temperature_k: float,
    relative_err: float | None = None,
    propagated: tuple[str, ...] = (),
) -> CrossSection:
    """Assemble a `CrossSection` from a mobility and the conditions it was measured at.

    `relative_err` is the fractional uncertainty on the mobility, which passes straight
    through to the cross section because eqn (4) has K in the denominator and nothing
    else uncertain in this call. `propagated` names what that fraction covers, and is
    what makes an error bar readable six months later.
    """
    gas_entry = constants.gas(gas)
    omega = constants.ccs_from_mobility(
        mobility_cm2_v_s,
        charge=charge,
        ion_mass_da=ion_mass_da,
        gas_mass_da=gas_entry.mass_da,
        pressure_torr=pressure_torr,
        temperature_k=temperature_k,
    )
    return CrossSection(
        omega_a2=omega,
        omega_a2_err=None if relative_err is None else abs(relative_err) * omega,
        mobility_cm2_v_s=mobility_cm2_v_s,
        reduced_mobility_cm2_v_s=constants.reduce_mobility(
            mobility_cm2_v_s, pressure_torr, temperature_k
        ),
        charge=abs(int(charge)),
        ion_mass_da=ion_mass_da,
        gas=gas_entry.name,
        pressure_torr=pressure_torr,
        temperature_k=temperature_k,
        propagated=propagated,
    )


def regress(
    points: Sequence[DriftPoint],
    profile: InstrumentProfile,
    *,
    exclude: Sequence[str] = (),
    weighted: bool = False,
) -> Regression:
    """Least-squares fit of eqn (3) over a drift-voltage series.

    Arrival time in milliseconds against reciprocal drift voltage, giving a slope of
    `L^2 / K` in ms V and an intercept of `t_0` in ms. The standard errors follow
    Excel's `LINEST`: the residual sum of squares over `n - 2` degrees of freedom,
    divided by the spread of the abscissa. That convention is deliberate, and is what
    lets a 2013 spreadsheet and this function be compared digit for digit rather than
    approximately.

    `exclude` names acquisitions to leave out, and they are recorded on the result. A
    drift-voltage series is routinely fitted over a subset -- the 2013 polyalanine
    analysis of record kept ten of fourteen acquisitions, dropping the four lowest
    voltages -- and the thing that analysis did not do was write down which. Prefer
    `use = false` in the conditions table, which keeps the reason with the exclusion;
    this argument is for the case where the subset is being varied deliberately.

    `weighted` fits with `1 / sigma^2` weights from each point's `drift_time_ms_err`,
    and needs every point to carry one. It is off by default, and the default is not
    laziness: those sigmas are the formal errors of a single-Gaussian fit to an
    arrival-time distribution that is not quite Gaussian, so they measure how tightly
    the wrong model was determined rather than how well the centroid is known. Across
    the ten acquisitions of one 2013 regression they span a factor of 3.2 at the median,
    and weighting by them moves a cross section by 0.3 % on average and 1.7 % at worst
    with no improvement in R^2 to show for it (lab record, task 06).
    Weighting is here for series whose centroid errors are trustworthy; unweighted is
    what the workbooks did, and what reproduces them.

    `profile` is required because a drift-voltage series is only meaningful against the
    instrument it was measured on: the voltages in `points` came from this profile's
    formula, and the mobility that comes out of the slope will use this profile's drift
    length. The line fit itself needs nothing more from it, so what it does here is
    check the series -- enough points to have a residual, and no two acquisitions
    sitting at the same drift voltage, which is the signature of a conditions table
    that was copied rather than filled in.
    """
    dropped = tuple(exclude)
    kept = [point for point in points if point.acquisition not in set(dropped)]

    if len(kept) < 3:
        raise ValueError(
            "a drift-time regression needs at least 3 acquisitions to have any "
            f"residual to speak of; got {len(kept)}"
            + (f" after excluding {', '.join(dropped)}" if dropped else "")
        )

    seen: dict[float, str] = {}
    for point in kept:
        clash = seen.get(point.v_drift_v)
        if clash is not None:
            raise ValueError(
                f"{point.acquisition} and {clash} are both at {point.v_drift_v} V on "
                f"{profile.name}; a drift-voltage series must step the drift voltage"
            )
        seen[point.v_drift_v] = point.acquisition

    xs = [point.inverse_v for point in kept]
    ys = [point.drift_time_ms for point in kept]

    if weighted:
        missing = [point.acquisition for point in kept if not point.drift_time_ms_err]
        if missing:
            raise ValueError(
                "a weighted fit needs a drift_time_ms_err on every point; "
                f"{len(missing)} lack one, starting with {missing[0]}"
            )
        weights = [1.0 / float(point.drift_time_ms_err) ** 2 for point in kept]  # type: ignore[arg-type]
    else:
        weights = [1.0] * len(kept)

    n = len(kept)
    w_sum = math.fsum(weights)
    x_bar = math.fsum(w * x for w, x in zip(weights, xs)) / w_sum
    y_bar = math.fsum(w * y for w, y in zip(weights, ys)) / w_sum
    s_xx = math.fsum(w * (x - x_bar) ** 2 for w, x in zip(weights, xs))
    if s_xx <= 0.0:
        raise ValueError(
            "every acquisition is at the same drift voltage, so there is no line to fit"
        )
    s_xy = math.fsum(w * (x - x_bar) * (y - y_bar) for w, x, y in zip(weights, xs, ys))
    s_yy = math.fsum(w * (y - y_bar) ** 2 for w, y in zip(weights, ys))

    slope = s_xy / s_xx
    intercept = y_bar - slope * x_bar

    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    sse = math.fsum(w * r * r for w, r in zip(weights, residuals))
    # LINEST's convention: the residual variance over n - 2 degrees of freedom. In a
    # weighted fit this scales the errors by the reduced chi-square, so that
    # `slope_ms_v_err` keeps one meaning -- the standard error of the fitted line --
    # in both modes.
    variance = sse / (n - 2)
    slope_err = math.sqrt(variance / s_xx)
    intercept_err = math.sqrt(variance * (1.0 / w_sum + x_bar**2 / s_xx))
    r_squared = 1.0 if s_yy <= 0.0 else max(0.0, 1.0 - sse / s_yy)

    return Regression(
        slope_ms_v=slope,
        slope_ms_v_err=slope_err,
        intercept_ms=intercept,
        intercept_ms_err=intercept_err,
        r_squared=r_squared,
        n=n,
        excluded=dropped,
        rms_residual_ms=math.sqrt(math.fsum(r * r for r in residuals) / n),
        weighted=weighted,
    )


def cross_section_from_regression(
    regression: Regression,
    profile: InstrumentProfile,
    *,
    charge: int,
    ion_mass_da: float,
    gas: str,
    pressure_torr: float,
    temperature_k: float,
    pressure_torr_err: float | None = None,
    temperature_k_err: float | None = None,
) -> CrossSection:
    """The whole chain: a fitted line, plus conditions, to a cross section.

    `K = L^2 / slope` with `L` from `profile`, then eqn (4) at `pressure_torr` and
    `temperature_k`, then the uncertainty. This is the function analysis code calls;
    `mobility_from_slope` and `cross_section` are its two halves, kept separate so
    that either can be used alone.

    Uncertainties are optional one at a time. Whatever is supplied is combined in
    quadrature over the logarithmic derivatives of
    `Omega ~ slope . L^-2 . P^-1 . T^+1/2`, and `CrossSection.propagated` comes back
    naming exactly the terms that went in:

    * ``slope`` -- always, from `regression.slope_ms_v_err`.
    * ``drift_length`` -- when the profile carries a `drift_length_cm_err`, at twice
      its relative size, because eqn (3)'s slope goes as `L^2`.
    * ``pressure``, ``temperature`` -- when the conditions table supplies them.

    On the 2013 series only ``slope`` is available, and that is the honest answer for
    it: the pressure was read off a manometer's front panel and the temperature has no
    recorded provenance whatever, so a bar claiming to include them would be an
    invention. Naming the terms is what lets a series measured properly get a wider
    and truer bar without the format or this function changing.
    """
    mobility = mobility_from_slope(regression.slope_ms_v, profile.drift_length_cm)

    terms = ["slope"]
    variance = (regression.slope_ms_v_err / regression.slope_ms_v) ** 2

    if profile.drift_length_cm_err:
        terms.append("drift_length")
        variance += (2.0 * profile.drift_length_cm_err / profile.drift_length_cm) ** 2
    if pressure_torr_err:
        terms.append("pressure")
        variance += (pressure_torr_err / pressure_torr) ** 2
    if temperature_k_err:
        terms.append("temperature")
        variance += (0.5 * temperature_k_err / temperature_k) ** 2

    return cross_section(
        mobility,
        charge=charge,
        ion_mass_da=ion_mass_da,
        gas=gas,
        pressure_torr=pressure_torr,
        temperature_k=temperature_k,
        relative_err=math.sqrt(variance),
        propagated=tuple(terms),
    )

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

The regression itself is task 06. What is implemented here is the closed-form step
either side of it -- slope to mobility, mobility to cross section -- and the object
that carries a cross section, which exists to enforce one rule: **a cross section is
never quoted without its gas, its charge state, and the pressure and temperature it
was reduced with.** `CrossSection` cannot be constructed without all four.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import constants
from .profiles import InstrumentProfile

__all__ = [
    "CrossSection",
    "DriftPoint",
    "Regression",
    "cross_section",
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
    `t_0`. The errors are the standard errors of the fit -- and are the *only*
    uncertainty the legacy analysis had, propagated onto the cross section as a pure
    relative error with nothing from P, T or L. Where the conditions table supplies
    uncertainties for those, task 06 propagates them too; where it does not, this is
    what a cross section's error bar means, and `CrossSection.propagated` says so.
    """

    slope_ms_v: float
    slope_ms_v_err: float
    intercept_ms: float
    intercept_ms_err: float
    r_squared: float
    n: int
    excluded: tuple[str, ...] = ()


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


def regress(points: list[DriftPoint], profile: InstrumentProfile) -> Regression:
    """Least-squares fit of eqn (3) over a drift-voltage series.

    Not implemented: task 06, which also owns the uncertainty propagation and the
    reproduction of the legacy workbooks. The signature and `Regression` are fixed now.

    `profile` is taken here rather than only at the mobility step because a series is
    only meaningful against the instrument it was measured on, and because task 06
    weighs whether to weight the fit -- which needs the profile's view of what varies
    between acquisitions.
    """
    raise NotImplementedError(
        "the drift-time regression is task 06; mobility_from_slope and cross_section "
        "work today"
    )

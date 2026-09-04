"""Physical constants, drift gases, and the closed-form conversions built from them.

Nothing here is typed as a collected magic number. The Mason-Schamp prefactor, the
standard number density and the Torr/Pascal ratio are all evaluated from
`scipy.constants` at import, because the 2013 workflow's worst single systematic error
was a hand-typed prefactor that had drifted 7.6e-4 low (lab record, task 01).

Equation numbers are those of Allen, Giles, Gilbert & Bush, *Analyst* 2016, **141**,
884 (`docs/2016-allen-rf-confining-drift-cell.pdf`, CC-BY). The paper writes eqn (4)
with N and K at the experimental pressure and temperature rather than with N0 and K0.
Both routes are implemented and they agree to double precision, because N K = N0 K0
at fixed T -- but only at fixed T. Eqn (4)'s temperature is the drift-gas temperature
in either route; reducing a mobility to standard conditions does not reduce the
temperature in the square root along with it.

Standard conditions are 273.15 K and 760 Torr, the mobility convention rather than any
of the later IUPAC ones, and they are named `T_STANDARD` / `P_STANDARD_TORR` so that
no caller has to guess which convention a K0 was reduced with.

No SDK, no data, no I/O: pure arithmetic, and the core of the public self-check.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import scipy.constants as _sc

__all__ = [
    "ATOMIC_MASS",
    "BOLTZMANN",
    "DriftGas",
    "ELEMENTARY_CHARGE",
    "GASES",
    "N_STANDARD",
    "PA_PER_TORR",
    "P_STANDARD_PA",
    "P_STANDARD_TORR",
    "T_STANDARD",
    "ccs_from_mobility",
    "ccs_from_reduced_mobility",
    "gas",
    "mason_schamp_prefactor",
    "mobility_from_ccs",
    "number_density",
    "reduce_mobility",
    "reduced_mass",
    "unreduce_mobility",
]

# --- constants ------------------------------------------------------------------

ELEMENTARY_CHARGE = _sc.e  # C
BOLTZMANN = _sc.k  # J/K
ATOMIC_MASS = _sc.atomic_mass  # kg per Da

T_STANDARD = 273.15  # K
P_STANDARD_TORR = 760.0
P_STANDARD_PA = _sc.atm  # 101325 Pa exactly: 760 Torr is one standard atmosphere
PA_PER_TORR = P_STANDARD_PA / P_STANDARD_TORR

# Loschmidt: the number density of an ideal gas at T_STANDARD and P_STANDARD.
N_STANDARD = _sc.value("Loschmidt constant (273.15 K, 101.325 kPa)")  # m^-3


@dataclass(frozen=True)
class DriftGas:
    """A drift gas, by the only property Mason-Schamp needs of it.

    `mass_da` is the mass of the neutral collision partner, so nitrogen is 28.0134 Da,
    the whole molecule, not 14.0067. Getting that wrong moves the reduced mass and
    therefore every cross section in the experiment.
    """

    name: str
    mass_da: float
    formula: str
    aliases: tuple[str, ...] = ()


GASES: dict[str, DriftGas] = {
    entry.name: entry
    for entry in (
        DriftGas("helium", 4.002602, "He"),
        DriftGas("nitrogen", 28.0134, "N2", aliases=("dinitrogen",)),
    )
}


def gas(name: str) -> DriftGas:
    """Look a drift gas up by name, formula or alias, case-insensitively.

    Raises `KeyError` naming what is supported. Helium and nitrogen are the two the
    project covers (decision of record); a third needs its mass added to `GASES` and
    nothing else.
    """
    wanted = name.strip().lower()
    for entry in GASES.values():
        candidates = (entry.name, entry.formula.lower(), *(a.lower() for a in entry.aliases))
        if wanted in candidates:
            return entry
    known = ", ".join(f"{e.name} ({e.formula})" for e in GASES.values())
    raise KeyError(f"unknown drift gas {name!r}; schamp knows {known}")


# --- conversions ----------------------------------------------------------------


def _require_conditions(pressure_torr: float, temperature_k: float) -> None:
    if not (pressure_torr > 0.0) or not (temperature_k > 0.0):
        raise ValueError(
            "pressure and temperature must both be positive; got "
            f"{pressure_torr!r} Torr, {temperature_k!r} K"
        )


def number_density(pressure_torr: float, temperature_k: float) -> float:
    """Ideal-gas number density N, in m^-3, from pressure in Torr and temperature in K.

    N = P / (k_B T). At the couple of Torr this cell runs at, in helium or nitrogen,
    the ideal-gas law is far inside the experiment's own uncertainty; what actually
    limits N is how P and T were measured.
    """
    _require_conditions(pressure_torr, temperature_k)
    return pressure_torr * PA_PER_TORR / (BOLTZMANN * temperature_k)


def reduced_mass(ion_mass_da: float, gas_mass_da: float) -> float:
    """Ion-neutral reduced mass mu = m M / (m + M). Both arguments and the result in Da."""
    if not (ion_mass_da > 0.0) or not (gas_mass_da > 0.0):
        raise ValueError(
            f"masses must be positive; got {ion_mass_da!r} Da ion, {gas_mass_da!r} Da gas"
        )
    return ion_mass_da * gas_mass_da / (ion_mass_da + gas_mass_da)


def reduce_mobility(
    mobility_cm2_v_s: float, pressure_torr: float, temperature_k: float
) -> float:
    """K -> K0: reduce a mobility to 273.15 K and 760 Torr.

    K0 = K (T_STANDARD / T) (P / P_STANDARD_TORR). K is what a drift-time regression
    measures, at the experiment's own P and T; K0 is what gets reported and compared
    between instruments.
    """
    _require_conditions(pressure_torr, temperature_k)
    return mobility_cm2_v_s * (T_STANDARD / temperature_k) * (pressure_torr / P_STANDARD_TORR)


def unreduce_mobility(
    reduced_mobility_cm2_v_s: float, pressure_torr: float, temperature_k: float
) -> float:
    """K0 -> K at the given P and T. The inverse of `reduce_mobility`."""
    _require_conditions(pressure_torr, temperature_k)
    return (
        reduced_mobility_cm2_v_s
        * (temperature_k / T_STANDARD)
        * (P_STANDARD_TORR / pressure_torr)
    )


def mason_schamp_prefactor() -> float:
    """The collected coefficient C of the K0 form of eqn (4), evaluated from constants.

    With Omega in A^2, K0 in cm^2 V^-1 s^-1, mu in Da and the drift-gas temperature T
    in K, eqn (4) rewritten over K0 collapses to

        Omega = (z / K0) mu^-1/2 T^-1/2 C,      C = (3 e / 16 N0) (2 pi / k_B)^1/2

    in those units. This exists so that C can be compared against a value someone once
    typed into a spreadsheet, which is how the legacy analysis was found to run low.
    Analysis code calls `ccs_from_mobility`, which needs no C at all.
    """
    si = (3.0 * ELEMENTARY_CHARGE / (16.0 * N_STANDARD)) * math.sqrt(
        2.0 * math.pi / (ATOMIC_MASS * BOLTZMANN)
    )
    # K0 in cm^2 V^-1 s^-1 -> m^2 V^-1 s^-1 is 1e-4; Omega in m^2 -> A^2 is 1e20.
    return si * 1e20 / 1e-4


def ccs_from_mobility(
    mobility_cm2_v_s: float,
    *,
    charge: int,
    ion_mass_da: float,
    gas_mass_da: float,
    pressure_torr: float,
    temperature_k: float,
) -> float:
    """Eqn (4): the Mason-Schamp collision cross section, in A^2.

        Omega = (3 z e / 16 N) (2 pi / (mu k_B T))^1/2 (1 / K)

    N and K are both at the experiment's own pressure and temperature, as the paper
    prints it, and `temperature_k` is that same temperature in both places.

    The K0 route is `ccs_from_reduced_mobility`, and the two agree exactly because
    N K = N0 K0 *at fixed T*. Note what that does not say: T appears explicitly in
    eqn (4) and is always the drift-gas temperature, never 273.15 K. Reducing the
    mobility to standard conditions does not reduce the temperature in the square
    root with it, and treating it as though it did is a ~5 % error on this cell.

    `charge` is the charge number z and is taken as its absolute value, so a doubly
    deprotonated ion is 2 rather than -2. `ion_mass_da` is the mass of the ion,
    `gas_mass_da` that of the neutral collision partner (`gas(name).mass_da`), and
    `temperature_k` is the drift-gas temperature.

    A cross section means nothing without its gas, its charge state and the P and T it
    was reduced with. This function is the arithmetic alone; `mobility.CrossSection`
    is what carries the number together with all of those.
    """
    z = abs(int(charge))
    if z < 1:
        raise ValueError(f"charge must be a nonzero integer; got {charge!r}")
    if not (mobility_cm2_v_s > 0.0):
        raise ValueError(f"mobility must be positive; got {mobility_cm2_v_s!r}")
    n_m3 = number_density(pressure_torr, temperature_k)
    mu_kg = reduced_mass(ion_mass_da, gas_mass_da) * ATOMIC_MASS
    k_si = mobility_cm2_v_s * 1e-4  # cm^2 V^-1 s^-1 -> m^2 V^-1 s^-1
    omega_m2 = (
        (3.0 * z * ELEMENTARY_CHARGE / (16.0 * n_m3))
        * math.sqrt(2.0 * math.pi / (mu_kg * BOLTZMANN * temperature_k))
        / k_si
    )
    return omega_m2 * 1e20


def ccs_from_reduced_mobility(
    reduced_mobility_cm2_v_s: float,
    *,
    charge: int,
    ion_mass_da: float,
    gas_mass_da: float,
    temperature_k: float,
) -> float:
    """Eqn (4) by way of K0: the same cross section, in A^2, from a reduced mobility.

        Omega = (3 z e / 16 N0) (2 pi / (mu k_B T))^1/2 (1 / K0)

    which in the units this module uses is `(z / K0) mu^-1/2 T^-1/2 C` with C from
    `mason_schamp_prefactor`. This is the route the legacy workbooks took.

    There is deliberately no pressure argument. K0 has already absorbed the pressure,
    but `temperature_k` has not been absorbed by anything: it is still the drift-gas
    temperature, exactly as in `ccs_from_mobility`. The two functions agree to full
    double precision for any pressure, which is what the public self-check asserts.
    """
    z = abs(int(charge))
    if z < 1:
        raise ValueError(f"charge must be a nonzero integer; got {charge!r}")
    if not (reduced_mobility_cm2_v_s > 0.0):
        raise ValueError(f"mobility must be positive; got {reduced_mobility_cm2_v_s!r}")
    if not (temperature_k > 0.0):
        raise ValueError(f"temperature must be positive; got {temperature_k!r}")
    mu_da = reduced_mass(ion_mass_da, gas_mass_da)
    return (
        z
        / reduced_mobility_cm2_v_s
        * mu_da**-0.5
        * temperature_k**-0.5
        * mason_schamp_prefactor()
    )


def mobility_from_ccs(
    ccs_a2: float,
    *,
    charge: int,
    ion_mass_da: float,
    gas_mass_da: float,
    pressure_torr: float,
    temperature_k: float,
) -> float:
    """Eqn (4) solved for K, in cm^2 V^-1 s^-1. The exact inverse of `ccs_from_mobility`.

    Two uses: proving the closed form round trips, which the public self-check does,
    and placing a literature cross section onto a measured drift-time plot.
    """
    z = abs(int(charge))
    if z < 1:
        raise ValueError(f"charge must be a nonzero integer; got {charge!r}")
    if not (ccs_a2 > 0.0):
        raise ValueError(f"cross section must be positive; got {ccs_a2!r}")
    n_m3 = number_density(pressure_torr, temperature_k)
    mu_kg = reduced_mass(ion_mass_da, gas_mass_da) * ATOMIC_MASS
    k_si = (
        (3.0 * z * ELEMENTARY_CHARGE / (16.0 * n_m3))
        * math.sqrt(2.0 * math.pi / (mu_kg * BOLTZMANN * temperature_k))
        / (ccs_a2 * 1e-20)
    )
    return k_si * 1e4

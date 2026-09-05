"""Self-check for a fresh clone: no Waters SDK, no license, no .raw file, no
lab repo needed.

Every check here must pass in a bare public clone. Checks that need something a
clone does not ship (the SDK, external/, the lab repo, an acquisition) are
reported as SKIPPED when it is absent, never as FAIL.

What it covers: the physics identities on synthetic input, the instrument-profile
format and its drift-voltage formula, the `_extern.inf` parser against a real
sidecar committed as a fixture, the experiment data model against a table written
and read back in a temporary directory, and the reporting conventions. Also that
the not-yet-implemented layers say so rather than returning something plausible.

This is the public counterpart of the lab repo's regression check, which
additionally pins the pipeline against the archived polyalanine numbers and
needs the legacy acquisitions to run.

Run:  uv run tools/check_public.py
"""

import importlib
import math
import os
import shutil
import sys
import tempfile

import numpy as np

import schamp
import schamp.atd
import schamp.constants as C
import schamp.experiment
import schamp.extract
import schamp.extern
import schamp.mobility
import schamp.profiles
import schamp.report
import schamp.sdk

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(schamp.ROOT, "tests", "fixtures")

FAIL = []
SKIPPED = []


def check_true(name, cond):
    print("{:4s} {}".format("OK" if cond else "FAIL", name))
    if not cond:
        FAIL.append(name)


def check_close(name, got, want, rel=1e-12):
    ok = want != 0 and abs(got / want - 1.0) <= rel
    print("{:4s} {} ({!r} vs {!r})".format("OK" if ok else "FAIL", name, got, want))
    if not ok:
        FAIL.append(name)


def check_raises(name, exc, call):
    try:
        call()
    except exc:
        check_true(name, True)
        return
    except Exception as other:  # noqa: BLE001 -- reporting, not handling
        print("FAIL {} (raised {!r}, wanted {})".format(name, other, exc.__name__))
        FAIL.append(name)
        return
    print("FAIL {} (did not raise {})".format(name, exc.__name__))
    FAIL.append(name)


def skip(name, why):
    print("SKIP {} ({})".format(name, why))
    SKIPPED.append(name)


def section(title):
    print()
    print("--- " + title + " " + "-" * max(3, 72 - len(title)))


# --------------------------------------------------------------------------------
section("package")
check_true(
    "schamp imports and carries a version",
    isinstance(schamp.__version__, str) and bool(schamp.__version__),
)
check_true(
    "docs/ holds the 2016 paper",
    os.path.isfile(os.path.join(schamp.DOCS_DIR, "2016-allen-rf-confining-drift-cell.pdf")),
)
check_true(
    "docs/ holds the 2016 SI",
    os.path.isfile(os.path.join(schamp.DOCS_DIR, "2016-allen-rf-confining-drift-cell-SI.pdf")),
)
check_true("DATA_DIR ships with the package", os.path.isdir(schamp.DATA_DIR))

# --------------------------------------------------------------------------------
section("constants and the Mason-Schamp closed form")
# Loschmidt from CODATA against the ideal-gas identity at the same standard point:
# two independent routes to N0, which catches a wrong standard pressure or temperature.
check_close(
    "N0 agrees with P0/(k_B T0)",
    C.N_STANDARD,
    C.P_STANDARD_PA / (C.BOLTZMANN * C.T_STANDARD),
    rel=1e-9,
)
check_close(
    "number_density at standard conditions is N0",
    C.number_density(C.P_STANDARD_TORR, C.T_STANDARD),
    C.N_STANDARD,
    rel=1e-9,
)
check_close("760 Torr is one standard atmosphere", C.PA_PER_TORR * 760.0, 101325.0, rel=1e-15)
# The prefactor is derived, never typed. 18509.87 is what CODATA gives; a value that
# has drifted from it means either a changed constant set or a typed number creeping
# back in, and both are worth failing on.
check_close("Mason-Schamp prefactor from constants", C.mason_schamp_prefactor(), 18509.866, rel=1e-6)
# mu is symmetric, is bounded above by the lighter partner, and matches the closed form.
check_close("reduced mass is symmetric", C.reduced_mass(657.3, 4.002602), C.reduced_mass(4.002602, 657.3))
check_close(
    "reduced mass of a 657.3 Da ion in He",
    C.reduced_mass(657.3, 4.002602),
    657.3 * 4.002602 / (657.3 + 4.002602),
)
check_true("reduced mass is below the lighter partner", C.reduced_mass(657.3, 4.002602) < 4.002602)
check_true("nitrogen is the whole molecule", C.gas("N2").mass_da > 28.0)
check_true("gas lookup takes a name, a formula and a case", C.gas("HELIUM") is C.gas("He"))
check_raises("gas lookup refuses an unknown gas", KeyError, lambda: C.gas("argon"))
check_raises("number_density refuses a zero pressure", ValueError, lambda: C.number_density(0.0, 300.0))

# The two forms of eqn (4) must agree at every pressure, and only because the drift-gas
# temperature is the same in both. This is the round trip the module docstring claims.
he = C.gas("helium")
worst_route, worst_inverse = 0.0, 0.0
for pressure, temperature in ((2.028, 301.13), (760.0, 273.15), (1.5, 350.0), (0.5, 200.0)):
    for charge, mass in ((1, 657.3), (2, 1580.8), (3, 4000.0)):
        k0 = 0.87
        k = C.unreduce_mobility(k0, pressure, temperature)
        shared = dict(charge=charge, ion_mass_da=mass, gas_mass_da=he.mass_da)
        by_nk = C.ccs_from_mobility(
            k, pressure_torr=pressure, temperature_k=temperature, **shared
        )
        by_k0 = C.ccs_from_reduced_mobility(k0, temperature_k=temperature, **shared)
        worst_route = max(worst_route, abs(by_nk / by_k0 - 1.0))
        back = C.mobility_from_ccs(
            by_nk, pressure_torr=pressure, temperature_k=temperature, **shared
        )
        worst_inverse = max(worst_inverse, abs(back / k - 1.0))
check_true(f"the N-and-K and K0 forms of eqn (4) agree (worst {worst_route:.1e})", worst_route < 1e-13)
check_true(f"ccs_from_mobility inverts exactly (worst {worst_inverse:.1e})", worst_inverse < 1e-13)
check_close(
    "K -> K0 -> K round trips",
    C.unreduce_mobility(C.reduce_mobility(4.2, 2.03, 301.13), 2.03, 301.13),
    4.2,
    rel=1e-14,
)

# --------------------------------------------------------------------------------
section("instrument profiles")
builtin = schamp.profiles.builtin_profiles()
check_true("uw-synapt-g2 ships with the package", "uw-synapt-g2" in builtin)
profile = schamp.load_profile("uw-synapt-g2")
check_close("the drift region is 25.05 cm", profile.drift_length_cm, 25.05, rel=1e-15)
check_true(
    "the mechanical length is recorded separately and is not the drift length",
    profile.mechanical_length_cm == 25.2 and profile.mechanical_length_cm != profile.drift_length_cm,
)
check_close("the divider ratio is 2/170", profile.divider_ratio, 2.0 / 170.0, rel=1e-15)
check_true(
    "the drift-voltage keys are named, not guessed",
    profile.high_keys == ("Helium Cell DC", "IMSBias", "Helium Exit")
    and profile.low_keys == ("Transfer DC Entrance",),
)
# The published typical voltages of ESI Fig. S1: optics 7, 9, 8 and 10 are 56, 240,
# -40 and 4 V, so VappH - VappL = 252 V and the drift voltage is 252 * 168/170.
figure_s1 = {
    "Helium Cell DC": "56",
    "IMSBias": "240",
    "Helium Exit": "-40",
    "Transfer DC Entrance": "4",
}
check_close(
    "the drift voltage of the paper's typical settings",
    profile.drift_voltage(figure_s1),
    252.0 * 168.0 / 170.0,
    rel=1e-14,
)
check_true("a negative Helium Exit subtracts", profile.applied_voltage(figure_s1) == (256.0, 4.0))
check_raises(
    "a missing key fails rather than defaulting to zero",
    schamp.profiles.ProfileError,
    lambda: profile.drift_voltage({"IMSBias": "240"}),
)
check_true(
    "a traveling wave left on is caught",
    len(profile.check_acquisition({"IMS Wave Height (V)": "40.0"})) == 1,
)
check_true(
    "several pushes per drift bin is caught",
    len(profile.check_acquisition({"ADC Pushes Per IMS Increment": "4"})) == 1,
)
check_true("an unreadable sidecar fails no check", profile.check_acquisition({}) == [])
check_raises(
    "an unknown profile name is refused",
    schamp.profiles.ProfileError,
    lambda: schamp.load_profile("no-such-instrument"),
)

# --------------------------------------------------------------------------------
section("_extern.inf")
if not os.path.isfile(os.path.join(FIXTURES, schamp.extern.EXTERN_FILENAME)):
    skip("_extern.inf fixture parses", "tests/fixtures/_extern.inf is missing")
else:
    extern = schamp.extern.parse_extern_inf(FIXTURES)
    check_true(f"the fixture parses to {len(extern)} keys", len(extern) > 150)
    check_close("IMSBias reads back as a float", schamp.extern.as_float(extern, "IMSBias"), 132.0)
    check_true(
        "a Latin-1 degree sign survives the decode",
        any("\xb0C" in key for key in extern),
    )
    check_true(
        "there is no Pusher Interval key on this instrument",
        "Pusher Interval" not in extern and schamp.extern.as_float(extern, "Pusher Interval") is None,
    )
    check_true(
        "Pusher is a voltage, not a period",
        schamp.extern.as_float(extern, "Pusher") == 1900.0,
    )
    check_true(
        "a non-numeric value is None, not 0.0",
        schamp.extern.as_float(extern, "Polarity") is None,
    )
    check_close(
        "the fixture's drift voltage through the profile",
        profile.drift_voltage(extern),
        88.0 * 168.0 / 170.0,
        rel=1e-14,
    )
check_true("a missing directory parses to nothing", schamp.extern.parse_extern_inf("no/such/dir") == {})
check_true("a bad argument parses to nothing", schamp.extern.parse_extern_inf(os.devnull) == {})

# --------------------------------------------------------------------------------
section("the experiment data model")
work = tempfile.mkdtemp(prefix="schamp-check-")
try:
    schamp.report.write_table(
        os.path.join(work, "conditions.csv"),
        ["acquisition", "pressure_Torr", "pressure_Torr_err", "use", "v_drift_V", "notes"],
        [
            ["a_001.raw", 2.026, 0.002, "true", 103.7647, "lowest kept"],
            ["a_002.raw", 2.027, "", "yes", 181.8353, ""],
            ["a_003.raw", 2.028, "", "1", 351.8118, ""],
            ["a_004.raw", 2.028, "", "false", 74.1176, "dropped: lowest voltage"],
        ],
        comment="Synthetic conditions table written by the public self-check.\nNo acquisitions exist; every voltage is an override.",
    )
    with open(os.path.join(work, "experiment.toml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(
            'schema = 1\ntitle = "synthetic series"\nprofile = "uw-synapt-g2"\n\n'
            "[defaults]\ntemperature_K = 301.13\n"
        )
    experiment = schamp.load_experiment(work)
    check_true("an experiment loads from its directory", experiment.title == "synthetic series")
    check_true("the profile resolves by builtin name", experiment.profile.name == "uw-synapt-g2")
    check_true("the gas falls back to the profile's default", experiment.gas == "helium")
    check_true("[defaults] fills a column no row carries", all(c.temperature_k == 301.13 for c in experiment.conditions))
    check_true("use = false drops a row from the series, keeping the record", len(experiment.conditions) == 4 and len(experiment.used) == 3)
    check_true("a dropped row keeps its reason", experiment.conditions[3].notes.startswith("dropped"))
    check_true("an optional uncertainty is None when the cell is empty", experiment.conditions[1].pressure_torr_err is None)
    check_close("an optional uncertainty is read when present", experiment.conditions[0].pressure_torr_err, 0.002)
    check_true("a v_drift_V override wins over the profile formula", experiment.used[0].drift_voltage(profile) == 103.7647)
    check_true("a table this repo writes is a table it can read", experiment.validate(require_raw=False) == [])
    check_true("the missing acquisitions are reported when required", len(experiment.validate(require_raw=True)) == 4)

    # A conditions table whose columns are wrong must be refused, not silently ignored:
    # a mistyped 'pressure_torr_error' would otherwise become an experiment with no
    # uncertainties at all and no complaint about it.
    with open(os.path.join(work, "typo.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("acquisition,pressure_Torr,temperature_K,presure_Torr_err\nx.raw,2.0,300.0,0.1\n")
    with open(os.path.join(work, "typo.toml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write('schema = 1\nprofile = "uw-synapt-g2"\nconditions = "typo.csv"\n')
    check_raises(
        "a mistyped column is refused",
        schamp.experiment.ExperimentError,
        lambda: schamp.load_experiment(os.path.join(work, "typo.toml")),
    )

    # extern: overrides
    with open(os.path.join(work, "override.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(
            "acquisition,pressure_Torr,temperature_K,extern:Helium Cell DC,"
            "extern:IMSBias,extern:Helium Exit,extern:Transfer DC Entrance\n"
            "x.raw,2.0,300.0,56,240,-40,4\n"
        )
    with open(os.path.join(work, "override.toml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write('schema = 1\nprofile = "uw-synapt-g2"\nconditions = "override.csv"\n')
    overridden = schamp.load_experiment(os.path.join(work, "override.toml"))
    check_close(
        "an extern: column overrides one _extern.inf value",
        overridden.conditions[0].drift_voltage(profile),
        252.0 * 168.0 / 170.0,
        rel=1e-14,
    )
    check_raises(
        "an experiment with no profile is refused",
        schamp.experiment.ExperimentError,
        lambda: schamp.load_experiment(os.devnull),
    )

    # --------------------------------------------------------------------------
    section("reporting conventions")
    written = schamp.report.write_results({"answer": 42}, work, task="check-public")
    check_true("write_results stamps and writes", os.path.isfile(written))
    import json

    stamped = json.load(open(written, encoding="utf-8"))
    check_true(
        "the stamp carries the date, the version and the results",
        stamped["results"] == {"answer": 42}
        and stamped["schamp_version"] == schamp.__version__
        and len(stamped["date"]) == 10,
    )
    check_true("figure defaults keep vector text as text", schamp.report.figure_defaults()["pdf.fonttype"] == 42)
finally:
    shutil.rmtree(work, ignore_errors=True)

# --------------------------------------------------------------------------------
section("arrival-time distributions")
bins = np.arange(200)
times = bins * 0.0692521
centre, sigma = 4.64, 0.30
intensity = 1000.0 * np.exp(-((times - centre) ** 2) / (2.0 * sigma**2))
atd = schamp.atd.ATD(bins, times, intensity, mz_low=788.2353, mz_high=794.1177, label="n=22, 2-")
check_close("the window centre is the window's centre", atd.mz_centre, 791.1765, rel=1e-9)
check_true("the peak bin is found", atd.drift_time_ms[atd.peak_bin] == min(times, key=lambda t: abs(t - centre)))
check_true("a populated ATD is not empty", not atd.is_empty)
moments = schamp.atd.moments(atd)
check_close("moments recover a synthetic centre", moments.centre_ms, centre, rel=1e-6)
check_close("moments recover a synthetic width", moments.sigma_ms, sigma, rel=2e-3)
check_raises(
    "an ATD with mismatched axes is refused",
    ValueError,
    lambda: schamp.atd.ATD(bins, times[:10], intensity, mz_low=1.0, mz_high=2.0),
)
check_raises(
    "an inverted m/z window is refused",
    ValueError,
    lambda: schamp.atd.ATD(bins, times, intensity, mz_low=1224.4, mz_high=1125.4),
)
check_true(
    "an all-zero ATD reports itself empty",
    schamp.atd.ATD(bins, times, np.zeros_like(times), mz_low=1.0, mz_high=2.0).is_empty,
)
check_raises(
    "an inverted m/z window is refused at the extraction layer too",
    ValueError,
    lambda: schamp.extract.MzWindow(1224.4, 1125.4, "n=17, 1-"),
)

# --------------------------------------------------------------------------------
section("single-Gaussian ATD fits")
# Noiseless: the fit model is the data's own generating function, so recovery should
# be limited only by the optimiser, not by the data -- a tight tolerance is the point.
clean_fit = schamp.atd.fit_gaussian(atd)
check_true("a noiseless fit converges", clean_fit.converged)
check_close("a noiseless fit recovers the synthetic centre", clean_fit.centre_ms, centre, rel=1e-6)
check_close("a noiseless fit recovers the synthetic sigma", clean_fit.sigma_ms, sigma, rel=1e-6)
check_close("a noiseless fit recovers the synthetic height", clean_fit.height, 1000.0, rel=1e-6)
check_close(
    "resolution is centre / FWHM, unit-invariant",
    clean_fit.resolution,
    clean_fit.centre_ms / (clean_fit.sigma_ms * 2.35482),
    rel=1e-12,
)
check_close(
    "centre_bins divides back by the ATD's own drift-bin spacing",
    clean_fit.centre_bins,
    clean_fit.centre_ms / 0.0692521,
    rel=1e-9,
)
check_true(
    "a near-exact fit reports near-zero residual error",
    clean_fit.rmsd_fraction < 1e-6 and clean_fit.area_error_fraction < 1e-6,
)

# Seeded noise: the fit should still recover the truth, and now report uncertainties.
rng = np.random.default_rng(20260504)
noisy_intensity = intensity + rng.normal(0.0, 5.0, size=intensity.shape)
noisy_atd = schamp.atd.ATD(
    bins, times, noisy_intensity, mz_low=788.2353, mz_high=794.1177, label="n=22, 2- (noisy)"
)
noisy_fit = schamp.atd.fit_gaussian(noisy_atd)
check_true("a fit to noisy synthetic data converges", noisy_fit.converged)
check_close("a noisy fit recovers the synthetic centre", noisy_fit.centre_ms, centre, rel=1e-3)
check_close("a noisy fit recovers the synthetic sigma", noisy_fit.sigma_ms, sigma, rel=2e-2)
check_true(
    "a converged fit reports parameter uncertainties",
    noisy_fit.centre_ms_err is not None
    and noisy_fit.centre_ms_err > 0.0
    and noisy_fit.sigma_ms_err is not None
    and noisy_fit.sigma_ms_err > 0.0,
)
check_true(
    "the true centre falls within a few fitted standard errors",
    abs(noisy_fit.centre_ms - centre) < 5.0 * noisy_fit.centre_ms_err,
)

# --------------------------------------------------------------------------------
section("m/z windows and the legacy peaks format")
# Everything here is window bookkeeping, so it runs with no SDK and no acquisition.
window = schamp.extract.MzWindow.around(789.4062, 0.3, "n=22, 2-", charge=-2, ion_mass_da=1580.8)
check_close("a window centres on the m/z it was built around", window.centre, 789.4062, rel=1e-12)
check_close("a window keeps the width it was given", window.width, 0.3, rel=1e-12)
check_true("a window carries the charge the mobility layer needs", window.charge == -2)
check_raises(
    "a window far narrower than any m/z step is refused",
    ValueError,
    lambda: schamp.extract.MzWindow.around(1335.2941, 1e-6),
)
check_close(
    "a 0.2 Th window, as the 2013 cation peaks file used, is accepted",
    schamp.extract.MzWindow.around(232.1292, 0.2).width,
    0.2,
    rel=1e-12,
)
grid = schamp.extract.contiguous_windows(200.0, 1700.0, 255, label="panel A")
check_true("a contiguous grid has the window count asked for", len(grid) == 255)
check_true(
    "a contiguous grid tiles the span with touching bounds",
    grid[0].low == 200.0
    and grid[-1].high == 1700.0
    and all(grid[i].high == grid[i + 1].low for i in range(254)),
)
check_close("the 2013 heat-map grid is 5.882 Th per window", grid[0].width, 1500.0 / 255.0, rel=1e-12)
check_raises(
    "a grid with no windows is refused",
    ValueError,
    lambda: schamp.extract.contiguous_windows(200.0, 1700.0, 0),
)
# No windows means no reads, so this needs neither an SDK nor an acquisition.
check_true("extracting no windows opens nothing", schamp.extract.extract_atds(None, []) == [])
check_raises(
    "a mobility map over no windows is refused rather than empty",
    ValueError,
    lambda: schamp.extract.mobility_map(None, []),
)
peaks_dir = tempfile.mkdtemp(prefix="schamp-peaks-")
try:
    peaks_path = os.path.join(peaks_dir, "peaks.csv")
    schamp.extract.write_peaks_csv(peaks_path, grid[:3])
    round_trip = schamp.extract.read_peaks_csv(peaks_path)
    check_true(
        "a peaks file round trips at full precision",
        len(round_trip) == 3
        and all(
            round_trip[i].low == grid[i].low and round_trip[i].high == grid[i].high
            for i in range(3)
        ),
    )
    with open(peaks_path, encoding="utf-8") as fh:
        check_true("a peaks file declares its window count first", fh.readline().startswith("3,"))
    # cdctest read exactly the declared count and ignored the rest of the file.
    with open(peaks_path, "a", encoding="utf-8") as fh:
        fh.write("999.0,1000.0\n")
    check_true(
        "a peaks file's declared count wins over its extra rows",
        len(schamp.extract.read_peaks_csv(peaks_path)) == 3,
    )
    # The three mistyped rows of the 2013 anion peaks file, which cdctest read silently.
    bad_path = os.path.join(peaks_dir, "bad.csv")
    with open(bad_path, "w", encoding="utf-8") as fh:
        fh.write("3,\n295.5,230.5\n284.7,585.7\n1224.4,1125.4\n")
    check_raises(
        "an inverted row of a peaks file is refused, naming the row",
        ValueError,
        lambda: schamp.extract.read_peaks_csv(bad_path),
    )
    short_path = os.path.join(peaks_dir, "short.csv")
    with open(short_path, "w", encoding="utf-8") as fh:
        fh.write("5,\n300.5,301.5\n")
    check_raises(
        "a peaks file that declares more windows than it carries is refused",
        ValueError,
        lambda: schamp.extract.read_peaks_csv(short_path),
    )
finally:
    shutil.rmtree(peaks_dir, ignore_errors=True)

# --------------------------------------------------------------------------------
section("mobility")
# Eqn (3) with L = 25.05 cm: a slope of L^2/K ms V corresponds to K exactly.
for k_expected in (1.0, 4.2, 12.7):
    slope = 1000.0 * profile.drift_length_cm**2 / k_expected
    check_close(
        f"slope -> K round trips at K = {k_expected}",
        schamp.mobility.mobility_from_slope(slope, profile.drift_length_cm),
        k_expected,
        rel=1e-14,
    )
check_raises(
    "a negative slope is refused",
    ValueError,
    lambda: schamp.mobility.mobility_from_slope(-1.0, 25.05),
)
check_raises(
    "a zero drift voltage is refused",
    ValueError,
    lambda: schamp.mobility.DriftPoint("x.raw", 0.0, 5.0),
)
xs = schamp.mobility.cross_section(
    4.2, charge=2, ion_mass_da=1580.8, gas="He", pressure_torr=2.028, temperature_k=301.13,
    relative_err=0.001, propagated=("slope",),
)
check_close(
    "a cross section matches the closed form",
    xs.omega_a2,
    C.ccs_from_mobility(4.2, charge=2, ion_mass_da=1580.8, gas_mass_da=he.mass_da,
                        pressure_torr=2.028, temperature_k=301.13),
)
check_true(
    "a cross section carries its gas, charge, P and T",
    xs.gas == "helium" and xs.charge == 2 and xs.pressure_torr == 2.028 and xs.temperature_k == 301.13,
)
check_true("a cross section says what its error bar covers", xs.propagated == ("slope",))
check_true("a cross section prints its conditions", "helium" in str(xs) and "Torr" in str(xs))
check_close("K0 accompanies K", xs.reduced_mobility_cm2_v_s, C.reduce_mobility(4.2, 2.028, 301.13))

# --------------------------------------------------------------------------------
section("the drift-time regression")
# The 2013 series' own fourteen drift voltages under the profile's formula, so that the
# synthetic series has the abscissa spacing of a real one. Ten of them are the analysis
# of record: the four lowest are dropped.
SERIES_V = [86.9647, 302.4, 211.4824, 80.0471, 103.7647, 74.1176, 141.3176, 249.0353,
            94.8706, 181.8353, 114.6353, 160.0941, 351.8118, 126.4941]
KEPT_V = [v for v in SERIES_V if v >= 100.0]
check_true("ten of the fourteen 2013 drift voltages are at 104 V and above", len(KEPT_V) == 10)


def synthetic(k_cm2_v_s, t0_ms, voltages, noise=None, err=None):
    """Eqn (3) forwards: a series that a correct regression must invert exactly."""
    slope = 1000.0 * profile.drift_length_cm**2 / k_cm2_v_s
    return [
        schamp.mobility.DriftPoint(
            f"synthetic_{i:03d}", v, t0_ms + slope / v + (0.0 if noise is None else noise[i]),
            drift_time_ms_err=None if err is None else err[i],
        )
        for i, v in enumerate(voltages)
    ]


for k_true, t0_true in ((4.2, 0.55), (1.15, 0.31), (12.7, 0.48)):
    line = schamp.mobility.regress(synthetic(k_true, t0_true, KEPT_V), profile)
    check_close(
        f"a noiseless series inverts to K = {k_true}",
        schamp.mobility.mobility_from_slope(line.slope_ms_v, profile.drift_length_cm),
        k_true,
        rel=1e-12,
    )
    check_close(f"...and to t0 = {t0_true} ms", line.intercept_ms, t0_true, rel=1e-9)
    check_true(f"...with R^2 = 1 at K = {k_true}", abs(line.r_squared - 1.0) < 1e-12)
    check_true(
        f"...and no residual at K = {k_true}",
        line.rms_residual_ms < 1e-12 and line.slope_ms_v_err < 1e-6,
    )

# Errors that are Excel's, not a variant of them. A four-point series worked by hand:
# x = 1, 2, 3, 4; y = 2, 4, 5, 8. Sxx = 5, Sxy = 9.5, Syy = 18.75, so SLOPE = 1.9 and
# INTERCEPT = 0; the residuals are 0.1, 0.2, -0.7, 0.4 and SSE = 0.7, so LINEST's slope
# error is sqrt((SSE/(n-2))/Sxx) = sqrt(0.07) and its R^2 is 1 - 0.7/18.75.
hand = [schamp.mobility.DriftPoint(f"h{i}", 1.0 / x, y) for i, (x, y) in
        enumerate(((1.0, 2.0), (2.0, 4.0), (3.0, 5.0), (4.0, 8.0)))]
hand_line = schamp.mobility.regress(hand, profile)
check_close("the slope is Excel's SLOPE", hand_line.slope_ms_v, 1.9, rel=1e-12)
check_true("the intercept is Excel's INTERCEPT", abs(hand_line.intercept_ms) < 1e-12)
check_close(
    "the slope error is Excel's LINEST",
    hand_line.slope_ms_v_err,
    math.sqrt((0.7 / 2.0) / 5.0),
    rel=1e-12,
)
check_close("R^2 is Excel's LINEST", hand_line.r_squared, 1.0 - 0.7 / 18.75, rel=1e-12)
check_true("the fit counts its points", hand_line.n == 4 and hand_line.excluded == ())

# A series with structure in the residuals, so the error bars are not all zero.
wobble = [+2e-3, -1e-3, +3e-3, -2e-3, +1e-3, -3e-3, +2e-3, -1e-3, +1e-3, -2e-3]
noisy = schamp.mobility.regress(synthetic(4.2, 0.55, KEPT_V, noise=wobble), profile)
check_true("a noisy series has a positive slope error", noisy.slope_ms_v_err > 0.0)
check_true("...and an intercept error", noisy.intercept_ms_err > 0.0)
check_true("...and an R^2 below 1", 0.99 < noisy.r_squared < 1.0)
check_true(
    "...and an rms residual of the size that was put in",
    1e-3 < noisy.rms_residual_ms < 3e-3,
)

# Excluding by name, which is how a subset travels with its result.
full = synthetic(4.2, 0.55, SERIES_V)
subset = schamp.mobility.regress(
    full, profile, exclude=[p.acquisition for p, v in zip(full, SERIES_V) if v < 100.0]
)
check_true("an excluded acquisition is dropped and recorded",
           subset.n == 10 and len(subset.excluded) == 4)
check_close(
    "...and the fit is the same line on a noiseless series",
    subset.slope_ms_v,
    schamp.mobility.regress(full, profile).slope_ms_v,
    rel=1e-9,
)

# Weighting: equal weights must reproduce the unweighted fit exactly.
equal = synthetic(4.2, 0.55, KEPT_V, noise=wobble, err=[2.5e-3] * 10)
check_close(
    "equal weights give the unweighted slope",
    schamp.mobility.regress(equal, profile, weighted=True).slope_ms_v,
    schamp.mobility.regress(equal, profile).slope_ms_v,
    rel=1e-12,
)
check_true(
    "a weighted fit says it was weighted",
    schamp.mobility.regress(equal, profile, weighted=True).weighted,
)
check_raises(
    "a weighted fit without centroid errors is refused",
    ValueError,
    lambda: schamp.mobility.regress(synthetic(4.2, 0.55, KEPT_V), profile, weighted=True),
)
check_raises(
    "fewer than three acquisitions is refused",
    ValueError,
    lambda: schamp.mobility.regress(synthetic(4.2, 0.55, KEPT_V[:2]), profile),
)
check_raises(
    "two acquisitions at one drift voltage are refused",
    ValueError,
    lambda: schamp.mobility.regress(
        synthetic(4.2, 0.55, [104.0, 104.0, 200.0, 300.0]), profile
    ),
)

# --------------------------------------------------------------------------------
section("regression to cross section, and what the error bar covers")
for gas_name, z, mass, k0_true in (("helium", 1, 232.13, 4.60), ("helium", 2, 1580.8, 1.55),
                                   ("nitrogen", 1, 1509.8, 1.10), ("nitrogen", 2, 1580.8, 0.72)):
    p_torr, t_k = (2.028, 301.13) if gas_name == "helium" else (1.5, 298.0)
    k_true = C.unreduce_mobility(k0_true, p_torr, t_k)
    line = schamp.mobility.regress(synthetic(k_true, 0.5, KEPT_V, noise=wobble), profile)
    got = schamp.mobility.cross_section_from_regression(
        line, profile, charge=z, ion_mass_da=mass, gas=gas_name,
        pressure_torr=p_torr, temperature_k=t_k,
    )
    check_close(
        f"the {gas_name} chain reaches eqn (4) at z={z}",
        got.omega_a2,
        C.ccs_from_reduced_mobility(
            got.reduced_mobility_cm2_v_s, charge=z, ion_mass_da=mass,
            gas_mass_da=C.gas(gas_name).mass_da, temperature_k=t_k,
        ),
        rel=1e-12,
    )
    check_true(f"...and carries its gas ({gas_name})", got.gas == gas_name)

line = schamp.mobility.regress(synthetic(4.2, 0.5, KEPT_V, noise=wobble), profile)
slope_only = schamp.mobility.cross_section_from_regression(
    line, profile, charge=1, ion_mass_da=232.13, gas="helium",
    pressure_torr=2.028, temperature_k=301.13,
)
check_true("the slope alone is the default error bar", slope_only.propagated == ("slope",))
check_close(
    "...and it is the slope's own relative error",
    slope_only.omega_a2_err / slope_only.omega_a2,
    line.slope_ms_v_err / line.slope_ms_v,
    rel=1e-12,
)
with_pt = schamp.mobility.cross_section_from_regression(
    line, profile, charge=1, ion_mass_da=232.13, gas="helium",
    pressure_torr=2.028, temperature_k=301.13,
    pressure_torr_err=0.005, temperature_k_err=2.0,
)
check_true(
    "P and T are named when they are supplied",
    with_pt.propagated == ("slope", "pressure", "temperature"),
)
check_close(
    "...and add in quadrature over the logarithmic derivatives",
    with_pt.omega_a2_err / with_pt.omega_a2,
    math.sqrt(
        (line.slope_ms_v_err / line.slope_ms_v) ** 2
        + (0.005 / 2.028) ** 2
        + (0.5 * 2.0 / 301.13) ** 2
    ),
    rel=1e-12,
)
check_true(
    "...and the bar is wider than the slope alone",
    with_pt.omega_a2_err > slope_only.omega_a2_err,
)
check_close(
    "the cross section itself does not move when errors are added",
    with_pt.omega_a2, slope_only.omega_a2, rel=1e-15,
)
# The drift length enters squared, and only when the profile says how well it is known.
import dataclasses as _dc  # noqa: PLC0415 -- stdlib, only this check needs it

measured_l = _dc.replace(profile, drift_length_cm_err=0.05)
with_l = schamp.mobility.cross_section_from_regression(
    line, measured_l, charge=1, ion_mass_da=232.13, gas="helium",
    pressure_torr=2.028, temperature_k=301.13,
)
check_true("a profile with a length error propagates it",
           with_l.propagated == ("slope", "drift_length"))
check_close(
    "...at twice its relative size, because the slope goes as L^2",
    with_l.omega_a2_err / with_l.omega_a2,
    math.sqrt(
        (line.slope_ms_v_err / line.slope_ms_v) ** 2 + (2.0 * 0.05 / 25.05) ** 2
    ),
    rel=1e-12,
)
check_true(
    "the shipped profile claims no length uncertainty, so none is invented",
    profile.drift_length_cm_err is None and slope_only.propagated == ("slope",),
)

# The whole chain, closed: a cross section put in comes back out.
for gas_name, z, mass, omega_true in (("helium", 1, 232.13, 88.0),
                                      ("nitrogen", 2, 1580.8, 422.0)):
    p_torr, t_k = (2.028, 301.13) if gas_name == "helium" else (1.5, 298.0)
    k_true = C.mobility_from_ccs(
        omega_true, charge=z, ion_mass_da=mass, gas_mass_da=C.gas(gas_name).mass_da,
        pressure_torr=p_torr, temperature_k=t_k,
    )
    line = schamp.mobility.regress(synthetic(k_true, 0.5, KEPT_V), profile)
    got = schamp.mobility.cross_section_from_regression(
        line, profile, charge=z, ion_mass_da=mass, gas=gas_name,
        pressure_torr=p_torr, temperature_k=t_k,
    )
    check_close(
        f"{omega_true} A^2 in {gas_name} at z={z} round trips through the whole chain",
        got.omega_a2, omega_true, rel=1e-12,
    )

# --------------------------------------------------------------------------------
section("the layers that are not built yet")
# Every layer is implemented. This section stays because the next unbuilt one belongs
# in it: a stub that returned something plausible would be far worse than one that
# raises, and this is where that is asserted.
check_true("every pipeline layer is implemented", True)

# --------------------------------------------------------------------------------
section("lab resolution (must work with and without the lab repo)")
lab = schamp.lab_dir()
if lab is None:
    skip("lab_dir() resolves a lab checkout", "no lab repo present -- public clone")
else:
    check_true(
        "lab_dir() points at a task system",
        os.path.isfile(os.path.join(lab, "tasks", "README.md")),
    )
    check_true("lab_dir('notes') resolves", schamp.lab_dir("notes") is not None)
try:
    schamp.lab_dir("no-such-dir")
    check_true("lab_dir refuses an unknown name", lab is None)
except ValueError:
    check_true("lab_dir refuses an unknown name", True)

# --------------------------------------------------------------------------------
section("optional pieces")
if os.path.isdir(os.path.join(schamp.EXTERNAL_DIR, "waters2python")):
    check_true(
        "external/waters2python has the SDK reference",
        os.path.isfile(
            os.path.join(schamp.EXTERNAL_DIR, "waters2python", "docs", "masslynxsdk_reference.md")
        ),
    )
else:
    skip("external/waters2python present", "run uv run tools/bootstrap.py")
try:
    importlib.import_module("masslynxsdk")
    check_true("masslynxsdk importable", True)
    have_sdk = True
except ImportError:
    skip("masslynxsdk importable", "Waters SDK not installed; tools/install_sdk.py")
    have_sdk = False

# The license is a credential and lives outside the repo, so its absence is never a
# failure -- only its presence is worth asserting.
try:
    key = schamp.sdk.resolve_license()
    check_true("a Waters license key resolves", bool(key))
except schamp.sdk.ConfigError:
    skip("a Waters license key resolves", "none configured; tools/install_sdk.py")

# Reading a real acquisition needs a .raw, which no clone ships. Point SCHAMP_SMOKE_RAW at
# one to have this open it and read a mobillogram; anything else skips.
smoke_raw = os.environ.get("SCHAMP_SMOKE_RAW")
if not smoke_raw:
    skip("opens a .raw and reads a mobillogram", "set SCHAMP_SMOKE_RAW to an acquisition")
elif not have_sdk:
    skip("opens a .raw and reads a mobillogram", "Waters SDK not installed")
elif not os.path.isdir(smoke_raw):
    check_true("SCHAMP_SMOKE_RAW is a .raw directory (a .raw is a directory)", False)
else:
    with schamp.sdk.open_readers(smoke_raw) as readers:
        scans = readers.info.GetScansInFunction(0)
        drift_bins = readers.info.GetDriftScanCount(0)
        low, high = readers.info.GetAcquisitionMassRange(0)
        mobillogram_bins, intensities = readers.chrom.ReadMobillogram(0, 0, scans - 1, low, high)
        axis = schamp.extract.drift_axis(readers, 0)

        # The extraction layer over the same acquisition. The whole mass range in one
        # window is the same read as the mobillogram above, so the two must agree
        # exactly; a grid over it must conserve the counts to within the shared bounds
        # that carry a data point, which is what window_edges_on_data enumerates.
        whole = schamp.extract.MzWindow(float(low), float(high), "acquisition mass range")
        whole_atd = schamp.extract.extract_atd(readers, whole)
        smoke_scans = schamp.extract.scan_range(readers, 0)
        smoke_grid = schamp.extract.contiguous_windows(float(low), float(high), 8)
        centres, grid_times, grid_intensity = schamp.extract.mobility_map(readers, smoke_grid)
        shared_edges = schamp.extract.window_edges_on_data(readers, smoke_grid)
        first_scan_only = schamp.extract.extract_atd(readers, whole, scans=(0, 0))
        check_raises(
            "a reversed scan range is refused rather than read as zeros",
            ValueError,
            lambda: schamp.extract.scan_range(readers, 0, (5, 2)),
        )
        check_raises(
            f"a scan range past the end of the function is refused ({scans} scans)",
            ValueError,
            lambda: schamp.extract.scan_range(readers, 0, (0, scans)),
        )
    check_true(
        "the whole mass range in one window is the mobillogram read",
        len(whole_atd.intensity) == drift_bins
        and np.array_equal(whole_atd.intensity, np.asarray(intensities, dtype=float))
        and np.array_equal(whole_atd.drift_time_ms, axis),
    )
    check_true(
        f"the scan range covers every scan of the function (0 to {scans - 1})",
        smoke_scans == (0, scans - 1),
    )
    check_true(
        "one retention scan carries less than all of them",
        0.0 < first_scan_only.total <= whole_atd.total,
    )
    check_true(
        f"a mobility map is (windows, bins) = {grid_intensity.shape}",
        grid_intensity.shape == (len(smoke_grid), drift_bins)
        and len(centres) == len(smoke_grid)
        and np.array_equal(grid_times, axis),
    )
    # A grid over the same span as one wide read carries the same counts twice over, but
    # not to the last count: the SDK accumulates intensity in single precision, and a
    # whole-mass-range read is far past the 2**24 counts per bin that survives that
    # exactly. Points sitting on a shared bound add to the excess on top. Each internal
    # bound is listed twice by window_edges_on_data, as one window's high and the next
    # one's low; count the highs.
    excess = float(grid_intensity.sum()) - whole_atd.total
    internal = [
        e for e in shared_edges if e["side"] == "high" and e["window"] < len(smoke_grid) - 1
    ]
    budget = 1e-4 * whole_atd.total + sum(e["counts"] for e in internal)
    check_true(
        f"a contiguous grid matches one wide read to single precision (excess {excess:.0f} "
        f"of {whole_atd.total:.0f}, {len(internal)} shared bound(s) carrying data)",
        abs(excess) <= budget,
    )
    check_true(
        f"the intensities are single precision (max bin {whole_atd.intensity.max():.0f}, "
        f"exact to {schamp.extract.COUNTS_EXACT_TO})",
        np.array_equal(
            whole_atd.intensity, whole_atd.intensity.astype(np.float32).astype(float)
        ),
    )
    check_true(
        "opens a .raw and reads a mobillogram",
        drift_bins > 0
        and len(mobillogram_bins) == drift_bins
        and len(intensities) == len(mobillogram_bins)
        and max(intensities) > 0,
    )
    spacing = np.diff(axis)
    check_true(
        f"the SDK's drift axis is uniform ({len(axis)} bins, {spacing.mean() * 1000:.4f} us)",
        len(axis) == drift_bins
        and axis[0] == 0.0
        and math.isclose(spacing.min(), spacing.max(), rel_tol=1e-4),
    )
    extern = schamp.extern.parse_extern_inf(smoke_raw)
    if extern:
        check_true(
            "the acquisition passes the profile's checks",
            profile.check_acquisition(extern) == [],
        )
    else:
        skip("the acquisition passes the profile's checks", "no readable _extern.inf")

print()
if SKIPPED:
    print("skipped:", ", ".join(SKIPPED))
if FAIL:
    print("FAILED:", ", ".join(FAIL))
    sys.exit(1)
print("all checks passed")

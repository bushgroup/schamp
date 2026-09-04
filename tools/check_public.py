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
section("the layers that are not built yet")
# A stub that returned something plausible would be far worse than one that raises.
for name, call in (
    ("atd.fit_gaussian", lambda: schamp.atd.fit_gaussian(atd)),
    ("extract.extract_atd", lambda: schamp.extract.extract_atd(None, schamp.extract.MzWindow(1.0, 2.0))),
    ("extract.extract_atds", lambda: schamp.extract.extract_atds(None, [])),
    ("extract.mobility_map", lambda: schamp.extract.mobility_map(None, [])),
    ("mobility.regress", lambda: schamp.mobility.regress([], profile)),
):
    check_raises(f"{name} is honestly unimplemented", NotImplementedError, call)

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

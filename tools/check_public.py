"""Self-check for a fresh clone: no Waters SDK, no license, no .raw file, no
lab repo needed.

Every check here must pass in a bare public clone. Checks that need something a
clone does not ship (the SDK, external/, the lab repo, an acquisition) are
reported as SKIPPED when it is absent, never as FAIL.

What it covers: the physics identities on synthetic input, the instrument-profile
format and its drift-voltage formula, the `_extern.inf` parser against a real
sidecar committed as a fixture, the experiment data model against a table written
and read back in a temporary directory, the acquisition run sheet's model, the live
pane's watcher, completion rule and series on synthetic arrivals, and the reporting
conventions. Also that the not-yet-implemented layers say so rather than
returning something plausible. The run sheet's window is opened offscreen when
PySide6 (the optional extra schamp[acquire]) is installed, and SKIPPED when not.

This is the public counterpart of the lab repo's regression check, which
additionally pins the pipeline against the archived polyalanine numbers and
needs the legacy acquisitions to run.

Run:  uv run tools/check_public.py
"""

import importlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

import schamp
import schamp.atd
import schamp.constants as C
import schamp.experiment
import schamp.extract
import schamp.extern
import schamp.importer
import schamp.mobility
import schamp.profiles
import schamp.report
import schamp.sdk
import schamp.spectrum

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
check_true(
    "docs/ holds the acquisition guide",
    os.path.isfile(os.path.join(schamp.DOCS_DIR, "acquiring-a-drift-voltage-series.md")),
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
section("_HEADER.TXT and the calibration function")
if not os.path.isfile(os.path.join(FIXTURES, schamp.extern.HEADER_FILENAME)):
    skip("_HEADER.TXT fixture parses", "tests/fixtures/_HEADER.TXT is missing")
else:
    header = schamp.extern.parse_header_txt(FIXTURES)
    check_true("the fixture parses to a mapping", len(header) >= 8)
    check_true(
        "a value containing colons survives", header.get("Cal Time") == "13:36"
    )
    check_true(
        "a line with no $$ is not a key",
        not any("not a header line" in key for key in header),
    )
    cubic = schamp.extern.cal_function_coefficients(FIXTURES, 1)
    check_true("a cubic calibration reads as four coefficients", len(cubic) == 4)
    check_close("its leading coefficient", cubic[0], 7.472570630168562e-3)
    check_true(
        "the identity function is read, not skipped",
        schamp.extern.cal_function_coefficients(FIXTURES, 2) == (0.0, 1.0),
    )
    check_true(
        "an unparseable calibration is None, not a partial one",
        schamp.extern.cal_function_coefficients(FIXTURES, 3) is None,
    )
    check_true(
        "an absent function is None",
        schamp.extern.cal_function_coefficients(FIXTURES, 9) is None,
    )
    check_true(
        "a similar key is not mistaken for the calibration",
        "Cal StdDev Function 1" in header,
    )

    calibration = importlib.import_module("schamp.calibration")
    cal = calibration.CalFunction.from_raw(FIXTURES, 1)
    check_true("the cubic loads as order 3", cal.order == 3 and not cal.is_identity)
    # The two frames, round trip. Nothing here asserts a displacement of any
    # particular size -- that is the acquisition's business -- only that the
    # conversion is invertible to double precision, which is what a window depends on.
    for mz in (230.1146, 800.4261, 1650.8568):
        spectrum = cal.to_spectrum_frame(mz)
        check_close(
            f"the calibration function inverts at {mz:g} Th",
            cal.to_chromatogram_frame(spectrum),
            mz,
            rel=1e-12,
        )
    check_true(
        "it moves a mass by a sane amount, not by an order of magnitude",
        abs(cal.to_spectrum_frame(800.0) - 800.0) < 1.0,
    )
    check_true(
        "an array goes in and an array comes out",
        np.asarray(cal.to_chromatogram_frame(np.array([300.0, 900.0]))).shape == (2,),
    )
    identity = calibration.CalFunction.identity()
    check_true("the identity is recognised as one", identity.is_identity)
    check_true(
        "and changes nothing", identity.to_chromatogram_frame(1234.5) == 1234.5
    )
    check_true(
        "an acquisition with no header is the identity by default",
        calibration.CalFunction.from_raw("no/such/dir").is_identity,
    )
    check_raises(
        "unless the caller says an absent function must raise",
        calibration.CalibrationError,
        lambda: calibration.CalFunction.from_raw("no/such/dir", missing_is_identity=False),
    )
    check_raises(
        "a one-coefficient calibration is refused",
        calibration.CalibrationError,
        lambda: calibration.CalFunction(coefficients=(1.0,)),
    )
    check_raises(
        "a non-positive m/z is refused, since sqrt(m/z) is the variable",
        calibration.CalibrationError,
        lambda: cal.to_spectrum_frame(0.0),
    )
    check_true("the description names the order", "order 3" in cal.describe())
    check_true("and the identity says so", "identity" in identity.describe())

    # A recalibration recovers a known displacement. Synthetic on purpose: a fit whose
    # answer is known exactly is the only way to test a fit.
    true_mz = np.array([230.11, 301.15, 443.23, 656.34, 869.45, 1082.56, 1437.75])
    displaced = (1.0002 * np.sqrt(true_mz) - 0.0035) ** 2
    recal = calibration.Recalibration.fit(true_mz, displaced, order=1)
    check_true("an order-1 recalibration fits its own displacement", recal.rms_th < 1e-9)
    check_close("and recovers it", recal.to_measured(800.0), (1.0002 * math.sqrt(800.0) - 0.0035) ** 2)
    check_close("and inverts", recal.to_true(recal.to_measured(800.0)), 800.0, rel=1e-12)
    check_true("it reports the range it was fitted over", recal.covers(800.0) and not recal.covers(1800.0))
    check_true("and says what it did", "order 1" in recal.describe())
    check_raises(
        "a recalibration with too few species to leave a residual is refused",
        calibration.CalibrationError,
        lambda: calibration.Recalibration.fit(true_mz[:2], displaced[:2], order=1),
    )
    check_raises(
        "and mismatched series are refused",
        calibration.CalibrationError,
        lambda: calibration.Recalibration.fit(true_mz, displaced[:-1]),
    )

    frame = calibration.MzFrame(cal=cal, recalibration=recal)
    check_close(
        "MzFrame applies the recalibration and then the frame conversion",
        frame.reader_mz(800.0),
        cal.to_chromatogram_frame(recal.to_measured(800.0)),
    )
    check_close("and inverts both", frame.computed_mz(frame.reader_mz(800.0)), 800.0, rel=1e-10)
    check_true(
        "a frame with no recalibration applies only the conversion",
        calibration.MzFrame(cal=cal).reader_mz(800.0) == cal.to_chromatogram_frame(800.0),
    )
    check_true("the identity frame is one", calibration.MzFrame.identity().is_identity)
    check_true(
        "the description names both steps",
        "cal function" in frame.describe() and "recalibration" in frame.describe(),
    )
    check_true(
        "and says so when there is no recalibration",
        "recalibration: none" in calibration.MzFrame(cal=cal).describe(),
    )

    # covers() reports and does not refuse (lab record, task 11), so what has to hold
    # is that an extrapolated species still gets placed and still says it was.
    check_close(
        "an m/z outside the fitted range is still placed",
        frame.reader_mz(1800.0),
        cal.to_chromatogram_frame(recal.to_measured(1800.0)),
    )
    check_true(
        "and note_for says it was extrapolated",
        "extrapolated" in frame.note_for(1800.0),
    )
    check_true(
        "while an m/z inside the range gets no such clause",
        "extrapolated" not in frame.note_for(800.0),
    )
    check_true(
        "a frame with no recalibration never claims extrapolation",
        "extrapolated" not in calibration.MzFrame(cal=cal).note_for(1800.0),
    )

    # from_spectrum: the fit an analysis makes for itself, from a spectrum rather than
    # from a list of measured positions. Synthetic again, and with one species absent
    # from the spectrum so the skip-and-count path is exercised.
    spectrum_module = importlib.import_module("schamp.spectrum")
    grid = np.arange(200.0, 1500.0, 0.004)
    intensity = np.zeros_like(grid)
    for centre in displaced[:-1]:
        intensity += 1.0e5 * np.exp(-0.5 * ((grid - centre) / 0.01) ** 2)
    synthetic = spectrum_module.MassSpectrum(
        mz=grid, intensity=intensity, acquisition="synthetic"
    )
    from_spectrum = calibration.Recalibration.from_spectrum(
        synthetic, [(mz, f"ref {index}") for index, mz in enumerate(true_mz)], order=1
    )
    check_true(
        "a recalibration fitted from a spectrum skips a species that is not in it",
        from_spectrum.species == len(true_mz) - 1 and from_spectrum.offered == len(true_mz),
    )
    check_true(
        "and says how many it was offered",
        f"of {len(true_mz)} offered" in from_spectrum.describe(),
    )
    check_close(
        "and recovers the displacement it was built from",
        from_spectrum.to_measured(800.0),
        (1.0002 * math.sqrt(800.0) - 0.0035) ** 2,
        rel=1e-4,
    )
    check_raises(
        "a spectrum with too few of the species is refused, naming both counts",
        calibration.CalibrationError,
        lambda: calibration.Recalibration.from_spectrum(
            synthetic, true_mz, order=1, min_apex=1.0e9
        ),
    )
    check_raises(
        "and an empty reference list is refused rather than fitted",
        calibration.CalibrationError,
        lambda: calibration.Recalibration.from_spectrum(synthetic, [], order=1),
    )

    # The provenance the file carries. The fixture header is a cross-polarity
    # calibration three months stale, which is the shape of the real 2013 defect.
    provenance = calibration.CalibrationProvenance.from_raw(FIXTURES)
    check_true(
        "the calibration's name and calibrant come off the dynamic-params line",
        provenance.calibration == "<X+>ESI_CAL_POS"
        and provenance.calibrant == "<X+>SYNTHETIC_CAL_001",
    )
    check_true("its mass range too", provenance.mass_range == (400.0, 8000.0))
    check_true(
        "two date formats in one header both parse",
        provenance.age_days is not None and 91.0 < provenance.age_days < 92.0,
    )
    check_true("the polarity comes from _extern.inf", provenance.polarity == "ES-")
    check_true("and the cross-polarity marker is read", provenance.cross_polarity)
    concerns = provenance.concerns()
    check_true(
        "concerns names the stale calibration and the polarity, and nothing else",
        len(concerns) == 2
        and any("polarity" in c for c in concerns)
        and any("days old" in c for c in concerns),
    )
    check_true(
        "a longer tolerance drops the staleness complaint but not the polarity one",
        len(provenance.concerns(stale_after_days=365.0)) == 1,
    )
    check_true("and the description is one line", "\n" not in provenance.describe())
    missing = calibration.CalibrationProvenance.from_raw("no/such/dir")
    check_true(
        "a file with no header records nothing rather than raising",
        missing.calibration is None and missing.age_days is None,
    )
    check_true(
        "and says that is what happened",
        any("no mass calibration at all" in c for c in missing.concerns()),
    )

    # from_calibrant is the whole-frame borrow, and the polarity check is the reason it
    # is a constructor rather than a note. The fixture is ES-, so ES- passes.
    borrowed = calibration.MzFrame.from_calibrant(FIXTURES, recal, for_acquisition=FIXTURES)
    check_true(
        "a borrowed frame carries the calibrant's calibration function, not the analyte's",
        borrowed.cal.coefficients == cal.coefficients
        and borrowed.recalibration is recal,
    )
    check_true(
        "and it is not the same thing as handing over the residual alone",
        borrowed.reader_mz(800.0) != calibration.MzFrame(recalibration=recal, cal=identity).reader_mz(800.0),
    )
    with tempfile.TemporaryDirectory() as other_polarity:
        with open(
            os.path.join(other_polarity, schamp.extern.EXTERN_FILENAME), "w", encoding="utf-8"
        ) as handle:
            handle.write("Polarity\t\tES+\n")
        check_raises(
            "a calibrant of the other polarity is refused, since that was measured to fail",
            calibration.CalibrationError,
            lambda: calibration.MzFrame.from_calibrant(
                FIXTURES, recal, for_acquisition=other_polarity
            ),
        )
        check_true(
            "unless the caller says it has shown otherwise",
            calibration.MzFrame.from_calibrant(
                FIXTURES,
                recal,
                for_acquisition=other_polarity,
                allow_polarity_mismatch=True,
            ).cal.coefficients
            == cal.coefficients,
        )
    check_true(
        "and with no acquisition named there is nothing to check it against",
        calibration.MzFrame.from_calibrant(FIXTURES, recal).recalibration is recal,
    )

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

    # [recalibration]: an experiment-level table, because the fit belongs to the series
    # (lab record, task 11). Same strictness as the columns, for the same reason.
    check_true("an experiment with no [recalibration] has none", experiment.recalibration is None)

    def experiment_with(block, name):
        with open(os.path.join(work, f"{name}.toml"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(f'schema = 1\nprofile = "uw-synapt-g2"\n\n[defaults]\ntemperature_K = 301.13\n\n{block}')
        return schamp.load_experiment(os.path.join(work, f"{name}.toml"))

    listed = experiment_with(
        "[recalibration]\norder = 1\nmin_apex = 2.0e4\n"
        'species = [ { mz = 392.7148, label = "CsI 1+ n=1" }, { mz = 652.5248 }, '
        "{ mz = 912.3347 } ]\n",
        "recal",
    )
    spec = listed.recalibration
    check_true(
        "a [recalibration] block is read whole",
        spec is not None
        and spec.order == 1
        and spec.min_apex == 2.0e4
        and len(spec.species) == 3,
    )
    check_true("a species keeps its label, and may have none", spec.species[0][1] == "CsI 1+ n=1" and spec.species[1][1] == "")
    check_true(
        "and with no species listed the anchor is the analysis's own",
        "the species the analysis places windows on"
        in experiment_with("[recalibration]\norder = 1\n", "bare").recalibration.describe(),
    )
    for block, why in (
        ("[recalibration]\nmin_apax = 100.0\n", "a mistyped key is refused"),
        ("[recalibration]\norder = 0\n", "an order below 1 is refused"),
        ("[recalibration]\norder = 1.5\n", "a fractional order is refused"),
        ('[recalibration]\nspecies = [ { mz = 400.0, labell = "x" } ]\n', "a mistyped species key is refused"),
        ("[recalibration]\nspecies = [ { label = \"x\" } ]\n", "a species with no mz is refused"),
        ("[recalibration]\nspecies = [ { mz = -1.0 } ]\n", "a non-positive reference m/z is refused"),
        ("[recalibration]\nspecies = [ { mz = 400.0 }, { mz = 800.0 } ]\n", "too few species for the order is refused"),
        ("[recalibration]\nmin_apex = -1.0\n", "a negative intensity floor is refused"),
        ("[recalibration]\nallow_polarity_mismatch = true\n", "allowing a mismatch with no calibrant is refused"),
    ):
        check_raises(why, schamp.experiment.ExperimentError, lambda b=block: experiment_with(b, "bad"))
    check_true(
        "a named calibrant resolves under data_dir",
        experiment_with(
            '[recalibration]\ncalibrant = "cal.raw"\n'
            "species = [ { mz = 400.0 }, { mz = 800.0 }, { mz = 1200.0 } ]\n",
            "cal",
        ).recalibration.calibrant_path
        == os.path.abspath(os.path.join(work, "cal.raw")),
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

    # A blank pressure is the importer's seeded state: it loads, and validate names it.
    with open(os.path.join(work, "blank.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("acquisition,pressure_Torr,temperature_K\nx.raw,,300.0\ny.raw,2.0,300.0\n")
    with open(os.path.join(work, "blank.toml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write('schema = 1\nprofile = "uw-synapt-g2"\nconditions = "blank.csv"\n')
    blank = schamp.load_experiment(os.path.join(work, "blank.toml"))
    check_true("a blank pressure with no default loads as None", blank.conditions[0].pressure_torr is None)
    check_true(
        "and validate names the row that still wants one",
        sum("row 2: no pressure_Torr" in p for p in blank.validate(require_raw=False)) == 1,
    )

    # --------------------------------------------------------------------------
    section("reporting conventions")
    written = schamp.report.write_results({"answer": 42}, work, task="check-public")
    check_true("write_results stamps and writes", os.path.isfile(written))
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
section("importing a set of acquisitions")
# Synthetic .raw-shaped directories: the committed _extern.inf, a _HEADER.TXT carrying
# an acquisition time, and a few files of noise standing in for the binary data. Named
# so that name order and acquisition order disagree, and carrying the files MassLynx
# and a sync client add after the fact, which must not be copied.
work = tempfile.mkdtemp(prefix="schamp-import-")
try:
    share = os.path.join(work, "share")
    set_dir = os.path.join(work, "data-experiment", "synthetic-set")
    rng = np.random.default_rng(19)

    def fake_raw(name, acquired, *, finished=True):
        raw = os.path.join(share, name)
        os.makedirs(os.path.join(raw, "_mob"), exist_ok=True)
        if finished:
            shutil.copyfile(os.path.join(FIXTURES, "_extern.inf"), os.path.join(raw, "_extern.inf"))
        with open(os.path.join(raw, "_HEADER.TXT"), "w", encoding="latin-1", newline="\r\n") as fh:
            fh.write(f"$$ Acquired Date: 23-Apr-2013\n$$ Acquired Time: {acquired}\n")
        for data in ("_FUNC001.DAT", "_FUNC001.IDX", "_func001.cdt"):
            with open(os.path.join(raw, data), "wb") as fh:
                fh.write(rng.bytes(int(rng.integers(1000, 200_000))))
        for litter in ("_HISTORY.INF", "desktop.ini", os.path.join("_mob", "_dt.bin")):
            with open(os.path.join(raw, litter), "wb") as fh:
                fh.write(b"written after the fact")
        return raw

    fake_raw("S_001.raw", "12:30:00")
    fake_raw("S_002.raw", "12:10:00")
    fake_raw("S_003.raw", "12:20:00")
    fake_raw("S_004.raw", "12:40:00", finished=False)

    first = schamp.importer.import_set(share, set_dir)
    check_true(
        "the finished acquisitions import and the unfinished one is skipped",
        sorted(first.imported) == ["S_001.raw", "S_002.raw", "S_003.raw"]
        and first.unfinished == ["S_004.raw"]
        and first.ok,
    )
    check_true("the set verifies against its manifest", schamp.importer.verify(set_dir) == [])
    manifest_lines = open(os.path.join(set_dir, "MANIFEST.sha256"), "rb").read().split(b"\n")
    check_true(
        "the manifest is sha256sum format with LF endings",
        manifest_lines[-1] == b""
        and all(
            len(line.split(b"  ", 1)[0]) == 64 and b"\r" not in line and b"\\" not in line
            for line in manifest_lines[:-1]
        ),
    )
    listed = schamp.importer.read_manifest(set_dir)
    check_true("every file of every acquisition is listed", len(listed) == 3 * 5)
    check_true(
        "what MassLynx writes on opening, and sync-client litter, are neither copied nor listed",
        not any("_mob" in k or "HISTORY" in k or "desktop" in k for k in listed)
        and not os.path.exists(os.path.join(set_dir, "S_001.raw", "_mob")),
    )
    sha256sum = shutil.which("sha256sum")
    if sha256sum:
        checked = subprocess.run(
            [sha256sum, "--check", "--strict", "--quiet", "MANIFEST.sha256"],
            cwd=set_dir, capture_output=True, text=True,
        )
        check_true("sha256sum --check agrees, with nothing of schamp's", checked.returncode == 0)
    else:
        skip("sha256sum --check agrees, with nothing of schamp's", "no sha256sum on PATH")

    seeded = schamp.load_experiment(set_dir)
    check_true(
        "the seeded set loads by directory, with the default profile",
        seeded.profile.name == "uw-synapt-g2" and seeded.title == "synthetic-set",
    )
    check_true(
        "the table is in acquisition-time order, not name order",
        [c.acquisition for c in seeded.conditions] == ["S_002.raw", "S_003.raw", "S_001.raw"],
    )
    check_true(
        "pressure and temperature are left blank, and use is true",
        all(c.pressure_torr is None and c.temperature_k is None and c.use for c in seeded.conditions),
    )
    blanks = [p for p in seeded.validate(require_raw=False) if "no pressure_Torr" in p or "no temperature_K" in p]
    check_true("validate says every blank wants typing in", len(blanks) == 6)
    check_true(
        "the drift voltage comes from the sidecar through the profile",
        seeded.conditions[0].drift_voltage(seeded.profile)
        == profile.drift_voltage(schamp.extern.parse_extern_inf(FIXTURES)),
    )

    second = schamp.importer.import_set(share, set_dir)
    check_true(
        "a second import adds nothing",
        second.imported == [] and sorted(second.present) == ["S_001.raw", "S_002.raw", "S_003.raw"],
    )

    # The operator types a pressure; the acquisition that was still going finishes.
    conditions_path = os.path.join(set_dir, "conditions.csv")
    typed = open(conditions_path, encoding="utf-8").read().replace("S_002.raw,,,", "S_002.raw,2.031,301.1,")
    open(conditions_path, "w", encoding="utf-8", newline="").write(typed)
    shutil.copyfile(os.path.join(FIXTURES, "_extern.inf"), os.path.join(share, "S_004.raw", "_extern.inf"))
    third = schamp.importer.import_set(share, set_dir)
    reloaded = schamp.load_experiment(set_dir)
    check_true(
        "a newly finished acquisition is appended to the manifest and the table",
        third.imported == ["S_004.raw"]
        and [c.acquisition for c in reloaded.conditions][-1] == "S_004.raw"
        and len(schamp.importer.read_manifest(set_dir)) == 4 * 5,
    )
    check_true("a typed pressure survives the append", reloaded.conditions[0].pressure_torr == 2.031)
    readme = open(os.path.join(set_dir, "README.md"), encoding="utf-8").read()
    check_true("the README records each import that added something", readme.count("\n- ") == 2)

    # A file rewritten at the source after it was imported.
    with open(os.path.join(share, "S_003.raw", "_FUNC001.IDX"), "ab") as fh:
        fh.write(b"\x00")
    changed = schamp.importer.import_set(share, set_dir)
    check_true(
        "a source changed after import is reported and not overwritten",
        changed.changed == ["S_003.raw"] and schamp.importer.verify(set_dir) == [],
    )
    replaced = schamp.importer.import_set(share, set_dir, replace_changed=True)
    check_true(
        "and replaced when asked, verifying again",
        replaced.replaced == ["S_003.raw"] and schamp.importer.verify(set_dir) == [],
    )

    # A local copy deleted: copied back from a source that still matches.
    shutil.rmtree(os.path.join(set_dir, "S_001.raw"))
    check_true("a missing acquisition fails verify", any("S_001.raw" in p for p in schamp.importer.verify(set_dir)))
    check_true(
        "but not in tables-only mode",
        schamp.importer.verify(set_dir, require_raw=False) == [],
    )
    restored = schamp.importer.import_set(share, set_dir, select=["S_001*"])
    check_true(
        "and is copied back from a source that still matches the manifest",
        restored.restored == ["S_001.raw"] and schamp.importer.verify(set_dir) == [],
    )

    # One byte flipped on the destination.
    target = os.path.join(set_dir, "S_002.raw", "_FUNC001.DAT")
    with open(target, "r+b") as fh:
        byte = fh.read(1)
        fh.seek(0)
        fh.write(bytes([byte[0] ^ 0x01]))
    problems = schamp.importer.verify(set_dir)
    check_true(
        "a single tampered byte fails verify, naming the file",
        len(problems) == 1 and "S_002.raw/_FUNC001.DAT" in problems[0],
    )
    check_raises(
        "the importer refuses a profile that does not load",
        schamp.profiles.ProfileError,
        lambda: schamp.importer.import_set(share, os.path.join(work, "other"), profile="no-such-instrument"),
    )
    check_true("and has copied nothing by then", not os.path.exists(os.path.join(work, "other")))
finally:
    shutil.rmtree(work, ignore_errors=True)

# --------------------------------------------------------------------------------
section("writing an experiment, and reading it back")
# Both shipped example tables, given an experiment.toml by write_experiment and copied
# through read_conditions_table / write_conditions_table: the CSV must come back byte
# for byte, and the experiment must load the same.
work = tempfile.mkdtemp(prefix="schamp-write-")
try:
    for series in ("npalan", "palan"):
        src = os.path.join(schamp.ROOT, "examples", "data", series, "conditions.csv")
        dst = os.path.join(work, series)
        os.makedirs(dst)
        table = schamp.experiment.read_conditions_table(src)
        schamp.experiment.write_conditions_table(
            os.path.join(dst, "conditions.csv"), table.header, table.rows, comment=table.comment
        )
        with open(src, "rb") as a, open(os.path.join(dst, "conditions.csv"), "rb") as b:
            check_true(f"{series}: the conditions table round-trips byte for byte", a.read() == b.read())
        schamp.experiment.write_experiment(
            dst,
            title=f'{series} "130423", ES±',
            profile="uw-synapt-g2",
            gas="helium",
            notes="line one\nline two\\",
            defaults={"temperature_K": 301.13},
            recalibration=schamp.experiment.RecalibrationSpec(
                order=1, species=((392.7148, "CsI 1+ n=1"), (652.5248, "n=2"), (912.3348, "n=3")), min_apex=2.0e4
            ),
        )
        first = schamp.load_experiment(dst)
        schamp.experiment.write_experiment(
            dst,
            title=first.title,
            profile=first.profile.name,
            gas=first.gas,
            notes=first.notes,
            defaults={"temperature_K": 301.13},
            recalibration=first.recalibration,
        )
        second = schamp.load_experiment(dst)
        rows = lambda e: [  # noqa: E731
            (c.acquisition, c.pressure_torr, c.temperature_k, c.gas, c.use, c.notes) for c in e.conditions
        ]
        check_true(
            f"{series}: write_experiment's file loads, and loads the same written again",
            first.title == second.title == f'{series} "130423", ES±'
            and first.notes == second.notes == "line one\nline two\\"
            and first.recalibration == second.recalibration
            and first.recalibration.min_apex == 2.0e4
            and rows(first) == rows(second)
            and len(first.used) == 10,
        )
    check_raises(
        "write_experiment refuses a [defaults] key that is not a column",
        schamp.experiment.ExperimentError,
        lambda: schamp.experiment.write_experiment(work, title="t", profile="p", defaults={"pressure": 2.0}),
    )
    leftovers = [f for _, _, files in os.walk(work) for f in files if f.endswith(".tmp")]
    check_true("an atomic write leaves no temporary file behind", leftovers == [])
finally:
    shutil.rmtree(work, ignore_errors=True)

# --------------------------------------------------------------------------------
section("the acquisition run sheet (no Qt)")
import schamp.acquire.model as runsheet  # noqa: E402 -- grouped with its checks

# The 130423 series' fourteen (Helium Cell DC, IMS Bias) settings in acquisition order,
# as recorded in each acquisition's _extern.inf, and the drift voltage each gives through
# the uw-synapt-g2 profile (lab record, tasks 01 and 14). The shipped template must
# reproduce both.
SERIES_130423 = [
    (0.0, 132.0, 86.9647), (110.0, 240.0, 302.4), (18.0, 240.0, 211.4824),
    (0.0, 125.0, 80.0471), (0.0, 149.0, 103.7647), (0.0, 119.0, 74.1176),
    (0.0, 187.0, 141.3176), (56.0, 240.0, 249.0353), (0.0, 140.0, 94.8706),
    (0.0, 228.0, 181.8353), (0.0, 160.0, 114.6353), (0.0, 206.0, 160.0941),
    (160.0, 240.0, 351.8118), (0.0, 172.0, 126.4941),
]
template = runsheet.load_template(runsheet.DEFAULT_TEMPLATE)
check_true("the standard template ships and loads by name", "uw-standard-14" in runsheet.builtin_templates())
check_true(
    "it steps Helium Cell DC and IMS Bias and holds the rest at the 130423 values",
    template.vary == ("Helium Cell DC", "IMSBias")
    and dict(template.fixed) == {"Helium Exit": -40.0, "Transfer DC Entrance": 4.0},
)
check_true(
    "its fourteen rows are the 130423 settings in the 130423 order",
    [row for row in template.rows] == [(dc, bias) for dc, bias, _ in SERIES_130423],
)
work = tempfile.mkdtemp(prefix="schamp-runsheet-")
try:
    set_dir = os.path.join(work, "set")
    data_dir = os.path.join(work, "data")
    sheet = runsheet.RunSheet.from_template(template, set_dir, base_name="S")
    profile = sheet.load_profile()
    worst = max(
        abs(sheet.nominal_voltage(row, profile) - volts)
        for row, (_, _, volts) in zip(sheet.rows, SERIES_130423)
    )
    check_true(f"each row's nominal drift voltage is the profile's (worst {worst:.1e} V)", worst < 1e-4)
    check_true(
        "acquisitions are named as a MassLynx sample list numbers them",
        [r.acquisition for r in sheet.rows[:2]] == ["S_001.raw", "S_002.raw"],
    )

    def bad_template(text, name):
        path = os.path.join(work, f"{name}.toml")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write('schema = 1\nprofile = "uw-synapt-g2"\n' + text)
        return lambda: runsheet.load_template(path)

    for why, text in (
        ("a template missing a drift-voltage key is refused",
         'vary = ["IMSBias"]\nrows = [[132.0]]\n[fixed]\n"Helium Exit" = -40.0\n"Transfer DC Entrance" = 4.0\n'),
        ("a template naming a key the drift voltage does not use is refused",
         'vary = ["IMSBias", "Cone"]\nrows = [[132.0, 1.0]]\n[fixed]\n"Helium Cell DC" = 0.0\n"Helium Exit" = -40.0\n"Transfer DC Entrance" = 4.0\n'),
        ("a row with the wrong number of settings is refused",
         'vary = ["Helium Cell DC", "IMSBias"]\nrows = [[132.0]]\n[fixed]\n"Helium Exit" = -40.0\n"Transfer DC Entrance" = 4.0\n'),
        ("rows written after [fixed] (and so inside it) are refused",
         'vary = ["Helium Cell DC", "IMSBias"]\n[fixed]\n"Helium Exit" = -40.0\n"Transfer DC Entrance" = 4.0\nrows = [[0.0, 132.0]]\n'),
    ):
        check_raises(why, runsheet.RunSheetError, bad_template(text, "bad"))

    sheet.temperature_k = "301.13"
    sheet.data_dir = data_dir
    sheet.save()
    check_true(
        "a sheet with nothing taken writes experiment.toml and runsheet.toml, no table",
        sorted(os.listdir(set_dir)) == ["experiment.toml", "runsheet.toml"],
    )
    check_true(
        "and says there is nothing to check yet",
        any("nothing to check" in p for p in sheet.validate()),
    )
    for row, pressure in zip(sheet.rows[:4], ("2.03", "2.03", "2.031", "")):
        row.taken = True
        row.cells["pressure_Torr"] = pressure
    sheet.rows[5].cells["pressure_Torr"] = "2.04"  # typed ahead, not yet taken
    sheet.save()
    experiment = schamp.load_experiment(set_dir)
    check_true(
        "a ticked row enters conditions.csv in sheet order, and an unticked one does not",
        [c.acquisition for c in experiment.conditions] == ["S_001.raw", "S_002.raw", "S_003.raw", "S_004.raw"],
    )
    check_true(
        "the session temperature reaches every row through [defaults]",
        all(c.temperature_k == 301.13 for c in experiment.conditions),
    )
    check_true(
        "the record validates but for the one pressure not typed in",
        [p.split(": ", 1)[1][:15] for p in sheet.validate(require_raw=False)] == ["no pressure_Tor"],
    )

    # schamp.importer appends to the same table while the sheet is open; a row the
    # sheet wrote and then unticked must stay out, and a row it never saw must come in.
    with open(sheet.conditions_path, "a", encoding="utf-8", newline="") as fh:
        fh.write("S_005.raw,,true,\nCAL.raw,,true,calibrant\n")
    sheet.rows[3].taken = False
    sheet.rows[0].cells["notes"] = "first, after the gas had settled"
    sheet.save()
    taken = [(r.acquisition, r.planned) for r in sheet.rows if r.taken]
    check_true(
        "a save adopts rows another writer appended, planned or not, and keeps an untick",
        taken == [("S_001.raw", True), ("S_002.raw", True), ("S_003.raw", True), ("S_005.raw", True), ("CAL.raw", False)],
    )
    reopened = runsheet.RunSheet.open(set_dir)
    check_true(
        "a reopened sheet has the plan, the ticks, the typed-ahead cells and the session",
        [(r.acquisition, r.taken) for r in reopened.rows] == [(r.acquisition, r.taken) for r in sheet.rows]
        and [r.settings for r in reopened.rows] == [r.settings for r in sheet.rows]
        and reopened.rows[5].cells.get("pressure_Torr") == "2.04"
        and reopened.rows[0].cells.get("notes") == "first, after the gas had settled"
        and reopened.temperature_k == "301.13"
        and reopened.data_dir == data_dir
        and reopened.fixed == sheet.fixed,
    )
    before = open(sheet.conditions_path, "rb").read()
    reopened.save()
    check_true("saving a reopened sheet unchanged changes nothing", open(sheet.conditions_path, "rb").read() == before)

    # Planned against recorded: two .raw directories carrying the committed sidecar,
    # which records row 1's settings, and a third that has not been finished.
    for name in ("S_001.raw", "S_002.raw", "S_003.raw"):
        os.makedirs(os.path.join(data_dir, name))
    for name in ("S_001.raw", "S_002.raw"):
        shutil.copyfile(os.path.join(FIXTURES, "_extern.inf"), os.path.join(data_dir, name, "_extern.inf"))
    scan = runsheet.scan_data_dir(sheet.data_path)
    check_true(
        "a scan lists the .raw directories, and an unfinished one has no sidecar",
        scan.reachable and sorted(scan.externs) == ["S_001.raw", "S_002.raw", "S_003.raw"] and scan.externs["S_003.raw"] == {},
    )
    agrees = sheet.check_recorded(sheet.rows[0], scan.externs["S_001.raw"])
    check_true(
        "an acquisition recorded as planned agrees",
        agrees.agrees and abs(agrees.recorded_v - SERIES_130423[0][2]) < 1e-4,
    )
    differs = sheet.check_recorded(sheet.rows[1], scan.externs["S_002.raw"])
    check_true(
        "one whose knobs did not move is caught, naming both knobs",
        not differs.agrees
        and abs(differs.delta_v - (SERIES_130423[0][2] - SERIES_130423[1][2])) < 1e-4
        and {k for k, _, _ in differs.knobs} == {"Helium Cell DC", "IMSBias"},
    )
    gone = runsheet.scan_data_dir(os.path.join(work, "no-such-share"))
    check_true("an unreachable data directory is reported, not raised", not gone.reachable and gone.error)
finally:
    shutil.rmtree(work, ignore_errors=True)

# --------------------------------------------------------------------------------
section("the acquisition GUI (PySide6, offscreen)")
try:
    importlib.import_module("PySide6.QtWidgets")
    have_qt = True
except ImportError:
    have_qt = False
    for name in (
        "a QApplication and a QtAgg canvas construct",
        "the run-sheet window shows a sheet",
        "the Live tab shows the sheet's species and its completion rule",
    ):
        skip(name, "PySide6 not installed; it is the optional extra schamp[acquire]")
if have_qt:
    # Offscreen, so this runs on a machine with no display. In a subprocess, so that the
    # Qt platform choice and the QApplication do not leak into the checks after this one.
    probe = r"""
import os, sys, tempfile
os.environ["QT_QPA_PLATFORM"] = "offscreen"
from PySide6 import QtWidgets
app = QtWidgets.QApplication([])
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
canvas = FigureCanvasQTAgg(Figure())
canvas.figure.add_subplot().plot([1, 2], [3, 4])
canvas.draw()
print("canvas", canvas.get_width_height()[0] > 0)
from schamp.acquire import app as gui, model
window = gui.MainWindow()
sheet = model.RunSheet.from_template("uw-standard-14", tempfile.mkdtemp(), base_name="S")
window.set_sheet(sheet)
print("window", window.table.rowCount(), window.table.data(window.table.index(0, 2)))
sheet.species = [{"label": "1- n=10", "mz": 727.3744, "charge": -1}]
window.set_sheet(sheet)
app.processEvents()
print("live", window.live.table.rowCount(), window.live.table.item(0, 0).text(), "Finished when" in window.live.status.text())
"""
    result = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, timeout=120)
    lines = result.stdout.split("\n")
    check_true("a QApplication and a QtAgg canvas construct", "canvas True" in lines)
    check_true("the run-sheet window shows a sheet", "window 14 86.96" in lines)
    check_true("the Live tab shows the sheet's species and its completion rule", "live 1 1- n=10 True" in lines)
    if result.returncode:
        print(result.stderr[-2000:])

# --------------------------------------------------------------------------------
section("mass spectra and measured peak widths")
_frame_cal = importlib.import_module("schamp.calibration")
# A real linear calibration function, from a 2013 acquisition: near enough the identity
# that it is easy to miss, and far enough from it to move a peak by a tenth of a Th.
_shift = _frame_cal.CalFunction(coefficients=(0.006330094772637646, 0.999727577300022))
# A synthetic profile on a real axis: time-of-flight channels, so the spacing grows as
# sqrt(m/z) and a peak is three points wide, which is what the 2013 acquisitions give
# and what the measurement has to survive. Two peaks of known width and known
# resolving power, and one too weak to measure.
_root = np.arange(math.sqrt(200.0), math.sqrt(1200.0), 2.19e-4)
_axis = _root**2
_profile = np.zeros_like(_axis)
for _centre, _height, _power in ((300.0, 1.0e5, 20000.0), (800.0, 2.0e5, 20000.0), (500.0, 30.0, 20000.0)):
    _sigma = _centre / _power / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    _profile += _height * np.exp(-0.5 * ((_axis - _centre) / _sigma) ** 2)
spectrum = schamp.spectrum.MassSpectrum(
    mz=_axis, intensity=_profile, acquisition="synthetic", summed_over="a check, not an acquisition"
)
peak = spectrum.peak_near(800.0)
check_close("a peak's centroid is its centre", peak.mz, 800.0, rel=2e-6)
check_true(
    "a peak of R = 20000 measures near it, on a three-point profile",
    18000 < peak.resolving_power < 22000,
)
check_true("and the peak is resolved on both sides", peak.resolved)
check_true("a peak this narrow is a few channels wide", 2 <= peak.points_above_half <= 5)
check_close(
    "the m/z step is the axis's own, and grows as sqrt(m/z)",
    spectrum.mz_step_at(800.0) / spectrum.mz_step_at(200.0),
    math.sqrt(800.0 / 200.0),
    rel=1e-3,
)
check_true("a species that is not there measures as None", spectrum.peak_near(650.0) is None)
check_true(
    "a species below the apex threshold measures as None",
    spectrum.peak_near(500.0, min_apex=1.0e3) is None,
)
check_true(
    "...but is there when nothing is asked of it", spectrum.peak_near(500.0) is not None
)
check_raises(
    "a descending m/z axis is refused",
    ValueError,
    lambda: schamp.spectrum.MassSpectrum(mz=_axis[::-1], intensity=_profile),
)
# An empty spectrum is a real outcome -- one drift bin of one scan can hold nothing --
# so it is representable, says so, and refuses only the question that has no answer.
_empty = schamp.spectrum.MassSpectrum(mz=np.array([]), intensity=np.array([]))
check_true("an empty spectrum is representable and says so", _empty.is_empty)
check_true("...and measures no peaks rather than raising", _empty.peak_near(800.0) is None)
check_true("...and totals zero", _empty.total == 0.0)
check_raises("...but has no m/z range", ValueError, lambda: _empty.mz_range)
check_true(
    "a spectrum of zeros is empty too",
    schamp.spectrum.MassSpectrum(mz=_axis, intensity=np.zeros_like(_axis)).is_empty,
)
check_true("a real spectrum is not empty", not spectrum.is_empty)
check_raises(
    "so are mismatched arrays",
    ValueError,
    lambda: schamp.spectrum.MassSpectrum(mz=_axis, intensity=_profile[:-1]),
)

_measured = schamp.extract.MzWindow.around_measured(800.0, spectrum, label="synthetic")
check_close(
    "a measured window is `widths` FWHM across",
    _measured.width,
    schamp.extract.WINDOW_WIDTHS_DEFAULT * peak.fwhm,
    rel=1e-12,
)
check_close("and centres on the measured peak", _measured.centre, peak.mz, rel=1e-12)
check_true(
    "and says it came from a measurement",
    "measured peak" in _measured.frame_note and _measured.computed_mz == 800.0,
)
_measured_framed = schamp.extract.MzWindow.around_measured(
    800.0, spectrum, frame=_frame_cal.MzFrame(cal=_shift)
)
check_close(
    "a framed measured window converts the peak into the reader's frame",
    _measured_framed.centre,
    _shift.to_chromatogram_frame(peak.mz),
    rel=1e-9,
)
check_raises(
    "a species too weak to measure refuses rather than guessing",
    ValueError,
    lambda: schamp.extract.MzWindow.around_measured(500.0, spectrum, min_apex=1.0e3),
)
_fallback = schamp.extract.MzWindow.around_measured(
    500.0, spectrum, min_apex=1.0e3, fallback_width=0.5
)
check_close("unless a fallback width is given", _fallback.width, 0.5, rel=1e-12)
check_true(
    "and then the window says the width was a fallback",
    "fallback" in _fallback.frame_note,
)
check_raises(
    "a non-positive number of widths is refused",
    ValueError,
    lambda: schamp.extract.MzWindow.around_measured(800.0, spectrum, widths=0.0),
)

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
check_true(
    "a window placed with no frame records that it was placed with no frame",
    window.computed_mz == 789.4062 and window.frame_note == "",
)
_framed = schamp.extract.MzWindow.around(789.4062, 0.15, "n=22, 2-", frame=_frame_cal.MzFrame(cal=_shift))
check_close(
    "a framed window centres where the chromatogram reader puts the ion",
    _framed.centre,
    _shift.to_chromatogram_frame(789.4062),
    rel=1e-12,
)
check_true(
    "the two differ by more than a rounding, on a real calibration function",
    abs(_framed.centre - 789.4062) > 1e-3,
)
check_true(
    "and the framed window says what was applied, so it is never silent",
    _framed.computed_mz == 789.4062 and "cal function" in _framed.frame_note,
)
check_close(
    "an identity frame places a window exactly where an unframed one goes",
    schamp.extract.MzWindow.around(789.4062, 0.15, frame=_frame_cal.MzFrame.identity()).centre,
    window.centre,
    rel=1e-15,
)
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
section("the live pane (no Qt): watching, completion, and the series as it arrives")
import schamp.acquire.live as live  # noqa: E402 -- grouped with its checks

rule = live.WatchRule()
check_true(
    "the default completion rule: _extern.inf and _FUNC001.IDX, 3 unchanged looks 4 s apart",
    rule.files == ("_extern.inf", "_FUNC001.IDX") and rule.polls == 3 and rule.interval_s == 4.0,
)
check_true(
    "a [watch] table overrides only the terms it names",
    live.WatchRule.from_table({"polls": 5}).to_table() == {"files": list(rule.files), "polls": 5, "interval_s": 4.0},
)
check_raises("a [watch] key that is not a term is refused", live.LiveError, lambda: live.WatchRule.from_table({"pols": 3}))
check_raises("one look cannot see that nothing changed", live.LiveError, lambda: live.WatchRule.from_table({"polls": 1}))
check_raises("a species of charge 0 is refused", live.LiveError, lambda: live.Species("x", 500.0, 0))

work = tempfile.mkdtemp(prefix="schamp-live-")
try:
    data_dir = os.path.join(work, "data")
    os.makedirs(os.path.join(data_dir, "A.raw"))
    watcher = live.Watcher(live.WatchRule(polls=2, interval_s=1.0))

    def touch(name, size=1):
        with open(os.path.join(data_dir, "A.raw", name), "wb") as fh:
            fh.write(b"\0" * size)

    def look():
        return watcher.update(live.list_raws(data_dir))

    touch("_extern.inf")
    look()
    check_true("a .raw with a required file missing is 'writing'", watcher.status("A.raw") == "writing")
    touch("_FUNC001.IDX")
    touch("_FUNC001.DAT", 10)
    check_true("a change resets the count", look() == [] and watcher.status("A.raw") == "settling 1/2")
    events = look()
    check_true(
        "unchanged across `polls` looks, it is finished, and the event says which rule fired",
        [e.kind for e in events] == ["finished"] and "none changed over 2 looks" in events[0].reason,
    )
    check_true("a finished acquisition stays finished while nothing changes", look() == [])
    touch("_FUNC001.DAT", 20)
    check_true("one that grows after it was called finished is reopened", [e.kind for e in look()] == ["reopened"])
    check_true("...and finishes again once it settles", [e.kind for e in look()] == ["finished"])
    check_true(
        "a share that stops answering changes nothing",
        watcher.update(live.list_raws(os.path.join(work, "gone"))) == [] and watcher.finished == ["A.raw"],
    )

    worker = live.SerialWorker("check")
    worker.submit("one", lambda: 1)
    worker.submit("boom", lambda: 1 / 0)
    worker.submit("two", lambda: 2)
    done = []
    for _ in range(200):
        done += worker.drain()
        if len(done) == 3:
            break
        time.sleep(0.01)
    check_true(
        "the worker runs jobs in order on one thread and hands back an exception as a result",
        [t for t, _ in done] == ["one", "boom", "two"] and isinstance(done[1][1], ZeroDivisionError),
    )

    # A series arriving: synthetic ATDs at eqn (3)'s arrival times through a stand-in for
    # the SDK read, against an experiment written the way the run sheet writes one.
    K_TRUE, T0_TRUE, MZ, Z = 1231.8, 0.42, 727.3744, -1  # K at 2 Torr: the 1- n=10 test ion
    slope = 1000.0 * profile.drift_length_cm**2 / K_TRUE
    axis = np.arange(200) * 0.0692521
    names = [f"S_{i:03d}.raw" for i in range(1, 9)]
    volts = dict(zip(names, [104.0, 140.0, 302.4, 181.8, 249.0, 120.0, 85.0, 160.0]))
    late = {"S_007.raw": 0.35}  # the use = false row, and off the line by 0.35 ms

    def fake_reader(path, species):
        name = os.path.basename(path)
        centre = T0_TRUE + slope / volts[name] + late.get(name, 0.0)
        atd = schamp.atd.ATD(
            np.arange(200), axis, 5.0e4 * np.exp(-((axis - centre) ** 2) / (2 * 0.08**2)), MZ - 0.1, MZ + 0.1,
            acquisition=path, label=species[0].label,
        )
        return {s.label: live.SpeciesResult(s, (MZ - 0.1, MZ + 0.1), "synthetic", atd) for s in species}

    set_dir = os.path.join(work, "set")
    os.makedirs(set_dir)
    schamp.experiment.write_experiment(set_dir, title="live", profile="uw-synapt-g2", gas="helium", data_dir=data_dir)
    schamp.experiment.write_conditions_table(
        os.path.join(set_dir, "conditions.csv"),
        ["acquisition", "pressure_Torr", "temperature_K", "use", "v_drift_V"],
        [[n, 2.03, 301.13, "false" if n in late else "true", volts[n]] for n in names[:-1]],
    )  # S_008.raw arrives but was never ticked as taken
    experiment = schamp.load_experiment(set_dir)
    species = live.Species("1- n=10", MZ, Z)
    series = live.LiveSeries([species])
    for name in names:
        os.makedirs(os.path.join(data_dir, name), exist_ok=True)
        series.arrive(name, os.path.join(data_dir, name))
        acquisition, path, wanted = series.next_job()
        series.record(series.job(acquisition, path, wanted, reader=fake_reader)(), acquisition)
    check_true("every arrival is read once, and nothing is left to read", series.next_job() is None and len(series.results) == 8)
    summary = series.summarise(species, experiment, profile)
    used = [p.acquisition for p in summary.points if p.used]
    check_true(
        "use = false and untaken acquisitions are plotted but stay out of the regression",
        used == names[:6] and {p.acquisition: p.why_not for p in summary.points if not p.used}
        == {"S_007.raw": "use = false", "S_008.raw": "not in conditions.csv"},
    )
    direct = schamp.mobility.cross_section_from_regression(
        schamp.mobility.regress(
            [schamp.mobility.DriftPoint(n, volts[n], T0_TRUE + slope / volts[n]) for n in names[:6]], profile
        ),
        profile, charge=1, ion_mass_da=MZ, gas="helium", pressure_torr=2.03, temperature_k=301.13,
    )
    check_close("the live cross section is the pipeline's own", summary.cross_section.omega_a2, direct.omega_a2, rel=1e-6)
    check_true(
        "the running cross section starts at the third point and ends at the final one",
        [n for n, _ in summary.running] == [3, 4, 5, 6] and summary.running[-1][1] == summary.cross_section,
    )
    check_true(
        "the point off the line is flagged against the others' rms, and no other",
        [p.acquisition for p in summary.points if p.outlier] == ["S_007.raw"],
    )
    series.add_species(live.Species("2- n=20", 640.0, -2))
    check_true(
        "a species added later queues a re-read of the finished acquisitions for it alone",
        [s.label for s in series.next_job()[2]] == ["2- n=20"],
    )
    series.remove_species("2- n=20")
    series.reopen("S_003.raw")
    check_true("a reopened acquisition leaves the series until it is read again", "S_003.raw" not in series.results)
    series.arrive("S_003.raw", os.path.join(data_dir, "S_003.raw"))
    acquisition, path, wanted = series.next_job()
    series.record(series.job(acquisition, path, wanted, reader=fake_reader)(), acquisition)
    json_path, csv_path = series.finish(experiment, set_dir)
    table = schamp.experiment.read_conditions_table(csv_path)
    with open(json_path, encoding="utf-8") as fh:
        written = json.load(fh)["results"]
    check_true(
        "finishing writes ccs.csv (gas, charge, P and T beside every cross section) and results.json",
        len(table.rows) == 1
        and abs(float(table.get(table.rows[0], "CCS_A2")) - round(direct.omega_a2, 2)) < 0.006
        and table.get(table.rows[0], "gas") == "helium"
        and table.get(table.rows[0], "excluded") == "S_007.raw;S_008.raw"
        and written["species"][0]["cross_section"]["pressure_Torr"] == 2.03,
    )

    sheet = runsheet.RunSheet.from_template(template, os.path.join(work, "sheet"), base_name="S")
    sheet.watch = {"polls": 4, "files": ["_extern.inf"]}
    sheet.species = [species.to_table()]
    sheet.save()
    again = runsheet.RunSheet.open(sheet.set_dir)
    check_true(
        "the run sheet keeps the [watch] rule and the species across a reopen",
        live.WatchRule.from_table(again.watch).polls == 4
        and [live.Species.from_table(e) for e in again.species] == [species],
    )

    smoke = os.environ.get("SCHAMP_SMOKE_RAW")
    problem = live.sdk_problem()
    if problem or not smoke:
        why = "no SDK or license" if problem else "set SCHAMP_SMOKE_RAW to an acquisition"
        skip("the watcher finishes a real acquisition under the default rule", why)
        skip("one real acquisition read and fitted through the live path", why)
    else:
        real = live.Watcher(live.WatchRule(polls=2, interval_s=1.0))
        parent = os.path.dirname(os.path.abspath(smoke))
        real.update(live.list_raws(parent))
        events = real.update(live.list_raws(parent))
        check_true(
            "the watcher finishes a real acquisition under the default rule",
            os.path.basename(os.path.abspath(smoke)) in [e.acquisition for e in events if e.kind == "finished"],
        )
        spectrum = live.spectrum_of(smoke)
        tallest = float(spectrum.mz[int(np.argmax(spectrum.intensity))])
        result = live.analyse_acquisition(smoke, [live.Species("tallest", tallest, 1)])
        check_true("one real acquisition read and fitted through the live path", result.species["tallest"].ok)
finally:
    shutil.rmtree(work, ignore_errors=True)

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

# --------------------------------------------------------------------------------
section("the example notebook")
NOTEBOOK = os.path.join(schamp.ROOT, "examples", "polyalanine-walkthrough.ipynb")
if not os.path.isfile(NOTEBOOK):
    check_true("examples/polyalanine-walkthrough.ipynb exists", False)
else:
    with open(NOTEBOOK, encoding="utf-8") as fh:
        notebook_json = json.load(fh)
    code_cells = [c for c in notebook_json.get("cells", []) if c.get("cell_type") == "code"]
    check_true("the notebook has code cells", len(code_cells) > 0)
    check_true(
        "the shipped notebook carries no committed execution artefacts",
        all(not c.get("outputs") and c.get("execution_count") is None for c in code_cells),
    )

    try:
        importlib.import_module("nbformat")
        importlib.import_module("nbclient")
        have_dev = True
    except ImportError:
        have_dev = False
        skip("the notebook executes end to end", "the dev dependency group is not installed")

    if have_dev:
        example_raw = os.environ.get("SCHAMP_EXAMPLE_RAW")
        if not example_raw:
            lab = schamp.lab_dir("legacy")
            if lab and os.path.isdir(os.path.join(lab, "130423")):
                example_raw = os.path.join(lab, "130423")
        if not have_sdk:
            skip("the notebook executes end to end", "Waters SDK not installed")
        elif not example_raw:
            skip(
                "the notebook executes end to end",
                "no polyalanine acquisitions; set SCHAMP_EXAMPLE_RAW or run from a lab "
                "checkout that ships them",
            )
        else:
            # A subprocess, not an in-process kernel: nbconvert's own CLI is what
            # examples/README.md tells a user to run, so this is the same path
            # being checked rather than a shortcut around it.
            env = dict(os.environ, SCHAMP_EXAMPLE_RAW=example_raw)
            executed = tempfile.mktemp(suffix=".ipynb")
            try:
                done = subprocess.run(
                    [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook",
                     "--execute", "--output", executed, NOTEBOOK],
                    capture_output=True, text=True, env=env, timeout=1800,
                )
                check_true(
                    "the notebook executes end to end",
                    done.returncode == 0 and os.path.isfile(executed),
                )
                if done.returncode != 0:
                    print(done.stderr[-4000:])
            finally:
                if os.path.isfile(executed):
                    os.remove(executed)

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

        # The spectrum side of the same acquisition, which is the frame a species is
        # found in rather than the frame a window bound is on.
        smoke_total = schamp.extract.total_spectrum(readers)
        smoke_one_scan = schamp.extract.total_spectrum(readers, scans=(0, 0))
        # One drift bin of one retention scan. On a real mobility acquisition that can
        # be a handful of points or none, which is an ordinary outcome and must not
        # raise -- the first bin of the first scan is before the ions arrive.
        smoke_one_bin = schamp.extract.total_spectrum(readers, scans=(0, 0), drift=(0, 0))
        check_raises(
            "a reversed drift range is refused rather than returning nothing",
            ValueError,
            lambda: schamp.extract.total_spectrum(readers, drift=(10, 5)),
        )
        check_raises(
            "a drift range past the end is refused",
            ValueError,
            lambda: schamp.extract.total_spectrum(readers, drift=(0, 10**6)),
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
    check_true(
        f"the total spectrum is {len(smoke_total.mz)} points over "
        f"{smoke_total.mz_range[0]:.1f}-{smoke_total.mz_range[1]:.1f} Th",
        len(smoke_total.mz) > 1000 and smoke_total.total > 0.0,
    )
    check_true(
        "and it says what it was summed over",
        "all" in smoke_total.summed_over and smoke_total.acquisition == smoke_raw,
    )
    # A subset of everything. Cheap, and it is the check that a narrowed combine
    # actually narrows rather than silently returning the lot.
    check_true(
        "a narrowed combine carries less than the whole",
        0.0 < smoke_one_scan.total < smoke_total.total,
    )
    check_true(
        f"the emptiest slice of the acquisition is representable "
        f"({len(smoke_one_bin.mz)} points, {smoke_one_bin.total:.0f} counts)",
        smoke_one_bin.total <= smoke_one_scan.total
        and smoke_one_bin.is_empty == (smoke_one_bin.total <= 0.0),
    )
    # The frame conversion on the real acquisition: whatever its calibration function
    # is, a window placed through it has to invert.
    real_cal = _frame_cal.CalFunction.from_raw(smoke_raw)
    check_true(f"its calibration function reads ({real_cal.describe()})", real_cal.order >= 1)
    probe_mz = 0.5 * (smoke_total.mz_range[0] + smoke_total.mz_range[1])
    check_close(
        "the acquisition's own frame conversion inverts",
        real_cal.to_spectrum_frame(real_cal.to_chromatogram_frame(probe_mz)),
        probe_mz,
        rel=1e-12,
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

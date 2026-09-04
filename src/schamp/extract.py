"""Extraction: a `.raw` plus m/z windows, to arrival-time distributions.

The only layer that needs the Waters SDK, and the only one that cannot run in a public
clone. It is thin on purpose. Mobillograms, drift-scan reads, scan combining and the
drift-time axis are the SDK's own and are called directly; what is added here is the
window bookkeeping, the summing convention, and packing the result into `atd.ATD`.
Nothing the SDK does faithfully is reimplemented.

This replaces the 2013 workflow's `cdctest.exe`, which took a CSV of m/z windows and a
`.raw` and wrote one summed ATD per window. The equivalent read is
`ReadMobillogram(function, first_scan, last_scan, low_mass, high_mass)` over every
retention scan of the function, and given the same windows it reproduces that tool bin
for bin on 1457 of the 1473 windows of the 2013 acquisitions -- every exception being a
window bound that lands on signal, which the next section is about (lab record,
task 04).

`masslynxsdk` is imported inside the functions that need it, so this module imports on
a machine with no SDK and only the reads fail.

Vendor conventions pass through untranslated: functions, retention scans and drift bins
are 0-indexed, and a scan range is inclusive at both ends.

### A window bound is snapped to the acquisition's own m/z steps

Measured, not assumed (lab record, task 04). The reader does not compare each data
point's m/z against the bounds it was given: it snaps each bound to the nearest step of
the acquisition's m/z axis, which is the time-of-flight channel spacing and therefore
grows as the square root of m/z -- about 0.012 Th at 800 Th on the 2013 acquisitions.
So a window's bounds are honoured to one step and no better. Two things follow, and the
first is easy to trip over:

- **A window narrower than one step is not narrower than one step.** Asked for a
  0.001 Th window it returns a whole step's counts; asked for a 1e-6 Th one it returns
  nothing. Widths below one step do not mean what they say.
- **Where a bound falls near a step edge, two tools can snap it either way.** That is
  the whole of the disagreement with `cdctest.exe`: given identical bounds, the two
  agree bin for bin on 1457 of the 1473 windows of the 2013 acquisitions, and every one
  of the 16 exceptions is a bound whose step carries signal, where the SDK takes the
  step in and cdctest leaves it out. The cost is one step's counts -- 0.04 to 0.12 % of
  the window -- and the peak bin is unchanged in all 1473.

Both are avoided the same way: **put a bound in a valley, not on a peak's flank.** A
step that holds nothing cannot be won or lost. `window_edges_on_data` is the check --
it reports, per bound, the counts in the step the bound lands in, so an empty list means
every bound of every window is in clear air.

A grid of touching windows nearly always partitions the counts: of the 254 shared bounds
of the 2013 heat map's grid, five land in a step that carries signal, and at four of
those the step goes to exactly one of the two neighbours. At the fifth it goes to both.
So a grid is a partition with rare exceptions rather than by construction, and the
exceptions are the same bounds `window_edges_on_data` names.

### The mass range is not the calibrated m/z axis

Third measured property, and the one to keep in mind when choosing a width. The mass
range this reader is given is not the calibrated m/z axis. On the 2013 acquisitions a
peak of known m/z appears low at low mass and high at high mass -- from about -0.08 Th at
300 Th to +0.16 Th at 1156 Th, crossing zero near 540 Th, which is the shape of a
calibration difference rather than a fixed shift. The same peak read through the
drift-scan reader sits within 1 ppm of its computed m/z, and the two readers agree
exactly on total counts over the whole mass range, so nothing is lost; it is redistributed
in m/z.

The 0.2 Th windows the 2013 analysis typed by hand sit clear of the peak at the ends of
its mass range for that reason. It cost that analysis nothing, because every isotopologue
of an ion has the same mobility and so sampling a peak's shoulder costs counts and not
drift time -- but a window narrow enough to clip is a window whose contents depend on a
calibration offset nobody wrote down.

**Give a window room.** One Th is generous at these masses and costs nothing when the
neighbouring species are 35 Th away. Where crowding does not allow it, measure where the
peak actually lands before trusting a narrow window (lab record, task 04, which measures
it for one instrument and leaves to the layer that places analysis windows the question
of whether to correct for it).

### Intensities are single precision

Also measured (lab record, task 04). The intensities come back float32-exact, so counts
in one drift bin are integers only up to `COUNTS_EXACT_TO`, and the SDK's accumulation
into that precision loses a little on the way: summing a grid of eight windows over a
whole acquisition mass range exceeds one read over the same span by 1e-5 to 5e-5 of the
total, entirely from rounding rather than from any shared bound. Nothing here corrects
it, because the loss is inside the SDK's own accumulation.

It does not touch the numbers this package exists to produce -- an analysis window holds
a single oligomer and peaks four orders of magnitude below `COUNTS_EXACT_TO`, and a
centroid is a ratio -- but it does mean that a total over a wide window is good to about
five figures and not more, and that summing sub-windows is not a way to check a wide
read.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass

import numpy as np

from .atd import ATD
from .sdk import Readers

__all__ = [
    "COUNTS_EXACT_TO",
    "MZ_BOUND_PROBE_RELATIVE",
    "MZ_WINDOW_FLOOR_RELATIVE",
    "MzWindow",
    "contiguous_windows",
    "drift_axis",
    "extract_atd",
    "extract_atds",
    "mobility_map",
    "read_peaks_csv",
    "scan_range",
    "window_edges_on_data",
    "write_peaks_csv",
]

# The half-width, relative to m/z, of the nominally zero-width probe that
# `window_edges_on_data` puts on a bound. It is deliberately far narrower than any
# instrument's m/z step: the reader widens a request to the step it lands in, so a probe
# this narrow reads exactly the counts of the one step containing the bound, which is
# the quantity that decides whether two tools can disagree there. Not a tolerance, and
# never something to subtract off a bound.
MZ_BOUND_PROBE_RELATIVE = 2.0e-7

# The narrowest window worth asking for, relative to m/z. An absurdity floor rather than
# an instrument constant: it is more than an order of magnitude below the m/z step of any
# time-of-flight acquisition this package is for, so a window narrower than this cannot
# mean what it says under any calibration. Choosing a width that suits the peak is the
# caller's business, and the instrument's own step is the floor that actually matters --
# `window_edges_on_data` is how to find out where it falls.
MZ_WINDOW_FLOOR_RELATIVE = 1.0e-6

# The largest count in one drift bin that the SDK's single-precision intensities still
# represent exactly, 2**24. Above it the returned value is quantised in steps of 2, then
# 4, and so on. The 2013 heat map's brightest bin reached 6.4e6, so nothing in that
# analysis is near the limit; a read over a whole acquisition mass range is well past
# it. See this module's docstring.
COUNTS_EXACT_TO = 2**24


@dataclass(frozen=True)
class MzWindow:
    """One m/z window, and what it is meant to contain.

    `low` and `high` are in Th, and each is honoured to one m/z step of the acquisition
    rather than exactly; see this module's docstring. `label` names the species --
    `"n=22, 2-"` -- because a
    window is placed around an expected m/z and the label is the only record of what
    was expected. `charge` and `ion_mass_da` are carried through to the mobility layer,
    which cannot compute a cross section without them.

    Windows must be built from the m/z series of the species under study, not read from
    a legacy peaks file: the surviving 2013 file has three of its nineteen rows
    mistyped -- two inverted pairs, which return all zeros, and one 300 Th wide instead
    of 1 Th, which returns a blend of several oligomers -- and is demonstrably not the
    file the published analysis used. `MzWindow.around` is the constructor for that.
    """

    low: float
    high: float
    label: str = ""
    charge: int | None = None
    ion_mass_da: float | None = None

    def __post_init__(self) -> None:
        if not (self.high > self.low):
            raise ValueError(
                f"m/z window {self.label or '(unlabelled)'} has low={self.low}, "
                f"high={self.high}: the bounds are inverted or zero-width. cdctest "
                "echoed such a window back as a column of zeros; schamp refuses it."
            )
        if self.low <= 0.0:
            raise ValueError(f"m/z window {self.label or ''} has a non-positive low bound")

    @classmethod
    def around(
        cls,
        mz: float,
        width: float,
        label: str = "",
        charge: int | None = None,
        ion_mass_da: float | None = None,
    ) -> MzWindow:
        """A window of full width `width` Th centred on `mz`.

        The way to build a window: from a species' expected m/z and a width chosen for
        its peak, so that both are visible in the calling code rather than typed as a
        pair of bounds. Choose the width for the peak and for where its bounds land --
        this module's docstring says why a bound in a valley is worth more than a narrow
        window, and `window_edges_on_data` says where the bounds actually fell. Only an
        absurd width is refused here, `MZ_WINDOW_FLOOR_RELATIVE`.
        """
        if width <= 0.0:
            raise ValueError(f"a window width must be positive; got {width}")
        floor = abs(float(mz)) * MZ_WINDOW_FLOOR_RELATIVE
        if width < floor:
            raise ValueError(
                f"a {width:g} Th window at {mz:g} Th is orders of magnitude below the "
                f"m/z step of any time-of-flight acquisition ({floor:g} Th is already "
                "far below it); the reader would widen it to one step and the width "
                "would mean nothing"
            )
        return cls(
            low=float(mz) - 0.5 * float(width),
            high=float(mz) + 0.5 * float(width),
            label=label,
            charge=charge,
            ion_mass_da=ion_mass_da,
        )

    @property
    def centre(self) -> float:
        """The window's centre in Th."""
        return 0.5 * (self.low + self.high)

    @property
    def width(self) -> float:
        """The window's full width in Th."""
        return self.high - self.low


def contiguous_windows(
    low: float, high: float, count: int, *, label: str = ""
) -> list[MzWindow]:
    """`count` touching windows tiling `low` to `high`, lowest first.

    Adjacent windows share a bound exactly, which is what the 2013 heat map's
    255-window grid did over 200 to 1700 Th, and it is how that grid is reproduced.
    Read this module's docstring before summing such a grid and expecting the total of
    a single wide read: a point sitting on a shared bound is counted in both
    neighbours, and `window_edges_on_data` says which bounds those are.
    """
    if count < 1:
        raise ValueError(f"a grid needs at least one window; got {count}")
    if not (high > low):
        raise ValueError(f"a grid needs high > low; got {low} to {high}")
    edges = np.linspace(float(low), float(high), count + 1)
    return [
        MzWindow(
            low=float(edges[i]),
            high=float(edges[i + 1]),
            label=f"{label} {i}".strip() if label else "",
        )
        for i in range(count)
    ]


def read_peaks_csv(path: str | os.PathLike[str]) -> list[MzWindow]:
    """The m/z windows in a legacy `cdctest.exe` peaks file.

    The format, for anyone migrating off that tool: the first field of the first row is
    the number of windows, and each row after it is `low,high` in Th. cdctest read
    exactly the declared count and ignored the rest of the file; so does this.

    An inverted or zero-width pair raises rather than being read, naming the row.
    cdctest echoed such a row back as a column of zeros, and three of the nineteen rows
    of the one 2013 anion peaks file that survives are wrong that way. Windows worth
    trusting are built with `MzWindow.around` from an m/z series; this reader exists to
    read what a lab already has, and to fail loudly on the rows a silent tool did not.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [row for row in csv.reader(fh) if any(field.strip() for field in row)]
    if not rows:
        raise ValueError(f"{path} is empty; a peaks file starts with a window count")
    try:
        declared = int(float(rows[0][0]))
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"{path} does not start with a window count; its first row is {rows[0]!r}"
        ) from exc
    body = rows[1:]
    if len(body) < declared:
        raise ValueError(
            f"{path} declares {declared} windows and carries {len(body)} rows after the "
            "count"
        )
    windows = []
    for number, row in enumerate(body[:declared], start=1):
        try:
            low, high = float(row[0]), float(row[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(
                f"{path} row {number + 1} is not a low,high pair: {row!r}"
            ) from exc
        try:
            windows.append(MzWindow(low=low, high=high, label=f"peaks row {number}"))
        except ValueError as exc:
            raise ValueError(f"{path} row {number + 1}: {exc}") from exc
    return windows


def write_peaks_csv(path: str | os.PathLike[str], windows: list[MzWindow]) -> None:
    """Write windows in the legacy `cdctest.exe` peaks format.

    The counterpart of `read_peaks_csv`, and the reason it is here: it is what lets a
    lab put one set of windows through both tools on its own acquisitions and compare
    the two extractions, which is how this module was validated. Full precision is
    written; cdctest holds a bound to single precision and echoes it back rounded to
    four decimals, which is one of the two reasons the tools can disagree at a shared
    bound.
    """
    if not windows:
        raise ValueError("a peaks file with no windows would extract nothing")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([len(windows), ""])
        for window in windows:
            writer.writerow([repr(window.low), repr(window.high)])


def scan_range(
    readers: Readers, function: int = 0, scans: tuple[int, int] | None = None
) -> tuple[int, int]:
    """The inclusive, 0-indexed retention-scan range a read will cover.

    `None` means every scan of the function, which is what `cdctest.exe` did
    unconditionally and what reproducing the 2013 numbers therefore requires.

    Validated here rather than left to the SDK, which raises on a scan number past the
    end but accepts a reversed or negative range and returns a column of zeros for it.
    A silently empty ATD is exactly the class of failure this project is replacing.
    """
    count = int(readers.info.GetScansInFunction(function))
    if scans is None:
        return 0, count - 1
    first, last = int(scans[0]), int(scans[1])
    if first < 0 or last < 0:
        raise ValueError(f"retention scans are 0-indexed and non-negative; got {scans}")
    if last < first:
        raise ValueError(
            f"the scan range {scans} is reversed; the SDK returns a column of zeros for "
            "one of those rather than refusing it"
        )
    if last >= count:
        raise ValueError(
            f"the scan range {scans} runs past the end of function {function}, which has "
            f"{count} scans (0 to {count - 1}; MassLynx displays 1 to {count})"
        )
    return first, last


def drift_axis(readers: Readers, function: int = 0) -> np.ndarray:
    """The drift-time axis of one function, in milliseconds, from the SDK.

    This is `GetDriftTime` over every drift bin, and it is the only source of the drift
    time that schamp uses. It is not a bin index times a pusher period: `_extern.inf`
    carries no pusher-interval key on this instrument, and the drift-bin cache that
    does carry one is written by MassLynx only on copies of an acquisition that have
    been opened in it.

    Functions and drift bins are 0-indexed, as the SDK numbers them; MassLynx displays
    both from 1.
    """
    count = int(readers.info.GetDriftScanCount(function))
    return np.array(
        [float(readers.info.GetDriftTime(bin_index)) for bin_index in range(count)],
        dtype=float,
    )


def _read_window(
    readers: Readers,
    window: MzWindow,
    function: int,
    first: int,
    last: int,
    axis: np.ndarray,
) -> ATD:
    """One `ReadMobillogram` call packed into an `ATD`. The whole of the extraction."""
    bins, intensities = readers.chrom.ReadMobillogram(
        int(function), int(first), int(last), float(window.low), float(window.high)
    )
    drift_bin = np.asarray(bins, dtype=int)
    intensity = np.asarray(intensities, dtype=float)
    if len(drift_bin) != len(axis):
        raise ValueError(
            f"the mobillogram of function {function} has {len(drift_bin)} bins and its "
            f"drift-time axis has {len(axis)}; the acquisition is not what it says it is"
        )
    return ATD(
        drift_bin=drift_bin,
        drift_time_ms=axis,
        intensity=intensity,
        mz_low=window.low,
        mz_high=window.high,
        acquisition=readers.path,
        label=window.label,
    )


def extract_atd(
    readers: Readers,
    window: MzWindow,
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> ATD:
    """One window's arrival-time distribution, summed over retention-time scans.

    `scans` is an inclusive `(first, last)` pair of 0-indexed retention scans, or None
    for the whole function -- which is what `cdctest.exe` did unconditionally, and what
    reproducing the legacy numbers therefore requires.

    One SDK call per window, and the summing over retention scans is the SDK's own. Use
    `extract_atds` for more than one window on the same acquisition: it reads the
    drift-time axis once instead of once per window.
    """
    first, last = scan_range(readers, function, scans)
    return _read_window(readers, window, function, first, last, drift_axis(readers, function))


def extract_atds(
    readers: Readers,
    windows: list[MzWindow],
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> list[ATD]:
    """`extract_atd` over a list of windows, in order.

    Separate from `extract_atd` because the native library is not assumed thread-safe
    and a series of windows on one open acquisition is the case worth optimising; a
    caller should never be looping over `open_readers`. The scan range and the
    drift-time axis are resolved once and shared, and the returned ATDs therefore share
    one axis array -- safe because `ATD` is frozen and nothing downstream writes to it.
    """
    if not windows:
        return []
    first, last = scan_range(readers, function, scans)
    axis = drift_axis(readers, function)
    return [_read_window(readers, w, function, first, last, axis) for w in windows]


def mobility_map(
    readers: Readers,
    windows: list[MzWindow],
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`(mz_centres, drift_time_ms, intensity)` over a grid of windows.

    `intensity` is shaped `(len(windows), n_bins)`. This is the drift-time-against-m/z
    heat map of the published figure, which the 2013 workflow built by running its
    extraction over 255 contiguous 5.882 Th windows tiling 200 to 1700 Th; the modern
    one is the same read over the same windows, from `contiguous_windows`.

    No normalisation, no interpolation and no smoothing: the counts as read, so that
    the figure's square-root intensity scaling stays a choice made when plotting.
    """
    if not windows:
        raise ValueError("a mobility map needs at least one window")
    atds = extract_atds(readers, windows, function=function, scans=scans)
    return (
        np.array([atd.mz_centre for atd in atds], dtype=float),
        atds[0].drift_time_ms,
        np.vstack([atd.intensity for atd in atds]),
    )


def window_edges_on_data(
    readers: Readers,
    windows: list[MzWindow],
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> list[dict[str, float]]:
    """Which window bounds land in an m/z step that carries counts.

    One entry per such bound, carrying the window's index in `windows`, which bound it
    is, the bound's m/z, the half-width of the probe, and the counts the step holds. An
    empty list means every bound of every window sits in clear air, which is the state
    an analysis wants: a step holding nothing cannot be claimed by one tool and not
    another, and cannot be double counted across a shared bound.

    One read per bound, of a probe `MZ_BOUND_PROBE_RELATIVE` wide, which the reader
    widens to the step the bound falls in -- so this measures the step's contents rather
    than applying a rule. A non-empty result is a warning and not a verdict: at four of
    the five bounds it flags on the 2013 heat-map grid, the SDK and `cdctest.exe` snap
    the same way and agree exactly anyway. It is the check the 2013 workflow could not
    make; the one bound of that grid where they differ is invisible in its output and
    cost the window 0.12 % of its counts.
    """
    first, last = scan_range(readers, function, scans)
    found: list[dict[str, float]] = []
    for index, window in enumerate(windows):
        for side, bound in (("low", window.low), ("high", window.high)):
            band = abs(bound) * MZ_BOUND_PROBE_RELATIVE
            _, intensities = readers.chrom.ReadMobillogram(
                int(function), first, last, float(bound - band), float(bound + band)
            )
            counts = float(np.sum(intensities))
            if counts > 0.0:
                found.append(
                    {
                        "window": index,
                        "side": side,
                        "mz": bound,
                        "band_th": band,
                        "counts": counts,
                    }
                )
    return found

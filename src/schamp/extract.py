"""Extraction: a `.raw` plus m/z windows, to arrival-time distributions.

The only layer that needs the Waters SDK, and the only one that cannot run in a public
clone. It is thin on purpose. Mobillograms, drift-scan reads, scan combining and the
drift-time axis are the SDK's own and are called directly; what is added here is the
window bookkeeping, the summing convention, and packing the result into `atd.ATD`.
Nothing the SDK does faithfully is reimplemented.

This replaces the 2013 workflow's `cdctest.exe`, which took a CSV of m/z windows and a
`.raw` and wrote one summed ATD per window. Task 02 established that
`ReadMobillogram(function, first_scan, last_scan, low_mass, high_mass)` reproduces it:
on the 255-window grid of the published heat map, summed over all retention scans, 253
of 255 windows agree bin for bin and all 255 peak in the same bin. The two that differ
are the pair either side of one boundary, where the SDK counts the same ions in both
windows -- so the open question is the boundary convention, not the read, and it
belongs to task 04.

`masslynxsdk` is imported inside the functions that need it, so this module imports on
a machine with no SDK and only the reads fail.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .atd import ATD
from .sdk import Readers

__all__ = ["MzWindow", "drift_axis", "extract_atd", "extract_atds", "mobility_map"]


@dataclass(frozen=True)
class MzWindow:
    """One m/z window, and what it is meant to contain.

    `low` and `high` are in Th. `label` names the species -- `"n=22, 2-"` -- because a
    window is placed around an expected m/z and the label is the only record of what
    was expected. `charge` and `ion_mass_da` are carried through to the mobility layer,
    which cannot compute a cross section without them.

    Windows must be built from the m/z series of the species under study, not read from
    a legacy peaks file: the surviving 2013 file has three of its nineteen rows
    mistyped -- two inverted pairs, which return all zeros, and one 300 Th wide instead
    of 1 Th, which returns a blend of several oligomers -- and is demonstrably not the
    file the published analysis used.
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

    @property
    def centre(self) -> float:
        """The window's centre in Th."""
        return 0.5 * (self.low + self.high)


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


def extract_atd(
    readers: Readers,
    window: MzWindow,
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> ATD:
    """One window's arrival-time distribution, summed over retention-time scans.

    Not implemented: task 04, which owns the window-boundary convention and validates
    the output bin for bin against the legacy `cdctest.exe` CSVs. The signature is
    fixed now.

    `scans` is an inclusive `(first, last)` pair of 0-indexed retention scans, or None
    for the whole function -- which is what `cdctest.exe` did unconditionally, and what
    reproducing the legacy numbers therefore requires.
    """
    raise NotImplementedError("ATD extraction through the SDK is task 04")


def extract_atds(
    readers: Readers,
    windows: list[MzWindow],
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> list[ATD]:
    """`extract_atd` over a list of windows, in order. Not implemented: task 04.

    Separate from `extract_atd` because the native library is not assumed thread-safe
    and a series of windows on one open acquisition is the case worth optimising; a
    caller should never be looping over `open_readers`.
    """
    raise NotImplementedError("ATD extraction through the SDK is task 04")


def mobility_map(
    readers: Readers,
    windows: list[MzWindow],
    *,
    function: int = 0,
    scans: tuple[int, int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`(mz_centres, drift_time_ms, intensity)` over a contiguous grid of windows.

    Not implemented: task 04. `intensity` is shaped `(len(windows), n_bins)`. This is
    the drift-time-against-m/z heat map of the published figure, which the 2013
    workflow built by running its extraction over 255 contiguous 5.882 Th windows
    tiling 200 to 1700 Th; the modern one is the same read with the same windows.
    """
    raise NotImplementedError("the mobility map is task 04")

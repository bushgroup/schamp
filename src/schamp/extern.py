"""The two sidecars schamp reads for itself: `_extern.inf` and one line of `_HEADER.TXT`.

Everything else about a `.raw` comes from the Waters SDK's own readers. Two things
are not there, and each is here because it is not there:

* **The acquisition parameters.** The SDK exposes no accessor for the DC optic
  voltages, so the drift voltage of a linear-field measurement is unreadable through
  it (decision of record).
* **`$$ Cal Function N` in `_HEADER.TXT`**, the per-function mass-calibration
  polynomial. `MassLynxHeaderItem` carries `CAL_MS1_STATIC`, `CAL_MS1_DYNAMIC_PARAMS`
  and the rest, but has no member for this line, and without it there is no way to
  convert between the m/z axis the spectrum readers return and the frame
  `ReadMobillogram` is addressed in -- `calibration` is what that conversion is for
  (decision of record, lab record task 10).

`parse_header_txt` returns the whole file because a one-key parser is a worse thing
to maintain than a flat one, but **read a header item through the SDK if the SDK has
it**: `GetHeaderItemValue` is the supported accessor and this is a fallback for the
one key it lacks. A *third* non-SDK read needs its own decision; two are not a
precedent either.

Contract: a `.raw` directory in, a flat `dict[str, str]` out, values verbatim apart
from surrounding whitespace. **It never raises.** A missing directory, a missing
file, an empty one, junk, or an encoding surprise all return `{}`, and the caller
leaves whatever it wanted unset. Parsing this sidecar must never be what stops an
acquisition from opening.

The file is undocumented and vendor-controlled, so key by name and never by line
number, and never assume a key exists. Two facts about it that cost time to learn
(lab record, task 01):

* **There is no pusher-interval key.** `Pusher` is the pusher pulse voltage in volts
  and `Pusher Offset` is a voltage too. The 2013 analysis grepped for a
  `Pusher Interval` that this instrument has never written, silently fell back to
  1.0, and so produced every one of its centroids in drift bins rather than
  milliseconds. The drift-time axis comes from the SDK's `GetDriftTime`.
* **Key naming differs between instruments.** Which keys make up the drift voltage is
  therefore a field of the instrument profile, named explicitly, rather than anything
  this module guesses.

Observed layout, from the 2013 acquisitions (192 CRLF lines, Latin-1 -- the degree
sign in `Source Temperature (\xb0C)` is a raw 0xb0 byte):

    Instrument Configuration:
    Lteff                   <TAB><TAB>1800.0
    ADC Pushes Per IMS Increment    <TAB>1

    Experimental Instrument Parameters
    Polarity                <TAB><TAB>ES-
    IMSBias                 <TAB><TAB>132.0

A line with no tab is a blank line or a bare section header; neither carries a value
and neither appears in the result. Keys are not guaranteed unique across sections, so
the last occurrence wins -- section-qualify here if that ever matters, rather than
changing the return type.
"""

from __future__ import annotations

import os
from typing import Mapping

__all__ = [
    "EXTERN_FILENAME",
    "HEADER_FILENAME",
    "as_float",
    "cal_function_coefficients",
    "parse_extern_inf",
    "parse_header_txt",
]

EXTERN_FILENAME = "_extern.inf"
HEADER_FILENAME = "_HEADER.TXT"

# Tried in order. MassLynx is a Windows application: the 2013 files are Latin-1 and
# cp1252 is its superset, so utf-8 first (in case a newer MassLynx writes it) and
# cp1252 second decodes everything seen so far. cp1252 has a few undefined bytes;
# latin-1 has none, so it is the last resort that cannot fail.
_ENCODINGS = ("utf-8", "cp1252", "latin-1")


def _read_text(path: str) -> str | None:
    """The file's text, or None if it is not usefully decodable.

    A UTF-16 file decodes fine under latin-1 -- into mojibake studded with NULs, which
    would become junk keys. Text carrying NULs is rejected as not our format, because
    the contract promises an empty mapping for an encoding surprise rather than a
    dictionary of nonsense.
    """
    for encoding in _ENCODINGS:
        try:
            with open(path, encoding=encoding) as handle:
                text = handle.read()
        except UnicodeDecodeError:
            continue
        except OSError:
            return None
        return None if "\x00" in text else text
    return None


def parse_extern_inf(raw_dir: str | os.PathLike[str]) -> dict[str, str]:
    """Parse `<raw_dir>/_extern.inf` into a flat `{key: value}` mapping.

    `raw_dir` is the `.raw` **directory**, not the sidecar itself: a Waters
    acquisition is a directory of binary files.

    Returns `{}` for every failure mode, and does not raise.
    """
    try:
        path = os.path.join(os.fspath(raw_dir), EXTERN_FILENAME)
    except (TypeError, ValueError):
        return {}
    try:
        if not os.path.isfile(path):
            return {}
    except OSError:
        return {}

    text = _read_text(path)
    if not text:
        return {}

    parsed: dict[str, str] = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue  # blank line, or a bare section header
        key, _, value = line.partition("\t")
        key = key.strip()
        if key:
            parsed[key] = value.strip()
    return parsed


def parse_header_txt(raw_dir: str | os.PathLike[str]) -> dict[str, str]:
    """Parse `<raw_dir>/_HEADER.TXT` into a flat `{key: value}` mapping.

    The format is one `$$ Key: value` line per item, CRLF, in MassLynx's own
    encoding. The `$$ ` prefix and the first colon are the delimiters; everything
    after that colon is the value, verbatim apart from surrounding whitespace,
    because several values are comma-separated lists that contain no colon and one
    (`Cal Time`) contains two.

    Same contract as `parse_extern_inf`: `{}` for every failure mode, never raises.
    Prefer the SDK's `GetHeaderItemValue` for anything `MassLynxHeaderItem` names;
    this exists for `Cal Function N`, which it does not.
    """
    try:
        path = os.path.join(os.fspath(raw_dir), HEADER_FILENAME)
    except (TypeError, ValueError):
        return {}
    try:
        if not os.path.isfile(path):
            return {}
    except OSError:
        return {}

    text = _read_text(path)
    if not text:
        return {}

    parsed: dict[str, str] = {}
    for line in text.splitlines():
        if not line.startswith("$$"):
            continue
        key, separator, value = line[2:].partition(":")
        key = key.strip()
        if key and separator:
            parsed[key] = value.strip()
    return parsed


def cal_function_coefficients(
    raw_dir: str | os.PathLike[str], function: int = 1
) -> tuple[float, ...] | None:
    """The `$$ Cal Function N` coefficients, lowest order first, or None.

    The line is `c0,c1,...,Tn`: the coefficients, then a trailing type token that is
    dropped. Functions are numbered from **1** here, as `_HEADER.TXT` numbers them,
    and not from 0 as the SDK numbers its own -- the caller converts.

    `(0.0, 1.0)` is the identity, which is what MassLynx writes for a function that
    carries no calibration of its own, and is a different fact from the line being
    absent, which is `None`. A line that cannot be read as numbers is also `None`:
    a calibration that is half-parsed is worse than one that is missing.

    What the numbers mean is `calibration.CalFunction`, which is what to use.
    """
    value = parse_header_txt(raw_dir).get(f"Cal Function {int(function)}")
    if not value:
        return None
    fields = [f.strip() for f in value.split(",")]
    coefficients = [f for f in fields if f and not f.upper().startswith("T")]
    if len(coefficients) < 2:
        return None
    try:
        return tuple(float(f) for f in coefficients)
    except ValueError:
        return None


def as_float(extern: Mapping[str, str], key: str, default: float | None = None) -> float | None:
    """One `_extern.inf` value as a float, or `default` if it is missing or not a number.

    Never raises and never returns 0.0 for an absent key -- a voltage that was not
    recorded is `None`, and code that needs it says so. The 2013 workflow's silent
    numeric fallback is exactly what this signature is shaped to prevent.
    """
    raw = extern.get(key)
    if raw is None:
        return default
    try:
        return float(raw.strip())
    except (AttributeError, ValueError):
        return default

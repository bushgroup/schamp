"""The `_extern.inf` sidecar: the one read schamp does for itself.

Everything else about a `.raw` comes from the Waters SDK's own readers. The
acquisition parameters do not: the SDK exposes no accessor for the DC optic
voltages, so the drift voltage of a linear-field measurement is unreadable through
it, and that one gap is filled here (decision of record). Any *further* non-SDK read
needs its own decision; this is not a precedent.

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

__all__ = ["EXTERN_FILENAME", "as_float", "parse_extern_inf"]

EXTERN_FILENAME = "_extern.inf"

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

"""Getting a licensed Waters SDK reader open on a `.raw`, and nothing more.

Two things live here and both are boundary concerns rather than physics: finding the
license key the SDK insists on, and opening the reader trio that every other read goes
through. The reads themselves -- mobillograms, drift scans, spectra, combining, the CCS
conversions -- are the SDK's own and are called directly at the point of use; schamp does
not wrap them (lab record, task 02).

`masslynxsdk` is imported inside the functions that need it, so this module imports on a
machine that has no SDK and every non-extraction code path keeps working. Install the SDK
with `uv run tools/install_sdk.py <your Waters download>`.

The license key resolves in this order, first hit winning:

1. `license_key=` argument
2. `license_path=` argument
3. `$SCHAMP_LICENSE_KEY`
4. `$SCHAMP_LICENSE_PATH`
5. `[license] key` then `[license] path` in `schamp.ini`, looked for at `config_path=`,
   then `$SCHAMP_CONFIG`, then `./schamp.ini`, then `~/.schamp.ini`
6. `./license.key`
7. `external/masslynxsdk/license.key`, where `tools/install_sdk.py` puts it

1-6 are waters2python's convention with schamp's names (lab record, task 02); 7 is added
so that an install through `install_sdk.py` needs no configuration at all. The key is a
credential: it is never logged, never put in an exception message, and never written out.
"""

from __future__ import annotations

import configparser
import os
from dataclasses import dataclass, field

from . import EXTERNAL_DIR

__all__ = ["ConfigError", "Readers", "open_readers", "resolve_license"]

ENV_KEY = "SCHAMP_LICENSE_KEY"
ENV_PATH = "SCHAMP_LICENSE_PATH"
ENV_CONFIG = "SCHAMP_CONFIG"
CONFIG_NAME = "schamp.ini"
INSTALLED_KEY = os.path.join(EXTERNAL_DIR, "masslynxsdk", "license.key")


class ConfigError(Exception):
    """No usable license key, or one that could not be read."""


def _read_key_file(path: str) -> str:
    """The single-line key in a `license.key` file."""
    try:
        with open(path, encoding="utf-8-sig") as fh:
            key = fh.read().strip()
    except OSError as exc:
        raise ConfigError(f"cannot read the license key at {path}: {exc.strerror}") from exc
    if not key:
        raise ConfigError(f"the license key at {path} is empty")
    return key


def _config_candidates(config_path: str | None) -> list[str]:
    return [
        p
        for p in (
            config_path,
            os.environ.get(ENV_CONFIG),
            os.path.join(os.getcwd(), CONFIG_NAME),
            os.path.join(os.path.expanduser("~"), "." + CONFIG_NAME),
        )
        if p
    ]


def _from_config(config_path: str | None) -> tuple[str | None, str | None]:
    """`(key, source)` from the first `schamp.ini` that has one, else `(None, None)`.

    Read as utf-8-sig, not utf-8: Windows editors add byte-order marks freely, and under
    strict utf-8 the mark lands inside the `[license]` header and hides the whole section
    (waters2python found this the hard way, 2026-07-29).
    """
    for path in _config_candidates(config_path):
        if not os.path.isfile(path):
            continue
        parser = configparser.ConfigParser()
        try:
            parser.read(path, encoding="utf-8-sig")
        except configparser.Error as exc:
            raise ConfigError(f"{path} is not readable as an INI file: {exc}") from exc
        if not parser.has_section("license"):
            continue
        key = parser.get("license", "key", fallback="").strip()
        if key:
            return key, f"[license] key in {path}"
        key_path = parser.get("license", "path", fallback="").strip()
        if key_path:
            return _read_key_file(key_path), f"[license] path in {path}"
    return None, None


def resolve_license(
    *,
    license_key: str | None = None,
    license_path: str | None = None,
    config_path: str | None = None,
) -> str:
    """The Waters SDK license key, by the precedence in this module's docstring.

    Raises `ConfigError`, naming every place it looked, when nothing resolves.
    """
    if license_key:
        return license_key.strip()
    if license_path:
        return _read_key_file(license_path)
    env_key = os.environ.get(ENV_KEY, "").strip()
    if env_key:
        return env_key
    env_path = os.environ.get(ENV_PATH, "").strip()
    if env_path:
        return _read_key_file(env_path)
    from_config, _ = _from_config(config_path)
    if from_config:
        return from_config
    for path in (os.path.join(os.getcwd(), "license.key"), INSTALLED_KEY):
        if os.path.isfile(path):
            return _read_key_file(path)
    raise ConfigError(
        "no Waters SDK license key. Looked at: the license_key and license_path "
        f"arguments, ${ENV_KEY}, ${ENV_PATH}, {CONFIG_NAME} "
        f"({', '.join(_config_candidates(config_path))}), ./license.key, and "
        f"{INSTALLED_KEY}. The key ships with your Waters SDK download; "
        "`uv run tools/install_sdk.py <download>` puts it in the last of those."
    )


@dataclass
class Readers:
    """The three `Ex` readers for one open acquisition, sharing a single native handle.

    Readers open from a path or from another reader, and only the first takes the license
    key; the rest inherit its handle. The `Ex` variants are a superset of the plain ones
    and are the only ones that can read mobility data, so schamp opens those throughout.

    The native library is not thread-safe. Serialise calls on one instance, or open one
    per process.
    """

    path: str
    info: object = field(repr=False)
    scan: object = field(repr=False)
    chrom: object = field(repr=False)

    def close(self) -> None:
        """Drop the readers. Idempotent; the SDK has no explicit close."""
        self.info = self.scan = self.chrom = None

    def __enter__(self) -> Readers:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()


def open_readers(
    path: str | os.PathLike[str],
    *,
    license_key: str | None = None,
    license_path: str | None = None,
    config_path: str | None = None,
) -> Readers:
    """Open the `Ex` reader trio on a `.raw` acquisition directory.

    A Waters `.raw` is a directory, not a file; `path` is that directory. Functions and
    scans are 0-indexed here and 1-indexed in the MassLynx display, and retention times
    come back in minutes.
    """
    path = os.fspath(path)
    if not os.path.isdir(path):
        raise FileNotFoundError(
            f"{path} is not a directory. A Waters .raw acquisition is a directory of "
            "binary files, not a single file."
        )
    key = resolve_license(
        license_key=license_key, license_path=license_path, config_path=config_path
    )
    # One import site for every SDK name this module binds, and it is deferred so that
    # importing schamp.sdk works with no SDK installed.
    from masslynxsdk import (  # noqa: PLC0415
        MassLynxRawChromatogramReaderEx,
        MassLynxRawInfoReaderEx,
        MassLynxRawScanReaderEx,
    )

    info = MassLynxRawInfoReaderEx(path, key)
    return Readers(
        path=path,
        info=info,
        scan=MassLynxRawScanReaderEx(info),
        chrom=MassLynxRawChromatogramReaderEx(info),
    )

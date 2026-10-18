"""Installing the Waters MassLynx SDK from the user's own download, as a library.

Two callers share this code: `tools/install_sdk.py` in a checkout, and the installed
`schamp-acquire`'s first-launch wizard (`acquire/sdkwizard.py`), so that the frozen
application verifies and places the download exactly as the development tool does.
The SDK is licensed to the user, not to this project, and its EULA forbids
redistribution and alteration (lab record, task 01): every file is placed as shipped and
nothing is edited. Nothing here prints; progress goes to an optional `log` callable.

Layout under the install directory, `schamp.SDK_DIR` by default (gitignored
`external/masslynxsdk/` in a checkout, `%LOCALAPPDATA%\\schamp\\masslynxsdk\\` for the
installed application):

    pkg/            the unpacked wheel, and the only copy Python imports
    wheel/          the wheel as shipped, kept so the install can be reproduced
    dll/            MassLynxRaw.dll and .lib for C++ callers; the wheel carries its own copy
    headers/        the C++ headers, MassLynxRawDefs.h being the enum reference
    help/           Waters' Doxygen help, with the Python set expanded into help/python/
                    (a checkout only: the installed application does not read it)
    license.key     the per-licensee key, which schamp.sdk finds here by default
    EULA.pdf        the licence the above is used under
    INSTALLED.json  what was installed, from which download, when

This module and `confirm` in particular import `masslynxsdk`, which makes it the third
module beside `extract` and `sdk` that touches the SDK. Everything else in it -- reading
the download, unpacking, writing the key's path -- runs without one, and is what
`tools/check_public.py` exercises on a synthetic download.
"""

from __future__ import annotations

import configparser
import hashlib
import importlib
import json
import os
import shutil
import sys
import sysconfig
import tempfile
import zipfile
from datetime import datetime, timezone
from typing import Callable

from . import FROZEN, SDK_DIR

__all__ = [
    "InstallError",
    "activate",
    "confirm",
    "install",
    "probe",
    "read_manifest",
    "register",
    "set_license_path",
    "unpack",
    "verify",
]

# Members of the inner download, by the name they carry there. Everything not named here is
# left in the download: the Linux .so build (Windows only) and the C# nuget package.
WHEEL_ZIP = "MassLynxSDK wheel.zip"
DLL_ZIP = "MassLynx windows dll.zip"
HEADERS_ZIP = "MassLynx c++ headers.zip"
HELP_ZIP = "MasslynxSDK help zip/help.zip"
LICENSE_KEY = "license.key"
MANIFEST = "INSTALLED.json"
PTH_NAME = "_schamp_masslynxsdk.pth"

Log = Callable[[str], None]


class InstallError(Exception):
    """Anything that stops the install, reported without a traceback."""


def _quiet(_line: str) -> None:
    pass


# -- reading the download ----------------------------------------------------------------


def _members(zf: zipfile.ZipFile) -> dict[str, str]:
    """Map each file in the archive to its name with any leading directory stripped.

    Waters wraps both zips in a folder named for the download, and the MD5 list names files
    without it, so matching on the trailing component is what lets the two line up.
    """
    out = {}
    for name in zf.namelist():
        if name.endswith("/"):
            continue
        out.setdefault(name.rsplit("/", 1)[-1], name)
    return out


def _md5(zf: zipfile.ZipFile, name: str) -> str:
    digest = hashlib.md5()
    with zf.open(name) as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_md5summer(text: str) -> dict[str, str]:
    """Read an MD5summer list: `<hex> *<filename>`, with `#` comments."""
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, _, name = line.partition(" ")
        name = name.lstrip("*").strip()
        if len(digest) == 32 and name:
            out[name.rsplit("\\", 1)[-1].rsplit("/", 1)[-1]] = digest.lower()
    return out


def verify(zf: zipfile.ZipFile) -> list[str]:
    """Check the download against the MD5 list it carries. Returns the lines to report."""
    members = _members(zf)
    listing = next((m for m in members if m.lower().endswith(".md5")), None)
    if listing is None:
        raise InstallError(
            "this archive carries no .md5 list, so it is not the Waters download this "
            "installs (expected 667007504DDRevA.zip)"
        )
    expected = _parse_md5summer(zf.read(members[listing]).decode("utf-8", "replace"))
    if not expected:
        raise InstallError(f"the MD5 list {listing} is empty or unreadable")

    lines, bad = [], []
    for name, want in sorted(expected.items()):
        if name not in members:
            bad.append(f"{name}: named in {listing} but missing from the archive")
            continue
        got = _md5(zf, members[name])
        lines.append(f"  {'ok  ' if got == want else 'BAD '}{name}  {got}")
        if got != want:
            bad.append(f"{name}: expected {want}, got {got}")
    if bad:
        raise InstallError(
            "the download does not match its own checksums, so it is damaged or "
            "incomplete. Download it again from Waters.\n  " + "\n  ".join(bad)
        )
    return lines


def _extract_zip(zf: zipfile.ZipFile, member: str, dest: str) -> None:
    """Extract a zip that is itself a member of an open archive, flat, into dest."""
    os.makedirs(dest, exist_ok=True)
    with zf.open(member) as fh, tempfile.TemporaryFile() as tmp:
        shutil.copyfileobj(fh, tmp)
        tmp.seek(0)
        with zipfile.ZipFile(tmp) as inner:
            inner.extractall(dest)


# -- writing the install -----------------------------------------------------------------


def _fresh(path: str) -> str:
    """An empty directory at path, replacing whatever was there. Keeps installs idempotent."""
    if os.path.isdir(path):
        try:
            shutil.rmtree(path)
        except OSError as exc:
            raise InstallError(
                f"cannot replace {path}: {exc.strerror}. If a program has the SDK loaded "
                "(schamp-acquire, a notebook), close it and install again."
            ) from exc
    os.makedirs(path)
    return path


def unpack(
    outer_path: str, target: str = SDK_DIR, log: Log = _quiet, *, with_help: bool = True
) -> dict:
    """Verify the Waters download and unpack it into `target`. Returns the manifest.

    Does not import the SDK and does not write `INSTALLED.json`; `install` does both.
    `with_help=False` leaves Waters' developer help in the download: the installed
    application never reads it, and its deepest file names are long enough to pass
    Windows' 260-character path limit under a long profile path.
    """
    try:
        outer = zipfile.ZipFile(outer_path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise InstallError(f"{outer_path} is not a readable zip: {exc}") from exc
    with outer:
        log(f"verifying {os.path.basename(outer_path)} against its own MD5 list")
        for line in verify(outer):
            log(line)

        outer_members = _members(outer)
        inner_name = next(
            (m for m in outer_members if m.startswith("MassLynxSDKDownload") and m.endswith(".zip")),
            None,
        )
        if inner_name is None:
            raise InstallError("no MassLynxSDKDownload_*.zip inside the download")
        eula = next((m for m in outer_members if m.lower().endswith("eula.pdf")), None)

        log(f"unpacking into {target}")
        try:
            os.makedirs(target, exist_ok=True)
            with tempfile.TemporaryDirectory() as tmp:
                outer.extract(outer_members[inner_name], tmp)
                inner_path = os.path.join(tmp, outer_members[inner_name])
                manifest = _unpack_inner(inner_path, inner_name, target, with_help)
                if eula:
                    with outer.open(outer_members[eula]) as fh, \
                         open(os.path.join(target, "EULA.pdf"), "wb") as out:
                        shutil.copyfileobj(fh, out)
        except OSError as exc:
            where = exc.filename or target
            too_long = len(os.path.abspath(where)) >= 260
            raise InstallError(
                f"cannot write {where}: {exc.strerror or exc}"
                + (". The path is longer than Windows allows (260 characters); install "
                   "under a shorter one" if too_long else "")
            ) from exc

    manifest["source_zip"] = os.path.abspath(outer_path)
    manifest["installed_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return manifest


def _unpack_inner(inner_path: str, inner_name: str, target: str, with_help: bool) -> dict:
    with zipfile.ZipFile(inner_path) as inner:
        members = _members(inner)
        for required in (WHEEL_ZIP, LICENSE_KEY):
            if os.path.basename(required) not in members:
                raise InstallError(f"the download is missing {required}")

        # The wheel, kept as shipped and also unpacked. A py3-none-any wheel with
        # Root-Is-Purelib and no .data directory installs by extraction and nothing else,
        # which is what keeps this a placement rather than an alteration (EULA section 6).
        wheel_dir = _fresh(os.path.join(target, "wheel"))
        _extract_zip(inner, members[WHEEL_ZIP], wheel_dir)
        wheels = [f for f in os.listdir(wheel_dir) if f.endswith(".whl")]
        if len(wheels) != 1:
            raise InstallError(f"expected one wheel in {WHEEL_ZIP}, found {wheels}")
        wheel = wheels[0]

        pkg_dir = _fresh(os.path.join(target, "pkg"))
        with zipfile.ZipFile(os.path.join(wheel_dir, wheel)) as whl:
            if any(".data/" in n for n in whl.namelist()):
                raise InstallError(
                    f"{wheel} carries a .data directory, so extracting it is not a "
                    "complete install; schamp needs updating for this SDK version"
                )
            whl.extractall(pkg_dir)

        for member, name in ((DLL_ZIP, "dll"), (HEADERS_ZIP, "headers")):
            if os.path.basename(member) in members:
                _extract_zip(inner, members[os.path.basename(member)],
                             _fresh(os.path.join(target, name)))

        if with_help:
            help_dir = _fresh(os.path.join(target, "help"))
        if with_help and os.path.basename(HELP_ZIP) in members:
            _extract_zip(inner, members[os.path.basename(HELP_ZIP)], help_dir)
            python_help = os.path.join(help_dir, "help_python.zip")
            if os.path.isfile(python_help):
                with zipfile.ZipFile(python_help) as ph:
                    ph.extractall(os.path.join(help_dir, "python"))

        key_path = os.path.join(target, LICENSE_KEY)
        with inner.open(members[LICENSE_KEY]) as fh, open(key_path, "wb") as out:
            shutil.copyfileobj(fh, out)

    return {
        "download": inner_name,
        "wheel": wheel,
        "version": wheel.split("-")[1] if "-" in wheel else None,
        "package_dir": os.path.join(target, "pkg"),
        "license_key": key_path,
        "help": with_help,
    }


def register(pkg_dir: str) -> str:
    """Put pkg_dir on a checkout environment's import path with a .pth file.

    Why a .pth file and not `uv pip install` (lab record, task 02): `uv sync` uninstalls
    anything uv.lock does not name, and declaring masslynxsdk in pyproject.toml would make
    the lock unresolvable without the SDK. A .pth file carries no dist-info, so uv neither
    manages nor removes it. The installed application has no environment to register in;
    `activate` is its equivalent, run at every start by `schamp/__init__.py`.
    """
    purelib = sysconfig.get_paths()["purelib"]
    if not os.path.isdir(purelib):
        raise InstallError(f"no site-packages at {purelib}")
    pth = os.path.join(purelib, PTH_NAME)
    with open(pth, "w", encoding="utf-8") as fh:
        fh.write(os.path.abspath(pkg_dir) + "\n")
    return pth


def activate(target: str = SDK_DIR) -> bool:
    """Put an installed SDK's `pkg/` on this process's import path. Returns whether one was."""
    pkg_dir = os.path.abspath(os.path.join(target, "pkg"))
    if not os.path.isdir(os.path.join(pkg_dir, "masslynxsdk")):
        return False
    if pkg_dir not in sys.path:
        sys.path.append(pkg_dir)
    importlib.invalidate_caches()
    return True


def confirm(pkg_dir: str) -> tuple[str, str]:
    """Import the SDK from pkg_dir here and now. Returns its version and its DLL path.

    This is the check that matters: it loads the native MassLynxRaw.dll through ctypes, so
    it fails here rather than at the first read if the Visual C++ runtime is missing.
    """
    pkg_dir = os.path.abspath(pkg_dir)
    if pkg_dir not in sys.path:
        sys.path.insert(0, pkg_dir)
    return probe()


def probe() -> tuple[str, str]:
    """Import whichever SDK this process can see. Returns its version and its DLL path.

    Raises `InstallError` saying why not. `confirm` is this with the path put first.
    """
    importlib.invalidate_caches()
    try:
        import masslynxsdk  # noqa: PLC0415
        from masslynxsdk.Providers.MassLynxProvider import MassLynxProvider  # noqa: PLC0415
    except OSError as exc:
        raise InstallError(
            f"the SDK's native library would not load: {exc}\n"
            "MassLynxRaw.dll needs the Microsoft Visual C++ redistributable."
        ) from exc
    except ImportError as exc:
        raise InstallError(f"the SDK would not import: {exc}") from exc
    version = getattr(masslynxsdk, "__version__", None)
    return version or "5.0.0", MassLynxProvider.MassLynxPath


def read_manifest(target: str = SDK_DIR) -> dict | None:
    """The `INSTALLED.json` of an install, or None if there is none."""
    try:
        with open(os.path.join(target, MANIFEST), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _write_manifest(target: str, manifest: dict) -> None:
    with open(os.path.join(target, MANIFEST), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")


def set_license_path(key_path: str, ini_path: str | None = None) -> str:
    """Write `[license] path = key_path` into `~/.schamp.ini` (or `ini_path`). Returns it.

    For a key that Waters sent separately from the download; a download's own key is
    found where `unpack` put it without this. Other sections and keys of an existing file
    are kept (its comments are not: configparser does not carry them). Only the key's
    path is written, never the key.
    """
    from . import sdk  # noqa: PLC0415 -- the reader of what this writes

    key_path = os.path.abspath(key_path)
    try:
        sdk._read_key_file(key_path)
    except sdk.ConfigError as exc:
        raise InstallError(str(exc)) from exc
    ini_path = ini_path or os.path.join(os.path.expanduser("~"), "." + sdk.CONFIG_NAME)
    parser = configparser.ConfigParser()
    if os.path.isfile(ini_path):
        try:
            parser.read(ini_path, encoding="utf-8-sig")
        except configparser.Error as exc:
            raise InstallError(f"{ini_path} is not readable as an INI file: {exc}") from exc
    if not parser.has_section("license"):
        parser.add_section("license")
    parser.set("license", "path", key_path)
    with open(ini_path, "w", encoding="utf-8") as fh:
        parser.write(fh)
    return ini_path


def install(
    outer_path: str,
    *,
    target: str = SDK_DIR,
    license_file: str | None = None,
    link_environment: bool = False,
    with_help: bool | None = None,
    log: Log = _quiet,
) -> dict:
    """Verify, unpack, import and record the SDK. Returns the manifest written.

    `license_file`, when given, is a key Waters sent apart from the download, and its path
    goes into `~/.schamp.ini`. `link_environment` also writes the checkout's `.pth` file
    (`register`); the installed application leaves it False and relies on `activate`.
    `with_help` defaults to a checkout's yes and the installed application's no (`unpack`).
    """
    help_wanted = not FROZEN if with_help is None else with_help
    manifest = unpack(outer_path, target, log, with_help=help_wanted)
    if link_environment:
        manifest["pth"] = register(manifest["package_dir"])
    log("importing the SDK and loading its native library")
    version, dll = confirm(manifest["package_dir"])
    manifest["version"] = version
    manifest["dll"] = dll
    if license_file:
        manifest["license_ini"] = set_license_path(license_file)
        log(f"key path written to {manifest['license_ini']}")
    _write_manifest(target, manifest)
    activate(target)
    return manifest

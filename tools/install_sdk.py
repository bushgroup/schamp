"""Install the Waters MassLynx SDK into this project, from the user's own download.

    uv run tools/install_sdk.py <path-to-667007504DDRevA.zip>

Verifies the download against the MD5 list Waters ships inside it, unpacks it into
gitignored `external/masslynxsdk/`, and puts that tree on the project environment's import
path. Windows only. Idempotent: run it again to repair or upgrade an install.

Nothing this script writes is ever tracked by git, and nothing is redistributed: the SDK is
licensed to the user, not to this project (lab record, task 01 reads the EULA). Files are
placed, never edited.

Layout it creates under external/masslynxsdk/:

    pkg/            the unpacked wheel, and the only copy Python imports
    wheel/          the wheel as shipped, kept so the install can be reproduced
    dll/            MassLynxRaw.dll and .lib for C++ callers; the wheel carries its own copy
    headers/        the C++ headers, MassLynxRawDefs.h being the enum reference
    help/           Waters' Doxygen help, with the Python set expanded into help/python/
    license.key     the per-licensee key, which schamp.sdk finds here by default
    EULA.pdf        the licence the above is used under
    INSTALLED.json  what was installed, from which download, when

Why a .pth file and not `uv pip install` (lab record, task 02). `uv sync` makes the
environment match uv.lock exactly and uninstalls anything the lock does not name, so a wheel
installed with `uv pip install` is removed by the next sync; this was measured, not assumed.
Declaring masslynxsdk in pyproject.toml instead would break the other half of the
requirement, because then a clone without the SDK could not resolve the lock at all. A .pth
file carries no dist-info, so uv does not manage it and does not remove it: a clone with the
SDK keeps it across `uv sync`, and a clone without the SDK is untouched. The cost is that
`uv sync` after a Python version change rebuilds the environment and drops the .pth; rerun
this script when `import masslynxsdk` stops working.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import sysconfig
import tempfile
import zipfile
from datetime import datetime, timezone

import schamp

# Members of the inner download, by the name they carry there. Everything not named here is
# left in the download: the Linux .so build (Windows only) and the C# nuget package.
WHEEL_ZIP = "MassLynxSDK wheel.zip"
DLL_ZIP = "MassLynx windows dll.zip"
HEADERS_ZIP = "MassLynx c++ headers.zip"
HELP_ZIP = "MasslynxSDK help zip/help.zip"
LICENSE_KEY = "license.key"

TARGET = os.path.join(schamp.EXTERNAL_DIR, "masslynxsdk")
PTH_NAME = "_schamp_masslynxsdk.pth"


class InstallError(Exception):
    """Anything that stops the install, reported without a traceback."""


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
            "script installs (expected 667007504DDRevA.zip)"
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
        shutil.rmtree(path)
    os.makedirs(path)
    return path


def unpack(outer_path: str) -> dict:
    """Unpack the Waters download into external/masslynxsdk/. Returns the manifest."""
    with zipfile.ZipFile(outer_path) as outer:
        print(f"verifying {os.path.basename(outer_path)} against its own MD5 list")
        for line in verify(outer):
            print(line)

        outer_members = _members(outer)
        inner_name = next(
            (m for m in outer_members if m.startswith("MassLynxSDKDownload") and m.endswith(".zip")),
            None,
        )
        if inner_name is None:
            raise InstallError("no MassLynxSDKDownload_*.zip inside the download")
        eula = next((m for m in outer_members if m.lower().endswith("eula.pdf")), None)

        print(f"unpacking into {TARGET}")
        os.makedirs(TARGET, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            outer.extract(outer_members[inner_name], tmp)
            inner_path = os.path.join(tmp, outer_members[inner_name])
            manifest = _unpack_inner(inner_path, inner_name)
            if eula:
                with outer.open(outer_members[eula]) as fh, \
                     open(os.path.join(TARGET, "EULA.pdf"), "wb") as out:
                    shutil.copyfileobj(fh, out)

    manifest["source_zip"] = os.path.abspath(outer_path)
    manifest["installed_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return manifest


def _unpack_inner(inner_path: str, inner_name: str) -> dict:
    with zipfile.ZipFile(inner_path) as inner:
        members = _members(inner)
        for required in (WHEEL_ZIP, LICENSE_KEY):
            if os.path.basename(required) not in members:
                raise InstallError(f"the download is missing {required}")

        # The wheel, kept as shipped and also unpacked. A py3-none-any wheel with
        # Root-Is-Purelib and no .data directory installs by extraction and nothing else,
        # which is what keeps this a placement rather than an alteration (EULA section 6).
        wheel_dir = _fresh(os.path.join(TARGET, "wheel"))
        _extract_zip(inner, members[WHEEL_ZIP], wheel_dir)
        wheels = [f for f in os.listdir(wheel_dir) if f.endswith(".whl")]
        if len(wheels) != 1:
            raise InstallError(f"expected one wheel in {WHEEL_ZIP}, found {wheels}")
        wheel = wheels[0]

        pkg_dir = _fresh(os.path.join(TARGET, "pkg"))
        with zipfile.ZipFile(os.path.join(wheel_dir, wheel)) as whl:
            if any(".data/" in n for n in whl.namelist()):
                raise InstallError(
                    f"{wheel} carries a .data directory, so extracting it is not a "
                    "complete install; this script needs updating for this SDK version"
                )
            whl.extractall(pkg_dir)

        for member, name in ((DLL_ZIP, "dll"), (HEADERS_ZIP, "headers")):
            if os.path.basename(member) in members:
                _extract_zip(inner, members[os.path.basename(member)],
                             _fresh(os.path.join(TARGET, name)))

        help_dir = _fresh(os.path.join(TARGET, "help"))
        if os.path.basename(HELP_ZIP) in members:
            _extract_zip(inner, members[os.path.basename(HELP_ZIP)], help_dir)
            python_help = os.path.join(help_dir, "help_python.zip")
            if os.path.isfile(python_help):
                with zipfile.ZipFile(python_help) as ph:
                    ph.extractall(os.path.join(help_dir, "python"))

        key_path = os.path.join(TARGET, LICENSE_KEY)
        with inner.open(members[LICENSE_KEY]) as fh, open(key_path, "wb") as out:
            shutil.copyfileobj(fh, out)

    return {
        "download": inner_name,
        "wheel": wheel,
        "version": wheel.split("-")[1] if "-" in wheel else None,
        "package_dir": os.path.join(TARGET, "pkg"),
        "license_key": key_path,
    }


def register(pkg_dir: str) -> str:
    """Put pkg_dir on the environment's import path with a .pth file. Returns its path."""
    purelib = sysconfig.get_paths()["purelib"]
    if not os.path.isdir(purelib):
        raise InstallError(f"no site-packages at {purelib}")
    pth = os.path.join(purelib, PTH_NAME)
    with open(pth, "w", encoding="utf-8") as fh:
        fh.write(os.path.abspath(pkg_dir) + "\n")
    return pth


def confirm(pkg_dir: str) -> tuple[str, str]:
    """Import the SDK from pkg_dir here and now. Returns its version and its DLL path.

    This is the check that matters: it loads the native MassLynxRaw.dll through ctypes, so
    it fails here rather than at the first read if the Visual C++ runtime is missing.
    """
    sys.path.insert(0, os.path.abspath(pkg_dir))
    try:
        import masslynxsdk  # noqa: PLC0415
        from masslynxsdk.Providers.MassLynxProvider import MassLynxProvider  # noqa: PLC0415
    except OSError as exc:
        raise InstallError(
            f"the SDK unpacked but its native library would not load: {exc}\n"
            "MassLynxRaw.dll needs the Microsoft Visual C++ redistributable."
        ) from exc
    except ImportError as exc:
        raise InstallError(f"the SDK unpacked but would not import: {exc}") from exc
    version = getattr(masslynxsdk, "__version__", None)
    return version or "5.0.0", MassLynxProvider.MassLynxPath


# -- entry point -------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="uv run tools/install_sdk.py",
        description="Install the Waters MassLynx SDK from your own Waters download.",
        epilog="The SDK is licensed to you by Waters, not by this project. Nothing this "
               "script writes is tracked by git, and none of it may be redistributed.",
    )
    parser.add_argument(
        "download",
        nargs="?",
        help="the Waters download, e.g. 667007504DDRevA.zip",
    )
    parser.add_argument(
        "--relink",
        action="store_true",
        help="skip the download and put an SDK that is already unpacked back on this "
             "environment's import path. Use it for a second environment, or after a "
             "Python version change rebuilt the environment.",
    )
    args = parser.parse_args(argv)

    if sys.platform != "win32":
        print("This project reads Waters .raw files on Windows only.", file=sys.stderr)
        return 2

    try:
        import schamp as _  # noqa: F401,PLC0415
    except ImportError:  # pragma: no cover - only reachable outside the environment
        print("Run this with `uv run tools/install_sdk.py`, from the repository root.",
              file=sys.stderr)
        return 2

    if args.relink:
        try:
            with open(os.path.join(TARGET, "INSTALLED.json"), encoding="utf-8") as fh:
                manifest = json.load(fh)
        except OSError:
            print(f"Nothing is unpacked at {TARGET}. Run this with the Waters download "
                  "first.", file=sys.stderr)
            return 2
        return _report(manifest, relink_only=True)

    if not args.download:
        parser.print_usage(sys.stderr)
        print("\nGive it the SDK zip you downloaded from Waters (document 667007504).",
              file=sys.stderr)
        lab = schamp.lab_dir("sdk-download")
        if lab:
            for name in sorted(os.listdir(lab)):
                if name.lower().endswith(".zip"):
                    print(f"  there is one here: {os.path.join(lab, name)}", file=sys.stderr)
        return 2
    if not os.path.isfile(args.download):
        print(f"No such file: {args.download}", file=sys.stderr)
        return 2

    try:
        manifest = unpack(args.download)
    except InstallError as exc:
        print(f"\nInstall failed. {exc}", file=sys.stderr)
        return 1
    return _report(manifest)


def _report(manifest: dict, relink_only: bool = False) -> int:
    """Link the unpacked SDK into this environment, prove it imports, and say where."""
    try:
        pth = register(manifest["package_dir"])
        version, dll = confirm(manifest["package_dir"])
    except InstallError as exc:
        print(f"\nInstall failed. {exc}", file=sys.stderr)
        return 1
    manifest["pth"] = pth
    manifest["version"] = version
    manifest["dll"] = dll
    if not relink_only:
        with open(os.path.join(TARGET, "INSTALLED.json"), "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)
            fh.write("\n")

    print()
    print(f"masslynxsdk {version} is installed and imports in this environment.")
    print(f"  package      {manifest['package_dir']}")
    print(f"  native DLL   {dll}")
    print(f"  license key  {manifest['license_key']}")
    print(f"  import path  {manifest['pth']}")
    print()
    print("schamp finds that key on its own, so there is nothing else to set up. To use a")
    print("different key, set SCHAMP_LICENSE_PATH, or write schamp.ini in the working")
    print("directory:")
    print()
    print("    [license]")
    print("    path = C:\\path\\to\\license.key")
    print()
    print("Check the install any time with `uv run tools/check_public.py`.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

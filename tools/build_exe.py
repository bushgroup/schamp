"""Freeze `schamp-acquire` once and wrap it two ways: a per-user installer and a portable zip.

    uv run tools/build_exe.py                  build, screen, smoke-test, zip, installer
    uv run tools/build_exe.py --sdk-zip ZIP    also install the SDK from ZIP into a scratch
                                               profile and read $SCHAMP_SMOKE_RAW, frozen
    uv run tools/build_exe.py --no-installer   stop after the zip
    uv run tools/build_exe.py --no-prune       keep every Qt piece PyInstaller collects, to
                                               measure what the pruning saves

Windows 10/11 x64, from this project's own 3.12 environment with the `dev` group (which
carries PyInstaller and PySide6); there is no second requirements file. Outputs, all in
gitignored `dist/`, named with `schamp.__version__`:

    dist/schamp-acquire/                               the one-folder build
    dist/schamp-acquire-<version>-win64-portable.zip   that folder, zipped
    dist/schamp-acquire-<version>-win64-setup.exe      tools/installer/schamp-acquire.iss

PyInstaller's work, the generated spec, the build stamp and the cached Visual C++
redistributable live in gitignored `build/`. The installer needs Inno Setup 6's `ISCC.exe`
(a development-machine tool, not a project dependency); without it the build says so and
stops after the zip.

**Nothing SDK-derived is in any output.** The Waters MassLynx SDK's EULA forbids
redistribution (lab record, task 01): each lab installs it from its own download at first
launch, into its own profile (`schamp.sdk_install`). PyInstaller would otherwise find
`masslynxsdk` through this environment's `.pth` file, because it follows imports inside
functions too, so the spec excludes it and refuses to build if any of it is collected
anyway, and the finished folder is screened again by name before it is zipped.

Facts this inherits from the frozen-build spike (lab record, task 12): one folder, never
one file (a one-file build unpacks itself into `%TEMP%` on every start); `numpy.random`
reaches `secrets` from compiled Cython, invisible to the analysis; the SDK's own loader
finds its DLL beside its package, so the installed SDK's `pkg/` goes on `sys.path` at
start (`schamp/__init__.py`) and nothing of it is needed in the bundle -- only the few
standard-library modules its Python imports, which the application may not import itself.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone

import schamp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DIST = os.path.join(ROOT, "dist")
BUILD = os.path.join(ROOT, "build")
NAME = "schamp-acquire"
FOLDER = os.path.join(DIST, NAME)
ISS = os.path.join(HERE, "installer", "schamp-acquire.iss")

# The standard-library modules the SDK's own Python imports (`array`, `ctypes`, `enum`),
# which a frozen program contains only if something in the bundle imports them: the SDK
# is not in the bundle, so its imports are not seen. `secrets` is task 12's numpy fact.
HIDDEN_IMPORTS = ("array", "ctypes", "ctypes.wintypes", "enum", "secrets")

# Never collected. The SDK above all; then other Qt bindings and backends matplotlib
# would probe for, and the notebook stack the dev group carries.
EXCLUDES = (
    "masslynxsdk",
    "tkinter", "_tkinter", "PyQt5", "PyQt6", "PySide2",
    "IPython", "ipykernel", "jupyter", "jupyter_client", "jupyter_core", "nbclient",
    "nbformat", "notebook", "zmq", "tornado", "PyInstaller",
)

# The Qt the application does not use. PySide6-Essentials is already without WebEngine
# and Multimedia (those are Addons), and PyInstaller's hooks collect only the Qt modules
# something imports: QtCore, QtGui and QtWidgets from the application, QtSvg from
# matplotlib's Qt backend (unconditionally, so it stays), QtNetwork from PySide6 itself.
# Excluded on top, so that nothing reaches them by accident: QML/Quick, SQL, the designer
# and help tooling.
QT_EXCLUDES = (
    "PySide6.QtQml", "PySide6.QtQuick", "PySide6.QtQuickWidgets", "PySide6.QtQuickControls2",
    "PySide6.QtSql", "PySide6.QtTest", "PySide6.QtDesigner", "PySide6.QtUiTools",
    "PySide6.QtHelp", "PySide6.QtOpenGL", "PySide6.QtOpenGLWidgets",
    # Pulled in through optional imports (cffi's distutils shim, numpy's build-config
    # printer) and never reached by the application.
    "setuptools", "distutils", "_distutils_hack", "cffi", "pycparser", "yaml",
)
# Collected files to drop, by a pattern on their path in the bundle: Qt's software OpenGL
# (20 MB; the canvas is raster), the TLS, network-status, touch and alternative platform
# plugins, and every translation but Qt's own base strings for the standard dialogs.
QT_PRUNE = (
    r"PySide6[\\/]opengl32sw\.dll$",
    r"PySide6[\\/]plugins[\\/](tls|networkinformation|generic)[\\/]",
    r"PySide6[\\/]plugins[\\/]platforms[\\/]q(minimal|direct2d)\.dll$",
    r"PySide6[\\/]translations[\\/](?!qtbase_)",
)

# Microsoft's Visual C++ v14 redistributable (x64), by its versioned download URL and
# its SHA-256, which the URL's own path carries. Redistributable with the program by a
# licensed Visual Studio user under the Distributable Code terms, so an installer built
# for distribution is built on a machine with Visual Studio Build Tools or Community;
# tools/installer/schamp-acquire.iss cites the terms. When Microsoft retires this file,
# take the current one from https://aka.ms/vc14/vc_redist.x64.exe, read the terms again,
# and update both.
VC_REDIST_URL = (
    "https://download.visualstudio.microsoft.com/download/pr/"
    "ebdab8e5-1d7b-4d9f-a11b-cbb1720c3b12/"
    "843068991DAAA1F73AD9F6239BCE4D0F6A07A51F18C37EA2A867E9BECA71295C/VC_redist.x64.exe"
)
VC_REDIST_SHA256 = "843068991daaa1f73ad9f6239bce4d0f6a07a51f18c37ea2a867e9beca71295c"

# Names that must never appear in an output: the SDK, its DLL, its key, its wheel.
FORBIDDEN = re.compile(r"(masslynx|license\.key$|\.key$|\.whl$)", re.IGNORECASE)

SPEC = '''\
# Generated by tools/build_exe.py on every build; edit that, not this.
import re

FORBIDDEN = re.compile({forbidden!r}, re.IGNORECASE)
PRUNE = [re.compile(p) for p in {prune!r}]

a = Analysis(
    [{entry!r}],
    pathex=[{src!r}],
    datas=[({data!r}, "schamp/data"), ({stamp!r}, "schamp")],
    hiddenimports={hidden!r},
    excludes={excludes!r},
    hooksconfig={{"matplotlib": {{"backends": ["QtAgg"]}}}},
    noarchive=False,
)
leaked = sorted({{name for name, *_ in a.pure + a.binaries + a.datas if FORBIDDEN.search(name)}})
if leaked:
    raise SystemExit("refusing to build: SDK-derived names were collected: " + ", ".join(leaked))
before = len(a.binaries) + len(a.datas)
a.binaries = [t for t in a.binaries if not any(p.search(t[0]) for p in PRUNE)]
a.datas = [t for t in a.datas if not any(p.search(t[0]) for p in PRUNE)]
print("pruned", before - len(a.binaries) - len(a.datas), "collected files")

pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name={name!r},
    console=False,
    upx=False,
    version={version_file!r},
)
coll = COLLECT(exe, a.binaries, a.datas, name={name!r}, upx=False)
'''

VERSION_FILE = '''\
VSVersionInfo(
  ffi=FixedFileInfo(filevers={tup}, prodvers={tup}, mask=0x3f, flags=0x0, OS=0x40004,
                    fileType=0x1, subtype=0x0, date=(0, 0)),
  kids=[
    StringFileInfo([StringTable("040904B0", [
      StringStruct("CompanyName", "University of Washington"),
      StringStruct("FileDescription", "schamp-acquire: run sheet and live analysis for drift-voltage series"),
      StringStruct("FileVersion", "{version}"),
      StringStruct("InternalName", "schamp-acquire"),
      StringStruct("LegalCopyright", "Copyright (c) University of Washington. BSD 3-Clause."),
      StringStruct("OriginalFilename", "schamp-acquire.exe"),
      StringStruct("ProductName", "schamp"),
      StringStruct("ProductVersion", "{version}")])]),
    VarFileInfo([VarStruct("Translation", [1033, 1200])])
  ]
)
'''

PORTABLE_README = """\
schamp-acquire {version}, portable
==================================

Unzip anywhere you can write to and run schamp-acquire.exe; nothing is installed. To
remove it, delete the folder.

The Live tab reads acquisitions through the Waters MassLynx SDK, which Waters licenses
to each lab and which schamp does not include. At first start schamp-acquire offers to
install it from your own Waters download (document 667007504); it goes into

    %LOCALAPPDATA%\\schamp\\masslynxsdk\\

and stays there when this folder is deleted. Delete it by hand when no copy of schamp
needs it. A license key that Waters sent apart from the download is recorded in
%USERPROFILE%\\.schamp.ini as [license] path.

If Windows reports that MSVCP140.dll or VCRUNTIME140.dll is missing, install
Microsoft's Visual C++ Redistributable (x64): https://aka.ms/vc14/vc_redist.x64.exe

schamp is BSD 3-Clause licensed, copyright University of Washington (LICENSE.txt). The
components it is built with are under their own licences, in licenses\\.

Source and documentation: https://github.com/bushgroup/schamp
Build: {commit}, {built}
"""


class BuildError(Exception):
    pass


def say(text: str = "") -> None:
    print(text, flush=True)


# -- preconditions and the stamp -----------------------------------------------------------


def preconditions() -> None:
    if sys.platform != "win32" or platform.machine().lower() not in ("amd64", "x86_64"):
        raise BuildError("schamp-acquire is built on, and for, Windows 10/11 x64")
    if sys.version_info[:2] != (3, 12):
        raise BuildError("build from this project's 3.12 environment: uv run tools/build_exe.py")
    for module, why in (("PyInstaller", "the dev group"), ("PySide6", "the dev group or the acquire extra")):
        try:
            __import__(module)
        except ImportError:
            raise BuildError(f"{module} is not installed here; it comes with {why} (uv sync)") from None


def git(*args: str) -> str:
    try:
        return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def stamp() -> dict:
    commit = git("rev-parse", "--short", "HEAD") or "unknown"
    dirty = bool(git("status", "--porcelain", "--untracked-files=no"))
    import PyInstaller  # noqa: PLC0415
    import PySide6  # noqa: PLC0415

    return {
        "version": schamp.__version__,
        "commit": commit + ("+dirty" if dirty else ""),
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "pyinstaller": PyInstaller.__version__,
        "pyside6": PySide6.__version__,
    }


def version_tuple(version: str) -> tuple[int, int, int, int]:
    parts = [int(p) for p in re.findall(r"\d+", version)[:4]]
    return tuple(parts + [0] * (4 - len(parts)))  # type: ignore[return-value]


# -- the PyInstaller build ------------------------------------------------------------------


def freeze(build_stamp: dict, prune: bool) -> None:
    os.makedirs(BUILD, exist_ok=True)
    stamp_path = os.path.join(BUILD, "BUILD.json")
    with open(stamp_path, "w", encoding="utf-8") as fh:
        json.dump(build_stamp, fh, indent=2)
    entry = os.path.join(BUILD, "schamp_acquire_main.py")
    with open(entry, "w", encoding="utf-8") as fh:
        fh.write("from schamp.acquire.app import main\n\nraise SystemExit(main())\n")
    version_file = os.path.join(BUILD, "version_info.txt")
    with open(version_file, "w", encoding="utf-8") as fh:
        fh.write(VERSION_FILE.format(tup=version_tuple(build_stamp["version"]), version=build_stamp["version"]))
    spec = os.path.join(BUILD, f"{NAME}.spec")
    with open(spec, "w", encoding="utf-8") as fh:
        fh.write(SPEC.format(
            forbidden=FORBIDDEN.pattern,
            prune=list(QT_PRUNE) if prune else [],
            entry=entry,
            src=os.path.join(ROOT, "src"),
            data=schamp.DATA_DIR,
            stamp=stamp_path,
            hidden=list(HIDDEN_IMPORTS),
            excludes=list(EXCLUDES + (QT_EXCLUDES if prune else ())),
            name=NAME,
            version_file=version_file,
        ))
    if os.path.isdir(FOLDER):
        shutil.rmtree(FOLDER)
    args = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--log-level", "WARN",
            "--distpath", DIST, "--workpath", os.path.join(BUILD, "pyinstaller"), spec]
    say("$ " + " ".join(args))
    if subprocess.run(args).returncode:
        raise BuildError("PyInstaller failed; its output is above")


def add_notices(build_stamp: dict) -> list[str]:
    """LICENSE.txt, README.txt and licenses/<component>/ beside the executable."""
    shutil.copyfile(os.path.join(ROOT, "LICENSE"), os.path.join(FOLDER, "LICENSE.txt"))
    with open(os.path.join(FOLDER, "README.txt"), "w", encoding="utf-8", newline="\r\n") as fh:
        fh.write(PORTABLE_README.format(version=build_stamp["version"], commit=build_stamp["commit"],
                                        built=build_stamp["built_utc"]))
    out = os.path.join(FOLDER, "licenses")
    os.makedirs(out, exist_ok=True)
    missing = []
    # Python itself, and the PyInstaller bootloader that is the executable's front end.
    python_license = os.path.join(sys.base_prefix, "LICENSE.txt")
    if os.path.isfile(python_license):
        os.makedirs(os.path.join(out, "Python"), exist_ok=True)
        shutil.copyfile(python_license, os.path.join(out, "Python", "LICENSE.txt"))
    else:
        missing.append("Python")
    for dist in sorted(runtime_distributions() | {"pyinstaller"}):
        if not copy_license_files(dist, os.path.join(out, dist)):
            missing.append(dist)
    return missing


def runtime_distributions() -> set[str]:
    """Every distribution `schamp[acquire]` needs at run time, by walking its requirements."""
    from packaging.requirements import Requirement  # noqa: PLC0415 -- PyInstaller's own dependency

    seen: set[str] = set()
    todo = [("schamp", {"acquire"})]
    while todo:
        name, extras = todo.pop()
        key = name.lower().replace("_", "-")
        if key in seen:
            continue
        try:
            dist = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue
        seen.add(key)
        for line in dist.requires or ():
            req = Requirement(line)
            if req.marker and not any(req.marker.evaluate({"extra": e}) for e in (extras or {""})):
                continue
            todo.append((req.name, set(req.extras)))
    seen.discard("schamp")
    return seen


def copy_license_files(name: str, dest: str) -> bool:
    try:
        dist = metadata.distribution(name)
    except metadata.PackageNotFoundError:
        return False
    # The licence files a wheel carries in its .dist-info: under licenses/ (PEP 639), or
    # named LICENSE, COPYING, NOTICE or AUTHORS directly in it.
    found = []
    for f in dist.files or ():
        parts = str(f).replace("\\", "/").split("/")
        if not parts[0].endswith(".dist-info"):
            continue
        if (len(parts) > 2 and parts[1] == "licenses") or re.match(
            r"(LICEN[CS]E|COPYING|NOTICE|AUTHORS)", parts[-1], re.IGNORECASE
        ):
            if os.path.isfile(f.locate()):
                found.append(f)
    os.makedirs(dest, exist_ok=True)
    for f in found:
        shutil.copyfile(f.locate(), os.path.join(dest, os.path.basename(str(f))))
    licence = dist.metadata.get("License-Expression") or dist.metadata.get("License") or ""
    with open(os.path.join(dest, "METADATA.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"{dist.metadata['Name']} {dist.version}\n")
        if licence and len(licence) < 200:
            fh.write(f"License: {licence}\n")
        for line in dist.metadata.get_all("Classifier") or ():
            if line.startswith("License"):
                fh.write(f"Classifier: {line}\n")
        if dist.metadata.get("Home-page"):
            fh.write(f"Home-page: {dist.metadata['Home-page']}\n")
    return bool(found)


# -- screening the folder -------------------------------------------------------------------


def folder_size(folder: str) -> tuple[int, int]:
    total = count = 0
    for root, _dirs, files in os.walk(folder):
        for name in files:
            total += os.path.getsize(os.path.join(root, name))
            count += 1
    return total, count


def screen_names(folder: str) -> None:
    bad = []
    for root, dirs, files in os.walk(folder):
        for name in dirs + files:
            if FORBIDDEN.search(name):
                bad.append(os.path.relpath(os.path.join(root, name), folder))
    if bad:
        raise BuildError("SDK-derived names in the build, which must never ship: " + ", ".join(bad))


def screen_pe(folder: str) -> list[str]:
    """Every DLL each PE in the folder imports resolves, in the folder or in Windows.

    The spike's `pe_audit.py` walk without its Windows 7 question: a cheap static screen
    that the folder is complete. It does not see delay-loaded or `LoadLibrary` calls.
    """
    import pefile  # noqa: PLC0415 -- PyInstaller's own dependency on Windows

    present = {name.lower() for _root, _dirs, files in os.walk(folder) for name in files}
    system = os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32")
    unresolved: dict[str, set[str]] = {}
    for root, _dirs, files in os.walk(folder):
        for name in files:
            if os.path.splitext(name)[1].lower() not in (".exe", ".dll", ".pyd"):
                continue
            path = os.path.join(root, name)
            try:
                pe = pefile.PE(path, fast_load=True)
            except pefile.PEFormatError:
                continue
            try:
                pe.parse_data_directories(directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"]])
                for entry in getattr(pe, "DIRECTORY_ENTRY_IMPORT", []) or []:
                    dll = entry.dll.decode("ascii", "replace").lower()
                    if dll in present or dll.startswith(("api-ms-win-", "ext-ms-")):
                        continue
                    if os.path.isfile(os.path.join(system, dll)):
                        continue
                    unresolved.setdefault(dll, set()).add(os.path.relpath(path, folder))
            finally:
                pe.close()
    return [f"{dll} (imported by {', '.join(sorted(users)[:3])})" for dll, users in sorted(unresolved.items())]


def crt_inventory(folder: str) -> list[str]:
    """The Visual C++ runtime DLLs the build carries beside itself, with their versions."""
    out = []
    for root, _dirs, files in os.walk(folder):
        for name in sorted(files):
            if re.match(r"(msvcp140|vcruntime140|concrt140)\w*\.dll$", name, re.IGNORECASE):
                out.append(os.path.relpath(os.path.join(root, name), folder))
    return out


# -- running the frozen application ---------------------------------------------------------


def frozen_check(local: str, *, install_zip: str | None = None) -> dict:
    """Run the built executable with `--self-check` (after `--install-sdk`) in a scratch profile."""
    report = os.path.join(local, "check.json")
    args = [os.path.join(FOLDER, f"{NAME}.exe")]
    if install_zip:
        args += ["--install-sdk", install_zip]
    args += ["--self-check", report]
    env = dict(os.environ, LOCALAPPDATA=local, USERPROFILE=os.path.join(local, "home"))
    for name in ("SCHAMP_LICENSE_KEY", "SCHAMP_LICENSE_PATH", "SCHAMP_CONFIG", "SCHAMP_LAB"):
        env.pop(name, None)
    os.makedirs(env["USERPROFILE"], exist_ok=True)
    done = subprocess.run(args, env=env, cwd=local, timeout=600)
    try:
        with open(report, encoding="utf-8") as fh:
            result = json.load(fh)
    except (OSError, ValueError) as exc:
        raise BuildError(f"the frozen application wrote no report (exit {done.returncode}): {exc}") from exc
    result["exit"] = done.returncode
    return result


# -- the two wrappers -----------------------------------------------------------------------


def portable_zip(version: str) -> str:
    path = os.path.join(DIST, f"{NAME}-{version}-win64-portable.zip")
    top = f"{NAME}-{version}"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for root, _dirs, files in os.walk(FOLDER):
            for name in sorted(files):
                full = os.path.join(root, name)
                zf.write(full, os.path.join(top, os.path.relpath(full, FOLDER)))
    return path


def find_iscc() -> str | None:
    """ISCC on PATH, else where Inno Setup 6 installs it per user or per machine."""
    candidates = [shutil.which("ISCC.exe")]
    for base, sub in (("LOCALAPPDATA", "Programs"), ("ProgramFiles(x86)", ""), ("ProgramFiles", "")):
        if os.environ.get(base):
            candidates.append(os.path.join(os.environ[base], sub, "Inno Setup 6", "ISCC.exe"))
    return next((c for c in candidates if c and os.path.isfile(c)), None)


def vc_redist() -> str:
    """The pinned redistributable, downloaded once into build/redist/ and checked every time."""
    folder = os.path.join(BUILD, "redist")
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "vc_redist.x64.exe")
    if not os.path.isfile(path):
        say(f"downloading {VC_REDIST_URL}")
        with urllib.request.urlopen(VC_REDIST_URL, timeout=120) as response, open(path + ".part", "wb") as fh:
            shutil.copyfileobj(response, fh)
        os.replace(path + ".part", path)
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    if digest.hexdigest() != VC_REDIST_SHA256:
        os.remove(path)
        raise BuildError(
            f"vc_redist.x64.exe has SHA-256 {digest.hexdigest()}, not the pinned {VC_REDIST_SHA256}. "
            "Microsoft has replaced the file: read the redistribution terms again and update "
            "VC_REDIST_URL and VC_REDIST_SHA256 in this script."
        )
    return path


def file_version(path: str) -> tuple[int, int, int, int]:
    import pefile  # noqa: PLC0415

    pe = pefile.PE(path, fast_load=True)
    try:
        pe.parse_data_directories(directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_RESOURCE"]])
        info = pe.VS_FIXEDFILEINFO[0]
        ms, ls = info.FileVersionMS, info.FileVersionLS
        return ms >> 16, ms & 0xFFFF, ls >> 16, ls & 0xFFFF
    finally:
        pe.close()


def visual_studio_present() -> bool:
    """Whether vswhere finds any Visual Studio product, Build Tools included."""
    vswhere = os.path.join(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"),
                           "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if not os.path.isfile(vswhere):
        return False
    try:
        found = subprocess.run([vswhere, "-all", "-products", "*", "-property", "instanceId"],
                               capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.TimeoutExpired):
        return False
    return bool(found.strip())


def installer(iscc: str, version: str) -> str:
    redist = vc_redist()
    major, minor, bld, _ = file_version(redist)
    say(f"  carrying vc_redist.x64.exe {major}.{minor}.{bld}, run only where the installed runtime is older")
    if not visual_studio_present():
        say("  no Visual Studio or Build Tools on this machine: Microsoft licenses the redistributable "
            "to Visual Studio users, so install Build Tools before publishing this installer")
    base = f"{NAME}-{version}-win64-setup"
    args = [iscc, "/Q", f"/DAppVersion={version}", f"/DSourceDir={FOLDER}", f"/DOutputDir={DIST}",
            f"/DOutputBaseFilename={base}", f"/DVcRedist={redist}", f"/DVcRedistMajor={major}",
            f"/DVcRedistMinor={minor}", f"/DVcRedistBld={bld}", ISS]
    say("$ " + " ".join(args))
    if subprocess.run(args).returncode:
        raise BuildError("ISCC failed; its output is above")
    return os.path.join(DIST, base + ".exe")


# -- entry point ----------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="uv run tools/build_exe.py", description=__doc__.split("\n\n")[0])
    parser.add_argument("--sdk-zip", metavar="ZIP",
                        help="your own Waters SDK download: the frozen build installs it into a scratch "
                             "profile, which is deleted afterwards, and reads $SCHAMP_SMOKE_RAW if set")
    parser.add_argument("--no-installer", action="store_true", help="stop after the portable zip")
    parser.add_argument("--no-prune", action="store_true",
                        help="keep all the Qt PyInstaller collects (for measuring the pruning)")
    args = parser.parse_args(argv)

    try:
        preconditions()
        build_stamp = stamp()
        say(f"schamp-acquire {build_stamp['version']} from {build_stamp['commit']}, Python "
            f"{build_stamp['python']}, PyInstaller {build_stamp['pyinstaller']}, PySide6 {build_stamp['pyside6']}")
        if build_stamp["commit"].endswith("+dirty"):
            say("  the working tree has uncommitted changes: fine for a trial, not for a release")
        freeze(build_stamp, prune=not args.no_prune)
        missing = add_notices(build_stamp)
        size, count = folder_size(FOLDER)
        say(f"\n{FOLDER}: {size / 1048576:.1f} MB in {count} files"
            + (" (unpruned)" if args.no_prune else ""))
        if missing:
            say("  no licence file found for: " + ", ".join(missing) + " (see licenses/*/METADATA.txt)")

        screen_names(FOLDER)
        say("  no SDK-derived names in the folder")
        unresolved = screen_pe(FOLDER)
        if unresolved:
            raise BuildError("DLLs imported but neither in the folder nor in Windows: " + "; ".join(unresolved))
        say("  every DLL import resolves, in the folder or in System32")
        crt = crt_inventory(FOLDER)
        say("  Visual C++ runtime carried beside the program: " + (", ".join(crt) or "none"))

        with tempfile.TemporaryDirectory(prefix="schamp-frozen-") as local:
            result = frozen_check(local)
            if result["frozen"] is not True or not result["gui"]["ok"]:
                raise BuildError(f"the frozen application did not build its window: {json.dumps(result, indent=2)}")
            if result["import"]["ok"]:
                raise BuildError(f"the frozen application found an SDK in a fresh profile: {result['import']}")
            say(f"  frozen, fresh profile: the run sheet builds ({result['gui']['rows']} rows, Qt "
                f"{result['gui']['qt']}), build {result['build']['commit']}; no SDK, as it should be")
        if args.sdk_zip:
            with tempfile.TemporaryDirectory(prefix="schamp-frozen-sdk-") as local:
                result = frozen_check(local, install_zip=os.path.abspath(args.sdk_zip))
                if not result["ok"]:
                    raise BuildError(f"the frozen application did not install and use the SDK: "
                                     f"{json.dumps(result, indent=2)}")
                read = result["raw"]
                say(f"  frozen, SDK installed from the download: masslynxsdk {result['import']['version']} "
                    f"imports from the profile, key from {os.path.relpath(result['license']['source'], local)}"
                    + (f"; {os.path.basename(read['path'])}: ATD fitted at {read['centre_ms']:.3f} ms"
                       if read else "; set SCHAMP_SMOKE_RAW to read an acquisition too"))
            # The scratch profile, SDK and all, went with the temporary directory.

        archive = portable_zip(build_stamp["version"])
        say(f"\nportable  {archive}  ({os.path.getsize(archive) / 1048576:.1f} MB)")
        if args.no_installer:
            return 0
        iscc = find_iscc()
        if iscc is None:
            say("installer skipped: Inno Setup 6's ISCC.exe is not on PATH or in its usual places")
            return 0
        setup = installer(iscc, build_stamp["version"])
        say(f"installer {setup}  ({os.path.getsize(setup) / 1048576:.1f} MB)")
        screen_names(DIST)
    except BuildError as exc:
        say(f"\nBuild failed. {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

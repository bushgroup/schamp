"""Install the Waters MassLynx SDK into this project, from the user's own download.

    uv run tools/install_sdk.py <path-to-667007504DDRevA.zip>

Verifies the download against the MD5 list Waters ships inside it, unpacks it into
gitignored `external/masslynxsdk/`, and puts that tree on the project environment's import
path. Windows only. Idempotent: run it again to repair or upgrade an install.

Nothing this script writes is ever tracked by git, and nothing is redistributed: the SDK is
licensed to the user, not to this project (lab record, task 01 reads the EULA). Files are
placed, never edited.

The work is `schamp.sdk_install`, which the installed application's first-launch wizard
calls too; this is its command line, and its docstring has the layout it creates.

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
import os
import sys

import schamp
from schamp.sdk_install import InstallError, confirm, install, read_manifest, register

TARGET = schamp.SDK_DIR


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

    if args.relink:
        manifest = read_manifest(TARGET)
        if manifest is None:
            print(f"Nothing is unpacked at {TARGET}. Run this with the Waters download "
                  "first.", file=sys.stderr)
            return 2
        try:
            manifest["pth"] = register(manifest["package_dir"])
            manifest["version"], manifest["dll"] = confirm(manifest["package_dir"])
        except InstallError as exc:
            print(f"\nInstall failed. {exc}", file=sys.stderr)
            return 1
        return _report(manifest)

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
        manifest = install(args.download, target=TARGET, link_environment=True, log=print)
    except InstallError as exc:
        print(f"\nInstall failed. {exc}", file=sys.stderr)
        return 1
    return _report(manifest)


def _report(manifest: dict) -> int:
    """Say what is installed and where."""
    print()
    print(f"masslynxsdk {manifest['version']} is installed and imports in this environment.")
    print(f"  package      {manifest['package_dir']}")
    print(f"  native DLL   {manifest['dll']}")
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

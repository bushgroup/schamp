# schamp

schamp is software and documentation for acquiring and analyzing ion mobility mass spectrometry
data on Waters Synapt G2, G2-S, and G2-Si instruments in which the traveling-wave ion mobility
cell has been replaced by an RF-confining drift cell. The drift cell applies a uniform DC field
along its axis while RF confines ions radially, so drift times measured as a function of drift
voltage yield absolute mobilities and collision cross sections without calibration against
reference ions. The instrument and its first measurements are described in Allen, Giles, Gilbert,
and Bush, *Analyst* 2016, **141**, 884 (`docs/`, open access, CC-BY 3.0). The name is from the
Mason-Schamp equation, which converts a reduced mobility into a collision cross section.

## Status

This repository is at its initial commit. The package is a stub, the workflows below are planned,
and the development record that sequences them lives in a private companion repository. The first
release will reproduce the polyalanine measurements in Figure 1 of the 2016 paper from the
original acquisitions.

## What it will do

- Read arrival-time distributions for chosen m/z windows directly from Waters `.raw` acquisitions
  through the vendor MassLynx SDK, replacing the `cdctest.exe` step of the original workflow.
- Fit each arrival-time distribution and report its centroid, width, and fit quality.
- Regress drift time against reciprocal drift voltage across an acquisition series, then report
  mobility, reduced mobility, and collision cross section in helium or nitrogen with propagated
  uncertainties, using a per-instrument profile for the drift-voltage definition and cell length.
- Record the conditions the `.raw` file does not carry, i.e., pressure, temperature, and gas, in
  one conditions table per experiment.
- Document how to acquire a drift-voltage series on an RF-confining drift cell instrument.

## Requirements

- Windows. The instruments and MassLynx run Windows, and the SDK's native library is tested here on
  Windows only.
- [uv](https://docs.astral.sh/uv/). The interpreter (CPython 3.12) and every library version are
  pinned by `pyproject.toml`, `.python-version`, and the committed `uv.lock`.
- The Waters MassLynx SDK (v5.0.0), obtained from Waters under its EULA. It is not on PyPI and
  this repository cannot redistribute it. You need the download zip and the `license.key` that
  came with it. Everything downstream of data extraction runs without the SDK.

## Fresh clone

```
git clone https://github.com/bushgroup/schamp
cd schamp
git config core.hooksPath .githooks
uv sync                        # creates .venv/ from uv.lock
uv run tools/check_public.py   # SDK-free self-check
uv run tools/bootstrap.py      # clones the reference checkout into external/
```

`uv run <script>` needs no flags and no venv activation. Use `uv run` or activate `.venv/` rather
than a bare `python` on `PATH`.

`git config core.hooksPath .githooks` is per-clone and not tracked by git. Without it, commits do
not receive the `Assisted-by:` trailer rewrite, and nothing refuses an oversized file or a `.raw`
directory entering history.

## Installing the Waters SDK

To read `.raw` acquisitions, first obtain the MassLynx SDK from Waters. Waters distributes it as a
single zip, document number 667007504, containing the `masslynxsdk` wheel, the native
`MassLynxRaw.dll`, the C++ headers, the API help, and the `license.key` issued to you. Then point
the installer at that zip:

```
uv run tools/install_sdk.py path/to/667007504DDRevA.zip
```

The installer verifies the download against the MD5 list Waters ships inside it, unpacks it into
`external/masslynxsdk/`, puts the package on the environment's import path, and loads the native
library to confirm it works. It reports where the key landed. Run it again at any time to repair
or upgrade an install, and run it with `--relink` to attach an SDK that is already unpacked to a
second environment.

Nothing needs to be configured after that. schamp reads the key from
`external/masslynxsdk/license.key`, which is where the installer puts it. To use a key held
elsewhere, set `SCHAMP_LICENSE_PATH`, or write a `schamp.ini` in the working directory:

```ini
[license]
path = C:\path\to\license.key
```

The full precedence, first hit winning, is the `license_key` and `license_path` arguments, then
`SCHAMP_LICENSE_KEY` and `SCHAMP_LICENSE_PATH`, then the `[license]` section of `schamp.ini`, then
`./license.key`, then the installed copy. Verify an install with `uv run tools/check_public.py`,
which reports the SDK, the key, and, when `SCHAMP_SMOKE_RAW` names an acquisition, a real read
from it.

The SDK is licensed to you by Waters, not by this project, and its EULA does not permit
redistribution. Everything the installer writes goes into `external/`, which is gitignored;
`.gitignore` additionally excludes the wheel, the DLL, and every `*.key` file by name. Every lab
installs from its own download. Note that a `uv sync` following a Python version change rebuilds
`.venv/` and drops the import path entry. Rerun the installer with `--relink` to restore it.

## Layout

```
src/schamp/      the package; sdk.py opens licensed SDK readers on a .raw, the rest is a stub
tools/           install_sdk.py (the Waters SDK), bootstrap.py (fetches external/),
                 check_public.py (the self-check)
docs/            the 2016 Analyst paper and its supporting information, with an index
external/        the Waters SDK installed by install_sdk.py, and reference code fetched by
                 bootstrap.py and read in place, never imported or vendored
                 (not in this repository)
```

## Reference code

`bushgroup/waters2python` is the Bush lab's wrapper around the same SDK, and the most complete set
of worked examples of calling it from Python. `tools/bootstrap.py` clones it at a pinned commit
into `external/` for reading. schamp calls `masslynxsdk` directly and does not depend on it.

## Development record

The task history, working notes, exploration scripts, and the legacy analysis this project
modernizes live in a private companion repository. When it is present as a sibling checkout, or is
named by the `SCHAMP_LAB` environment variable, `schamp.lab_dir()` resolves it; without it the
package is fully functional.

This repository is private and carries no license file; licensing will be resolved before any
visibility change.

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

The pipeline runs end to end, from a set of `.raw` acquisitions to collision cross sections, and
it reproduces the polyalanine measurements in Figure 1 of the 2016 paper from the original
acquisitions. The worked example in [`examples/`](examples/README.md) is that reproduction,
figure and all. There is no tagged release yet, the API is still moving, and the development
record that sequences the work lives in a private companion repository.

## What it does

- Reads arrival-time distributions for chosen m/z windows directly from Waters `.raw`
  acquisitions through the vendor MassLynx SDK, replacing the `cdctest.exe` step of the original
  workflow.
- Places each window on the ion's own measured peak, in the m/z frame the acquisition's own
  calibration puts the mobility reader in, rather than on the m/z its formula predicts.
- Fits each arrival-time distribution and reports its centroid, width, and fit quality.
- Regresses drift time against reciprocal drift voltage across an acquisition series, then
  reports mobility, reduced mobility, and collision cross section in helium or nitrogen with
  propagated uncertainties, using a per-instrument profile for the drift-voltage definition and
  cell length.
- Records the conditions the `.raw` file does not carry, i.e., pressure, temperature, and gas, in
  one conditions table per experiment.
- Documents how to acquire a drift-voltage series on an RF-confining drift cell instrument, in
  [`docs/acquiring-a-drift-voltage-series.md`](docs/acquiring-a-drift-voltage-series.md).

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
```

That is the whole setup. The self-check passes with no Waters SDK, no license key, and no data:
the checks that need any of those report themselves skipped rather than failed.

`uv run <script>` needs no flags and no venv activation. Use `uv run` or activate `.venv/` rather
than a bare `python` on `PATH`.

`tools/bootstrap.py` is the one script that is not for you. It clones a private Bush lab
repository of MassLynx SDK examples into `external/`, as reading material for work on schamp
itself, and it fails with an authentication error for anyone without access. Nothing in the
package depends on it.

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
src/schamp/      the package: extraction through sdk.py, then fitting, regression,
                 cross sections, and the two config formats
tools/           install_sdk.py (the Waters SDK), bootstrap.py (fetches external/),
                 check_public.py (the self-check)
docs/            the acquisition guide, the 2016 Analyst paper and its supporting
                 information, with an index
external/        the Waters SDK installed by install_sdk.py, and reference code fetched by
                 bootstrap.py and read in place, never imported or vendored
                 (not in this repository)
```

## Development record

The task history, working notes, exploration scripts, and the legacy analysis this project
modernizes live in a private companion repository. When it is present as a sibling checkout, or is
named by the `SCHAMP_LAB` environment variable, `schamp.lab_dir()` resolves it; without it the
package is fully functional.

## License

BSD 3-Clause, in [`LICENSE`](LICENSE). The two papers in `docs/` are not covered by it: they are
open access under CC-BY 3.0 and carry their own terms, and `docs/README.md` says so. Neither
covers the Waters MassLynx SDK, which Waters licenses to you directly and which this repository
never redistributes.

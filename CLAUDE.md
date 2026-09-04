# schamp — public-role code repo

Shared-core package `schamp`: software and documentation for acquiring and analysing ion mobility
mass spectrometry data on Waters Synapt G2 / G2-S / G2-Si instruments fitted with an RF-confining
drift cell. Named for the Mason–Schamp equation. This is the **public-role half of a two-repo
arrangement**: it carries `src/schamp/`, `tools/`, `docs/` (the 2016 Analyst paper + SI, CC-BY),
`examples/` and the public self-check — nothing else. The development record lives in the private
sibling **`bushgroup/schamp-lab`** (tasks, notes, explorations, the legacy analysis, the SDK zip,
and the working CLAUDE.md with the full project state and never-do list).

**Sessions are based here.** A gitignored `CLAUDE.local.md` imports the lab repo's CLAUDE.md; if
you are reading this file *without* that import, you have a public clone — the package and
`tools/check_public.py` are fully usable, the lab material is not yours to see.

## Rules that live with the code

- **Run anything with `uv run <path/to/script.py>`** — no flags, no venv activation. Never rely on
  the `python` on PATH.
- **`uv run tools/check_public.py`** after touching `src/schamp/` or `examples/`. Data-free and
  SDK-free: a clone with no Waters SDK, no license and no `.raw` must pass it. Checks that need the
  SDK report SKIPPED, never FAIL, when it is absent.
- **The Waters MassLynx SDK is never committed** — not the wheel, not `MassLynxRaw.dll`, not
  `license.key`, not the zip, not the help files. Its EULA does not permit redistribution. Users
  install it from their own Waters download into gitignored `external/masslynxsdk/`
  (`tools/install_sdk.py`, lab record task 02); `.gitignore` also excludes the artefacts by name.
- **Never reimplement what the SDK does faithfully.** Mobillograms, drift-time reads, scan combining
  and CCS conversions come from `masslynxsdk`'s own readers; `schamp` adds the drift-cell physics
  and the workflow around them. Only what the SDK does not expose (the `_extern.inf` acquisition
  parameters) is parsed here, in one module.
- **`external/` is read-only reference material**: `bushgroup/waters2python` (fetched by
  `uv run tools/bootstrap.py` at a pinned commit) is the lab's own SDK wrapper and the best set of
  worked SDK examples; **it is never imported**, only read. Never vendor, edit or reformat anything
  under `external/`.
- **Windows only.** MassLynx and the instruments run Windows; the SDK's Linux build is untested here
  and not promised.
- **Lab-side paths resolve through `schamp.lab_dir()`**: `$SCHAMP_LAB`, then this root, then the
  sibling `../schamp-lab`. Code in `src/` refers to lab material only in opaque form ("lab record,
  task NN"), never by a path that only resolves lab-side.
- **Public commit messages are self-contained statements of the change.** Task IDs may appear as
  opaque references ("lab record, task 04") at most. Trailer is `Assisted-by: <model name>`, no
  email — never `Co-Authored-By:`. `.githooks/commit-msg` rewrites, `.githooks/pre-commit` rejects
  staged files over 5 MiB and any path inside a `.raw` directory; both need
  `git config core.hooksPath .githooks` once per clone.
- **`.gitattributes` pins `* text=auto eol=lf`.**
- **Outward-facing prose (README, `docs/`, `examples/`, user guides) follows the `manuscript-voice`
  skill** from the lab repo. Repo-internal prose (this file, docstrings, commit messages) does not.
- **No license while private**; licensing and the go-public date are open decisions recorded in the
  lab repo (task 09).

## Maintaining this file

This file stays lean: rules for working *in this repo*, nothing about the science or the project's
state. Those belong in the lab repo's CLAUDE.md and notes. Keep it under 60 lines.

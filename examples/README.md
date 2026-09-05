# Examples

`polyalanine-walkthrough.ipynb` runs `schamp`'s whole pipeline (extraction, fitting,
regression, cross sections) on the poly-DL-alanine drift-voltage series of Allen,
Giles, Gilbert & Bush, *Analyst* 2016, **141**, 884, and draws the paper's Figure 1: a
mobiligram heat map, drift time against reciprocal drift voltage per oligomer, and
collision cross section against m/z for all four ion series (1+, 2+, 1-, 2-).

## What is and is not here

`data/npalan/conditions.csv` and `data/palan/conditions.csv` carry the series' pressures
and temperatures, one row per acquisition, the same shape `schamp.experiment` reads
anywhere, and small enough to commit. To fill one in for a series of your own, start
from [`docs/acquiring-a-drift-voltage-series.md`](../docs/acquiring-a-drift-voltage-series.md),
which covers the acquisition end of this workflow: the instrument state, choosing the
drift voltages, and what to read off the manometer and the cell thermocouples while
the series runs. **The acquisitions themselves are not here.** The
series is two sets of 14 Waters `.raw` directories, about 5.4 GB, and this repository
never ships one: `docs/README.md` and the top-level `CLAUDE.md` say why.

Why: this is the public-role half of a two-repo project (top-level `CLAUDE.md`), and a
Waters `.raw` acquisition is a directory of binary files with nothing generic about it
to redistribute. These particular ones are 5.4 GB and belong to a private lab
checkout. To run the notebook against the real series, set `SCHAMP_EXAMPLE_RAW` to a
directory holding `130423_SJA_NPALAN_1_001.raw` .. `_014.raw` and
`130423_SJA_PALAN_1_001.raw` .. `_014.raw`:

```
SCHAMP_EXAMPLE_RAW=/path/to/130423 uv run --group dev jupyter notebook examples/polyalanine-walkthrough.ipynb
```

Without that series, which for most readers of this repository means without a Waters
license to run alongside it, the notebook is still worth reading: every cell is real
code against `schamp`'s public API, and the notebook's own last section says exactly
what four things to change to run it against a series of your own.

## Running it

Needs the `dev` dependency group (`uv sync --group dev`), which installs Jupyter and
nothing schamp itself depends on at runtime, and the Waters SDK
(`uv run tools/install_sdk.py`). `uv run tools/check_public.py` runs a lighter check
that needs neither: it confirms the shipped notebook carries no committed outputs, and
executes it in full whenever the `dev` group, the SDK and the acquisitions are all
present.

The notebook writes its tables and figure to `examples/out/`, which is gitignored:
nothing an execution produces is meant to be committed, including the notebook's own
outputs. If you have run it locally and want to hand the `.ipynb` back to the repository,
clear its outputs first (`jupyter nbconvert --clear-output --inplace`). A committed
execution artefact is exactly what `check_public.py`'s notebook check refuses.

## Where this differs from the published figure

This reproduces the published figure's substance, not its exact rendering. The heat map
omits the two dashed guide lines the original drew by eye, at slopes typed in by hand;
the colour map, the axis ranges and the marker choices are this notebook's own rather
than the original artwork's. The lab record has the cross sections themselves against
the published table: every ion series agrees to better than 0.6 % in the mean, and three
of the four to better than the rounding of the published integers.

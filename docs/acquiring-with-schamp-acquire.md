# Acquiring with schamp-acquire

To keep the record of a drift-voltage series while it is acquired, run `schamp-acquire`, the
run sheet, on an analysis PC beside the instrument. The run sheet lists the drift-voltage
settings of the series in the order they are to be acquired, takes the pressure for each
acquisition as you read it off the manometer, and writes the two files the analysis reads,
`experiment.toml` and `conditions.csv`, after every edit. When the acquisitions are visible to
it, the run sheet also checks each one against its plan and analyzes the series as it arrives.
This guide covers installing the run sheet on Windows, supplying the Waters MassLynx SDK,
pointing the run sheet at the data, a session from the first row to Finish series, and what to
do with the files it writes.

This is the software half of an acquisition. The instrument half, i.e., the instrument state a
linear-field measurement requires, why the voltages are stepped as they are, and what the
pressure and temperature readings must come from, is
[`acquiring-a-drift-voltage-series.md`](acquiring-a-drift-voltage-series.md), and this guide
does not repeat it. Read that guide first.

## Where it runs

The run sheet runs on a Windows 10 or Windows 11 PC (64-bit) that can see the acquisitions,
and never on the instrument's acquisition PC. The acquisition PC belongs to MassLynx and to
the instrument, and nothing from schamp is installed on it. The `.raw` acquisitions reach the
analysis PC in one of two ways, either over a network share while the series is acquired or by
copying them after it is finished, and the run sheet supports both (see "Where the acquisitions
come from" below). Making a folder on the acquisition PC visible over the network is a matter
for your lab's IT and is outside this guide.

The run sheet itself needs nothing from the instrument. Rows can be ticked and pressures typed
with no acquisition visible for the entire session, and the checks and the live analysis that
need the acquisitions start as soon as they appear.

## Installing

To install the run sheet, download one of the two Windows builds attached to the
[latest release](https://github.com/bushgroup/schamp/releases/latest):

- `schamp-acquire-<version>-win64-setup.exe`, the installer. It installs for the current user
  under `%LOCALAPPDATA%\Programs\schamp-acquire`, adds a Start-menu shortcut, and needs no
  administrator rights. On a PC whose Microsoft Visual C++ runtime is older than the one it
  carries, it offers to install that runtime, which does need an administrator; the option can
  be unticked, because the application carries the runtime libraries it uses in its own folder.
  Uninstall it from Windows' list of installed apps.
- `schamp-acquire-<version>-win64-portable.zip`, the same program unpacked. Unzip it into any
  folder you can write to and run `schamp-acquire.exe`. Nothing is installed, and deleting the
  folder removes it.

Both are the same build, so choose by how your PC is managed. The installer is the better
choice for a PC that one person uses, whereas the portable zip suits a PC on which you cannot
install software, or one where several versions should sit side by side.

To run the run sheet from a development checkout instead, see the [README](../README.md), which
uses `uv run schamp-acquire` and the same window.

## Supplying the Waters SDK

To read `.raw` acquisitions, the run sheet needs the Waters MassLynx SDK, which Waters licenses
to each lab under its own agreement and which schamp therefore does not include. Obtain it from
Waters: it arrives as one zip, Waters document number 667007504 (for example
`667007504DDRevA.zip`), together with a license key issued to your lab, which is either inside
the download or sent separately. The licence is between your lab and Waters, and its terms, not
schamp's, govern the SDK.

At the first start, the run sheet offers to install the SDK. Accept, or later choose
**Tools > Install the Waters SDK**, and the wizard asks for two things:

1. **The download**, i.e., the zip exactly as Waters supplied it. Do not unpack it first.
2. **The license key.** Leave the default if the key came inside the download. If Waters sent
   the key separately, choose that option and browse to the file.

The wizard then checks the download against the checksums Waters ships inside it, unpacks it
unaltered into `%LOCALAPPDATA%\schamp\masslynxsdk\`, and loads it once to prove it works. It
reports the SDK version and where the key it will use comes from. A key sent separately stays
where it is, and its location is recorded in `%USERPROFILE%\.schamp.ini`. Nothing is sent
anywhere, and the key's contents are never displayed.

The installed SDK belongs to your Waters licence rather than to the run sheet, so uninstalling
the run sheet or deleting the portable folder leaves it in place, and a later version of the
run sheet finds it without asking again. Delete that folder by hand when no copy of schamp on
the PC needs it.

To install the SDK without the window, for example on several PCs, run the program from its
folder (`%LOCALAPPDATA%\Programs\schamp-acquire\` for the installer) with

```
schamp-acquire --install-sdk path\to\667007504DDRevA.zip --self-check report.json
```

adding `--license-file path\to\license.key` for a key sent separately. The program has no
console, so the report is where the outcome goes: a JSON file saying whether the window builds,
whether the SDK loads, and where the key would come from (never the key itself).
`--self-check` alone asks the same of a copy that is already set up.

Until the SDK is installed, the Live tab says what it is waiting for and the rest of the run
sheet works as usual.

## Where the acquisitions come from

The **Data** field of the session panel names the folder in which the `.raw` acquisitions
appear. It can be a share (`\\acquisition-pc\data\my-series`), a mapped drive, or `.`, which
means the set's own folder. The run sheet looks in that folder every 4 s, and two routes lead
into it.

### Over a share, while the series is acquired

To follow a series as it is acquired, point **Data** at the share MassLynx writes into. Each
acquisition then appears in the folder while MassLynx is still writing it, and the run sheet has
to decide when it is finished before reading it. By default, an acquisition counts as finished
when `_extern.inf` and `_FUNC001.IDX` are both present and no file in its directory has been
added, removed, or changed size over three looks 4 s apart, i.e., after at least 8 s unchanged.
That rule is a default rather than a measurement of how MassLynx writes an acquisition, so it is
configuration: a `[watch]` table in the set's `runsheet.toml` changes any of its three terms,

```toml
[watch]
files = ["_extern.inf", "_FUNC001.IDX"]
polls = 3
interval_s = 4.0
```

and the Live tab states the rule in force. An acquisition that changes after it was called
finished is read again, so a rule that fires early costs a second read rather than a wrong
result. With the default rule, each of the fourteen 2013 anion acquisitions, copied one at a
time into a local folder standing in for a share, was called finished and read 10 to 13 s after
it landed.

A share that stops answering is shown as unreachable after 8 s, in the Checks panel and in the
`.raw` column, and the run sheet keeps working without it. Nothing typed into the sheet depends
on the share.

### By copying, after the series is acquired

To work without a share, keep the run sheet during the session exactly as you would with one,
ticking rows and typing pressures with **Data** left at `.`, and bring the acquisitions in
afterwards. The pressures exist only on the manometer's display while the series runs, so the
run sheet has to be open at the instrument either way; what the copy route gives up is the
planned-against-recorded check and the live analysis during the session.

To bring the acquisitions onto the analysis PC, use **File > Import acquisitions** rather than a
file manager, and choose the folder that holds the `.raw` directories, whether a share, a
mapped drive, or a USB stick. The import copies each acquisition into the set's folder, hashes
every file with SHA-256 as it is read and again after it is written, and moves an acquisition
into place only when every file agrees, so an interrupted copy never leaves half an acquisition
under its real name. It writes a `MANIFEST.sha256` that `sha256sum -c` checks without schamp,
adds a row to `conditions.csv` for each acquisition the table does not already hold, and
records the import in the set's `README.md`. An acquisition whose name matches a row on the
sheet joins that row, its typed cells untouched, and an acquisition matching no row is added at
the bottom of the sheet without a plan. A row you forgot to tick is therefore ticked by the
import, and its empty pressure is listed in the Checks panel until you type it in. On the 2013
anion series, the import copied and verified fourteen acquisitions (1.5 GB) between two local
disks in 11 s, and every typed pressure survived it.

Running the import again copies only what is new, skips an acquisition that is still being
written, and reports rather than overwrites one that has changed at the source since it was
copied, so it can also be run after each acquisition. If **Data** pointed elsewhere, the run
sheet then asks whether to read the acquisitions from the imported copies from now on. The
Live tab reads the copies as they land, once they pass the same completion rule as a share.

Importing into an empty folder also works, and gives a set with one taken row per acquisition,
no plan, and every pressure and the session temperature still to type. Use it for a series that
was acquired without the run sheet, with the readings from the lab notebook.

The same import is `uv run tools/import_set.py <source> <set>` in a development checkout.

## A session

### Starting a series

To start a series, choose **File > New run sheet**, pick an empty folder for the set, accept the
standard template, and give the base name your MassLynx sample list will use. The set's folder
is where the run sheet writes its files, and it is a folder of your own choosing on the analysis
PC, which need not be the folder the acquisitions appear in. The standard template is the lab's
fourteen (Helium Cell DC, IMS Bias) settings, in the order the 2013 polyalanine series was
acquired, which visits 74 to 352 V out of sequence so that a pressure drifting over a morning
does not track the drift voltage. A lab with an order of its own copies the shipped template
(`uw-standard-14.toml` in the package's `data/series/` folder) and chooses **Another file** in
the template dialog. Starting the program with a folder,
`schamp-acquire D:\data\my-series`, opens the set there or starts one from the standard template
if the folder holds none.

Then fill in the session panel on the left:

- **Title**, for `experiment.toml`, and **Operator**, kept in the run sheet only.
- **Profile**, the instrument profile whose drift-voltage formula and cell length the analysis
  uses. `uw-synapt-g2` is the instrument of the 2016 paper; for any other instrument, name your
  own profile file by its path relative to the set (the acquisition guide covers writing one).
- **T / K**, the session temperature, typed once.
- **Data**, the folder above.
- **Base name**. Acquisitions are named `<base>_NNN.raw`, numbered as a sample list numbers them,
  and **Rename rows to come** renames every row not yet ticked after a change.
- **Gas**, left at the profile's default unless the cell is running a different gas.
- **Held for the session**, the voltages the series does not step, as the template sets them.
  Change them to what your tune file holds.
- **Notes**, free text for `experiment.toml`.

Each edit is saved as it is made. The **Series** tab on the right plots the planned drift
voltage of every row, with open circles for rows still to come.

### Working down the sheet

Each row of the sheet is one acquisition: the settings to apply, the drift voltage they give
through the profile (**V nominal**), and the name MassLynx will give the file. For each row:

1. Set the instrument to the pair of voltages the row shows, and acquire.
2. Tick **Taken**.
3. Type the manometer reading into **P / Torr**.
4. Untick **Use** and say why in **Notes** for an acquisition that should be kept on record but
   left out of the regression, rather than removing the row.

Ticking **Taken** is what puts the acquisition into `conditions.csv`, and unticking it takes the
acquisition back out. A ticked row's settings and name are frozen, because they now describe a
file that exists. A ticked row whose pressure is still empty is shaded until it has one. Rows
can be inserted, removed, and reordered from the **Rows** menu, and a template can be laid over
the sheet from the same menu.

The **Checks** panel at the bottom runs the analysis's own checks on the record every time it
is saved, and states how many of the sheet's acquisitions are visible in the data folder.
Before any acquisition is visible, it checks the tables alone; after, it also opens every
visible acquisition's `_extern.inf` and checks it against the instrument profile, including the
mobility wave height and the pushes per drift bin. "No problems found, acquisitions included",
with no red cell in the **ΔV** column (below), is the state to leave the instrument in.

### Planned against recorded

Once an acquisition is visible, the run sheet reads the voltages it recorded and fills
**V recorded** and **ΔV**, the difference from the plan. MassLynx records each setpoint as typed,
so a knob that was set as planned gives a ΔV of zero, and any difference larger than 0.05 V
turns both cells red. The tooltip on a red cell names the knob, its planned value, and its
recorded value. The **Series** tab marks each recorded voltage with a cross beside its planned
circle.

A red row is caught while the instrument is still available. Note that the analysis always uses
the voltage the acquisition recorded, never the planned one, so a red row is a series that did
not visit the voltage you meant, not a wrong number in the analysis.

## The Live tab

To follow the analysis while the series arrives, open the **Live** tab and name the species to
watch. **Add** takes a label, an m/z, and a signed charge. **Pick from spectrum** opens the total
spectrum of the first finished acquisition; zoom with the toolbar, click a peak, set its charge,
and add it, and the dialog stays open for the next. Species are kept in `runsheet.toml`, so
reopening the set brings them back, and adding one later reads it out of every acquisition
already finished.

As each acquisition finishes, the tab fits the arrival-time distribution of every species in it
and updates the table of species, which gives each one's points in the regression, R², and
cross section. Selecting a species draws three plots:

- **The latest arrival-time distribution** with its single-Gaussian fit, titled with the
  acquisition and the fitted centroid.
- **Arrival time against reciprocal drift voltage**, one point per finished acquisition, with
  the fitted line once three points are in the regression, titled with R², the intercept t₀, and
  the number of points. Open circles are acquisitions read but left out of the regression,
  because **Use** is unticked or the acquisition is not yet ticked as taken.
- **The cross section as it ran**, recomputed as each point joined the regression, from the
  third point on, titled with its gas, charge, pressure, and temperature.

A point ringed in red on the middle plot is off the line, i.e., its residual from the fitted
line is more than three times the rms residual of the other points in the regression, and more
than 1 µs. Every point is judged this way, whether it is in the regression or not, once at least
three others are. Ringing a point changes nothing in the analysis. It flags an acquisition to
look at while the instrument is still available: a knob that did not move, a pressure that was
mistyped, or a voltage too low for the ion. On the 2013 anion series, with the published
analysis's four exclusions entered as **Use** unticked, every ringed point is one of those four
acquisitions, i.e., the four below 100 V. The 1− ions of 10 and 15 residues ring all four, the
1− ion of 6 residues rings one, the 2− ion of 20 residues rings none, and no acquisition the
published analysis kept is ringed for any of them. The decision to exclude one is yours, and it
is made by unticking **Use** with a reason in **Notes**, which moves the point to an open circle
and recomputes the line.

The status lines at the top of the tab report how many acquisitions are finished, settling, or
still being written, the completion rule in force, and what the first acquisition's own record
says about its mass calibration, i.e., which calibrant it used, how old the calibration is, and
whether it was carried across polarity (the acquisition guide explains why this matters). If
`experiment.toml` asks for a mass recalibration, the tab says that the live view does not apply
it.

## Finishing the series

To close the series, choose **Finish series** on the Live tab. It writes two files beside
`experiment.toml`:

- `ccs.csv`, one row per species with a cross section: charge, m/z, gas, the mean pressure and
  temperature of the acquisitions in its regression, the number of points, R², the slope and
  intercept, the mobility and reduced mobility, the cross section with its uncertainty, the
  terms that uncertainty propagates, and the acquisitions excluded.
- `results.json`, everything behind those rows: every point, fit, and regression, each
  acquisition's calibration record, the completion rule in force, and what was excluded and why.

Finish series can be chosen as often as you like, and each time it rewrites both files from the
record as it stands. On the 2013 anion series, the cross sections it wrote for four test ions
(1− ions of 6, 10, and 15 residues and the 2− ion of 20) are the worked example's, to the
0.01 Å² both print.

## Where the files go next

The set's folder is the experiment, and every file in it has a reader:

| File | Written by | Read by |
|---|---|---|
| `experiment.toml` | the run sheet, the import | the analysis (`schamp.load_experiment`) |
| `conditions.csv` | the run sheet, the import | the analysis |
| `runsheet.toml` | the run sheet | the run sheet only: the plan, the species, the `[watch]` rule |
| `MANIFEST.sha256`, `README.md` | the import | `sha256sum -c`, and you |
| `ccs.csv`, `results.json` | Finish series | you, a spreadsheet, and any later comparison |
| `*.raw` | MassLynx, the import | the run sheet and the analysis, through the Waters SDK |

`experiment.toml` and `conditions.csv` are the record of the experiment, and the analysis reads no
other record beside the acquisitions themselves. They are plain text, edited by hand where needed in the format the acquisition
guide describes, and they travel with the acquisitions. To reanalyze a series with other
species, other windows, or a mass recalibration, load the set with `schamp.load_experiment`
and follow the worked example in [`examples/`](../examples/README.md), whose calls are the same
ones the Live tab makes.

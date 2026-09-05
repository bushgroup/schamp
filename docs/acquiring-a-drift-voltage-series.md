# Acquiring a drift-voltage series on an RF-confining drift cell

To measure an absolute collision cross section on this instrument, acquire the same sample
several times over and change nothing between acquisitions but the drift voltage. Arrival time
is then linear in reciprocal drift voltage, the slope of that line gives the mobility, and its
intercept absorbs the time the ion spends between the drift cell exit and the time-of-flight
analyzer. Those acquisitions together are a drift-voltage series, and a series is the unit
`schamp` analyzes. This guide covers the instrument state a series requires, how to choose the
voltages, what has to be written down by hand because the `.raw` file does not record it, and
the failure modes that make a series unusable.

The instrument is a Waters Synapt G2, G2-S, or G2-Si in which the traveling-wave ion mobility
cell has been replaced by an RF-confining drift cell: 168 ring electrodes spanning a 25.05 cm
drift region, two resistor chains that establish a uniform DC field along the axis, and two
capacitor chains that couple RF to every electrode to confine ions radially. Allen, Giles,
Gilbert, and Bush, *Analyst* 2016, **141**, 884 describes the device and its first measurements,
and `docs/` carries that paper and its supporting information. Waters' own guide to the cell is
the authority for installing and servicing it, and this guide does not repeat it. Where a step
below is a service-level change to the instrument, get the procedure from Waters.

## Prerequisites

1. The cell installed, with its gas manifold and its own absolute-pressure capacitance
   manometer. The published measurements used an MKS Baratron Type 626B on the tee that admits
   the drift gas at the center of the cell (ESI Fig. S3), reading the pressure inside the cell
   directly. The instrument's own vacuum readbacks are not the drift-cell pressure and are not a
   substitute for this gauge.
2. A temperature sensor on the cell. See "Recording the conditions" below for what the published
   measurements used.
3. The mobility cell's traveling wave disabled, both halves of it. See "Instrument state" below.
4. An instrument profile for your instrument. `schamp` ships one, `uw-synapt-g2`, which is the
   instrument of the 2016 paper and the worked example of the format. It carries the drift
   length, the resistor-ladder counts that set the divider ratio, the `_extern.inf` key names
   that enter the drift voltage, and the checks run against every acquisition. Copy it, edit it,
   and load it by path:

   ```python
   profile = schamp.load_profile("my-lab-g2s.toml")
   ```

   Do this before acquiring rather than after. The key names MassLynx writes differ between
   instruments, and a profile whose `high_keys` do not exist in your `_extern.inf` raises rather
   than guessing, which is a five-minute fix at the desk and a repeated acquisition at the
   instrument.

## Instrument state

**Turn the mobility cell's traveling wave off.** This has two parts, and only one of them leaves
a trace in the data. Setting the IMS wave amplitude to 0 V on the tune page is the visible part,
and it is recorded as `IMS Wave Height (V)` in every acquisition; `schamp`'s profile check
`require_wave_off` refuses an acquisition whose height is nonzero, because a wave-on acquisition
is not a linear-field measurement at all. Disabling the cell's periodic flushing pulses is the
other part, it is a service-level change rather than a tune-page setting, and nothing in the
`.raw` file records whether it was done. Get that procedure from Waters, do it once when the
cell is installed, and confirm it before a campaign rather than trusting it. An otherwise
perfect series acquired with the flushing pulses still enabled looks exactly like a good series
and gives wrong mobilities.

**Set the RF to the amplitude and frequency Waters specifies for the cell.** The published
measurements used 100 V peak to peak at 2.8 MHz (paper, Instrumentation).

**Let the drift gas settle before the first acquisition.** Cross section goes as the reciprocal
of pressure, so a pressure that is still falling while the series runs biases every cross
section in it, and biases them unequally. Waters' set-up guide states a settling time; the
practical test is the manometer itself. Read it at intervals after the gas is on and start when
successive readings agree to the precision you intend to quote, which for the published series
was 0.001 Torr on a nominal 2.03 Torr.

**Minimize the voltages that can activate the ion.** An ion that unfolds in the source and an
ion that unfolds in the trap both give a clean straight line and the wrong cross section. The
published conditions kept the voltages as low as the signal allowed, and for proteins and
protein complexes above 60 kDa they raised the source backing pressure to 4 to 6 mbar, raised
the trap cell flow rate to 10 mL/min, and set the transfer cell exit to 0 V (paper, Experimental
methods).

**Leave the transfer cell's traveling wave on.** Only the mobility cell's wave has to be off.
The published conditions set the transfer cell wave height to 2 to 4 V and adjusted its velocity
between 60 and 200 m/s to maximize transmission while retaining the separation the drift cell
produced (paper, Experimental methods).

**Set one pusher per drift bin.** `ADC Pushes Per IMS Increment` must be 1, so that one drift
bin is exactly one pusher period; the profile's `expect_pushes_per_bin` checks it. If several
pushes are summed into a bin, the drift-time axis is no longer the pusher period and everything
downstream of it is misscaled.

**Check that the slowest ion fits inside the drift-time window.** The acquisition mass range
sets the pusher period, the pusher period times the number of drift bins is the widest arrival
time the acquisition can record, and the slowest ion at the *lowest* drift voltage of the series
is the one that has to fit. The published anion series ran 200 to 2000 Th, which gave a pusher
period of 0.0693 ms and, over 200 drift bins, a 13.9 ms window against a slowest arrival of
9.2 ms; the cation series ran 200 to 1200 Th, which gave 0.0543 ms and a 10.9 ms window against
a slowest arrival of 7.4 ms. Both had headroom, and neither had much. Compute this before
acquiring, because an ion whose arrival time exceeds the window does not fail loudly.

## Calibrating, so that a mass recalibration is possible

Calibrate the instrument in the session, in the polarity you are about to acquire. Nothing else
in this guide buys as much: the drift-cell measurement itself is calibration-free, but placing
an m/z window on the species you want is not, and a mass calibration carried over from another
day puts the species somewhere other than where its formula says it is. On the 2013 polyalanine
series both polarities carry a calibration fitted 92 days earlier, and the anions carry a
positive-mode calibration applied to negative-mode data. They read 120 to 180 ppm and −350 to
+155 ppm out. A CsI calibrant acquired that same morning reads within 4 ppm.

Then acquire a calibrant run in each polarity you will analyze, in the same session, over a mass
range that covers the analyte's. This is the fallback that works when the analyte has no species
of computable m/z to anchor on, and `schamp.calibration.MzFrame.from_calibrant` takes both the
calibrant's own calibration function and the residual fitted from it. Two points about it are
worth stating plainly. A calibrant does not cross polarity, and one that was asked to carry a
fit across polarities placed 8 % of the 2013 anion series outside their m/z windows at every
m/z. And a calibrant's useful mass range is only what it actually produces; the 2013 CsI file
was acquired over 200 to 1200 Th and gave four clusters, which is enough for the order-1
correction these instruments need but leaves nothing above 1172 Th.

Whatever you do about calibration, run the check that costs nothing:

```python
from schamp.calibration import CalibrationProvenance

prov = CalibrationProvenance.from_raw("130423_SJA_NPALAN_1_001.raw")
print(prov.describe())
for concern in prov.concerns():
    print(concern)
```

`CalibrationProvenance` reads `_HEADER.TXT` and `_extern.inf` alone, so it needs no SDK, no
license, and no scan read, and `concerns()` names what the file itself says is wrong with its
calibration: no calibration recorded, a calibration carried across polarity, or one older than
the age you asked with. An empty tuple means the file gives no reason to doubt its masses. Run
it on the first acquisition of a series before acquiring the other thirteen.

## Choosing the voltage series

Three properties of the series matter: how many voltages, over what range, and how they are
spaced in reciprocal drift voltage rather than in drift voltage.

**How many.** `schamp` requires three usable acquisitions and reports a series with fewer as
invalid, because two points determine a slope and an intercept with nothing left over to say
whether the relation is linear. The published analysis of record used ten. Ten is a good target:
it costs an hour at the instrument and it lets four points be dropped, as the published analysis
dropped four, without the fit becoming a two-point extrapolation.

**Over what range.** The published series used drift voltages from 104 to 354 V (paper,
Instrumentation), a factor of 3.4. The top of the range is set by what the ladder and its supply
will hold and by the low-field condition. The bottom is set by the drift-time window above, and
in practice by the data: on the 2013 series the four lowest voltages, all below 104 V, lie off
the line by more than the rest of the series scatters, and the analysis of record drops all
four. Acquire below the range you expect to keep, so that dropping a point is a decision about a
point you have rather than about a gap.

**How they are spaced.** The regression abscissa is 1/V, so voltages spaced evenly in V bunch
together at the top of the series and contribute little to the fit's leverage there. Space them
geometrically instead. The ten voltages of the published analysis step by 10 to 21 % each, which
puts them evenly in 1/V to within a factor of two, from 2.84 × 10⁻³ to 9.64 × 10⁻³ V⁻¹.

**Give each acquisition its own voltage.** `Experiment.validate` reports any two used
acquisitions that share a drift voltage. A repeat at the same voltage is a useful check on
reproducibility and it adds no leverage to the fit, so mark one of the pair `use = false` and
say why in the table.

**The intercept is why one voltage is not enough.** Fitted over the ten voltages of the
published series, the transport time from the cell exit to the analyzer runs from 0.13 to
0.62 ms across the 73 polyalanine species, and at the highest drift voltage of the series that
is 5 to 46 % of the whole arrival time, with a median of 21 %. A single-voltage measurement that
treats arrival time as drift time is wrong by that much, in that direction, and the series is
what separates the two.

**Step the voltage with the two knobs the profile reads.** On the reference instrument the drift
voltage is

```
V_drift = (VappH − VappL) · (1 − 2/170)
VappH   = Helium Cell DC + IMS Bias + Helium Exit        (Helium Exit is negative)
VappL   = Transfer DC Entrance
```

which is what ESI Fig. S2 prints. The 2/170 is there because neither plate that bounds the drift
region sits at an end of the ladder: the ladder is 170 equal steps, and one of them lies outside
each end of the drift region. The published series is a worked example of stepping it, raising
IMS Bias alone from 119 to 240 V and then holding the bias at 240 V while raising Helium Cell DC
through 18, 56, 110, and 160 V. Never type a drift voltage into your analysis. `schamp` reads
these keys out of each acquisition's own `_extern.inf` and applies the profile's formula, so the
voltage that is analyzed is the voltage the instrument recorded.

Note that a G2-S or G2-Si may differ in the divider ratio, in the drift length, or in both. They
are separate profile fields, and neither is derived from the other.

## Recording the conditions

Pressure, temperature, and gas identity are not in the `.raw` file, and cross section depends on
all three. Write them down at the instrument.

**Pressure, per acquisition, from the manometer's own controller.** Read the display on the
capacitance manometer's controller and type the value into the conditions table row for that
acquisition. This is not bookkeeping: over the 2013 series the reading moved from 2.026 to
2.029 Torr, which is small and is not zero, and Ω goes as 1/P.

**Temperature, once per session, from a sensor on the cell.** The published measurements used
two T-type thermocouples, one on the foot where the cell meets the bottom of the vacuum chamber
and one on a screw holding the electrode board to the aluminum top plate, and took the average
of the two, which read closely together. Mounting them on cell metal rather than in the gas
stream is the right choice here because the drift gas passes through the chamber and the foot
before it enters the cell, so it arrives thermalized to the metal those sensors measure. The
cell sits above ambient, and that difference is attributable to heat from the nearby
turbomolecular pumps rather than to anything the experiment does, which is why one reading per
session is enough. Quote about 1 K, the standard limit of error for a T-type thermocouple, and
note what that buys: Ω goes as T^½, so 1 K on 300 K is 0.17 % on every cross section in the
series.

**Gas.** Helium and nitrogen are what `schamp` supports. The gas is a field in the instrument
profile, with a per-acquisition override in the conditions table, and it must be set explicitly
for anything other than the gas the instrument normally runs. A cross section without its gas
named is not a number anyone can use.

## Naming the acquisitions

Acquire the whole series into one MassLynx sample list so that the files are numbered
sequentially, and put the sample and the polarity in the base name. The published series is
`130423_SJA_NPALAN_1_001.raw` through `_014.raw` for the anions and `130423_SJA_PALAN_1_001.raw`
through `_014.raw` for the cations: date, operator, sample, replicate, acquisition number.

Nothing in `schamp` parses these names. The acquisition column of the conditions table matches
the directory name literally, and the drift voltage comes from inside the file, so a name is only
ever a label for a human reading the table. Sequential numbering earns its keep in one place: it
makes the row order of the table the order the acquisitions were taken in, which is what lets a
drifting pressure be seen as a trend rather than as scatter.

## Filling the conditions table

An experiment is two files, `experiment.toml` for the scalars and `conditions.csv` for the
per-acquisition table:

```toml
schema = 1
title = "130423 poly-DL-alanine anions (ES-)"
profile = "uw-synapt-g2"
gas = "helium"
conditions = "conditions.csv"
data_dir = "."

[defaults]
temperature_K = 301.13
```

```
# Pressures are the MKS Baratron 626B readings taken at the time of acquisition.
# Temperature is the session value, in experiment.toml's [defaults].
acquisition,pressure_Torr,use,notes
130423_SJA_NPALAN_1_001.raw,2.026,false,dropped: one of the four lowest drift voltages
130423_SJA_NPALAN_1_002.raw,2.026,true,
130423_SJA_NPALAN_1_003.raw,2.026,true,
```

`examples/data/npalan/conditions.csv` and `examples/data/palan/conditions.csv` are the two
tables of the published series, complete, and they are the shape to copy. They carry the
temperature as a per-row column rather than in `[defaults]`, which reads the same and is the
older form; `[defaults]` is the better place for a value that is the same on every row, and one
session temperature is exactly that.

Four things about the format are worth knowing before you fill one in.

`#` comment lines are allowed in the CSV, and the top of the file is the right place to record
where a pressure came from. A table `schamp` writes is a table `schamp` reads, so a comment
block survives a round trip through `report.write_table`.

An unrecognized column is an error rather than a shrug. A mistyped `presure_Torr` would
otherwise be an experiment with no pressures and no complaint about it.

`pressure_Torr_err` and `temperature_K_err` are optional columns, and filling them widens the
error bar on every cross section. `CrossSection.propagated` comes back naming the terms that
actually went in, so an error bar says what it covers. The example tables leave both columns
out, which is the honest answer for a series whose temperature was read once, from sensors that
were never checked against a reference.

Exclude an acquisition by setting `use = false` and saying why in `notes`, rather than by
deleting the row or by filtering in code. The record of what was excluded then travels with the
experiment, which is exactly what the 2013 spreadsheets did not do.

## Checks before you leave the instrument

Load the tables and validate them while the instrument is still available:

```python
import schamp

experiment = schamp.load_experiment("experiment.toml")
for complaint in experiment.validate():
    print(complaint)
print(experiment.drift_voltages())
```

`Experiment.validate` returns a list of one-line complaints rather than raising on the first, so
one call reports every problem in the series: duplicate rows, an unknown gas, a non-positive
pressure, temperature, or uncertainty, a `.raw` path that is not a directory, an unreadable
`_extern.inf`, a missing drift-voltage key, every profile check including the wave height and the
pushes per bin, fewer than three used acquisitions, and any two used acquisitions sharing a drift
voltage. Pass `require_raw=False` to check the tables alone on a machine that does not hold the
data.

Then look at the printed drift voltages. They should be monotone in the order you acquired them,
they should span the range you intended, and they should be evenly spread once inverted. A
voltage that is not where you put it means a knob that did not move, and catching that now is a
repeated acquisition rather than a lost series.

## Failure modes

Ordered by how quietly each one fails.

1. **The flushing pulses were never disabled.** No trace in the data, no check that can find it,
   and the series looks fine. Confirm the service-level change before a campaign.
2. **The drift-cell pressure had not settled.** Ω goes as 1/P, so a series acquired on a falling
   pressure carries a bias that changes down the series. Per-acquisition manometer readings are
   what make this visible after the fact, and a settled gauge before the first acquisition is
   what makes it absent.
3. **The mass calibration came from another day or another polarity.** Windows land off the peak,
   and on the 2013 anion series that reached 350 ppm. `CalibrationProvenance.concerns()` reads
   this off the file before anything is extracted.
4. **Something other than the drift voltage changed between acquisitions.** The regression
   assumes one variable. Retuning the source, changing the trap conditions, or altering the mass
   range partway through gives a series that is two series with one table.
5. **The ion was activated before it entered the cell.** A straight line and a wrong cross
   section. The published conditions minimize this, and the check for it is an independent
   measurement rather than anything in the data.
6. **The mobility cell's traveling wave was left on, or several pushes were summed into one
   drift bin.** Both caught, by the profile's `require_wave_off` and `expect_pushes_per_bin`.
7. **The slowest ion fell outside the drift-time window.** Compute the window from the pusher
   period and the bin count before acquiring, and look at the arrival-time distribution of the
   slowest ion at the lowest voltage while you are still at the instrument.
8. **Fewer than three usable voltages, or two at the same voltage.** Caught, by
   `Experiment.validate`.

## What comes next

With the tables filled and validated, `examples/polyalanine-walkthrough.ipynb` runs the whole
analysis on the published series: arrival-time distributions per m/z window, single-Gaussian
fits, the regression of arrival time against reciprocal drift voltage, and cross sections in
helium with their propagated uncertainties. Its last section names the four things to change to
point it at a series of your own.

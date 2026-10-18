"""The run sheet's **Live** tab. Imports PySide6; `live` holds everything that does not.

`LivePane` is a view over `live.LiveSeries`. Two `live.SerialWorker`s do every piece of
I/O on the data directory, one listing it for the watcher and one owning every SDK read,
so neither a share that has stopped answering nor a slow read ever stalls the window or
holds up the other. A timer collects their results. The species, and the `[watch]` rule
if the operator wrote one, live in the run sheet's `runsheet.toml`.

Three plots on one QtAgg canvas, for the species selected in the table: the latest
arrival-time distribution with its single-Gaussian fit; arrival time against reciprocal
drift voltage, with the line once three points are in the regression and any point more
than `live.OUTLIER_FACTOR` times the others' rms residual off it ringed in red; and the
cross section as it ran while the series arrived, from the third point on.
"""

from __future__ import annotations

import os

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6 import QtCore, QtGui, QtWidgets

from ..experiment import ExperimentError, load_experiment
from . import live, model

Qt = QtCore.Qt

DRAIN_EVERY_MS = 250


class LivePane(QtWidgets.QWidget):
    """The live analysis of a series as it arrives. `species_changed` asks for a save."""

    species_changed = QtCore.Signal()

    def __init__(self, unreachable_after_s: float = 8.0, parent=None) -> None:
        super().__init__(parent)
        self.unreachable_after_s = unreachable_after_s
        self.sheet: model.RunSheet | None = None
        self.experiment = None
        self.series = live.LiveSeries()
        self.watcher = live.Watcher()
        self.looks = live.SerialWorker("schamp-acquire watch")
        self.reader = live.SerialWorker("schamp-acquire sdk")
        self.sdk_problem = live.sdk_problem()
        self.rule_problem = ""
        self.listing: live.Listing | None = None
        self.picking = False

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.status)

        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Species", "m/z", "z", "Points", "R²", "Ω / Å²"])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(20)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setMinimumHeight(110)
        self.table.setMaximumHeight(150)
        self.table.itemSelectionChanged.connect(self.redraw)
        layout.addWidget(self.table)

        buttons = QtWidgets.QHBoxLayout()
        for text, slot in (
            ("Add…", self.add_species),
            ("Pick from spectrum…", self.pick_species),
            ("Remove", self.remove_species),
        ):
            button = QtWidgets.QPushButton(text)
            button.clicked.connect(slot)
            buttons.addWidget(button)
        buttons.addStretch(1)
        self.finish_button = QtWidgets.QPushButton("Finish series…")
        self.finish_button.setToolTip("Write results.json and ccs.csv beside experiment.toml")
        self.finish_button.clicked.connect(self.finish)
        buttons.addWidget(self.finish_button)
        layout.addLayout(buttons)

        self.figure = Figure(figsize=(4.2, 7.0), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(480)
        self.axes = [self.figure.add_subplot(3, 1, i + 1) for i in range(3)]
        layout.addWidget(self.canvas, 1)

        self.detail = QtWidgets.QPlainTextEdit()
        self.detail.setReadOnly(True)
        self.detail.setMaximumHeight(64)
        self.detail.setFrameShape(QtWidgets.QFrame.NoFrame)
        layout.addWidget(self.detail)

        self.look_timer = QtCore.QTimer(self)
        self.look_timer.timeout.connect(self.look)
        self.drain_timer = QtCore.QTimer(self, interval=DRAIN_EVERY_MS)
        self.drain_timer.timeout.connect(self.drain)
        self.drain_timer.start()
        self.set_sheet(None)

    # -- the sheet ---------------------------------------------------------------------

    def set_sheet(self, sheet: model.RunSheet | None) -> None:
        self.sheet = sheet
        self.looks.abandon()
        self.reader.abandon()
        self.listing = None
        species, self.rule_problem = [], ""
        rule = live.WatchRule()
        if sheet is not None:
            try:
                rule = live.WatchRule.from_table(sheet.watch)
            except live.LiveError as exc:
                self.rule_problem = f"{exc}; using the default rule instead"
            for entry in sheet.species:
                try:
                    species.append(live.Species.from_table(entry))
                except live.LiveError as exc:
                    self.rule_problem += f" {exc}"
        self.series = live.LiveSeries(species)
        self.watcher.reset(rule)
        self.look_timer.setInterval(int(rule.interval_s * 1000))
        self.look_timer.start()
        self.refresh_record()
        self.fill_table()
        if species:
            self.table.selectRow(0)
        self.setEnabled(sheet is not None)
        self.look()

    def data_dir_changed(self) -> None:
        """The operator pointed the sheet somewhere else: start watching from scratch."""
        self.set_sheet(self.sheet)

    def refresh_record(self) -> None:
        """Reload the experiment the run sheet just saved. Reads only the set directory."""
        self.experiment = None
        if self.sheet is not None and os.path.isfile(self.sheet.conditions_path):
            try:
                self.experiment = load_experiment(self.sheet.set_dir)
            except (ExperimentError, OSError, ValueError):
                self.experiment = None
        self.redraw()

    def save_species(self) -> None:
        if self.sheet is not None:
            self.sheet.species = [s.to_table() for s in self.series.species]
            self.species_changed.emit()

    # -- the workers -------------------------------------------------------------------

    def look(self) -> None:
        if self.sheet is None or self.looks.busy:
            self.show_status()
            return
        path = self.sheet.data_path
        self.looks.submit("look", lambda: live.list_raws(path))

    def drain(self) -> None:
        changed = False
        for _, listing in self.looks.drain():
            if self.sheet is None or listing.path != self.sheet.data_path:
                continue
            self.listing = listing
            for event in self.watcher.update(listing):
                if event.kind == "finished":
                    self.series.arrive(event.acquisition, event.path)
                else:
                    self.series.reopen(event.acquisition)
                changed = True
        for tag, result in self.reader.drain():
            if tag == "spectrum":
                self.picking = False
                self.open_picker(result)
                continue
            self.series.record(result, tag)
            changed = True
        if not self.sdk_problem:
            job = self.series.next_job()
            while job is not None:
                acquisition, path, wanted = job
                self.reader.submit(acquisition, self.series.job(acquisition, path, wanted))
                job = self.series.next_job()
        self.show_status()
        if changed:
            self.redraw()

    def show_status(self) -> None:
        lines = []
        if self.sheet is None:
            self.status.setText("Open or start a run sheet.")
            return
        if self.sdk_problem:
            lines.append(f"<b>{self.sdk_problem}</b>")
        running = self.looks.running()
        if running is not None and running[1] > self.unreachable_after_s:
            lines.append(f"Data directory: no answer after {running[1]:.0f} s.")
        elif self.listing is not None and not self.listing.reachable:
            lines.append(f"Data directory unreachable: {self.listing.error}.")
        elif self.listing is not None:
            states = [self.watcher.status(n) for n in self.listing.raws]
            finished = states.count("finished")
            settling = sum(s.startswith("settling") for s in states)
            writing = len(states) - finished - settling
            read = len(self.series.results)
            lines.append(
                f"{len(states)} .raw: {finished} finished, {settling} settling, {writing} being written; "
                f"{read} read" + (f", {len(self.series.failures)} failed" if self.series.failures else "")
                + (f", reading {reading[0]}" if (reading := self.reader.running()) else "")
            )
        rule = self.watcher.rule
        lines.append(
            f"Finished when {' + '.join(rule.files) or 'nothing required'} present and unchanged "
            f"over {rule.polls} looks {rule.interval_s:g} s apart ({rule.source})."
        )
        tips = [f"Completion rule: {rule.describe()}."]
        if self.rule_problem:
            lines.append(f"<b>{self.rule_problem}</b>")
        first = self.series.first_read()
        if first is not None:
            if first.concerns:
                lines.append(
                    f"<b>Mass calibration: {len(first.concerns)} concern(s)</b> on {first.acquisition}; hover for them."
                )
            else:
                lines.append(f"Mass calibration: no concerns from {first.acquisition}'s own record.")
            tips.append(f"{first.acquisition}: {first.calibration}.")
            tips += [f"- {c}" for c in first.concerns]
        lines += self.series.notes(self.experiment)
        self.status.setText("<br>".join(lines))
        self.status.setToolTip("\n".join(tips))

    # -- species -----------------------------------------------------------------------

    def selected(self) -> live.Species | None:
        rows = self.table.selectionModel().selectedRows()
        if not rows or rows[0].row() >= len(self.series.species):
            return None
        return self.series.species[rows[0].row()]

    def fill_table(self) -> None:
        keep = self.selected()
        self.table.setRowCount(len(self.series.species))
        profile = self.profile()
        for row, species in enumerate(self.series.species):
            summary = self.series.summarise(species, self.experiment, profile) if profile else None
            used = sum(p.used for p in summary.points) if summary else 0
            r2 = f"{summary.regression.r_squared:.5f}" if summary and summary.regression else ""
            xs = summary.cross_section if summary else None
            omega = "" if xs is None else f"{xs.omega_a2:.1f}" + (f" ± {xs.omega_a2_err:.1f}" if xs.omega_a2_err else "")
            for column, text in enumerate(
                (species.label, f"{species.mz:.4f}", f"{species.charge:+d}", str(used), r2, omega)
            ):
                item = QtWidgets.QTableWidgetItem(text)
                if column:
                    item.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                if summary and any(p.outlier for p in summary.points):
                    item.setToolTip("a point sits off the line; see the middle plot")
                self.table.setItem(row, column, item)
        self.table.resizeColumnsToContents()
        if keep is not None and keep in self.series.species:
            self.table.blockSignals(True)
            self.table.selectRow(self.series.species.index(keep))
            self.table.blockSignals(False)

    def add_species(self, *, mz: float | None = None) -> None:
        dialog = SpeciesDialog(self, mz=mz)
        if dialog.exec() != QtWidgets.QDialog.Accepted:
            return
        self._add(dialog.species())

    def _add(self, species: live.Species | None) -> None:
        if species is None:
            return
        try:
            self.series.add_species(species)
        except live.LiveError as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))
            return
        self.save_species()
        self.fill_table()
        self.table.selectRow(len(self.series.species) - 1)
        self.drain()

    def remove_species(self) -> None:
        species = self.selected()
        if species is not None:
            self.series.remove_species(species.label)
            self.save_species()
            self.fill_table()
            self.redraw()

    def pick_species(self) -> None:
        if self.sdk_problem:
            QtWidgets.QMessageBox.information(self, "schamp-acquire", self.sdk_problem)
            return
        if not self.series.arrived:
            QtWidgets.QMessageBox.information(
                self, "schamp-acquire",
                "No acquisition has finished yet. The picker shows the total spectrum of the "
                "first one to finish; until then, type a species with Add….",
            )
            return
        if self.picking:
            return
        first = next(iter(self.series.arrived.values()))
        self.picking = True
        self.status.setText(f"Reading the total spectrum of {os.path.basename(first)}…")
        self.reader.submit("spectrum", lambda: live.spectrum_of(first))

    def open_picker(self, spectrum) -> None:
        if isinstance(spectrum, Exception):
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", f"Reading the spectrum failed: {spectrum}")
            return
        dialog = SpectrumPicker(spectrum, self)
        dialog.picked.connect(self._add)
        dialog.exec()

    # -- the end of the series -----------------------------------------------------------

    def finish(self) -> None:
        if self.sheet is None:
            return
        if self.experiment is None:
            QtWidgets.QMessageBox.information(
                self, "schamp-acquire", "Nothing is ticked as taken, so there is no record to finish."
            )
            return
        waiting = len(self.series.in_flight) + sum(
            1 for a in self.series.arrived if a not in self.series.results and a not in self.series.failures
        )
        if waiting and QtWidgets.QMessageBox.question(
            self, "schamp-acquire", f"{waiting} acquisition(s) are still being read. Finish anyway?"
        ) != QtWidgets.QMessageBox.Yes:
            return
        try:
            paths = self.series.finish(self.experiment, self.sheet.set_dir, rule=self.watcher.rule)
        except OSError as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", f"Not written: {exc}")
            return
        QtWidgets.QMessageBox.information(self, "schamp-acquire", "Wrote\n" + "\n".join(paths))

    # -- plots -------------------------------------------------------------------------

    def profile(self):
        if self.experiment is not None:
            return self.experiment.profile
        try:
            return self.sheet.load_profile() if self.sheet is not None else None
        except Exception:  # noqa: BLE001 -- the run sheet already reports a bad profile
            return None

    def redraw(self) -> None:
        if self.sender() is not self.table:
            self.fill_table()
        for ax in self.axes:
            ax.clear()
        atd_ax, line_ax, ccs_ax = self.axes
        atd_ax.set_xlabel("drift time / ms")
        atd_ax.set_ylabel("intensity")
        line_ax.set_xlabel("1 / V$_{drift}$ / V$^{-1}$")
        line_ax.set_ylabel("arrival time / ms")
        ccs_ax.set_xlabel("points in the regression")
        ccs_ax.set_ylabel("Ω / Å$^2$")
        species, profile = self.selected(), self.profile()
        detail = ""
        if species is not None and profile is not None:
            detail = self._draw(species, profile)
        self.detail.setPlainText(detail)
        self.canvas.draw_idle()

    def _draw(self, species: live.Species, profile) -> str:
        atd_ax, line_ax, ccs_ax = self.axes
        latest = self.series.latest(species.label)
        notes = []
        if latest is not None and latest.atd is not None:
            atd = latest.atd
            atd_ax.plot(atd.drift_time_ms, atd.intensity, "-", color="0.35", lw=0.9)
            if latest.fit is not None:
                fit = latest.fit
                t = np.asarray(atd.drift_time_ms)
                atd_ax.plot(t, fit.height * np.exp(-((t - fit.centre_ms) ** 2) / (2 * fit.sigma_ms**2)), "-", color="C3", lw=1.0)
                atd_ax.axvline(fit.centre_ms, color="C3", lw=0.5, ls=":")
            atd_ax.set_title(
                f"{species.label}, {os.path.basename(atd.acquisition)}"
                + (f": {latest.fit.centre_ms:.3f} ms" if latest.fit else ""),
                fontsize=8,
            )
        if latest is not None and latest.problem:
            notes.append(f"latest: {latest.problem}")

        summary = self.series.summarise(species, self.experiment, profile)
        used = [p for p in summary.points if p.used and p.v_drift_v]
        others = [p for p in summary.points if not p.used and p.v_drift_v]
        if used:
            line_ax.plot([1 / p.v_drift_v for p in used], [p.centre_ms for p in used], "o", color="C0", ms=4)
        if others:
            line_ax.plot(
                [1 / p.v_drift_v for p in others], [p.centre_ms for p in others], "o",
                mfc="none", mec="C0", ms=4, label="not in the regression",
            )
        off = [p for p in summary.points if p.outlier and p.v_drift_v]
        if off:
            line_ax.plot([1 / p.v_drift_v for p in off], [p.centre_ms for p in off], "o", mfc="none", mec="C3", ms=9, mew=1.2, label="off the line")
            notes.append("off the line: " + ", ".join(p.acquisition for p in off))
        line = summary.regression
        if line is not None:
            x = np.array([min(1 / p.v_drift_v for p in summary.points if p.v_drift_v), max(1 / p.v_drift_v for p in summary.points if p.v_drift_v)])
            line_ax.plot(x, line.intercept_ms + line.slope_ms_v * x, "-", color="C0", lw=0.8)
            line_ax.set_title(f"R² = {line.r_squared:.5f}, t₀ = {line.intercept_ms:.3f} ms, n = {line.n}", fontsize=8)
        if others or off:
            line_ax.legend(fontsize="x-small", frameon=False, loc="upper left")

        if summary.running:
            counts = [n for n, _ in summary.running]
            omegas = [xs.omega_a2 for _, xs in summary.running]
            errors = [xs.omega_a2_err or 0.0 for _, xs in summary.running]
            ccs_ax.errorbar(counts, omegas, yerr=errors, fmt="o-", color="C2", ms=3, lw=0.8, capsize=2)
            ccs_ax.set_xlim(2.5, max(counts) + 0.5)
            xs = summary.cross_section
            ccs_ax.set_title(
                f"{xs.omega_a2:.1f}" + (f" ± {xs.omega_a2_err:.1f}" if xs.omega_a2_err else "")
                + f" Å², z = {species.charge:+d}, {xs.gas}, {xs.pressure_torr:.4g} Torr, {xs.temperature_k:.5g} K",
                fontsize=8,
            )
        elif summary.why_no_ccs:
            notes.append(f"cross section: {summary.why_no_ccs}")
        excluded = [f"{p.acquisition} ({p.why_not})" for p in summary.points if not p.used]
        if excluded:
            notes.append("not in the regression: " + ", ".join(excluded))
        return "\n".join(notes)


class SpeciesDialog(QtWidgets.QDialog):
    """A species typed in: label, m/z and signed charge."""

    def __init__(self, parent=None, *, mz: float | None = None, label: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle("Species")
        form = QtWidgets.QFormLayout(self)
        self.label = QtWidgets.QLineEdit(label)
        self.mz = QtWidgets.QLineEdit("" if mz is None else f"{mz:.4f}")
        self.mz.setValidator(QtGui.QDoubleValidator(0.0, 1e6, 6, self.mz))
        self.charge = QtWidgets.QSpinBox()
        self.charge.setRange(-99, 99)
        self.charge.setValue(1)
        form.addRow("Label", self.label)
        form.addRow("m/z", self.mz)
        form.addRow("Charge (signed)", self.charge)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        form.addRow(box)

    def species(self) -> live.Species | None:
        try:
            mz = float(self.mz.text())
            label = self.label.text().strip() or f"m/z {mz:.3f} {self.charge.value():+d}"
            return live.Species(label, mz, self.charge.value())
        except (ValueError, live.LiveError) as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))
            return None


class SpectrumPicker(QtWidgets.QDialog):
    """The total spectrum of the first finished acquisition; a click picks the peak nearest.

    The peak is `MassSpectrum.peak_near` -- measured, not fitted -- searched within 1 % of
    the m/z span in view, so zooming in narrows the search. Each species added goes
    straight to the table through `picked`; the dialog stays open for the next.
    """

    picked = QtCore.Signal(object)

    def __init__(self, spectrum, parent=None) -> None:
        super().__init__(parent)
        self.spectrum = spectrum
        self.peak = None
        self.setWindowTitle(f"Pick species: {os.path.basename(spectrum.acquisition)}")
        self.resize(900, 560)
        layout = QtWidgets.QVBoxLayout(self)
        self.figure = Figure(layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, 1)
        self.ax = self.figure.add_subplot()
        self.ax.plot(spectrum.mz, spectrum.intensity, "-", color="0.25", lw=0.6)
        self.ax.set_xlabel("m/z")
        self.ax.set_ylabel("intensity (all scans, all drift bins)")
        self.marker = None
        self.canvas.mpl_connect("button_press_event", self.on_click)

        row = QtWidgets.QHBoxLayout()
        self.info = QtWidgets.QLabel("Zoom with the toolbar, then click a peak (with no toolbar mode active).")
        row.addWidget(self.info, 1)
        row.addWidget(QtWidgets.QLabel("z"))
        self.charge = QtWidgets.QSpinBox()
        self.charge.setRange(-99, 99)
        self.charge.setValue(1)
        row.addWidget(self.charge)
        self.label = QtWidgets.QLineEdit()
        self.label.setPlaceholderText("label")
        row.addWidget(self.label)
        self.add = QtWidgets.QPushButton("Add species")
        self.add.setEnabled(False)
        self.add.clicked.connect(self.on_add)
        row.addWidget(self.add)
        done = QtWidgets.QPushButton("Done")
        done.clicked.connect(self.accept)
        row.addWidget(done)
        layout.addLayout(row)

    def on_click(self, event) -> None:
        if event.inaxes is not self.ax or event.xdata is None or self.toolbar.mode:
            return
        low, high = self.ax.get_xlim()
        peak = self.spectrum.peak_near(float(event.xdata), halfwidth=0.01 * (high - low))
        self.peak = peak
        if self.marker is not None:
            self.marker.remove()
            self.marker = None
        if peak is None:
            self.info.setText(f"No peak near m/z {event.xdata:.3f}.")
            self.add.setEnabled(False)
        else:
            self.marker = self.ax.axvline(peak.mz, color="C3", lw=0.8)
            self.info.setText(
                f"Peak at m/z {peak.mz:.4f}, FWHM {peak.fwhm:.4f}, apex {peak.apex:.3g}"
                + ("" if peak.resolved else " (unresolved shoulder)")
            )
            self.add.setEnabled(True)
        self.canvas.draw_idle()

    def on_add(self) -> None:
        if self.peak is None or self.charge.value() == 0:
            return
        label = self.label.text().strip() or f"m/z {self.peak.mz:.3f} {self.charge.value():+d}"
        try:
            self.picked.emit(live.Species(label, self.peak.mz, self.charge.value()))
        except live.LiveError as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))
            return
        self.label.clear()

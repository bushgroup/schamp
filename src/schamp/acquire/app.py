"""`schamp-acquire`: the run sheet's main window, and the entry point.

    schamp-acquire                      start empty: File > New or File > Open
    schamp-acquire SET_DIR              open the set there, or start one from the
                                        standard template if the folder holds none
    schamp-acquire SET_DIR --template T start a new set from template T (name or path)

The window is a view over `model.RunSheet`. Every edit is saved at once: the record
on disk is always what the window shows, and closing the window loses nothing. Calls
that touch the data directory run on a background thread (`widgets.Background`); a
share that does not answer within `UNREACHABLE_AFTER_S` is shown as unreachable and
asked again once it does.
"""

from __future__ import annotations

import argparse
import copy
import os
import sys

from PySide6 import QtCore, QtGui, QtWidgets

from .. import __version__
from ..constants import GASES
from ..importer import import_set
from . import model
from .widgets import Background, PlanCanvas, RunSheetTable

Qt = QtCore.Qt

SCAN_EVERY_MS = 4000
UNREACHABLE_AFTER_S = 8.0


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.sheet: model.RunSheet | None = None
        self.unsaved = False
        self.importing = False
        self.revalidate = False
        self.scan: model.DataDirScan | None = None

        self.table = RunSheetTable(model.RunSheet(set_dir=os.getcwd()), self)
        self.table.edited.connect(self.on_edit)
        self.view = QtWidgets.QTableView()
        self.view.setModel(self.table)
        self.view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.view.setAlternatingRowColors(False)
        self.view.verticalHeader().setDefaultSectionSize(22)
        self.view.horizontalHeader().setStretchLastSection(False)
        self.setCentralWidget(self.view)

        self._build_session_dock()
        self._build_side_dock()
        self._build_checks_dock()
        self._build_menus()

        self.worker = Background(self)
        self.worker.finished.connect(self.on_background)
        self.importer = Background(self)
        self.importer.finished.connect(self.on_imported)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(SCAN_EVERY_MS)
        self.timer.timeout.connect(self.tick)
        self.timer.start()

        self.resize(1500, 820)
        self.resizeDocks([self.session_dock, self.plot_dock], [330, 400], Qt.Horizontal)
        self.resizeDocks([self.checks_dock], [150], Qt.Vertical)
        self.set_sheet(None)

    # -- layout --------------------------------------------------------------------

    def _build_session_dock(self) -> None:
        form = QtWidgets.QFormLayout()
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.set_label = QtWidgets.QLabel()
        self.set_label.setWordWrap(True)
        self.set_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        form.addRow("Set", self.set_label)
        self.fields: dict[str, QtWidgets.QLineEdit] = {}
        for key, label, tip in (
            ("title", "Title", "The experiment's title, for experiment.toml"),
            ("operator", "Operator", "Who is acquiring; kept in the run sheet only"),
            ("profile", "Profile", "Instrument profile: a builtin name, or a path relative to the set"),
            ("temperature_k", "T / K", "Session temperature: the average of the cell thermocouples"),
            ("data_dir", "Data", "Where the .raw acquisitions appear: . for the set itself, a drive, or a share"),
            ("base_name", "Base name", "Acquisitions are named <base>_NNN.raw, as a MassLynx sample list numbers them"),
        ):
            edit = QtWidgets.QLineEdit()
            edit.setToolTip(tip)
            edit.editingFinished.connect(lambda k=key: self.on_session_field(k))
            self.fields[key] = edit
            if key == "data_dir":
                row = QtWidgets.QHBoxLayout()
                row.addWidget(edit)
                browse = QtWidgets.QToolButton(text="…")
                browse.clicked.connect(self.browse_data_dir)
                row.addWidget(browse)
                form.addRow(label, row)
            elif key == "base_name":
                row = QtWidgets.QHBoxLayout()
                row.addWidget(edit)
                rename = QtWidgets.QPushButton("Rename rows to come")
                rename.clicked.connect(self.rename_rows)
                row.addWidget(rename)
                form.addRow(label, row)
            else:
                form.addRow(label, edit)
        self.gas = QtWidgets.QComboBox()
        self.gas.addItem("profile default", "")
        for name in GASES:
            self.gas.addItem(name, name)
        self.gas.activated.connect(self.on_gas)
        form.addRow("Gas", self.gas)

        self.fixed_box = QtWidgets.QGroupBox("Held for the session / V")
        self.fixed_form = QtWidgets.QFormLayout(self.fixed_box)
        form.addRow(self.fixed_box)
        self.notes = QtWidgets.QPlainTextEdit()
        self.notes.setMaximumHeight(90)
        self.notes.setPlaceholderText("Notes for experiment.toml")
        self.notes_timer = QtCore.QTimer(self, singleShot=True, interval=800)
        self.notes_timer.timeout.connect(lambda: self.on_session_field("notes"))
        self.notes.textChanged.connect(self.notes_timer.start)
        form.addRow("Notes", self.notes)
        self.template_label = QtWidgets.QLabel()
        self.template_label.setWordWrap(True)
        form.addRow("Template", self.template_label)

        body = QtWidgets.QWidget()
        body.setLayout(form)
        body.setMaximumWidth(420)
        self.session_dock = QtWidgets.QDockWidget("Session", self)
        self.session_dock.setObjectName("session")
        self.session_dock.setWidget(body)
        self.session_dock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.session_dock)

    def _build_side_dock(self) -> None:
        tabs = QtWidgets.QTabWidget()
        self.canvas = PlanCanvas()
        tabs.addTab(self.canvas, "Series")
        live = QtWidgets.QLabel("Live arrival-time view: not available in this version.")
        live.setAlignment(Qt.AlignCenter)
        live.setWordWrap(True)
        tabs.addTab(live, "Live")
        self.side_tabs = tabs
        dock = QtWidgets.QDockWidget("Plot", self)
        dock.setObjectName("plot")
        dock.setWidget(tabs)
        dock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.plot_dock = dock

    def _build_checks_dock(self) -> None:
        body = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(body)
        layout.setContentsMargins(4, 4, 4, 4)
        self.data_status = QtWidgets.QLabel()
        self.data_status.setWordWrap(True)
        layout.addWidget(self.data_status)
        self.problems = QtWidgets.QListWidget()
        self.problems.setWordWrap(True)
        layout.addWidget(self.problems)
        dock = QtWidgets.QDockWidget("Checks", self)
        dock.setObjectName("checks")
        dock.setWidget(body)
        dock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)
        self.checks_dock = dock

    def _build_menus(self) -> None:
        bar = self.menuBar()
        file = bar.addMenu("&File")
        file.addAction("&New run sheet…", QtGui.QKeySequence.New, self.new_sheet)
        file.addAction("&Open set…", QtGui.QKeySequence.Open, self.open_sheet)
        self.import_action = file.addAction("&Import acquisitions…", self.import_acquisitions)
        file.addSeparator()
        file.addAction("&Quit", QtGui.QKeySequence.Quit, self.close)
        rows = bar.addMenu("&Rows")
        self.row_actions = [
            rows.addAction("Apply a &template…", self.apply_template),
            rows.addAction("&Insert row above", lambda: self.insert_row(0)),
            rows.addAction("Insert row &below", lambda: self.insert_row(1)),
            rows.addAction("&Remove row", self.remove_row),
            rows.addAction("Move &up", lambda: self.move_row(-1)),
            rows.addAction("Move &down", lambda: self.move_row(1)),
        ]

    # -- the sheet -----------------------------------------------------------------

    def set_sheet(self, sheet: model.RunSheet | None) -> None:
        self.sheet = sheet
        self.scan = None
        if hasattr(self, "worker"):
            self.worker.abandon()
        enabled = sheet is not None
        for widget in (self.session_dock, self.view, self.import_action, *self.row_actions):
            widget.setEnabled(enabled)
        if sheet is None:
            self.table.set_sheet(model.RunSheet(set_dir=os.getcwd()))
            self.setWindowTitle("schamp-acquire")
            self.data_status.setText("File > New run sheet, or File > Open set.")
            self.set_label.setText("")
            self.canvas.draw_sheet(self.table)
            return
        self.table.set_sheet(sheet)
        self.load_session_fields()
        self.view.resizeColumnsToContents()
        # Filled by the first scan, after this sizing has run.
        for head, width in (("Notes", 220), (".raw", 80), ("V recorded", 80), ("ΔV", 60)):
            self.view.setColumnWidth(self.table.column_of(head), width)
        self.save()

    def load_session_fields(self) -> None:
        sheet = self.sheet
        self.set_label.setText(sheet.set_dir)
        for key, edit in self.fields.items():
            edit.setText(str(getattr(sheet, key)))
        index = self.gas.findData(sheet.gas)
        if index < 0:
            self.gas.addItem(sheet.gas, sheet.gas)
            index = self.gas.count() - 1
        self.gas.setCurrentIndex(index)
        self.notes.blockSignals(True)
        self.notes.setPlainText(sheet.notes)
        self.notes.blockSignals(False)
        self.template_label.setText(sheet.template or "none")
        while self.fixed_form.rowCount():
            self.fixed_form.removeRow(0)
        for key, value in sheet.fixed.items():
            edit = QtWidgets.QLineEdit(f"{value:g}")
            edit.setValidator(QtGui.QDoubleValidator(edit))
            edit.editingFinished.connect(lambda k=key, e=edit: self.on_fixed(k, e))
            self.fixed_form.addRow(key, edit)
        self.setWindowTitle(f"schamp-acquire — {sheet.title}")

    def on_session_field(self, key: str) -> None:
        if self.sheet is None:
            return
        value = self.notes.toPlainText() if key == "notes" else self.fields[key].text().strip()
        if getattr(self.sheet, key) == value:
            return
        setattr(self.sheet, key, value)
        if key == "profile":
            complaint = self.table.refresh_profile()
            if complaint:
                self.statusBar().showMessage(complaint, 15000)
        if key == "title":
            self.setWindowTitle(f"schamp-acquire — {value}")
        if key == "data_dir":
            self.worker.abandon()
            self.scan = None
            self.table.set_scan(None)
        self.on_edit()

    def on_gas(self) -> None:
        if self.sheet is not None:
            self.sheet.gas = self.gas.currentData() or ""
            self.on_edit()

    def on_fixed(self, key: str, edit: QtWidgets.QLineEdit) -> None:
        try:
            value = float(edit.text())
        except ValueError:
            edit.setText(f"{self.sheet.fixed[key]:g}")
            return
        if self.sheet.fixed.get(key) != value:
            self.sheet.fixed[key] = value
            self.on_edit()

    def on_edit(self) -> None:
        self.table.everything_changed()
        self.save()

    def save(self) -> None:
        """Write the record, check it as written, and redraw. Never on the data directory."""
        if self.sheet is None or self.importing:
            return
        try:
            self.sheet.save()
            self.unsaved = False
        except OSError as exc:
            self.unsaved = True
            self.statusBar().showMessage(
                f"NOT SAVED: {exc.filename or ''} {exc.strerror or exc}. If a spreadsheet has it "
                "open, close it; saving is retried every few seconds.",
            )
            return
        rows_before = self.table.rowCount()
        self.table.layoutChanged.emit()
        if self.table.rowCount() != rows_before:
            self.view.resizeColumnsToContents()
        self.statusBar().showMessage(f"saved {self.sheet.set_dir}", 3000)
        # The tables-only check reads the set directory and nothing else, so it is safe
        # here; the check that needs the acquisitions follows from the background scan.
        self.show_problems(self.sheet.validate(require_raw=False), with_raw=False)
        self.canvas.draw_sheet(self.table)
        self.revalidate = True
        self.tick()

    # -- the background: scanning and the full check -------------------------------

    def tick(self) -> None:
        if self.sheet is None:
            return
        if self.unsaved:
            self.save()
            return
        if self.worker.busy:
            waited = self.worker.busy_for()
            if waited > UNREACHABLE_AFTER_S:
                self.table.set_scan(self.scan, unreachable=True)
                self.data_status.setText(
                    f"Data directory {self.sheet.data_path}: unreachable, no answer after "
                    f"{waited:.0f} s. The run sheet keeps working; nothing here needs it."
                )
            return
        snapshot = copy.deepcopy(self.sheet)
        self.revalidate = False

        def job():
            scan = model.scan_data_dir(snapshot.data_path)
            problems = None
            if scan.reachable and any(scan.externs):
                problems = snapshot.validate(require_raw=True)
            return snapshot.set_dir, scan, problems

        self.worker.submit("scan", job)

    def on_background(self, tag: str, result) -> None:
        if self.sheet is None or isinstance(result, Exception):
            if isinstance(result, Exception):
                self.data_status.setText(f"Checking the data directory failed: {result}")
            return
        set_dir, scan, problems = result
        if set_dir != self.sheet.set_dir:
            return
        self.scan = scan
        self.table.set_scan(scan)
        visible = sum(1 for row in self.sheet.rows if row.acquisition in scan.externs)
        if not scan.reachable:
            self.data_status.setText(f"Data directory {scan.path}: unreachable ({scan.error}).")
        else:
            others = len(set(scan.externs) - {row.acquisition for row in self.sheet.rows})
            self.data_status.setText(
                f"Data directory {scan.path}: {len(scan.externs)} .raw visible, "
                f"{visible} of {len(self.sheet.rows)} on the sheet"
                + (f", {others} not on the sheet" if others else "")
                + ("." if problems is not None else "; checked without the acquisitions.")
            )
        if problems is not None and not self.revalidate:
            self.show_problems(problems, with_raw=True)
        self.canvas.draw_sheet(self.table)
        if self.revalidate:
            self.tick()

    def show_problems(self, problems: list[str], *, with_raw: bool) -> None:
        self.problems.clear()
        if not problems:
            item = QtWidgets.QListWidgetItem(
                "No problems found" + (", acquisitions included." if with_raw else " in the tables.")
            )
            item.setForeground(QtGui.QColor(30, 120, 40))
            self.problems.addItem(item)
            return
        set_dir = self.sheet.set_dir + os.sep if self.sheet else ""
        for problem in problems:
            self.problems.addItem(problem.replace(set_dir, ""))

    # -- file actions --------------------------------------------------------------

    def new_sheet(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Folder for the new set (experiment.toml, conditions.csv)"
        )
        if not folder:
            return
        if any(os.path.isfile(os.path.join(folder, f)) for f in ("experiment.toml", model.RUNSHEET_FILENAME)):
            answer = QtWidgets.QMessageBox.question(
                self, "schamp-acquire", f"{folder} already holds a set. Open it instead?"
            )
            if answer == QtWidgets.QMessageBox.Yes:
                self.open_path(folder)
            return
        template = self.choose_template()
        if template is None:
            return
        base, ok = QtWidgets.QInputDialog.getText(
            self, "Acquisition names", "Base name; acquisitions are <base>_NNN.raw:",
            text=os.path.basename(folder),
        )
        if not ok:
            return
        try:
            sheet = model.RunSheet.from_template(template, folder, base_name=base.strip())
        except model.RunSheetError as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))
            return
        self.set_sheet(sheet)

    def choose_template(self) -> str | None:
        builtin = list(model.builtin_templates())
        choices = builtin + ["Another file…"]
        default = builtin.index(model.DEFAULT_TEMPLATE) if model.DEFAULT_TEMPLATE in builtin else 0
        choice, ok = QtWidgets.QInputDialog.getItem(
            self, "Series template", "Plan the series from:", choices, default, False
        )
        if not ok:
            return None
        if choice != "Another file…":
            return choice
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Series template", "", "TOML (*.toml)")
        return path or None

    def open_sheet(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Open a set")
        if folder:
            self.open_path(folder)

    def open_path(self, folder: str) -> None:
        try:
            self.set_sheet(model.RunSheet.open(folder))
        except (model.RunSheetError, OSError, ValueError) as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))

    def browse_data_dir(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Where the .raw acquisitions appear")
        if folder:
            self.fields["data_dir"].setText(os.path.normpath(folder))
            self.on_session_field("data_dir")

    def import_acquisitions(self) -> None:
        """Copy `.raw` acquisitions into the set, verified (`schamp.importer`)."""
        start = self.sheet.data_path if self.sheet.data_path != self.sheet.set_dir else ""
        source = QtWidgets.QFileDialog.getExistingDirectory(self, "Import .raw acquisitions from", start)
        if not source:
            return
        self.importing = True
        for widget in (self.session_dock, self.view, self.import_action, *self.row_actions):
            widget.setEnabled(False)
        self.data_status.setText(f"Importing from {source}: copying and verifying…")
        sheet = self.sheet
        self.importer.submit(
            "import",
            lambda: import_set(source, sheet.set_dir, profile=sheet.profile, title=sheet.title),
        )

    def on_imported(self, tag: str, report) -> None:
        self.importing = False
        for widget in (self.session_dock, self.view, self.import_action, *self.row_actions):
            widget.setEnabled(True)
        if isinstance(report, Exception):
            QtWidgets.QMessageBox.warning(self, "Import failed", str(report))
            return
        if self.sheet.data_dir not in ("", ".") and (report.imported or report.present):
            answer = QtWidgets.QMessageBox.question(
                self, "Import", report.describe()
                + "\n\nRead the acquisitions from the imported copies from now on?",
            )
            if answer == QtWidgets.QMessageBox.Yes:
                self.sheet.data_dir = "."
                self.fields["data_dir"].setText(".")
        else:
            QtWidgets.QMessageBox.information(self, "Import", report.describe())
        self.save()

    # -- rows ----------------------------------------------------------------------

    def current_row(self) -> int:
        index = self.view.currentIndex()
        return index.row() if index.isValid() else -1

    def apply_template(self) -> None:
        template = self.choose_template()
        if template is None:
            return
        try:
            self.table.beginResetModel()
            self.sheet.apply_template(template)
        except model.RunSheetError as exc:
            QtWidgets.QMessageBox.warning(self, "schamp-acquire", str(exc))
        finally:
            self.table.endResetModel()
        self.table.refresh_profile()
        self.load_session_fields()
        self.save()

    def insert_row(self, offset: int) -> None:
        if not self.sheet.vary:
            QtWidgets.QMessageBox.information(self, "schamp-acquire", "Apply a template first: it says which voltages a row sets.")
            return
        here = self.current_row()
        at = len(self.sheet.rows) if here < 0 else here + offset
        like = self.sheet.rows[here].settings if here >= 0 else None
        self.table.beginInsertRows(QtCore.QModelIndex(), at, at)
        self.sheet.insert_row(at, like)
        self.table.endInsertRows()
        self.view.selectRow(at)
        self.save()

    def remove_row(self) -> None:
        here = self.current_row()
        if here < 0:
            return
        row = self.sheet.rows[here]
        if row.taken:
            QtWidgets.QMessageBox.information(
                self, "schamp-acquire", f"{row.acquisition} is ticked as taken. Untick it first if it should leave the record."
            )
            return
        self.table.beginRemoveRows(QtCore.QModelIndex(), here, here)
        del self.sheet.rows[here]
        self.table.endRemoveRows()
        self.save()

    def move_row(self, step: int) -> None:
        here = self.current_row()
        there = here + step
        if here < 0 or not 0 <= there < len(self.sheet.rows):
            return
        rows = self.sheet.rows
        self.table.layoutAboutToBeChanged.emit()
        rows[here], rows[there] = rows[there], rows[here]
        self.table.layoutChanged.emit()
        self.view.selectRow(there)
        self.save()

    def rename_rows(self) -> None:
        self.on_session_field("base_name")
        self.sheet.propose_names()
        self.on_edit()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="schamp-acquire",
        description="The run sheet for acquiring a drift-voltage series: the planned settings, "
        "the conditions record, and the checks on both.",
    )
    parser.add_argument("set_dir", nargs="?", help="the set to open, or to start")
    parser.add_argument("--template", help=f"series template for a new set (default {model.DEFAULT_TEMPLATE})")
    parser.add_argument("--version", action="version", version=f"schamp {__version__}")
    args = parser.parse_args(argv)

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("schamp-acquire")
    window = MainWindow()
    if args.set_dir:
        folder = os.path.abspath(args.set_dir)
        holds_set = any(
            os.path.isfile(os.path.join(folder, f)) for f in ("experiment.toml", model.RUNSHEET_FILENAME)
        )
        try:
            if holds_set and not args.template:
                window.set_sheet(model.RunSheet.open(folder))
            elif holds_set:
                parser.error(f"{folder} already holds a set; open it without --template")
            else:
                window.set_sheet(
                    model.RunSheet.from_template(args.template or model.DEFAULT_TEMPLATE, folder)
                )
        except model.RunSheetError as exc:
            parser.error(str(exc))
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())

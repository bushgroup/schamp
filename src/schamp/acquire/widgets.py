"""The run sheet's widgets. Imports PySide6; `model` holds everything that does not.

Three pieces: `RunSheetTable`, the Qt item model over `model.RunSheet.rows`; `PlanCanvas`,
the planned series on a matplotlib QtAgg canvas; and `Background`, which runs the calls
that touch the data directory on a daemon thread and hands the answer back to the GUI
thread on a timer, so that a share which has gone away shows as unreachable rather than
freezing the window.
"""

from __future__ import annotations

import threading
import time
from typing import Any, Callable

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6 import QtCore, QtGui

from . import model

Qt = QtCore.Qt


class Background(QtCore.QObject):
    """One job at a time on a daemon thread, its result delivered on the GUI thread.

    A daemon thread rather than a pool: a stat on a dead share can block for as long as
    Windows takes to give up, and a pool's worker would hold the application open on
    exit for that long. `busy_for()` says how long the current job has been running, so
    the caller can say "no answer after N s" instead of waiting.
    """

    finished = QtCore.Signal(object, object)
    """(tag, result) -- or (tag, exception) if the job raised."""

    def __init__(self, parent: QtCore.QObject | None = None) -> None:
        super().__init__(parent)
        self._lock = threading.Lock()
        self._job: dict[str, Any] | None = None
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(100)
        self._timer.timeout.connect(self._poll)

    @property
    def busy(self) -> bool:
        return self._job is not None

    def busy_for(self) -> float:
        return time.monotonic() - self._job["started"] if self._job else 0.0

    def abandon(self) -> None:
        """Stop waiting for the current job. Its thread runs on and its answer is dropped;
        for a scan of a directory the user has since changed."""
        with self._lock:
            self._job = None
        self._timer.stop()

    def submit(self, tag: str, call: Callable[[], Any]) -> bool:
        """Start `call()` unless a job is already running. Returns whether it started."""
        if self._job is not None:
            return False
        job: dict[str, Any] = {"tag": tag, "started": time.monotonic(), "done": False}

        def run() -> None:
            try:
                result = call()
            except Exception as exc:  # noqa: BLE001 -- handed to the GUI to report
                result = exc
            with self._lock:
                job["result"], job["done"] = result, True

        self._job = job
        threading.Thread(target=run, name=f"schamp-acquire {tag}", daemon=True).start()
        self._timer.start()
        return True

    def _poll(self) -> None:
        with self._lock:
            job = self._job
            if job is None or not job["done"]:
                return
            self._job = None
        self._timer.stop()
        self.finished.emit(job["tag"], job["result"])


class RunSheetTable(QtCore.QAbstractTableModel):
    """The sheet's rows as a table: planned settings, the record, and what was recorded.

    Planned settings and the name are editable until a row is ticked as taken, and
    frozen after: they describe what the operator was asked to set, and the `.raw`
    under that name is the acquisition. `edited` fires after every change the user
    makes, which is the window's cue to save.
    """

    edited = QtCore.Signal()

    FIXED_HEADS = ("V nominal", "Acquisition", "Taken", "P / Torr", "Use", "Notes", ".raw", "V recorded", "ΔV")

    def __init__(self, sheet: model.RunSheet, parent: QtCore.QObject | None = None) -> None:
        super().__init__(parent)
        self.sheet = sheet
        self.scan: model.DataDirScan | None = None
        self.unreachable = False
        self._profile = None
        self.refresh_profile()

    # -- state from outside --------------------------------------------------------

    def set_sheet(self, sheet: model.RunSheet) -> None:
        self.beginResetModel()
        self.sheet = sheet
        self.scan = None
        self.refresh_profile()
        self.endResetModel()

    def refresh_profile(self) -> str:
        """Reload the profile; returns the complaint if it will not load."""
        try:
            self._profile = self.sheet.load_profile()
            return ""
        except Exception as exc:  # noqa: BLE001 -- shown in the window
            self._profile = None
            return str(exc)

    def set_scan(self, scan: model.DataDirScan | None, unreachable: bool = False) -> None:
        self.scan, self.unreachable = scan, unreachable
        self.everything_changed()

    def everything_changed(self) -> None:
        if self.rowCount():
            self.dataChanged.emit(self.index(0, 0), self.index(self.rowCount() - 1, self.columnCount() - 1))
        self.headerDataChanged.emit(Qt.Horizontal, 0, max(self.columnCount() - 1, 0))

    # -- columns -------------------------------------------------------------------

    def heads(self) -> list[str]:
        return [f"{key} / V" for key in self.sheet.vary] + list(self.FIXED_HEADS)

    def column_of(self, head: str) -> int:
        return len(self.sheet.vary) + self.FIXED_HEADS.index(head)

    def rowCount(self, parent=QtCore.QModelIndex()) -> int:  # noqa: N802 -- Qt's name
        return 0 if parent.isValid() else len(self.sheet.rows)

    def columnCount(self, parent=QtCore.QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self.heads())

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            heads = self.heads()
            return heads[section] if section < len(heads) else None
        return str(section + 1)

    # -- cells ---------------------------------------------------------------------

    def _head(self, column: int) -> tuple[str, str | None]:
        """`("vary", key)` for a settings column, else `(head, None)`."""
        if column < len(self.sheet.vary):
            return "vary", self.sheet.vary[column]
        return self.FIXED_HEADS[column - len(self.sheet.vary)], None

    def _extern(self, row: model.PlannedRow) -> dict[str, str] | None:
        if self.scan is None or not self.scan.reachable:
            return None
        return self.scan.externs.get(row.acquisition)

    def _check(self, row: model.PlannedRow) -> model.VoltageCheck | None:
        extern = self._extern(row)
        if not extern or self._profile is None:
            return None
        return self.sheet.check_recorded(row, extern, self._profile)

    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags
        row = self.sheet.rows[index.row()]
        head, _ = self._head(index.column())
        flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if head in ("Taken", "Use"):
            return flags | Qt.ItemIsUserCheckable
        if head in ("P / Torr", "Notes"):
            return flags | Qt.ItemIsEditable
        if head == "Acquisition" and not row.taken:
            return flags | Qt.ItemIsEditable
        if head == "vary" and row.planned and not row.taken:
            return flags | Qt.ItemIsEditable
        return flags

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        row = self.sheet.rows[index.row()]
        head, key = self._head(index.column())

        if role == Qt.CheckStateRole:
            if head == "Taken":
                return Qt.Checked if row.taken else Qt.Unchecked
            if head == "Use":
                return Qt.Checked if row.use else Qt.Unchecked
            return None

        if role in (Qt.DisplayRole, Qt.EditRole):
            if head == "vary":
                return f"{row.settings[key]:g}" if row.planned else ""
            if head == "V nominal":
                if not row.planned or self._profile is None:
                    return ""
                return f"{self.sheet.nominal_voltage(row, self._profile):.2f}"
            if head == "Acquisition":
                return row.acquisition
            if head == "P / Torr":
                return row.cells.get("pressure_Torr", "")
            if head == "Notes":
                return row.cells.get("notes", "")
            if head == ".raw":
                return self._raw_status(row) if role == Qt.DisplayRole else None
            check = self._check(row)
            if head == "V recorded" and check is not None and check.recorded_v is not None:
                return f"{check.recorded_v:.2f}"
            if head == "ΔV" and check is not None and check.delta_v is not None:
                return f"{check.delta_v:+.2f}"
            return ""

        if role == Qt.TextAlignmentRole and head in ("vary", "V nominal", "P / Torr", "V recorded", "ΔV"):
            return int(Qt.AlignRight | Qt.AlignVCenter)

        if role == Qt.BackgroundRole:
            if head in ("V recorded", "ΔV"):
                check = self._check(row)
                if check is not None and row.planned and not check.agrees:
                    return QtGui.QBrush(QtGui.QColor(255, 205, 200))
            if head == "P / Torr" and row.taken and not row.cells.get("pressure_Torr", "").strip():
                return QtGui.QBrush(QtGui.QColor(255, 235, 190))
            if row.taken:
                return QtGui.QBrush(QtGui.QColor(232, 244, 232))
            return None

        if role == Qt.ToolTipRole:
            check = self._check(row)
            if head in ("V recorded", "ΔV") and check is not None:
                if check.error:
                    return check.error
                if check.knobs:
                    return "\n".join(
                        f"{k}: planned {want:g} V, recorded "
                        + ("nothing" if got is None else f"{got:g} V")
                        for k, want, got in check.knobs
                    )
                return "recorded as planned" if check.agrees else None
            if head == ".raw" and self.unreachable:
                return self.scan.error if self.scan else "no answer from the data directory"
            return None
        return None

    def _raw_status(self, row: model.PlannedRow) -> str:
        if self.unreachable:
            return "unreachable"
        if self.scan is None:
            return "…"
        if not self.scan.reachable:
            return "unreachable"
        if row.acquisition not in self.scan.externs:
            return "not yet" if not row.taken else "absent"
        return "present" if self.scan.externs[row.acquisition] else "writing"

    def setData(self, index, value, role=Qt.EditRole):  # noqa: N802
        if not index.isValid():
            return False
        row = self.sheet.rows[index.row()]
        head, key = self._head(index.column())
        if role == Qt.CheckStateRole:
            checked = Qt.CheckState(value) == Qt.Checked
            if head == "Taken":
                row.taken = checked
            elif head == "Use":
                row.cells["use"] = "true" if checked else "false"
            else:
                return False
        elif role == Qt.EditRole:
            text = str(value).strip()
            if head == "vary":
                try:
                    row.settings[key] = float(text)
                except ValueError:
                    return False
            elif head == "Acquisition":
                if not text:
                    return False
                row.acquisition = text if text.lower().endswith(".raw") else f"{text}.raw"
            elif head == "P / Torr":
                row.cells["pressure_Torr"] = text
            elif head == "Notes":
                row.cells["notes"] = text
            else:
                return False
        else:
            return False
        self.dataChanged.emit(self.index(index.row(), 0), self.index(index.row(), self.columnCount() - 1))
        self.edited.emit()
        return True


class PlanCanvas(FigureCanvasQTAgg):
    """The series as planned, in acquisition order: nominal drift voltage per row.

    Filled markers are rows taken, open ones rows to come, and a cross is the voltage
    an acquisition recorded where its `.raw` is visible. The pseudo-random order shows
    as the scatter it is meant to be.
    """

    def __init__(self, parent=None) -> None:
        self.figure = Figure(figsize=(4.2, 3.2), layout="constrained")
        super().__init__(self.figure)
        self.setParent(parent)
        self.axes = self.figure.add_subplot()

    def draw_sheet(self, table: RunSheetTable) -> None:
        sheet, profile = table.sheet, table._profile
        ax = self.axes
        ax.clear()
        ax.set_xlabel("acquisition")
        ax.set_ylabel("drift voltage / V")
        if profile is not None:
            planned = [
                (n, sheet.nominal_voltage(r, profile), r.taken)
                for n, r in enumerate(sheet.rows, start=1)
                if r.planned
            ]
            for taken, style in ((True, {"color": "C0"}), (False, {"facecolor": "none", "edgecolor": "C0"})):
                points = [(n, v) for n, v, t in planned if t == taken]
                if points:
                    ax.scatter(*zip(*points), s=28, label="taken" if taken else "to come", **style)
            recorded = []
            for n, row in enumerate(sheet.rows, start=1):
                check = table._check(row)
                if check is not None and check.recorded_v is not None:
                    recorded.append((n, check.recorded_v))
            if recorded:
                ax.scatter(*zip(*recorded), marker="x", color="C3", s=36, label="recorded")
            if planned or recorded:
                ax.legend(loc="best", fontsize="small", frameon=False)
        ax.set_xlim(0.5, max(len(sheet.rows), 1) + 0.5)
        self.draw_idle()

"""The first-launch wizard: the Waters SDK from the lab's own download. Imports PySide6.

The installed `schamp-acquire` ships without the Waters MassLynx SDK, whose EULA forbids
redistribution (lab record, task 01): each lab supplies its own download and key, and
this wizard installs them under the user's profile through `schamp.sdk_install` -- the
same verify, unpack and import check as `tools/install_sdk.py` in a checkout. All of the
work is in that module; this is only its pages. The run sheet works before, during and
without it, and the Live tab says what it is waiting for until it has run.

Four pages: what is about to happen, the download, the key (the download's own, or one
Waters sent separately), and the install itself on a background thread with its log.
The key's contents are never shown, only where it was found.
"""

from __future__ import annotations

import os
import queue

from PySide6 import QtCore, QtWidgets

from .. import FROZEN, SDK_DIR, sdk, sdk_install
from .widgets import Background

TITLE = "Install the Waters MassLynx SDK"


class IntroPage(QtWidgets.QWizardPage):
    def __init__(self) -> None:
        super().__init__()
        self.setTitle(TITLE)
        text = QtWidgets.QLabel(
            "<p>The Live tab reads acquisitions through the Waters MassLynx SDK. Waters "
            "licenses it to each lab, not to schamp, and does not allow it to be passed "
            "on, so schamp-acquire does not include it: you install it once from your own "
            "Waters download (document 667007504, a zip named like "
            "<tt>667007504DDRevA.zip</tt>).</p>"
            "<p>The download is checked against the checksums Waters ships inside it, "
            "unpacked unaltered into</p>"
            f"<p><tt>{SDK_DIR}</tt></p>"
            "<p>and loaded once to prove it works. Nothing is sent anywhere. Uninstalling "
            "schamp-acquire leaves that folder in place; delete it by hand when no copy of "
            "schamp needs it.</p>"
            "<p>The run sheet does not need the SDK, and keeps working whether or not you "
            "install it now.</p>"
        )
        text.setWordWrap(True)
        text.setTextFormat(QtCore.Qt.RichText)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(text)


def _path_row(edit: QtWidgets.QLineEdit, browse) -> QtWidgets.QHBoxLayout:
    row = QtWidgets.QHBoxLayout()
    row.addWidget(edit)
    button = QtWidgets.QPushButton("Browse…")
    button.clicked.connect(browse)
    row.addWidget(button)
    return row


class DownloadPage(QtWidgets.QWizardPage):
    def __init__(self) -> None:
        super().__init__()
        self.setTitle("Your Waters download")
        self.setSubTitle("The zip as Waters supplied it. Do not unpack it first.")
        self.edit = QtWidgets.QLineEdit()
        self.edit.setPlaceholderText("…\\667007504DDRevA.zip")
        self.registerField("download*", self.edit)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(_path_row(self.edit, self.browse))
        layout.addStretch()

    def browse(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "The Waters SDK download", self.edit.text(), "Zip archives (*.zip)"
        )
        if path:
            self.edit.setText(os.path.normpath(path))

    def validatePage(self) -> bool:
        path = self.edit.text().strip()
        if not os.path.isfile(path):
            QtWidgets.QMessageBox.warning(self, TITLE, f"There is no file at {path}.")
            return False
        return True


class KeyPage(QtWidgets.QWizardPage):
    def __init__(self) -> None:
        super().__init__()
        self.setTitle("The license key")
        self.setSubTitle(
            "The SDK will not open an acquisition without the key Waters issued to your lab."
        )
        self.bundled = QtWidgets.QRadioButton("Use the license.key that came in the download")
        self.separate = QtWidgets.QRadioButton("Use a license.key file Waters sent separately:")
        self.bundled.setChecked(True)
        self.edit = QtWidgets.QLineEdit()
        browse_row = _path_row(self.edit, self.browse)
        self.separate.toggled.connect(lambda on: self.edit.setEnabled(on))
        self.edit.setEnabled(False)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.bundled)
        layout.addWidget(self.separate)
        layout.addLayout(browse_row)
        note = QtWidgets.QLabel(
            "A separate key stays where it is; its location is written to "
            f"<tt>{os.path.join(os.path.expanduser('~'), '.' + sdk.CONFIG_NAME)}</tt> as "
            "<tt>[license] path</tt>. Its contents are never displayed or copied."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        layout.addStretch()

    def browse(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "The license key", self.edit.text(), "License key (*.key);;All files (*)"
        )
        if path:
            self.edit.setText(os.path.normpath(path))

    def license_file(self) -> str | None:
        return self.edit.text().strip() if self.separate.isChecked() else None

    def validatePage(self) -> bool:
        path = self.license_file()
        if path is None:
            return True
        try:
            sdk._read_key_file(path)
        except sdk.ConfigError as exc:
            QtWidgets.QMessageBox.warning(self, TITLE, str(exc))
            return False
        return True


class InstallPage(QtWidgets.QWizardPage):
    """Runs `sdk_install.install` off the GUI thread and shows its log as it goes."""

    def __init__(self) -> None:
        super().__init__()
        self.setTitle("Installing")
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.log = QtWidgets.QPlainTextEdit(readOnly=True)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.log)
        self.lines: queue.Queue[str] = queue.Queue()
        self.worker = Background(self)
        self.worker.finished.connect(self.on_finished)
        self.timer = QtCore.QTimer(self, interval=100)
        self.timer.timeout.connect(self.drain)
        self.manifest: dict | None = None

    def initializePage(self) -> None:
        self.manifest = None
        self.log.clear()
        self.status.setText("Checking the download against its own checksums, then unpacking…")
        wizard = self.wizard()
        download = self.field("download").strip()
        license_file = wizard.key_page.license_file()
        self.timer.start()
        self.worker.submit(
            "install",
            lambda: sdk_install.install(
                download,
                license_file=license_file,
                link_environment=not FROZEN,
                log=self.lines.put,
            ),
        )
        self.completeChanged.emit()
        # No going back while it runs: a second install into the same folder would race
        # the first. After QWizard's own button update, which would re-enable it.
        QtCore.QTimer.singleShot(0, lambda: self._back(False))

    def _back(self, enabled: bool) -> None:
        button = self.wizard().button(QtWidgets.QWizard.BackButton)
        if button is not None:
            button.setEnabled(enabled)

    def drain(self) -> None:
        while True:
            try:
                self.log.appendPlainText(self.lines.get_nowait())
            except queue.Empty:
                return

    def on_finished(self, _tag: str, result) -> None:
        self.timer.stop()
        self.drain()
        if isinstance(result, Exception):
            self.status.setText(f"<b>Not installed.</b> {result}")
            self.status.setTextFormat(QtCore.Qt.RichText)
            self.completeChanged.emit()
            QtCore.QTimer.singleShot(0, lambda: self._back(True))
            return
        self.manifest = result
        try:
            source = sdk.license_source()
        except sdk.ConfigError as exc:
            source = f"none found ({exc})"
        wanted = self.wizard().key_page.license_file()
        warning = ""
        if wanted and os.path.abspath(wanted) not in source:
            warning = (
                "<p><b>Note:</b> a setting that comes first in schamp's order is chosen "
                "over the key you picked; remove it if that is not what you want.</p>"
            )
        self.status.setTextFormat(QtCore.Qt.RichText)
        self.status.setText(
            f"<p><b>masslynxsdk {result['version']} is installed</b> in "
            f"<tt>{SDK_DIR}</tt> and loads.</p>"
            f"<p>The license key schamp will use: {source}.</p>{warning}"
            "<p>The Live tab can read acquisitions now.</p>"
        )
        self.completeChanged.emit()
        # Installed and loaded: a second install now would find the SDK in use.
        QtCore.QTimer.singleShot(0, lambda: self._back(False))

    def isComplete(self) -> bool:
        return self.manifest is not None

    def cleanupPage(self) -> None:
        # Back from a failed install: forget it, so a second attempt starts clean.
        self.worker.abandon()
        self.timer.stop()


class SdkWizard(QtWidgets.QWizard):
    """`exec()` returns Accepted once the SDK is installed and imports."""

    def __init__(self, parent: QtWidgets.QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(TITLE)
        self.setWizardStyle(QtWidgets.QWizard.ModernStyle)
        self.key_page = KeyPage()
        self.install_page = InstallPage()
        for page in (IntroPage(), DownloadPage(), self.key_page, self.install_page):
            self.addPage(page)
        self.resize(640, 480)

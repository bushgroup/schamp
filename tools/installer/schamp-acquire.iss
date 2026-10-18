; schamp-acquire: the per-user installer for Windows 10/11 x64. Inno Setup 6.
;
; Compiled by tools/build_exe.py, never by hand: it builds the folder this packs and passes
;   /DAppVersion /DSourceDir /DOutputDir /DOutputBaseFilename
;   /DVcRedist /DVcRedistMajor /DVcRedistMinor /DVcRedistBld
;
; No administrator rights: PrivilegesRequired=lowest installs under
; %LOCALAPPDATA%\Programs\schamp-acquire and writes only to the current user's profile and
; registry hive. The Waters MassLynx SDK is not in this installer and may not be: its EULA
; forbids redistribution, and each lab installs it from its own Waters download when the
; application first starts (into %LOCALAPPDATA%\schamp\masslynxsdk). Uninstalling leaves
; that folder, and says so.
;
; The Microsoft Visual C++ v14 Redistributable (x64), vc_redist.x64.exe, is carried as
; Distributable Code under the Microsoft Software License Terms for Visual Studio 2026:
; "Distributable Code Files for Visual Studio 2026 and later", the REDIST list that the
; terms' Distributable Code section references (Community, Professional, Enterprise and
; Build Tools), https://learn.microsoft.com/en-us/visualstudio/releases/2026/redistribution
; (2025-11-11), which permits a validly licensed user to distribute the unmodified
; redistributable with their program, subject to those terms; and "Redistribute Visual C++
; Files", https://learn.microsoft.com/en-us/cpp/windows/redistributing-visual-cpp-files
; (2026-04-13), which states that distribution of the package downloaded from Microsoft is
; limited to licensed Visual Studio users. So the machine that builds a release installer
; for distribution carries a Visual Studio licence -- Visual Studio Build Tools or
; Community, whose terms cover academic research and open-source development -- and the
; file is shipped unmodified, inside this installer, never on its own. It is run only when
; the registry reports no x64 v14 runtime at least as new (Microsoft's own advice), and it
; asks for administrator rights itself; declined, the application still has the runtime
; that Qt and Python ship beside it.

#ifndef AppVersion
  #error Compile with tools/build_exe.py, which defines AppVersion, SourceDir and the rest
#endif

[Setup]
AppId={{18EEE4BA-B049-4581-9F8B-84BCFFF09FB9}
AppName=schamp-acquire
AppVersion={#AppVersion}
AppVerName=schamp-acquire {#AppVersion}
AppPublisher=Bush lab, University of Washington
AppPublisherURL=https://github.com/bushgroup/schamp
AppSupportURL=https://github.com/bushgroup/schamp/issues
AppCopyright=Copyright (c) University of Washington
VersionInfoVersion={#AppVersion}
DefaultDirName={autopf}\schamp-acquire
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
LicenseFile={#SourceDir}\LICENSE.txt
OutputDir={#OutputDir}
OutputBaseFilename={#OutputBaseFilename}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
UninstallDisplayName=schamp-acquire {#AppVersion}
UninstallDisplayIcon={app}\schamp-acquire.exe
SetupLogging=yes

[Tasks]
Name: desktopicon; Description: "Put a shortcut on the desktop"; Flags: unchecked
; Offered only where the PC's runtime is missing or older. It needs an administrator, so a
; user without one can untick it and rely on the runtime the application carries.
Name: vcredist; Description: "Install the Microsoft Visual C++ runtime (asks for administrator rights)"; Check: VCRedistNeeded

[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "{#VcRedist}"; DestDir: "{tmp}"; Flags: deleteafterinstall; Tasks: vcredist

[Icons]
Name: "{autoprograms}\schamp-acquire"; Filename: "{app}\schamp-acquire.exe"
Name: "{autodesktop}\schamp-acquire"; Filename: "{app}\schamp-acquire.exe"; Tasks: desktopicon

[Run]
Filename: "{tmp}\vc_redist.x64.exe"; Parameters: "/install /passive /norestart"; \
  StatusMsg: "Installing the Microsoft Visual C++ runtime..."; Tasks: vcredist; \
  Flags: waituntilterminated
Filename: "{app}\schamp-acquire.exe"; Description: "Start schamp-acquire"; \
  Flags: nowait postinstall skipifsilent

[UninstallDelete]
; The application's own cache (matplotlib's font list); never the SDK beside it.
Type: filesandordirs; Name: "{localappdata}\schamp\cache"

[Code]
const
  RuntimeKey = 'SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64';

function VCRedistNeeded: Boolean;
var
  Installed, Major, Minor, Bld: Cardinal;
begin
  Result := True;
  if RegQueryDWordValue(HKLM64, RuntimeKey, 'Installed', Installed) and (Installed = 1)
     and RegQueryDWordValue(HKLM64, RuntimeKey, 'Major', Major)
     and RegQueryDWordValue(HKLM64, RuntimeKey, 'Minor', Minor)
     and RegQueryDWordValue(HKLM64, RuntimeKey, 'Bld', Bld) then
    Result := (Major < {#VcRedistMajor})
      or ((Major = {#VcRedistMajor}) and (Minor < {#VcRedistMinor}))
      or ((Major = {#VcRedistMajor}) and (Minor = {#VcRedistMinor}) and (Bld < {#VcRedistBld}));
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Sdk: String;
begin
  if CurUninstallStep <> usPostUninstall then
    Exit;
  Sdk := ExpandConstant('{localappdata}\schamp\masslynxsdk');
  if DirExists(Sdk) and not UninstallSilent then
    MsgBox('schamp-acquire is removed.' + #13#10#13#10 +
      'The Waters MassLynx SDK you installed for it is still in' + #13#10 + Sdk + #13#10 +
      'and the path of any separately supplied license key is still in ' +
      ExpandConstant('{%USERPROFILE}') + '\.schamp.ini. They belong to your lab''s Waters ' +
      'licence, not to schamp: delete them by hand when no copy of schamp needs them.',
      mbInformation, MB_OK);
end;

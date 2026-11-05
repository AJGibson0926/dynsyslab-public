# Install DynSysLab

Repository: **AJGibson0926/dynsyslab-public**. The native desktop release supports
Windows x64. It includes Solver/Post-processor, Viewer, the scientific CLI,
numerical libraries, and the C++ extension SDK. No Python, compiler, Git,
administrator account, or scientific datasets are required to run it.

## Install the latest prebuilt release

Download [install.ps1](https://github.com/AJGibson0926/dynsyslab-public/blob/main/install.ps1)
from this repository and run it in PowerShell. Equivalently:

```powershell
Invoke-WebRequest https://raw.githubusercontent.com/AJGibson0926/dynsyslab-public/main/install.ps1 -OutFile "$env:TEMP\Install-DynSysLab.ps1"
& "$env:TEMP\Install-DynSysLab.ps1"
```

The installer downloads the latest Windows ZIP and its SHA-256 manifest from
[GitHub Releases](https://github.com/AJGibson0926/dynsyslab-public/releases/latest),
checks the ZIP before extraction, installs under
`%LOCALAPPDATA%\Programs\DynSysLab`, and creates Start-menu entries for both
applications. Close the applications before running the installer to update.

New installations store data under `Documents\DynSysLab\archives`. Existing
archive settings are preserved. `-DataRoot C:\MyScience\DynSysLab` chooses
a different data workspace; `-InstallRoot` chooses the application location.
The installer does not delete scientific data. Verified downloads are retained
in the temporary directory reported at completion.

Pin a version with `-Version v0.5.1`. For offline installation, download the
Windows ZIP and `SHA256SUMS.txt` from the same release, then run:

```powershell
.\install.ps1 -PackagePath .\DynSysLab-0.5.1-windows-x64.zip -ChecksumPath .\SHA256SUMS.txt
```

Alternatively extract the release ZIP and launch `bin\DynSysLab.exe` directly.
Both the application and CLI use the same native computational modules.

## Ask Codex to set it up

Give another Codex installation this instruction:

> Install AJGibson0926/dynsyslab-public from its latest official Windows x64 GitHub
> release. Read its README and docs/INSTALL.md, use install.ps1, verify the
> release checksum and native doctor, and run the bundled Duffing example in
> a new archive directory. Keep my existing scientific data.

Codex should use the prebuilt release for ordinary use. Clone the source only
when developing or inspecting code. A small checkout is sufficient:

```powershell
git clone --depth 1 https://github.com/AJGibson0926/dynsyslab-public.git
```

This project does not require the author's local Codex routing Hub or skills.
Read the repository's `AGENTS.md` for source-development conventions.

## Verify and run a first calculation

```powershell
$app = Join-Path $env:LOCALAPPDATA 'Programs\DynSysLab'
$data = Join-Path ([Environment]::GetFolderPath('MyDocuments')) 'DynSysLab\demo-archives'
& "$app\bin\dynsyslab-cli.exe" doctor --archive-root $data
& "$app\bin\dynsyslab-cli.exe" run "$app\resources\examples\duffing-free-solve.json" --archive-root $data
& "$app\bin\DynSysLab.exe"
```

The example generates fresh data. Large research datasets, case-study figures,
PDF exports, machine-specific audits, and personal configuration are not release
contents. Historical report links can point to local-only artifacts; consult
[local asset pointers](LOCAL_ASSET_POINTERS.json) and the report's job settings.
The public repository starts from a clean source snapshot. Historical scientific
media remain in the author's private repository and local archives.

## Develop and extend

See [native build instructions](BUILD_NATIVE.md), [pipeline boundaries](architecture/PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md),
and [eigenfunction/KROM APIs](architecture/EIGENFUNCTION_MODEL_API.md).
Windows desktop binaries are provided; optional scientific Python bindings have
their own platform build lane and are not a cross-platform desktop installer.
Third-party notices are bundled under `resources/third-party-licenses`.

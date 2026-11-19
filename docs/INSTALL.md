# Install DynSysLab

Repository: **AJGibson0926/dynsyslab-public**. The native desktop release supports
Windows x64. It includes Solver/Post-processor, Viewer, the scientific CLI,
numerical libraries, and the C++ extension SDK. No Python, compiler, Git,
administrator account, or scientific datasets are required to run it.

## Install and run (Windows x64, no admin, no Python)

```powershell
Invoke-WebRequest https://github.com/AJGibson0926/dynsyslab-public/releases/latest/download/install.ps1 -OutFile install.ps1
powershell -ExecutionPolicy Bypass -File .\install.ps1
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" quickstart
```

**Running DynSysLab means `dynsyslab-cli quickstart`.** It prints a welcome and overview, runs
the bundled examples (a Duffing double well and the chaotic Lorenz system), and reports the
results as Markdown: a table of the runs (samples, solver steps, dominant frequency, final
state), archive paths, and PNG figures rendered by the Viewer CLI. It opens no window and takes
a few seconds. `--format json` gives the same report as `dynsyslab.quickstart/v1` JSON.

The installer downloads the release ZIP and its SHA-256 manifest from
[GitHub Releases](https://github.com/AJGibson0926/dynsyslab-public/releases/latest),
checks the ZIP before extraction, installs under `%LOCALAPPDATA%\Programs\DynSysLab`, creates
Start-menu entries for both applications, and runs `quickstart` once at the end
(`-NoQuickstart` skips it). It does not start the desktop applications (Solver/Post-processor
and Viewer); open them from the Start menu when you want them. Close them before running the
installer to update. `release-install.json` in the install folder records the version, the
CLI, the welcome document and the quickstart command.

New installations store data under `Documents\DynSysLab\archives`, and every CLI command uses
that folder by default, wherever it is run from (`doctor` reports the folder and where the
setting came from). Existing archive settings are preserved. `-DataRoot C:\MyScience\DynSysLab`
chooses a different data workspace; `-InstallRoot` chooses the application location. The
installer does not delete scientific data. Verified downloads are retained in the temporary
directory reported at completion.

Pin a version with `-Version v0.6.2`. For offline installation, download the
Windows ZIP and `SHA256SUMS.txt` from the same release, then run:

```powershell
.\install.ps1 -PackagePath .\DynSysLab-0.6.2-windows-x64.zip -ChecksumPath .\SHA256SUMS.txt
```

Alternatively extract the release ZIP and run `bin\dynsyslab-cli.exe quickstart` directly.
The applications and the CLI use the same native computational modules.

## Python client from Git (any OS)

The optional, dependency-free `dynsyslab` Python client drives the compiled CLI and can be
installed straight from this repository, or from the wheel attached to each release:

```powershell
python -m pip install git+https://github.com/AJGibson0926/dynsyslab-public
```

It needs a compiled DynSysLab installation (above) to run computations.

## Install with an AI agent

In Claude Code, Codex or another agent, ask:

> Install and run DynSysLab at https://github.com/AJGibson0926/dynsyslab-public/releases/latest

The release notes and this page give the three commands above. The installer's own output
then contains the welcome, the example results and the figure paths.

Use the prebuilt release for ordinary use. Clone the source only when developing or inspecting
code. A small checkout is sufficient:

```powershell
git clone --depth 1 https://github.com/AJGibson0926/dynsyslab-public.git
```

This project does not require the author's local Codex routing Hub or skills.
Read the repository's `AGENTS.md` for source-development conventions.

## Your own calculations

```powershell
$cli = Join-Path $env:LOCALAPPDATA 'Programs\DynSysLab\bin\dynsyslab-cli.exe'
& $cli doctor
& $cli run "$env:LOCALAPPDATA\Programs\DynSysLab\resources\examples\duffing-free-solve.json"
```

`run` accepts any `dynsyslab.compute-job/v1` job; `--archive-root PATH` writes elsewhere. The
bundled jobs are in `resources\examples`; the quickstart's are in `resources\examples\quickstart`.

The examples generate fresh data. Large research datasets, case-study figures,
PDF exports, machine-specific audits, and personal configuration are not release
contents. Historical report links can point to local-only artifacts; consult
[local asset pointers](LOCAL_ASSET_POINTERS.json) and the report's job settings.
The public repository starts from a clean source snapshot. Historical scientific
media remain in the author's private repository and local archives.

## Develop and extend

See [native build instructions](BUILD_NATIVE.md), [pipeline boundaries](architecture/PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md),
and [eigenfunction/KROM APIs](architecture/EIGENFUNCTION_MODEL_API.md).

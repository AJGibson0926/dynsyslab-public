# Build, test, and install native DynSysLab

DynSysLab's canonical scientific engine, workflow service, CLI, Solver and
Post-processor, and Viewer are compiled C++. The installed application does
not contain, locate, or launch Python. Python bindings are optional developer
adapters and are disabled in the release preset.

## Prerequisites for developers

End users should use the [prebuilt release installer](INSTALL.md). The
requirements below apply to source builds and extension development.

- Visual Studio 2022 with Desktop development with C++
- CMake 3.24 or newer
- vcpkg, with `VCPKG_ROOT` set to its root directory

The checked-in presets set `VCPKG_MANIFEST_DIR` to `native/`, where the
compiled dependency manifest lives. No manual package list or Python runtime
is required.

These are build tools. End users receive the portable/install tree and need
none of them.

## One-command developer installation

From the repository root:

```powershell
.\tools\Install-DynSysLabNative.ps1 -PinToTaskbar
```

The script configures the unified CMake build, compiles every native module,
runs CTest and the isolated-package smoke solve, installs to
`%LOCALAPPDATA%\Programs\DynSysLab`, creates a Start-menu shortcut, and
optionally writes the taskbar shortcut. Rebuilding updates the executable and
DLLs at the same stable location, so the shortcut does not go stale. It also
writes the native shared-root configuration at
`%LOCALAPPDATA%\DynSysLab\runtime.json`, allowing no-argument Solver and Viewer
launches to reopen the repository's canonical archive. Pass `-ArchiveRoot` to
bind an installation to another archive directory.

To install and refresh shortcuts from an already validated build tree without
reconfiguring or compiling it again:

```powershell
.\tools\Install-DynSysLabNative.ps1 `
  -ExistingBuildRoot C:\path\to\build -SkipTests -PinToTaskbar
```

## Manual build

```powershell
cmake --preset windows-fast-release
cmake --build --preset windows-fast-release
cmake --build --preset windows-fast-release-gates
ctest --preset windows-fast-release
cmake --build build/native-suite-fast --target dynsyslab_package_smoke
cmake --install build/native-suite-fast --prefix "$env:LOCALAPPDATA\Programs\DynSysLab"
```

The public executables are:

- `bin\DynSysLab.exe` — Solver and Post-processor (primary app)
- `bin\DynSysLabViewer.exe` — Viewer
- `bin\dynsyslab-cli.exe` — headless workflow CLI

The coarse runtime layers are foundation, delay systems, lattice systems, PDE
systems, Koopman, control, operations, and workflow DLLs. Delay, lattice, and
PDE models live in `dynsyslab_delay.dll`, `dynsyslab_lattice.dll`, and
`dynsyslab_pde.dll`, so changes to those scientific modules do not compile
Eigen-heavy Koopman code. Third-party
numerical and compression DLLs and the Microsoft C/C++ runtime are installed
beside the executables. This is a normal app-local deployment, not a Python
environment.

The release build preset compiles only the three public applications and their
runtime dependencies. Test executables are an explicit `windows-release-gates`
target, so ordinary incremental builds do not walk every test project. Eigen's
large template implementation remains isolated behind compiled numerical DLLs;
an attempted PCH benchmark was rejected because it increased cold-build time.

New independent systems and methods should be separate extension DLL targets,
created with `dynsyslab_add_extension()` from
`cmake/DynSysLabExtensions.cmake`. That keeps ordinary scientific additions
from recompiling or relinking the established core.

Archive auditing and the durable job store live in `dynsyslab_operations.dll`.
They are not part of the scientific `dynsyslab_core` aggregate, so changes to
operational code do not recompile or relink Solver/Post-processor workflow or
Viewer services. The compiled CLI links operations explicitly.

The installed SDK contains the public headers, import libraries, extension
CMake helper, and a complete oscillator example under
`resources/sdk/examples/native-extension`. Copy that example outside the source
tree and point `DYNSYSLAB_SDK_DIR` at the installed application to verify an
independent one-DLL build.

## Optional Python front end

Python is not part of the installed application. Users who install the
dependency-free optional wheel can submit the same versioned JSON jobs to the
canonical compiled workflow:

```powershell
python -m dynsyslab_client plan job.json
python -m dynsyslab_client run job.json --archive-root C:\data\dynsyslab-archives
```

Programmatic use through `dynsyslab_client.NativeWorkflow` loads
`dynsyslab_workflow_c_api` with `ctypes`; it does not run a second Python solver.
The adapter exposes progress callbacks, cancellation, thread limits, and memory
limits. The wheel contains only `dynsyslab_client` and has no NumPy or SciPy
dependency.

`validate`, `plan`, and `run` use the compiled registries by default. Native
plans report the Batch/System/Cluster/Trajectory counts, sample count, coordinate
system, and peak/total solution-storage estimates without writing an archive.
This means new compiled controllers or post-processors are not blocked by an
older Python-side allow-list.

The source-checkout `python -m dynsyslab` front door delegates to that same
service. It is not part of the distributed Python client.

## Release gates

Native CTest covers systems, ODE, DDE, PDE, and lattice execution, archives, spectral
methods, KROMs, control, estimation, the durable job store, both C ABIs,
extension loading, workflow execution, and Viewer payloads.
`dynsyslab_package_smoke` installs an isolated tree, rejects any bundled Python
runtime, initializes both packaged desktop applications, starts the packaged CLI, performs real archived Duffing,
Mackey--Glass, Potts, and locally partitioned heat-equation solves, and
compiles and loads the installed SDK example.

## Continuous-integration lanes

The hosted checks deliberately separate the shipped product from migration
and adapter evidence:

- **Native product** configures the top-level CMake project once, builds both
  desktop applications and every runtime DLL, runs the full native CTest
  suite, enforces the machine-readable UI parity contract, exercises an
  isolated installation with no Python runtime, and publishes the portable
  Windows package. This is the release-authoritative lane.
- **Thin Python client** builds the dependency-free convenience wheel and
  verifies outside the source tree that it contains no historical solver,
  control, post-processing, ROM, or UI implementation.
- **Optional native bindings** builds the small cross-platform Python binding
  for developers who explicitly want direct Python access. It is not needed by
  either desktop application.
- **Reference parity** installs NumPy and SciPy only to preserve comparison
  evidence while historical Python implementations are retired. Passing this
  lane does not make Python part of the product architecture.

SUNDIALS behavior is tested by the canonical compiled integrator and package
smoke gates. CI must not reinstall the historical `sundials4py` execution path
as though it were part of the native product.

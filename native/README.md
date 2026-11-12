# DynSysLab native core

`dynsyslab_core` is the canonical, independently compilable scientific engine.
It has no Python headers, interpreter, callbacks, or runtime dependency. Native
desktop applications and the compiled CLI link this library directly.

Python is an optional client. The source-tree compute CLI uses `ctypes` and the
stable `dynsyslab_workflow_c_api` ABI for complete jobs. When
`DYNSYSLAB_BUILD_PYTHON=ON`, the private `_dynsyslab_native` extension also
provides low-level pybind11 kernel access. Algorithms do not live in either
binding layer.

## Python-free build

```powershell
cmake -S native -B build/native-core `
  -DDYNSYSLAB_BUILD_PYTHON=OFF -DDYNSYSLAB_BUILD_CLI=ON
cmake --build build/native-core --config Release
ctest --test-dir build/native-core -C Release --output-on-failure
build/native-core/Release/dynsyslab_native_cli.exe capabilities
```

The current native core includes all registered ODE systems, dimensional and
nondimensional model contracts, controlled ensemble integration, fixed-step
Euler/RK4, adaptive Dormand--Prince RK45, and optional statically linked CVODE
BDF for stiff systems. It also owns NPZ/LZMA archive compatibility, SVD/PCA/POD,
HAVOK, DMD/H-DMD/EDMD, polynomial SINDy eigenfunctions, Koopman models, exact and regressed
control maps, continuous-time LQR/KLQR, control-map error fields, classical
observability analysis, and continuous-time LQE/KLQE. The completed
migration plan records the implemented feature packages and verified parity gates.
PCA, POD, seeded bagged PCA, principal-component-pursuit robust PCA, and
Hankel POD share this compiled decomposition boundary and native archive path.

Enable all compiled scientific dependencies with:

```powershell
cmake -S native -B build/native-core `
  -DDYNSYSLAB_ENABLE_ARCHIVE_IO=ON `
  -DDYNSYSLAB_ENABLE_LINALG=ON `
  -DDYNSYSLAB_ENABLE_SUNDIALS=ON `
  -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows-static
```

## Optional low-level Python binding

```powershell
python -m build --wheel native
```

The binding wheel has its own Python build dependency, but the C++ library,
native CLI, Viewer, and Solver/Post-processor do not. Full workflow calls use
the dependency-free adapter in `compute_workbench.native_workflow`. The stable
C++ headers are installed from `native/include/dynsyslab/core`.

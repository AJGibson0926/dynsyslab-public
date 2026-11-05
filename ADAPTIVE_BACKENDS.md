# Adaptive integration in native DynSysLab

Adaptive ODE integration is a compiled engine capability. The shipped runtime
uses SUNDIALS/CVODE directly through `dynsyslab_foundation.dll`; it does not
load Python, SciPy, `sundials4py`, or a Python callback during a solve.

## Available methods

- Dormand--Prince RK45 is the general adaptive explicit method.
- CVODE BDF is the canonical stiff method and the closest supported counterpart
  to MATLAB `ode15s` for the retained comparison workflows.
- Fixed-step Euler and RK4 remain available when a fixed sample clock is part
  of the experiment definition.

Every method uses the same native workflow job, cancellation, progress,
archive, control-history, and provenance contracts. Output sampling is distinct
from internal adaptive steps. In particular, continuous-RHS LQR, SDRE, and KLQR
re-evaluate their state-dependent laws at internal RHS evaluations; sampled
control uses the configured zero-order-hold clock instead.

## Build and verification

SUNDIALS is declared in `native/vcpkg.json` and linked into the native runtime.
The normal top-level build and package smoke gate verify that its required DLLs
are installed beside the applications:

```powershell
cmake --preset windows-fast-release
cmake --build --preset windows-fast-release
cmake --build --preset windows-fast-release-gates
ctest --preset windows-fast-release
cmake --build build/native-suite-fast --target dynsyslab_package_smoke
```

Historical Python adaptive backends remain source-only comparison evidence
while retirement work proceeds. They are not selectable from the shipped
applications, not included in the optional thin-client wheel, and not a release
dependency.

Primary upstream references:

- https://sundials.readthedocs.io/
- https://github.com/LLNL/sundials

Admission to `auto` requires cross-platform wheel tests, numerical tests for
stiff and nonstiff systems, lifecycle parity, and a representative end-to-end
performance advantage over SciPy on the same runtime.

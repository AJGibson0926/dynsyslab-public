# Native extension example

This example adds a two-state harmonic oscillator and an energy observable to
an installed DynSysLab without rebuilding DynSysLab itself. It includes only
the stable extension ABI; Eigen and private core headers are not required.

```powershell
cmake -S . -B build -G Ninja `
  -DDYNSYSLAB_SDK_DIR="$env:LOCALAPPDATA\Programs\DynSysLab"
cmake --build build
```

The output is one independent DLL. Add its path to a compute job's
`extensions` array, then use `example_oscillator` as the system name. Rebuilding
this example does not relink any DynSysLab executable or runtime DLL.

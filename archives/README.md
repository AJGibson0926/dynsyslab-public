# archives

The default archive root of a source checkout: solver output, derived objects, hierarchies and
job records written by the native suite. Everything here except this file is local and
git-ignored; archives are regenerated from job specifications, never committed.

The directory must exist: the command-line tools recognize a source checkout by its `native/`,
`docs/` and `archives/` folders (`compute_workbench/native/src/CliRuntimePaths.cpp`). An
installed release keeps its archives under `Documents\DynSysLab\archives` instead.

# Optional Python client

DynSysLab does not require Python. The optional `dynsyslab` wheel contains one
dependency-free package, `dynsyslab_client`, which calls an already installed
native DynSysLab engine. It contains no equations, solver, archive engine,
post-processing algorithm, controller, estimator, or plotting implementation.

## Install and locate the engine

Install native DynSysLab first. Then install the optional Python wheel or an
editable checkout:

```powershell
python -m pip install dynsyslab
# source checkout:
python -m pip install -e .
python -m dynsyslab_client doctor
```

The client checks `DYNSYSLAB_INSTALL_DIR` first. On Windows it otherwise uses
`%LOCALAPPDATA%\Programs\DynSysLab`. `DYNSYSLAB_WORKFLOW_LIBRARY` may name an
exact workflow C-ABI library for development or embedding tests. It never
searches `PATH` for an arbitrary similarly named DLL.

## CLI

The `dynsyslab` command launches the native applications or passes commands to
`dynsyslab-cli`:

```powershell
dynsyslab doctor
dynsyslab ui
dynsyslab viewer
dynsyslab systems
dynsyslab plan job.json
dynsyslab run job.json --archive-root C:\data\dynsyslab-archives
dynsyslab archive audit --all-ages
dynsyslab queue list --archive-root C:\data\dynsyslab-archives
```

Exit status and standard output come from the compiled CLI. Scientific arrays
remain in native archives; the client does not deserialize and recompute them.

## Programmatic workflow binding

```python
from dynsyslab_client import NativeWorkflow

workflow = NativeWorkflow()
plan = workflow.plan(job)
result = workflow.run(job, r"C:\data\dynsyslab-archives")
```

`NativeWorkflow.run` also accepts progress and cancellation callbacks plus
thread and memory limits. These callbacks cross the versioned C ABI; Python is
the caller, never the numerical execution runtime.

## Source-reference code

The repository temporarily retains historical Python scientific packages for
parity evidence. They are excluded from the wheel. Maintainers working on an
explicit migration comparison may install their numerical dependencies with
`pip install -r requirements/reference-parity.txt`; new product behavior must
go into the native API and cannot be added to those reference implementations.

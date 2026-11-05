# DynSysLab tutorial: solve, analyze, and inspect

This tutorial follows the shipped, Python-free workflow. The Solver and
Post-processor, Viewer, native CLI, and optional thin Python client all use the
same compiled scientific services and archive tree.

## 1. Install and launch

From a developer checkout:

```powershell
.\tools\Install-DynSysLabNative.ps1 -PinToTaskbar
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLab.exe"
```

An end-user portable package starts `bin\DynSysLab.exe` directly. No Python,
compiler, source checkout, or environment activation is needed after the native
package has been built.

The main window has separate **Solver** and **Post-processor** panels over one
system-centric archive tree. **Open in Viewer** launches the compiled Viewer on
the selected object. **Refresh** re-indexes systems, presets, KROMs, and archive
relationships without changing any data.

## 2. Solve a free Duffing trajectory

In the Solver panel:

1. Select `duffing` and the `double_well_reference` preset.
2. Choose explicit initial conditions and enter `[-1.0, 0.1]`.
3. Set the time interval, output step, and an RK45 or fixed-step RK4 method.
4. Leave the controller at `none`.
5. Select **Plan / validate**. Review the resolved system, trajectory count,
   sample count, coordinates, and storage estimate.
6. Select **Run and archive** only after the plan is correct.

The solve publishes an immutable trajectory archive with the resolved preset,
parameter overrides, initial condition, solver tolerances, coordinate policy,
and data role. Viewer can open it immediately as a time series, phase portrait,
Fourier view, or scalar-observable view.

The archive stores only the time coordinate and authoritative sample values.
Select **Fourier** to compute a frequency-domain view on demand, or **Time +
frequency** to place each state's time-domain and Fourier plots in one display.
For controlled trajectories the control history remains first, followed by
state time series and then spectra. The native O(N log N) result is cached only
for the selected in-memory `TimeSeries`; it is not added to the archive.

The same operation is available without the UI:

```powershell
dynsyslab-cli fourier trajectory.npz --component 0 --window hann
```

## 3. Generate an ensemble or campaign

Initial-condition generation is part of the same solver definition:

- `default` uses the model default; a controlled system may define a distinct
  controlled default such as dimensionless RPE equilibrium `(1,0)`.
- `explicit` specifies every state or field directly.
- `uniform` and `normal` produce deterministic seeded samples from component
  ranges.

Set trajectory and cluster counts to form the
Batch/System/Cluster/Trajectory hierarchy. Use **Build campaign** when system,
controller, time, solver, or IC values should vary across a sweep. Campaign
axes change the configured values; they do not create a separate solver or
archive format.

Data roles such as free training, forced identification, control attempt, and
forcing campaign are metadata on the same resolved dynamical system. Distinct
physical parameter sets remain distinct system variants. Dimensional and
nondimensional records may be related only when their metadata proves they are
representations of the same physical dynamics.

## 4. Configure control correctly

The native controller selector includes prescribed waveforms, random pulses,
PID, LQR, state-dependent LQR (`sdre`), finite-horizon MPC, KLQR, and KMPC.
KLQR/KMPC additionally select a KROM and its exact, constant, or bilinear
lifted control map.

The control coordinate convention always follows the dynamics. A dimensionless
RPE solve receives dimensionless target, controller, and control-history
coordinates; a dimensional solve remains dimensional. The archive records the
convention and scale metadata so Viewer labels the control plot from evidence,
not mutable defaults.

Use **Controller evaluation** to distinguish:

- `continuous_rhs`: recompute a stateless state-dependent controller at every
  internal adaptive RHS evaluation; and
- `sampled`: update on an exact integer-tick clock and hold the value between
  updates.

Automatic evaluation keeps the established KLQR behavior continuous. PID,
SDRE, MPC, and KMPC use a sampled clock. Control output and target histories are
saved separately beside the solution; Viewer overlays the target on the state
plots and gives the actual driver its own time-series panel.

## 5. Reuse and queue validated jobs

**Save job…** writes the active panel as a versioned JSON job. **Open job…**
accepts that same document for read-only planning or explicit execution.

**Native jobs…** submits the current job to the compiled SQLite-backed queue,
lists its exact manifest and lifecycle, runs the next item asynchronously,
cancels queued or running work, and recovers a stale worker after an explicit
age. The footer reports throttled progress without copying numerical states or
querying the database at every RHS call.

**Archives…** audits duplicate candidates read-only by default. Applying a
cleanup requires confirmation and moves only byte-verified, unreferenced
duplicates to recoverable archive-local storage.

## 6. Post-process an existing solve

Select a trajectory in the archive tree and switch to **Post-processor**. A
post-processing job consumes the selected archive and never reruns its dynamics.

For the standard Koopman workflow:

1. Choose **Configure standard H-DMD → SINDy**.
2. Select any real scalar observable published by the system. The conventional
   default is the first state component, but radius, echo, energy, surface
   tension, or another registered observable is equally admissible.
3. Set H-DMD matrix size, sampling window, stride, SVD tolerance, rank policy,
   and energy sorting. Rank zero retains the full admissible spectrum.
4. Set SINDy observable coordinates, reconstruction observables, polynomial
   degree, sparsity threshold, and derivative policy.
5. Plan and run the selected stages.

The discrete H-DMD spectrum and continuous generator spectrum remain distinct;
the workflow applies `lambda = log(mu) / delta_t` at their explicit boundary.
The current native SINDy stage stores coefficient-first analytic
eigenfunctions. H-DMD itself stores its aligned trajectory coordinates as
scattered-interpolation eigenfunctions and its modal reconstruction row for
the selected scalar observable. That archive is therefore a directly
inspectable and extractable KROM; finite-difference gradients make its exact
`B_Psi` available to the standard control and regression services.

To validate an existing KROM, select a free trajectory, enable **Validate
KROM**, and open **Configure KROM validation…**. Choose a KROM for the same
system, mark the trajectory as `training` or `test`, and either leave the
semigroup horizons blank or enter explicit model-time values. This stage never
refits the KROM. It archives per-mode generator, one-step, anchored long-term,
and finite-horizon errors together with training-domain extrapolation coverage.
Viewer displays the comparisons and a mode-by-horizon error heat map.

To construct fitted-boundary eigenfunctions instead of SINDy polynomials,
enable **Boundary eigenfunctions** and choose its supplied spectrum. The first
selected trajectory is the primary source; additional same-system trajectories
may be listed in the structured dialog. Select an observable, and provide mode
groups only for a multi-output observable. This stage keeps the eigenvalues
fixed, fits only per-trajectory modal boundary amplitudes, and publishes an
independent scattered-interpolation KROM with its observable reconstruction
matrix.

The panel also exposes DMD, EDMD, HAVOK, SVD, PCA, POD, bagged PCA, robust PCA,
Hankel POD, KROM extraction, KROM validation, control-map regression/error
grids, observability, and KLQE. Every stage has independent settings and
publishes a new immutable derived object with direct lineage to its inputs.

## 7. Inspect results in Viewer

Select **Open in Viewer** from either panel, or launch
`bin\DynSysLabViewer.exe` directly. Viewer organizes records as:

```text
system → physical parameter variant → dataset/campaign → derived object
```

It exposes the relationships between training data, KROMs, extracted KROMs,
controlled regression data, fitted control maps, KLQR attempts, observability
analyses, and KLQE results. Relevant native views include:

- trajectory time series, phase portraits, control histories, targets, and
  scalar-observable time series or 2-D heat maps;
- discrete/continuous Koopman spectra and standard eigenfunction dashboards;
- KROM generator/one-step/long-term validation and finite-horizon error maps;
- constant/bilinear versus exact `B_Psi` error heat maps;
- decomposition spectra, modes, scores, low-rank/sparse traces;
- PDE space-time fields and lattice initial/final configurations; and
- estimator measurements, reconstructed state, and observability directions.

Use **Copy plot** or **Save PNG** for the active plot. Axis and square-aspect
controls are view-level operations and do not mutate the archive.

## 8. Use the native CLI for repeatable automation

The installed CLI consumes the same job files as the UI:

```powershell
$cli = "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe"
& $cli doctor
& $cli capabilities
& $cli systems
& $cli validate examples\native\duffing-free-solve.json
& $cli plan examples\native\duffing-free-solve.json
& $cli run examples\native\duffing-free-solve.json
& $cli queue list
& $cli archive audit --all-ages
```

Planning and validation are read-only. Running publishes through the same
workflow service as the desktop app. Unsupported family-specific operations
fail explicitly; the CLI never falls back to Python.

## 9. Optional Python access

The dependency-free `dynsyslab_client` wheel is a convenience adapter for an
already installed native suite:

```powershell
python -m pip install dynsyslab
dynsyslab doctor
dynsyslab plan job.json
dynsyslab run job.json --archive-root C:\data\dynsyslab-archives
```

Programmatic callers may use `dynsyslab_client.NativeWorkflow` through the
versioned C ABI. Python remains the caller; the compiled workflow owns solving,
post-processing, control, progress, cancellation, and archiving.

Historical Python scientific packages in the repository are migration-only
reference evidence. They are not installed by the wheel and are not a second
supported implementation. New systems and algorithms belong in separate native
extension DLLs using the installed SDK example.

## 10. Reproducibility checklist

Retain, and normally let the archive record:

- system, preset, full resolved parameters, and dimensionless groups;
- initial conditions, seed, batch/cluster membership, and campaign axes;
- controller, target, evaluation policy, KROM, and control-map representation;
- time span, output mesh, coordinate system, solver, and tolerances;
- post-processing algorithms, observables, dictionaries, ranks, and sorting;
- source archive identifiers and complete derived-object lineage; and
- software/schema versions and resource limits.

For a detailed field-by-field guide, see
[`SOLVER_POSTPROCESSOR.md`](SOLVER_POSTPROCESSOR.md). For application
architecture and the enforced UI coverage ledger, see
[`architecture/NATIVE_UI_PARITY.md`](architecture/NATIVE_UI_PARITY.md).

# MATLAB RPE → H-DMD → SINDy parity exporter

This external wrapper runs the existing MATLAB DynSysLab suite without
editing it. It writes its results outside the MATLAB suite as:

- an authoritative MATLAB `parity` struct in a v7.3 MAT file;
- a compact, readable JSON manifest beside it.

The default export contains the essentials needed to compare MATLAB with
the Python/C++ DynSysLab implementation:

- exact RPE physical parameters, derived dimensionless groups, and natural-period scaling;
- exact dimensionless initial condition, uniform time vector, state trajectory, and solver tolerances;
- the radius samples and exact Hankel construction used by H-DMD;
- H-DMD matrix fingerprints, singular values, retained rank, full discrete and continuous spectra, projected eigenfunctions, terminal amplitudes, and energy ordering;
- the SINDy state/derivative data, `Theta` and `Gamma`, ordered monomial exponent table, coefficients, fitted eigenfunction values, residuals, and spectrum conversion policy;
- MATLAB/suite provenance and wall-clock timings for trajectory integration, H-DMD, and SINDy.

## Canonical run

In MATLAB:

```matlab
addpath("scripts/matlab_parity")

[matFile, jsonFile] = export_dynsyslab_rpe_hdmd_sindy( ...
    OutputFile="<local-data>/matlab_rpe_hdmd_sindy_parity.mat");
```

The defaults are the canonical comparison configuration:

- initial state `[2.3333865986017348, 0]` (dimensionless);
- `tau = 0:0.001:25`;
- `ode15s` through `dynbatch.compute("soft")`;
- `RelTol = AbsTol = 1e-9`;
- H-DMD rows `20000`, columns `2000`, horizontal radius embedding;
- SVD absolute tolerance `1e-10`, full retained spectrum, energy sorting;
- SINDy degree 4, no constant term, no library normalization, one-pass `1e-6` coefficient threshold.

The suite's `Koopman.process` implementation contains code to handle a
`tol` option but does not declare that option in its MATLAB arguments block.
Consequently, its public H-DMD pipeline can only use the internal `1e-10`
default. The wrapper enforces that value and reports a clear error if a
different `HDMDTolerance` is requested.

This is intentionally a long calculation. `SaveHDMDBasis=false` avoids
duplicating the large Hankel `X/Y` and `U/V` arrays in the MAT file while
still saving enough boundary data and norms to verify their construction.
Set it to `true` only when a byte-level investigation needs those matrices.

## Export an already-computed MATLAB system

If `rpe` already contains the trajectory you want to compare:

```matlab
[matFile, jsonFile] = export_dynsyslab_rpe_hdmd_sindy( ...
    ExistingSystem=rpe, ...
    OutputFile="C:\comparison\existing_rpe_parity.mat");
```

The wrapper will not alter the trajectory or physical parameters in that
object. H-DMD and SINDy are still run to produce the parity bundle. For the
cleanest isolation, use the canonical form, which creates a fresh object.

## Small smoke test

This verifies connectivity and data extraction without running the full
20,000×2,000 decomposition:

```matlab
export_dynsyslab_rpe_hdmd_sindy( ...
    HorizonTau=0.1, DTau=0.001, ...
    HankelRows=50, HankelColumns=20, ...
    OutputFile=fullfile(tempdir,"dynsyslab_parity_smoke.mat"));
```

The wrapper refuses to place output inside the MATLAB suite root. Existing
MATLAB source and data files are never written or modified.

## MATLAB KLQR control reproduction

`run_matlab_rpe_klqr_control_experiment.m` is a separate diagnostic wrapper
that extends the experiment through the specific successful MATLAB KLQR
pathway. It is deliberately not presented as a desired Python/C++ design.

```matlab
[matFile,jsonFile,figureFile] = run_matlab_rpe_klqr_control_experiment( ...
    OutputFile="C:\comparison\matlab_rpe_klqr_control.mat");
```

It reproduces the pathway in `generate_training.m` and
`control_dynamics.m`: the first complex SINDy eigenfunction/conjugate pair,
the MATLAB tuning-table values, state-dependent `B_psi`, LQR recomputation
at every adaptive-solver controller call, equilibrium control IC `[1,0]`,
and the two-tone target centered at radius 4. The MAT bundle contains the
uniform solution, target, control history, errors, metrics, and raw
controller-call trace. A PNG provides an immediate visual check.

## Korda--Mezic Duffing controlled-data parity candidate

`export_korda_mezic_duffing_control_candidate.m` is a non-mutating wrapper
around the data-generation and identification pathway in the separate
Korda--Mezic `Duffing.m` example. It uses the supplied precomputed eigenvalues,
reproduces the published random stream, exports free and random-forced data,
fits the published discrete constant lifted-input matrix, and exports an exact
continuous `B_psi` field for independent Python comparison.

The wrapper intentionally reproduces two source-specific details and records
them in metadata: the fourth RK stage advances from `k1` rather than `k3`, and
the source input contributes `0.5*u` to acceleration. Neither convention is
applied to DynSysLab's standard Duffing solver.

```matlab
addpath("scripts/matlab_parity")

summary = export_korda_mezic_duffing_control_candidate( ...
    "C:\comparison\korda_mezic_duffing_random_forced_full.mat");
```

The source script's interactive `keyboard` statements, MPC solve, and source
files are not invoked or modified. Generated MAT/JSON experiment bundles live
under `archives/experiments/`, which is intentionally ignored by Git.

## Cleanup-suite Duffing `duffrom2` parity

`export_playground_duffing_lambdaopt_sindy_regression.m` reproduces the exact
30-trajectory `playground.m` path without invoking `data_node.load`, a file
picker, or the unconditional `keyboard` statement in `KROM.get_B_parts`.
The authoritative retained ROM is `duffrom2`; `rom2` is not part of this
comparison.

```matlab
addpath("scripts/matlab_parity")

export_playground_duffing_lambdaopt_sindy_regression( ...
    "archives/experiments/matlab_duffing_exact_parity/matlab_playground_duffing_lambdaopt_sindy_regression_suite_hard.mat")
```

The wrapper uses `dynbatch.compute("hard")`; this is essential because the
legacy lambda optimizer amplifies trajectory perturbations at roughly machine
precision into different local spectra. It exports the exact SINDy domain,
polynomial coefficients, `C`, controlled lifted arrays, and all three regressed
matrices. Source/input hashes are verified unchanged after execution.

Run the corresponding Python comparison with:

```powershell
python scripts/run_duffing_duffrom2_matlab_python_parity.py
```

That runner retains a native Viewer spectrum comparison and keeps MATLAB's
parameter-loss and mixed-conjugation behavior inside an opt-in candidate
compatibility branch.

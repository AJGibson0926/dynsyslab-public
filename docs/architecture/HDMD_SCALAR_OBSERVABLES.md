# H-DMD scalar-observable contract

The native pipeline also supports a shared multichannel decomposition; see
[Joint H-DMD and linear diagnostics](JOINT_HDMD_AND_LINEAR_ANALYSIS.md).
The scalar Python reference contract below remains unchanged.

H-DMD delay windows are built from a scalar signal, not implicitly from the
first state component. Radius/position remains the standard operating default,
but it is one selection from the same public contract as echo, energy, pressure,
surface tension, or another system-owned scalar observable.

## Supported inputs

`hankel_dmd(...)` and `hankel_dmd_ensemble(...)` accept exactly one of:

- `observable=<integer>` for a state component;
- `observable=<catalog key>` for a named observable owned by the archived or
  explicitly supplied dynamical system;
- `observable=<Observable>` for an explicit state-to-scalar map; or
- `observable_series=<array or TimeSeries>` for a pre-evaluated real scalar
  signal on the exact state-trajectory time grid.

The scalar contract is active with `matrix_layout="window-snapshots"`. Supplying
an observable to the all-state `delay-state` layout is an error rather than a
silently ignored option. Vector-valued and complex-valued observables are also
rejected explicitly by the current real H-DMD backend.

Every result records `metadata["preprocessing"]["observable_contract"]`. Named
state observables use `state-observable/v1`; recorded signals use
`sampled-scalar-observable/v1`. The contract includes name, value type, state and
output dimensions, physical descriptor/units, coordinate system, and whether a
state callable exists. A raw recorded signal does not claim that it can later be
evaluated away from its trajectory.

## Named observables and coordinates

Catalog names are resolved by the owning dynamical system, never by a global
string alias. Archive workflows reconstruct the system from immutable source
metadata and evaluate the observable in the same integration coordinates used
for H-DMD. Examples include:

- RPE: `radius`, `radial_velocity`, `far_field_acoustic_echo`,
  `inertial_kinetic_energy`, `internal_gas_pressure`, `laplace_pressure`,
  `viscous_damping`, and `momentum_balance`;
- Marmottant: radius/velocity plus `surface_tension`, `shell_viscosity`, energy,
  pressure, damping, and momentum-balance terms; and
- Duffing: position/velocity plus kinetic, potential, and total energy.

Marmottant surface tension deserves special care. A catalog evaluation is a
non-mutating pointwise representation. A trajectory that crosses its hysteretic
branches should record the actual surface-tension history during integration and
pass that `TimeSeries` through `observable_series=`. This preserves the branch
history rather than attempting to infer it from radius alone after the run.

## Spectrum comparison

`compare_hankel_observables(...)` runs one scalar H-DMD decomposition per
observable while holding the state training data and all numerical settings
fixed. It stores both unitless discrete-time eigenvalues and continuous-time
eigenvalues, the exact observable contract for every spectrum, and symmetric
nearest-spectrum RMS distances. `observable_spectrum_comparison_report(...)`
produces a normal Viewer report with one color per observable; archived
comparisons are recognized by the native Viewer catalog and workbench.

The installed CLI provides:

```text
dynsyslab-postprocess observables SOURCE
dynsyslab-postprocess compare-observables SOURCE --observables NAME [NAME ...] --execute
```

The comparison is descriptive. Differences between spectra are expected when
observables illuminate different Koopman-invariant subspaces; the tooling does
not imply that one observable is universally superior.

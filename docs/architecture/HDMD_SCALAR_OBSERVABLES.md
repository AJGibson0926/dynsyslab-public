# H-DMD scalar-observable contract

The native pipeline also supports a shared multichannel decomposition; see
[Joint H-DMD and linear diagnostics](JOINT_HDMD_AND_LINEAR_ANALYSIS.md).

H-DMD delay windows are built from a scalar signal, not implicitly from the
first state component. Radius/position remains the standard operating default,
but it is one selection from the same public contract as echo, energy, pressure,
surface tension, or another system-owned scalar observable.

## Supported inputs

The native `hdmd` post-processing stage selects its scalar signal with
`settings.observable`: either an integer state-component index or the name of
an observable owned by the archived dynamical system. Vector-valued and
complex-valued observables are rejected explicitly by the real H-DMD backend.
The result records the observable it used.

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
branches cannot have its branch history inferred from radius alone after the
run.

## Spectrum comparison

To compare spectra from several scalar observables, run separate immutable
`hdmd` jobs that hold the source archive and all matrix, window, rank,
sampling, sorting, and linear-algebra settings fixed, changing only
`settings.observable`, and compare them in Viewer.

The comparison is descriptive. Differences between spectra are expected when
observables illuminate different Koopman-invariant subspaces; the tooling does
not imply that one observable is universally superior.

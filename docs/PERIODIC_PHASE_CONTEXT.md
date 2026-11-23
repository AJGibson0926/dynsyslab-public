# Periodic phase-space context

The optional native sidecar supports a fixed-forcing-phase return map for systems with a known sinusoidal clock. It uses the actual nonautonomous RHS. No analytic solution, autonomous approximation, global basin certificate or Melnikov conclusion is substituted for numerical return-map evidence.

## API and clock

Use dynsyslab.phase-space-query/v1 with options.dynamics = periodic-forcing. system.driver has amplitude, frequency_hz and phase_radians in the native system's documented physical-input/physical-clock convention. For abstract models these are model input/time. The driver is amplitude*cos(2*pi*frequency*time+phase). The native phaseSpaceForcingPeriod API resolves a common period, including Reddy–Szeri's internal acoustic-translation clock.

- forcing_period: positive period in model time; omitted/zero infers a common period. An explicit period must contain integer cycles of every active clock. Automatic common-period search is bounded at 64 multiples of the longest clock.
- forcing_section_time: absolute model time of the section, default zero. Return integrations start exactly at this time and end after whole forcing periods; times are not selected from arbitrary solver-output samples.
- forcing_periods: positive integer, default one, to search for a subharmonic response. Proper divisors are checked before reporting a minimal response period within this requested clock count.
- periodic_seeds: states on this forcing section. The legacy period field may be omitted/zero; it does not override the forcing clock. Without seeds, at most 32 domain-spanning guesses are used.
- section_seeds: arbitrary states on this fixed forcing phase, with no periodicity requirement. With poincare_map enabled, return each state after forcing_periods clocks and its full variational derivative. This does not run periodic shooting.
- manifold_iterations, maximum_manifold_points, manifold_tolerance: finite continuation, point and domain-scaled interpolation limits. The regular integration tolerances, bounds and cancellation apply as well.

The query supplies period-scaled defaults for maximum_step, output_step and horizon when omitted. Explicit maximum_step must resolve the forcing period with at least eight steps. The result retains normalized options and a forcing_clock record with the section time/phase, units, physical conversion and full driver definition. The reference clock is theta=2*pi*t/T with time/phase origin zero; the driver has its separately recorded phase offset. Curve times remain absolute and unwrapped.

## Numerical meaning

Shooting solves P^m(x)-x=0 at a fixed forcing phase. Newton uses the integrated state-transition matrix, domain scaling and damped corrections. The variational equation differentiates the actual forced RHS by centered differences, so it includes state-dependent input gains even when a system's supplied Jacobian describes only the drift. A second integration with smaller step/tolerances checks the endpoint and closure.

All Floquet multipliers belong to the forcing-section map. No multiplier is discarded as an autonomous neutral phase direction. Results distinguish attracting, repelling, saddle and nonhyperbolic periodic responses. A periodic response is stored under periodic_orbits with kind forced_periodic_response and section_state; it is never stored as an equilibrium. A failed finite shooting search does not establish absence.

For each one-dimensional stable/unstable subspace of a saddle periodic response, the sidecar advances an initial fundamental segment under the forward/inverse return map. It refines mapped segments using midpoint defects and spacing, and repeats with half the local seed displacement and tighter integration/interpolation. Negative multipliers use the squared map to preserve branch orientation. Curves failing the geometric comparison are retained as unresolved metadata and omitted by Viewer. Boundary/budget truncation remains explicit. These are approximate invariant curves in one forcing section, not time trajectories: their time array is empty and parameterization is forcing-section curve. They are not complete stable/unstable flow surfaces. Higher-dimensional section surfaces still require a separate continuation method.

Viewer draws periodic responses and accepted section curves only on state phase plots, with muted dashed lines. Section curves are labeled as such. Sparse arrows show the instantaneous vector field at the recorded section time. They do not replace the main trajectory or appear in time series or observable/Koopman heat maps.

## Archive ownership and Viewer adapter

Opening a plot performs lookup only. Phase context computes missing results explicitly, in the background through the native workflow C ABI. Cache identity includes equations, parameters, forcing definition, section time, numerical options, guide metadata and method version 5. Bibliography edits do not invalidate this numerical cache. The .dslphase cache is system metadata and stays outside the solution data tree. Existing solve archives are never modified. Version 5 adds arbitrary fixed-phase state returns and their derivatives; it retains the first-return and continuation records introduced in version 4, described in [Poincare sections and bifurcation sidecar](POINCARE_AND_BIFURCATION_SIDECAR.md).

The C++ API is `stroboscopicReturn` in `core/PhaseSpace.hpp`. The thin Python
adapter `dynsyslab_client.phase_space.query(request, archive_root)` invokes the
existing compiled phase-space CLI and canonical cache. It does not integrate in
Python. For these fixed-time maps `return_derivative` and `reduced_derivative`
both contain the full row-major state derivative, and `tangent_basis` is the
identity. State-section transversality fields are null because there is no
state-defined crossing condition. Eigenvalues at an arbitrary seed are local
Jacobian eigenvalues; they are Floquet multipliers only at a verified periodic
point of the requested return map. Endpoint and derivative refinement failures
are explicit; these checks are numerical, not interval enclosures.

Compute's continuous single-tone waveform controller can be reconstructed directly from controller_configuration. The adapter validates the production controller settings, converts sine phase to the cosine driver, and respects separate physical/model amplitudes and frequencies. Sampled controllers, active windows, nonzero offsets and multiple tones are currently rejected instead of being silently reduced to an incorrect waveform. Feedback-control paths retain their existing explicitly labeled zero-input context.

## Evidence and consulted material

Native tests use an independently soluble, periodically translated curved saddle with time-varying stable rate. They verify the full orbit, exact forcing-phase times, multipliers exp(-0.7) and exp(0.4), curved section manifolds, phase-shift invariance, smaller-period detection, actual viscous RPE forcing and physical/native clock conversion. The end-to-end test solves a waveform-controlled harmonic oscillator, reconstructs its exact sinusoidal steady response from the archive, checks cache isolation by section phase, validates raw archive immutability and bibliography exclusion, and renders Viewer without a window.

Read-only Grothendieck Labs/Class Hub consultation: MTH.dynamical-systems.periodic-forcing-phase-lift-and-viewer-audit; CODE.dynamical-systems.dynamical-systems-model-library (periodic_forcing.py); CODE.dynamical-systems.basin-manifold-toolkit; MTH.dynamical-systems.invariant-manifold-parameterization-and-invariance-equation. The native implementation is a reformulation with explicit residual/refinement tests, not a copy of an unvalidated global-manifold claim.

## Configured RPE saddle check

A retained native-model check uses vapor_pressure_pa=105000, the remaining default RPE parameters, and a 20 Pa, 100 kHz cosine drive with phase 0.31. The declared domain is r in [0.8,4], v in [-0.5,0.5], at forcing section time zero. Two responses were verified: an attracting response near (1.00111,-0.00453) and a saddle near (2.81498,-0.000452). The saddle multipliers were approximately 0.45660 and 2.08206; return residual was 1.3e-15 in domain-scaled coordinates. Four section branches passed seed/interpolation refinement, with geometric differences below 4e-6 before domain truncation.

The native regression test additionally checks the independent Liouville determinant identity det(M)=exp(-viscosity*integral(dt/r^2)-3*log(r_end/r_start)), derived from the RPE field divergence. This is a configured native-model verification, not a reproduction of a historical paper's normalization or parameters.

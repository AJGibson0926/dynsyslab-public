# Controlled periodic-coordinate transport

This generic API advances a real two-dimensional periodic chart under a scalar
control-affine input. It is an approximate nonlinear coordinate model; it is not
a finite linear Koopman controller or a measurement-based state estimator.

For a chart \(x=X(a,\theta)\) of unforced periodic orbits, with consistent
phase origins and frequency \(\Omega(a)\), the chain rule gives
\[
\begin{pmatrix}\dot a\\\dot\theta\end{pmatrix}
=\begin{pmatrix}0\\\Omega(a)\end{pmatrix}
+[X_a\ X_\theta]^{-1}B(X(a,\theta))u.
\]
Both input components are retained. The chart and input map can be approximate;
local nonsingularity does not establish global injectivity or model accuracy.
Construction rejects a sampled orientation reversal on the connected chart:
continuity implies a singular Jacobian between the conflicting samples. Passing
this check does not exclude an undetected reversal between sampling nodes.

## Ownership and APIs

- `math/PeriodicChart.hpp`: real tensor cubic Hermite chart, frequency and first
  derivatives. Phase slopes use the supplied Fourier coefficients. Parameter
  slopes use nonuniform three-point differences, one-sided at endpoints (secants
  for a two-knot chart). Returned
  derivatives are of this interpolant. Existing linear atlas behavior is unchanged.
  Conjugate-symmetry defects below `1e-10 * (1 + coefficient l1 norm)` are projected
  to an exactly real Fourier series before differentiation; larger defects reject
  the chart. The maximum per-channel l1 correction is reported separately for
  chart and input coefficients. This prevents differentiated FFT roundoff from
  being mistaken for a complex physical coordinate.
- `pullbackPlanarField`: a scaled two-by-two solve with a normalized determinant
  guard. Normalization removes column units, but depends on physical state units.
- `core/PeriodicTransport.hpp`: precomputes both input components, Fourier
  projects them and discards the physical input callback. It implements the
  standard `ODESystem` interface and uses the existing native integrators.
- Workflow `periodic-analysis`, operation `controlled_predict`: composes a
  registered physical input map, the reduced model, and an optional reference solve.
- Python `dynsyslab_client.periodic.predict_controlled`: thin CLI adapter.

The control API requires two state channels in the same order and units as the
registered physical system. Other chart channels are passive readouts. Input maps
may depend on physical state; only the orbit parameter and phase are integrated.
The periodic chart needs at least two strictly ordered parameter knots and
resolved real Fourier series. Extrapolation is disabled. The equilibrium and
singular chart points must be excluded or treated by another chart.

## Query contract

Use `schema: dynsyslab.periodic-query/v1`, `operation: controlled_predict`:

- `system`: registered `name`, optional `preset` and `parameters`.
- `knots`: existing periodic atlas contract, with `parameter` and Fourier
  `channels`. `state_channels` defaults to `[0,1]`.
- `parameter`, `initial_phase`, `initial_time` (default 0), `final_time`,
  `output_step`. All times and controls use native model units.
- `table_samples` (8192), `input_maximum_harmonic` (2047),
  `minimum_normalized_determinant` (1e-10). Require twice the harmonic limit
  smaller than table samples. Input projection and chart interpolation both need
  convergence checks; phase oversampling does not fix missing parameter knots.
- `relative_tolerance` (1e-9), `absolute_tolerance` (1e-11), `maximum_step` (0).
- `input`: `bias`; `sines` with rows `[amplitude,angular_frequency,phase]`;
  `state_gain` and `state_reference`, each length two; `limit` (0 means no bound).
  The law is bias plus sines plus gain dot (physical state minus reference), then
  saturated to the symmetric limit. `control_sample_time` >0 selects held sampled
  feedback; default zero evaluates the law at each integrator RHS call.
- `frozen_frequency` >0 replaces only the autonomous phase speed, as a diagnostic
  baseline. It does not freeze the chart or input components.
- `include_reference` (true) performs a separate registered-system solve under
  the same feedback policy. Its default initial state is the decoded chart IC.
  Optional `physical_initial_state` selects an explicit reference IC and reports
  the initial decoding discrepancy.

Results include flat sample-major `coordinates` (two columns), `predictions`
(`channels` columns), `controls`, and optional `reference` (two columns) and
`reference_controls`. `time` is shared. Numerical failure throws; failed rollouts
are not silently replaced with a physical solution.

Diagnostics sample chart conditioning and input reconstruction at knot and
midpoint parameters with shifted phase nodes. A separate midpoint audit compares
the decoded autonomous drift against the physical drift. These are finite
diagnostics, not uniform error bounds. Raw conditioning and errors depend on units.
Readout channels supplied as unforced observables do **not** acquire an input
feedthrough term automatically. A measurement with direct input dependence must
be evaluated with that dependence explicitly by its owning observation model.

## Sources and verification

Coordinate transformation is the classical chain rule; the action-angle setting
is described in Arnold, *Mathematical Methods of Classical Mechanics*, chapter 10.
State-dependent control terms in Koopman coordinates are discussed by
[Kaiser, Kutz and Brunton](https://arxiv.org/abs/1707.01146).
Interpolation follows the cubic Hermite construction in de Boor, *A Practical
Guide to Splines* (Springer, 1978).

Native tests check Jacobians against finite differences, phase wrapping, singular
chart rejection, and forced/feedback harmonic-oscillator trajectories against
independent closed-form solutions. They also verify the physical input callback
is never used during reduced rollout. These tests do not certify an arbitrary
user-supplied chart.

## Past-window identification under a known input

`core/ShootingIdentification.hpp` owns generic initial-state single shooting for
any registered or caller-supplied `ODESystem`. Callers supply an initial-state
parameterization, bounded starting candidates, an input-aware measurement map,
and a scalar input as a function of time. The existing native integrator and
bounded damped Gauss–Newton solver perform the trajectory solves and fitting.
The objective is the sum of squared measurement residuals divided by their
channel variances (independent Gaussian observation errors). See Nocedal and
Wright, *Numerical Optimization*, second edition (2006), chapter 10.

The workflow operation `controlled_identify`, and the thin Python adapter
`dynsyslab_client.periodic.identify_controlled`, compare `model: atlas` (default)
and `model: physical`. Both use the same orbit-parameter/phase parameterization
of the initial state. The physical model decodes that initial state through the
chart and subsequently integrates the registered equations directly. It needs
`observables`, ordered to match all atlas readout channels. The atlas evolves
only parameter and phase, using its precomputed input map.

Required inputs are `system`, `knots`, `time`, flat sample-major `measurements`,
`measurement_channels`, positive `measurement_variances`, and
`initial_candidates: [[parameter,phase], ...]`. The measurement clock must be
uniform, strictly increasing, with at least two samples. Parameters stay inside
the atlas; each starting phase has a search interval of plus/minus pi. At most
256 starts are accepted. The chart, integrator and frequency options above apply.
`maximum_iterations` and `finite_difference_step` configure local optimization.

`input` accepts the analytic `bias`, `sines`, and `limit` above, or the recorded
zero-order-held schedule described below. State-feedback laws are rejected:
the estimator must receive the actual applied input. The C++ API also accepts
caller-defined continuous input functions.

Measurements may include an explicit input term:
\[
y_j=h_{c_j}(z)+d_j u\,s(z).
\]
`measurement_feedthrough` supplies the coefficients \(d_j\), default zero.
`feedthrough_scale_channel` selects \(s=h_k\); default -1 means \(s=1\).
The same observation formula is used by both compared dynamics models.
`predictions` contains the unforced readout channels; `predicted_measurements`
contains their input-aware measured subset. They are deliberately separate.

`forecast_final_time` defaults to the last supplied observation. `output_step`
defaults to the observation spacing. Extending either forecast setting cannot
alter the fitting objective: only supplied past-window data enter the fit.
Curves inside the fitted window are reconstructions using that whole window;
they must not be labeled online filtered estimates. Later curves are free
forecasts given the known future input. This is not recursive Bayesian filtering,
KLQE, a finite linear Koopman model, or a global observability certificate.

Results include the best fit, all starting-candidate diagnostics, trajectory
evaluation count, measurement cutoff, decoded predictions, and predicted
measurements. Integration failures are recorded. A finite initial candidate is
retained if its optimizer fails, explicitly marked as unconverged; a failed start
does not suppress other valid candidates. Independent exact forced-oscillator
tests check state recovery, measurement feedthrough and forecast independence.

## Repeated causal identification

`core/SequentialShooting.hpp` reuses the same arbitrary-dimensional C++ model,
initial-state parameterization, known input, observation map and optimizer.
It refits from expanding past prefixes, then predicts only to the next scheduled
update. This is deterministic sequential identification, not a Kalman filter,
finite-window moving-horizon estimator, or uncertainty calculation.

Workflow operation `controlled_sequential_identify` and Python adapter
`identify_controlled_sequential` accept the identification contract above, plus:

- `time` and `measurements`: the complete observation stream for offline replay.
  Each optimization receives only its permitted prefix.
- `update_indices`: strictly increasing indices into that stream. Each must be
  at least 1 and smaller than the stream length. At most 10000 updates are allowed.
- `prediction_substeps`: output intervals per measurement interval, default 1.
  The output begins at the first update and ends at the last observation time.
  `forecast_final_time` and `output_step` are rejected for this operation.
- `retain_original_starts`: false by default. After the first multi-start fit,
  the previous best parameter vector initializes the next fit. If true, original
  starts are also retained, subject to the 256-start limit.

The previous optimizer's box is retained with its warm start. The fitted initial
state is propagated to the current cutoff, then forward until the next update.
At an update timestamp the new estimate replaces the old prediction. Earlier
emitted samples are never overwritten. The result includes `information_time`
for every output sample, together with every update's fit, failure diagnostics
and trajectory evaluation count. Samples before the first update are omitted:
reconstructed past trajectories must not be mistaken for online estimates.

The existing physical/chart alternatives and measurement feedthrough apply.
The workflow accepts analytic or recorded known input schedules. It does not infer
controls from hidden physical states. Repeated optimization can be expensive;
the API makes no real-time performance or process-noise robustness claim.
Exact forced-oscillator tests verify the output clock and solution. Independent
future-data mutation tests require earlier outputs to remain bit-for-bit equal
while later estimates respond to changed observations.

## Recorded applied input and forward replay

`core/KnownInput.hpp` supplies `KnownScalarInput::recorded(time, values,
beforeFirst, measurementBeforeCommand)` and `integrateKnownInput`. Continuous
callables remain valid inputs to the shooting APIs. Recorded times must strictly
increase; at least two finite matching samples are required. Values are the
actual applied commands, already limited or gated. The value at each timestamp
is held until the next timestamp. No extrapolation beyond the record is allowed.

At a command timestamp, propagation uses the new command. A measurement sampled
before that command uses the preceding value (or `beforeFirst` at the first
sample). At other timestamps both use the current held value. The final recorded
value defines the final post-command sample, not an unbounded future hold.
The integrator stops and restarts at every recorded boundary, holding the old
interval input through its endpoint; it preserves the requested output grid.
This is the standard ZOH sampled-data convention (Franklin, Powell and Workman,
*Digital Control of Dynamic Systems*, 3rd ed., 1998), with explicit sensor timing.

Identification and replay JSON accept:

```json
"input": {
  "time": [0, 0.2, 0.5, 1.0],
  "values": [0.1, -0.2, 0.05, 0.0],
  "before_first": 0.0,
  "measurement_timing": "before_command"
}
```

`before_first` is required. `measurement_timing` defaults to `before_command`;
`after_command` is available for data explicitly sampled after each update.
Analytic modifiers (`bias`, `sines`, `limit`) cannot accompany a recorded input.
Use the actual sensor/command timestamps; nearby floating-point timestamps are
not silently treated as simultaneous. All forecast times must be covered.

Operation `controlled_replay` and Python `replay_controlled` forward-solve the
physical or orbit model without fitting. Supply `initial_parameters` as orbit
parameter and phase, `time` endpoints, and `output_step`, plus the same system,
atlas, ordered physical observables and measurement mapping as identification.
Returned `applied_input` and `sensor_input` expose the two timing conventions.
Supplying `measurements` is rejected to avoid confusing forward replay with an
estimation operation. This is also useful for preparing controlled benchmarks.

Exact harmonic tests include an off-grid pulse much shorter than the output
interval, windows inside a held interval, pre/post-command direct measurements,
and future-input mutation that cannot change a past-window fit.

## Causal online observation and feedback simulation

`core/OnlineShooting.hpp` supplies `OnlineShootingObserver`. Each `observe(t,y,u)`
call accepts one new measurement and the actual scalar input held over the
preceding sample interval. It stores only this history. Before the configured
`firstFitIndex`, the returned state is empty and `ready` is false. At the first
fit, fixed candidates or a caller-supplied past-prefix initializer start the
existing bounded single-shooting solver. Subsequent fits occur every `fitEvery`
samples, warm-started from the preceding fit. Between fits the model propagates
under the newly supplied applied input. No plant state enters this API.

The returned `informationTime` is the last assimilation cutoff; current inputs
can still advance the state between assimilations. `forecast(end,u,intervals)`
returns a prospective constant-input trajectory without modifying the observer.
It does not consume future measurements. A failed solve throws; a finite best
candidate retained by the nonlinear optimizer carries its convergence flag and
diagnostic reason. This is expanding-window deterministic identification, not
a stochastic filter, finite-window MHE, or a real-time algorithm.

Operation `controlled_feedback`, exposed as Python `simulate_feedback`, connects
this observer to a separate physical plant, simulated noisy input-aware sensor,
and sampled affine feedback law. The plant's `initial_parameters` are never
passed to the observer. Required additions to the measurement/model contract:

- `time`: initial and final endpoints, containing an integer number of samples.
- `sample_time`, `first_fit_index`, `fit_every`, and `prediction_substeps`.
  The dense plant/forecast output interval is `sample_time/prediction_substeps`.
- `controller`: two `state_gain` entries, optional two `state_reference` entries,
  optional `bias`, positive symmetric `limit`, optional `start_index` and
  `every` (command cadence in sensor samples). An optional `acquisition_input`
  is held before feedback starts. Feedback cannot precede the first fit.
- `measurement_noise`: the existing streaming-noise specification. Noise is
  added after direct input dependence. `measurement_noise_scale` optionally
  supplies per-channel scales; explicit absolute noise levels avoid ambiguity.
- Either `initial_candidates`, or `coarse_initialization` with a parameter grid,
  `phase_count`, and `starts`. Coarse initialization runs only on the acquired
  first prefix, uses the unforced periodic filter, and requires zero acquisition
  input. Its states are candidate seeds, not measurements of plant state.

`feedback_source` is `estimate` by default. `truth` explicitly requests a
full-state reference controller; `zero` requests the zero-input baseline.
Neither baseline is an estimator result. The simulator rejects externally
supplied `input`, `measurements`, and `output_step` because it generates them.
The same named observation formula is used for sensor and estimator in this
benchmark gate; standalone KLQE still supports independent sensor/model terms.

Outputs include dense physical truth, applied input, all sensor timestamps,
clean/noisy measurements and noise draws, estimates with a readiness mask,
measurement estimates, old/new commands, information cutoffs, every fit and its
attempt diagnostics, and estimation compute time. Dense model forecasts are
computed before generating the next measurement. Plant readout channels retain
their unforced definitions; input-aware measurement arrays are separate.
Integration or estimation failure retains the completed physical prefix and an
explicit termination reason. A returned object does not by itself imply success:
check `completed`. The simulator assumes instantaneous computation, exact known
actuation, no process noise and no sensor delay. It is not KLQR/KLQG.

Matched closed-loop comparisons share clocks, gain, limits, prior information
and noise draws; subsequent measurements may differ because each controller
changes its own plant. For an identical-data comparison, replay one saved
physical/sensor/input stream through `identify_controlled_sequential` instead.
Tests compare online estimates to such a replay, use independent exact harmonic
transitions, mutate future measurements, and verify acquired-data initialization,
sensor timing, noise matching, and retained failure diagnostics.

## Observable diagnostics on saved states

Compiled `observable` operation `evaluate_batch` and Python
`dynsyslab_client.observables.evaluate_batch` evaluate named real observables on
flat sample-major `states`. They neither solve an ODE nor refit a trajectory.
Results include `samples`, `state_dimension`, `output_dimension`, and flat
sample-major `values`. The gate rejects derivative requests; a pre-defined
derived observable may be evaluated instead.

RPE now exposes `mechanical_energy` (alias `energy`), `energy_dissipation`, and
`control_power_coefficient`, with analytic derivatives. For gas coefficient
\(G\), ambient-minus-vapor \(A\), capillarity \(S\), viscosity coefficient
\(\nu\), ambient ratio \(P\), and pressure input sign \(s\),
\[
H=\tfrac12r^3v^2+\frac{G}{3(\kappa-1)}
 (r^{3(1-\kappa)}-1)+\frac A3(r^3-1)+\frac S2(r^2-1),
\quad
\dot H=-\nu r v^2+sP r^2v u.
\]
At \(\kappa=1\) the gas term is \(-G\log r\). The additive constant sets
\(H(1,0)=0\); these are dimensionless energy units. This first integral follows
from the integrating factor \(r^3\) used in the RPE analytic reductions (see
[Kudryashov and Sinelshchikov](https://arxiv.org/abs/1409.6699)). Native tests
check the gradient/RHS balance at several states, inputs, viscosities and both
polytropic and logarithmic cases. Numerical work integrals and positivity near
equilibrium remain subject to floating-point and quadrature errors.

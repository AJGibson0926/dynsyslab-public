# RPE echo-only KLQR driver-rediscovery experiment

Date: 2026-09-14. Native DynSysLab. Source and model archives were read without modification; both new control outputs are immutable archives.

## Result

The echo-only KLQR did **not** rediscover the known sine-wave driver. It tracked the scalar echo extremely closely while following a physically different trajectory and generating a strongly biased, phase-shifted input. A matched state-coordinate KLQR baseline recovered the original pressure waveform accurately, so the negative echo-only result is not a failure of the target archive, waveform convention, solver, or comparison method.

![RPE driver-rediscovery comparison](../../viewer_outputs/rpe-echo-only-driver-rediscovery-comparison.png)

Metrics use the common post-transient interval `10 <= tau <= 30`.

| Quantity | Known driver | Echo-only eigenfunction KLQR | State-KROM baseline |
|---|---:|---:|---:|
| Mean input, model units | approximately 0 | -1.02188 | 0.0008645 |
| Fundamental amplitude | 0.2 | 0.20211 | 0.19788 |
| Amplitude ratio | 1 | 1.01056 | 0.98938 |
| Phase error | 0 deg | -93.42 deg | 0.0597 deg |
| Direct driver NRMSE | 0% | 744.02% | 2.609% |
| Correlation with driver | 1 | -0.0402 | 0.99971 |
| One-tone fit R-squared | 1 | 0.4703 | 0.99942 |
| Echo NRMSE | -- | 0.06346% | -- |
| Radius RMSE | -- | 15.5171 | 0.0004368 |
| Maximum radius | 1.09904 target maximum | 21.3329 | 1.09805 |

The echo-only input's mean is about -102.19 kPa because one model input unit is 100 kPa. Its initial full-run peak is 7.4567 model units. The close fundamental amplitude is therefore not recovery: the phase, bias, harmonics, direct waveform error, and hidden-state behavior all fail.

## Target selection

The selected source is:

`archives/solutions/rpe/rpe-periodic-100khz-20kpa-phase-portrait__20260911T214426Z.npz`

It is a native, dimensionless RPE solve for air in water, starting at `[1,0]`, driven by one physical waveform:

```text
frequency: 100 kHz
amplitude: 20 kPa
phase: pi/2
u*(tau) = 0.2 cos(1.8927362514011978 tau)
period: 3.3196306683132035 tau
duration: 300 tau = 90.37 cycles
```

It is a reasonable scalar-echo target: over `200 <= tau <= 300`, one sinusoid explains 99.2127% of echo variance; adjacent final cycles differ by only `9.76e-6` of the AC echo RMS; the steady radius remains in `[0.95280,1.06406]`. The echo range is `[-0.33981,0.44294] Pa m`, with fundamental amplitude `0.22024 Pa m`.

## KROM and echo-only boundary

The retained native clean-isolation reference KROM is:

`archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz`

Its existing echo-only refit is:

`archives/derived/koopman-with-sindy/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948-289776Z__koopman-with-sindy__20260911T110325Z.npz`

The refit retains the parent pair `lambda = -0.2106164458 +/- 4.8186345616 i`, but each polynomial eigenfunction is a function of only the instantaneous `far_field_acoustic_echo`. It has one coordinate variable and one nonredundant complex eigenfunction plus its conjugate. The parent KROM's original eigenfunction was not used for the echo-only trial because it depends on radius and velocity. Calling the parent the reference KROM identifies the validated native pathway; it is not a universal optimality claim across every historical KROM.

The exact lifted input map still evaluates `B_psi(x)=J_phi(h(x)) J_h(x) B(x)` from simulator state. Thus the eigenfunction is echo-only, but this experiment is not a proof of an echo-sensor-only implementation.

## Control experiment

Both arms use the same RPE parameters, target source, initial state `[1,0]`, `0 <= tau <= 30`, output step `0.005`, BDF tolerances `1e-9`, maximum step `0.001`, exact lifted input map, continuous-RHS control, `Q=5000 I`, `R=1e-4`, no cap, and no feedforward.

The echo-only arm sets `target.space="observable"` and maps the archived state trajectory through the exact echo expression retained by the echo KROM. The archived scalar target agrees with an independent evaluation of that expression to `2.50e-16` maximum absolute error.

The state baseline uses the reference parent KROM and the same archived state trajectory directly. Its 2.61% waveform NRMSE, 0.99 amplitude ratio, 0.06-degree phase error, and bounded state show that feedback-only KLQR can nearly reconstruct this modest driver when the state-distinguishing coordinates are retained.

## Interpretation

Instantaneous echo is non-injective: distinct RPE states can have essentially the same echo and different input authority. The echo-only controller minimizes its lifted scalar error, finds a different state path with almost identical echo, and has no state penalty telling it that the expanding radius is wrong. The result is more decisive than echo tracking alone: echo NRMSE is only 0.063%, yet the driver correlation is negative and the radius reaches 21.33.

This experiment therefore rejects the hypothesis that the current instantaneous echo-only eigenfunction can identify the generating pressure from this target. A next scientifically justified escalation would add information—echo delay coordinates, echo plus measured echo derivative, or a state/domain constraint—not retune this scalar model until one waveform happens to look favorable.

## Reproduction and artifacts

- Echo-only job: `examples/native/case-studies/rpe-echo-only-klqr-rediscover-100khz-20kpa.json`
- Matched baseline job: `examples/native/case-studies/rpe-state-klqr-rediscover-100khz-20kpa-baseline.json`
- Echo-only result: `archives/solutions/rpe/rpe-echo-only-klqr-rediscover-100khz-20kpa__20260914T160803Z.npz`
- State baseline: `archives/solutions/rpe/rpe-state-klqr-rediscover-100khz-20kpa-baseline__20260914T161049Z.npz`
- Machine-readable metrics: `viewer_outputs/rpe-echo-only-driver-rediscovery-comparison.json`
- Comparison figure: `viewer_outputs/rpe-echo-only-driver-rediscovery-comparison.png`

Both jobs passed native validation and planning before execution and completed with 6,001 published samples. SHA-256 identities and full numerical values are retained in the machine-readable report.

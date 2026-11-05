# RPE echo/echo-dot KLQR nonlinear sine-driver stress test

Date: 2026-09-14

## Executive result

The repaired echo/echo-dot eigenfunction can control a substantially nonlinear,
pure-sine-driven RPE orbit and recover the dominant sine driver, but it requires
a much softer KLQR gain than the earlier 20 kPa example and the recovered command
is not a perfectly clean sine.

The selected target is an air bubble in water with equilibrium radius
`R0 = 10 um`, driven by a `100 kHz`, `40 kPa` pure sinusoidal pressure. Its
settled radius spans `0.887216` to `1.171938 R0`; the second radius harmonic is
`30.65%` of the fundamental and the RMS of harmonics 2--5 is `40.38%` of the
fundamental. Those distortion measures are about 2.4 and 2.8 times the 20 kPa
case, respectively.

With the same archived echo/echo-dot KROM, the same selected Koopman mode, exact
state-dependent `B_psi`, and `Q = 5000`, changing only the input penalty from
`R = 1e-4` to `R = 0.006` produces a branch-stable long-horizon result:

| Metric, `tau = 100..300` | Tuned echo/echo-dot KLQR |
|---|---:|
| Direct driver NRMSE | `18.8599%` |
| Driver correlation | `0.982380` |
| Recovered fundamental amplitude | `39.4774 kPa` vs `40 kPa` |
| Fundamental phase error | `+0.3130 deg` |
| Fundamental fit R-squared | `0.96510` |
| Direct pressure RMSE | `5.3344 kPa` |
| Radius RMSE | `0.002014 R0` = `20.14 nm` |
| Radial-velocity RMSE | `0.019758 R0/T` = approximately `65.6 mm/s` |
| Echo NRMSE | `3.8217%` |
| Autonomous echo-dot NRMSE | `11.5698%` |
| Complete late cycles | `59` |
| Cycle driver NRMSE, min / median / max / last | `15.20 / 19.15 / 22.95 / 17.32%` |

The dominant sinusoid is therefore rediscovered accurately in amplitude and
phase and remains locked for the complete late-time interval. The remaining
`18.86%` waveform error comes principally from narrow, phase-locked correction
spikes: the command reaches `1.00365` in model units (`100.365 kPa`) even though
the desired sine amplitude is `0.4` (`40 kPa`). This is a successful nonlinear
orbit and fundamental-driver recovery, but not yet a laboratory-ready clean
waveform reconstruction.

![Nonlinear sine-driver stress-test summary](../../viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.png)

The exact machine-readable metrics, cycle rows, hashes, and lineage are in
`viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.json` and can be regenerated
with `python tools/analyze_rpe_echo_jet_nonlinear_sine.py`.

## 1. Archive search and target selection

The search inspected all `588` archived solution files: `21` native `rpe`
solutions and `567` legacy `rayleigh-plesset` solutions. Metadata identified
`418` pure-sine candidates.

The legacy collection contains `387` high-amplitude `155 kHz` sine trajectories,
mostly from `90` to `99.77 kPa`. Some are strongly cavitating; one `90 kPa`
example spans approximately `0.143` to `5.504 R0`. They were not used directly
as control targets because their plant parameters do not match the KROM:

| Parameter | Legacy sine sweep | Repaired-KROM plant |
|---|---:|---:|
| Ambient pressure | `101325 Pa` | `100000 Pa` |
| Polytropic exponent | `1.4` | `4/3` |
| Water density | `998.2 kg/m^3` | `998 kg/m^3` |
| Surface tension | `0.0728 N/m` | `0.0725 N/m` |
| Dynamic viscosity | `0.001002 Pa s` | `0.001 Pa s` |
| Vapor pressure | `2339 Pa` | `2330 Pa` |

Using such a trajectory would conflate controller error with plant mismatch and
would make recovery of its archived driver an ill-posed comparison. The legacy
sweep was instead used as a regime guide. Exact-parameter-matched pure-sine
targets were generated with the native RPE plant and screened at 100 and 155 kHz.

The 155 kHz targets became strongly inertial even at 40 kPa and lay far outside
the training domain. At 100 kHz, the useful progression was:

| Amplitude | Settled radius range | H2/H1 | RMS(H2:H5)/H1 | Control observation |
|---:|---:|---:|---:|---|
| `20 kPa` | `0.95276..1.06416` | `12.53%` | `14.46%` | Original success |
| `28 kPa` | `0.93103..1.10061` | `19.15%` | `23.55%` | Default gain succeeds long-term |
| `40 kPa` | `0.88722..1.17194` | `30.65%` | `40.38%` | Strong nonlinear target; gain retuning required |
| `45 kPa` | approximately `0.858..1.215` | `35.63%` | `48.22%` | Default-gain CARE loses a real solution early |
| `50 kPa` | approximately `0.827..1.257` | `40.53%` | `56.33%` | Default-gain CARE loses a real solution early |
| `55 kPa` | approximately `0.785..1.317` | `45.65%` | `65.75%` | Default-gain CARE loses a real solution early |

The 40 kPa target remains inside the axis-aligned state and echo-jet ranges of
the forced training archive, while 45 kPa and above increasingly leave them.
An axis-aligned range check is only a necessary condition--it does not prove
coverage of the learned manifold--but it makes 40 kPa the strongest defensible
same-KROM test before clear extrapolation.

## 2. The matched nonlinear target

The physical forcing is

```text
p_c(t) = 40 kPa sin(2 pi (100 kHz) t + pi/2).
```

With the archived RPE nondimensionalization this becomes

```text
u_*(tau) = 0.4 cos(1.8927362514011978 tau),
period = 3.3196306683132035 tau.
```

The target was first integrated to `tau = 300`. A phase-zero state on the
settled attractor was extracted at `tau = 295.447129`:

```text
[r, v] = [0.8875145511136061, -0.03863871558543898].
```

A second target solve was started from this state so that the KLQR test began on
the settled orbit instead of spending its first cycles following target spin-up.
Adjacent final cycles of the source target agree to approximately `1.36e-5`
NRMSE in radius and `3.65e-5` in velocity.

The frozen target archive is:

```text
archives/solutions/rpe/
rpe-periodic-100khz-40kpa-settled-target__20260914T183110Z.npz
```

## 3. Reused echo/echo-dot eigenfunction

No KROM was retrained for this amplitude. The controller reuses exactly:

```text
archives/derived/koopman-with-sindy/
rpe-echo-study-state-basis-baseline__20260911T104742Z__
koopman-with-sindy__20260914T165124Z.npz
```

Its observable coordinates are

```text
y = [far_field_acoustic_echo, far_field_acoustic_echo_dot]
  = [e, L_f0 e].
```

The eigenfunction is a normalized degree-four polynomial in this two-coordinate
echo jet. The selected complex mode has

```text
lambda = -0.21061644578341265 +/- 4.818634561556464 i.
```

The controller configuration held fixed across the amplitude and gain trials was:

| Setting | Value |
|---|---|
| Selected modes | `[0]` (one complex-conjugate pair) |
| Control representation | `exact` state-dependent `B_psi` |
| State cost | `Q = 5000` |
| Integrator | BDF |
| Relative / absolute tolerance | `1e-9 / 1e-9` |
| Maximum internal step | `0.001 tau` |
| Control evaluation | continuous RHS |

Only `R`, the input penalty, was varied.

## 4. Why the old gain failed

The 20 kPa value `R = 1e-4` is too aggressive on the 40 kPa orbit. Its
`tau = 100..300` driver NRMSE is `218.98%`; correlation is `0.2292`, median cycle
error is `254.88%`, and the controller repeatedly selects a competing feedback
branch. A phase-aligned periodic target does not remove this failure, so it is not
a target-startup artifact.

The state-coordinate KROM with the same one-mode KLQR setup and the original
`Q/R` succeeds on the same 40 kPa target:

| 40 kPa long-term metric | Echo jet, old gain | State-coordinate baseline |
|---|---:|---:|
| Driver NRMSE | `218.98%` | `1.9869%` |
| Driver correlation | `0.2292` | `0.999880` |
| Fundamental amplitude ratio | `0.2685` | `0.98816` |
| Radius RMSE | `0.13329 R0` | `0.001722 R0` |
| Velocity RMSE | `0.20021 R0/T` | `0.004275 R0/T` |

Thus the 40 kPa plant, the target, and a one-mode KLQR controller are not the
fundamental obstruction. The failure is specific to the gain interacting with
the nonlinear echo-coordinate eigenfunction and its state-dependent input map.

## 5. Gain search

The short-screen sweep used the same settled target and initial phase. Except for
the old-gain long run, scores below are over `tau = 60..100`:

| R | Driver NRMSE | Correlation | Result |
|---:|---:|---:|---|
| `0.0001` | `218.98%` | `0.2292` | long-run wrong branch |
| `0.0003` | `127.68%` | `0.5472` | wrong branch |
| `0.001` | `104.70%` | `0.6437` | wrong / mixed branch |
| `0.002` | `75.34%` | `0.7785` | mixed branch |
| `0.003` | `18.76%` | `0.9825` | good screen, later branch hopping |
| `0.004` | `75.36%` | `0.7751` | mixed branch |
| `0.006` | `20.68%` | `0.9788` | selected; robust long run |
| `0.01` | `21.06%` | `0.9780` | robust long-run comparator |
| `0.03` | `29.01%` | `0.9592` | stable but less accurate |
| `0.1` | `34.70%` | `0.9426` | stable but less accurate |

This response is not a smooth convex tuning curve. The nonmonotone jumps are
consistent with a folded echo-coordinate feedback geometry: small gain changes
can move the closed loop between different physical branches that look nearby in
the learned eigenfunction coordinate.

The apparently best short-screen value, `R = 0.003`, demonstrates why long tests
are essential. Over `tau = 100..300`, it intermittently jumps from cycles near
`15--20%` error to cycles near `255%`, producing `115.05%` aggregate NRMSE.

Both softer candidates remain branch-stable for all 59 late cycles:

| R | Long driver NRMSE | Correlation | Cycle min / median / max | Radius RMSE |
|---:|---:|---:|---:|---:|
| `0.006` | `18.8599%` | `0.982380` | `15.20 / 19.15 / 22.95%` | `0.002014 R0` |
| `0.01` | `21.4317%` | `0.977337` | `17.28 / 21.47 / 27.94%` | `0.002172 R0` |

## 6. What was recovered

For the selected `R = 0.006` run, a direct least-squares fit of the command to
the known 100 kHz frequency over `tau = 100..300` gives

```text
mean                  = -0.00339123       (-0.3391 kPa)
fundamental amplitude =  0.39477366       (39.4774 kPa)
phase error           = +0.312981 deg
sinusoidal R-squared  =  0.965099
```

The main sine is therefore recovered within `1.31%` in amplitude and `0.313 deg`
in phase. The direct NRMSE remains `18.86%` because of short correction features
near repeatable phases. These features do not cause secular radius growth or
late loss of lock, but they are undesirable control effort and should not be
hidden by quoting only the fundamental fit.

The selected archive is:

```text
archives/solutions/rpe/
rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300__
20260914T192234Z.npz
```

It required `1,882,131` accepted and `231,947` rejected BDF internal steps. The
large count is consistent with steep state-dependent feedback regions. The
robust `R = 0.01` comparator required `1,849,207` accepted and `226,623` rejected
steps.

## 7. Intermediate 28 kPa success

The same unchanged echo-jet KROM and the original `R = 1e-4` also succeed at
28 kPa for the full `tau = 300` horizon:

| Metric | 28 kPa result |
|---|---:|
| Driver NRMSE | `4.1849%` |
| Driver correlation | `0.999137` |
| Fundamental amplitude ratio | `1.00331` |
| Phase error | `-0.05254 deg` |
| Radius RMSE | `2.96e-5 R0` |
| Velocity RMSE | `8.50e-4 R0/T` |
| Cycle NRMSE min / median / max / last | `4.00 / 4.19 / 4.38 / 4.28%` |

This is a cleaner recovery than the 40 kPa case, but its target has only `19.15%`
second-harmonic content. It is best described as a moderate nonlinear extension,
not the final stress test.

## 8. Interpretation and limitations

The result changes the earlier diagnosis in an important way. At 40 kPa, the
echo/echo-dot eigenfunction is not categorically incapable of supplying a usable
control coordinate. A stable low-gain basin exists and recovers the correct
fundamental driver. However:

1. The basin is narrow and nonmonotone in `R`; short tests can select a controller
   that later jumps branches.
2. The recovered pressure contains phase-locked spikes and is materially worse
   than the state-coordinate baseline.
3. The `exact` `B_psi` implementation still evaluates physical RPE state. This
   experiment validates the eigenfunction coordinate and closed-loop geometry;
   it is not yet a deployable echo-sensor-only controller.
4. No measurement noise, echo-derivative estimation error, parameter mismatch,
   delay, bandwidth constraint, or actuator saturation is present.
5. The target remains within the training archive's axis-aligned bounds, but it
   is near a harder portion of the learned manifold. Stronger 45--55 kPa targets
   begin leaving that envelope and fail under the present controller construction.

The next technically meaningful improvement is not more manual scalar gain
tuning. It is a branch-aware or saturation-aware controller--for example KMPC,
gain scheduling using an observability/conditioning indicator, or an additional
physically accessible observable--combined with a genuinely sensor-side estimate
of `B_psi`. That should aim to remove the narrow correction spikes while retaining
the accurately recovered 40 kPa fundamental.

## 9. Reproduction map

Key jobs:

```text
examples/native/case-studies/rpe-periodic-100khz-40kpa-settled-target.json
examples/native/case-studies/rpe-echo-jet2-klqr-rediscover-100khz-40kpa-settled-tau300.json
examples/native/case-studies/rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300.json
examples/native/case-studies/rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r1e-2-tau300.json
examples/native/case-studies/rpe-state-klqr-rediscover-100khz-40kpa-settled-tau300.json
examples/native/case-studies/rpe-periodic-100khz-28kpa-nonlinear-tau300.json
examples/native/case-studies/rpe-echo-jet2-klqr-rediscover-100khz-28kpa-tau300.json
```

Analysis:

```text
python tools/analyze_rpe_echo_jet_nonlinear_sine.py
```

Outputs:

```text
viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.json
viewer_outputs/rpe-echo-jet2-nonlinear-sine-stress.png
```

The prior report
`docs/case-studies/RPE_ECHO_JET_SINE_DRIVER_REDISCOVERY_LONG_TERM_2026-09-14.md`
contains the full pipeline from training data through state-coordinate spectral
identification, forced-training generation, echo-jet eigenfunction refit, exact
`B_psi`, and the original 20 kPa KLQR result. This report changes only the target
amplitude study and gain selection; it deliberately reuses that same trained model.

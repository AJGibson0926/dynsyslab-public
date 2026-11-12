# Reproduction notes: RPE minimal-KROM KLQR stability case

Use these notes as the laboratory checklist. The complete analysis is in [the case-study report](../case-studies/RPE_MINIMAL_KROM_KLQR_STABILITY_CASE_STUDY_2026-08-26.md).

## Preserve first

- Do not overwrite the baseline solution or any upstream KROM/H-DMD/SINDy archive.
- The working tree was dirty when the case was documented; artifact and source-file hashes are more authoritative than Git commit `7037594f4db6874d3469a032b0cf245e446dac81` alone.
- Give every follow-up solve a distinct label.

## Baseline artifact

```text
archives/solutions/rayleigh-plesset/
rpe-equilibrium-minimal-krom-klqr-multisine-tau25__20260826T170659.584050Z.npz
SHA-256: 91a1b9255dcde946846712cb530ffd5ef80bc88877d410c5ed32e302aef43fdb
```

Viewer output:

```text
viewer_outputs/rpe-minimal-krom-klqr-multisine-controlled.html
```

## Essential configuration

```text
plant: Rayleigh--Plesset, air_water_20c
R0 = 1e-5 m
kappa = 1.3333
p0 = 100000 Pa
pV = 2330 Pa
rhoL = 998 kg/m^3
sigma = 0.0725 N/m
muL = 0.001 Pa s
pG0 = 112170 Pa

f0 = 331958.7786861145 Hz
T0 = 3.0124222168727627e-6 s
P0 = 10997.623748528325 Pa

initial state: (r,u) = (1,0)
target: rpe_multisine_equilibrium_small
r_d = 1 + 0.06 sin(4 tau) - 0.08 sin(3 tau)
u_d = 0.24 cos(4 tau) - 0.24 cos(3 tau)

integration: RK4, delta_tau = 0.001, tau = 0..25
nominal controller hold: delta_tau_c = 0.005
Q = I_2
R = [[1e-16]]
pressure limit: +/-200000 Pa
estimator: identity/full state
```

## Minimal KROM identity

```text
archive:
archives/derived/koopman-krom-extract/
rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__koopman-krom-extract__20260826T164339.807992Z__ebb5aedc.npz

SHA-256:
16dbd0675707021aa94c6452ff0d978a04a4f6c423e0a204993bbee628363712

reservation label:
rpe-natural-fundamental-single-eigenfunction

source full-spectrum indices:
433 and conjugate 434

mu:
0.9998540377155265 +/- 0.005598164269206836 i

lambda per tau:
-0.13029888671192555 +/- 5.598923003706494 i

representation:
one nonredundant complex eigenfunction = two real split coordinates
analytic state-dependent B_psi
no state reconstruction
```

Upstream lineage:

```text
free trajectory (r(0),u(0)) = (2.3334,0), tau=0..25, delta_tau=0.001
 -> H-DMD: 20000 x 2000 Hankel pair, full rank 2000, gesvd, tol 1e-10
 -> lambda = log(mu)/delta_tau
 -> SINDy: degree-4 polynomial, 15 features, threshold .05, max 20 iterations
 -> extract indices 433/434 with conjugate preservation
```

## Stability signature

Representative point:

```text
tau = 15.93
state (r,u) = (1.1043053887412562, 0.33476114547142466)
target       = (1.0959888348400604, 0.3400997180367939)
p_c          = 38154.56984837685 Pa
p_c*         = 3.469346717147202

continuous poles per tau:
-5.604306088257247
-342.5694922455953

exact held-input sampled poles at delta_tau_c=.005:
0.9723761687799314
-0.7362623165043962
```

Interpretation: the fast sampled pole is negative and reverses sign at every controller update. It is stable in magnitude at this point but creates alternating high-gain chatter. The largest observed pressure jump is `29825.4536 Pa = 2.71199 P0`. No pressure sample saturates.

Clock caveat: the original floating-point sample clock is not perfectly periodic. The call log contains 5,058 pressure changes: 4,969 intervals of `.005 tau`, 29 of `.0045`, 29 of `.0005`, 29 duplicate changes at identical solver times, and one roundoff variant. Around the principal event, updates occur near `.xxx5 tau`, one RK4 half-step after the nominal boundaries. Preserve this behavior only for baseline reproduction; fix it with an integer-tick clock in a separately labeled experiment.

## Baseline numerical metrics

```text
radius RMSE                 0.006360083992258022
velocity RMSE               0.007339480878497758
max |radius error|          0.011419017629996864
max |velocity error|        0.022717496907171505
pressure RMS                23775.90402176294 Pa
peak |pressure|             68040.66214461008 Pa
dimensionless pressure RMS  2.161912842757925
peak |dimensionless p|      6.1868512417252965
saturated samples           0
```

## Tuning warnings

- `Q` weights `[Re psi, Im psi]`, not radius and velocity directly.
- `R` weights dimensional Pa squared. Baseline equivalent dimensionless pressure weight is `R* = 1.209477281141942e-8`.
- Increasing `Q` or decreasing `R` increases aggressiveness and worsens the sampled mode.
- At the representative point, `Q=10I`, `R=1e-16` gives sampled fast pole `-4.49093`: genuinely unstable.
- `Q=I`, `R=1e-15` gives sampled fast pole `+0.45238`.
- Keeping baseline `Q,R` and reducing the hold to `delta_tau_c=0.001` removes sign reversal; the fast sampled pole is positive over the inspected `tau=14..17` window.

## First experiment matrix

Run and archive separately:

| Experiment | `Q` | `R` | `delta_tau_c` | Purpose |
|---|---:|---:|---:|---|
| baseline reproduction | `I` | `1e-16` | `.005` | verify environment |
| sampling only A | `I` | `1e-16` | `.002` | isolate hold effect |
| sampling only B | `I` | `1e-16` | `.001` | align control and solver grids |
| cost only A | `I` | `1e-15` | `.005` | remove negative fast pole |
| cost only B | `I` | `1e-14` | `.005` | smoother, slower bound |
| conservative combined | `I` | `1e-15` | `.001` | proposed safe starting point |

Then implement and compare an exact-ZOH discrete Riccati controller at `.005 tau` and a fully dimensionless cost/input formulation.

Also run a clock-isolation pair before tuning: original floating clock versus a monotone integer-tick clock with the same `Q`, `R`, nominal `.005 tau` interval, plant, KROM, and target.

## Completed tuning results

Two complete nonlinear follow-ups are preserved:

```text
R=1e-15, delta_tau_c=.001
archive SHA-256: 3d31831d6835bc84af24975f447fbe1b22f04a5bf2124c74fcf1eada8c8b6339

R=1e-16, delta_tau_c=.001
archive SHA-256: b94345fd39cba65f99c634937cc7093955b3724475982e0fd49bb64fb5af79ad
```

The `R=1e-15` run reduced the largest jump by about 98% but worsened radius RMSE by 81% and velocity RMSE by 187%. Do not adopt it as the accuracy-preserving default.

The selected result is:

```text
Q=I, R=1e-16, delta_tau_c=.001
radius RMSE:   0.006356763534513237  (0.052% better than baseline)
velocity RMSE: 0.007334030491399215  (0.074% better)
max jump:      2919.3398289241813 Pa (90.21% lower)
TV tau 14..17: 264171.594099237 Pa   (74.76% lower)
sampled poles near tau=15.93: 0.9944118352, 0.6573345727
```

This isolates the useful change as controller sampling, not a larger `R`. Peak pressure rose 1.11%, pressure RMS was essentially unchanged, and no sample saturated.

The `.001` floating-clock run still contains 148 duplicate same-time command changes. Correct that clock in a separately labeled experiment before declaring a production default.

## Alternate eigenfunction comparison

Full-spectrum index `431` plus conjugate `432` was extracted from the same original SINDy KROM and tested with all other selected-case settings unchanged.

```text
alternate lambda = -0.21123030699891182 + 4.822747859492991 i per tau
alternate KROM SHA-256 = f0cbfb21def5c4492b5c6a20b7f3a22ff224b41b3b80c57ed05645688d3d755c
alternate solution SHA-256 = 925619324dbab3dab8143e28bcf93bf1de049c891b4d6363158e8f43e1a9125a
```

Against the selected index-433 controller, index 431 produced:

```text
radius RMSE:       +3.84% worse
velocity RMSE:     +7.83% worse
pressure RMS:      -0.31%
peak pressure:     -1.06%
maximum jump:      -14.15%
TV over tau 14..17 -2.85%
sampled poles near tau=15.93: 0.9951816331, 0.6712951946
```

Conclusion: index 431 is smoother but less accurate. Retain index 433 as the preferred single-eigenfunction KROM. Frequency proximity to the target was not a sufficient selection criterion.

## Two-eigenfunction comparison

Full-spectrum indices `431` and `433`, together with conjugates `432` and `434`, were extracted from the original SINDy KROM without refitting. This gives two complex eigenfunctions represented by four real split coordinates. All selected-case settings were preserved except that the state cost became `Q=I4`.

```text
two-mode KROM SHA-256     = 29a2378d6716b948b00306b4a6994a20596896b51920a6498c7fbb4897320f92
two-mode solution SHA-256 = 1330ee074e51d73d7496c85f67e06d00bf689d8ae79484225aef361f2aec86a1
solve time                   110.777 s
```

Against the selected single index-433 controller, the two-function KROM produced:

```text
radius RMSE:          0.0051329694  (19.25% better)
velocity RMSE:        0.0055095604  (24.88% better)
max radius error:     0.0091523481  (19.78% better)
max velocity error:   0.0169266209  (25.37% better)
pressure RMS:         24062.3907 Pa (1.26% higher)
peak pressure:        68282.5938 Pa (0.74% lower)
maximum jump:         221.0613 Pa   (92.43% lower)
TV over tau 14..17:   253285.8414 Pa (4.12% lower)
sampled poles near tau=15.93:
  0.5237163908, 0.9948262850, 0.9995591297 +/- 0.0051929735 i
```

The dimensionless/model-scaled controllability matrix at `tau ~= 15.93` has full rank four but condition number `1019.55`. No pressure sample saturated. The known clock defect remains: 25,296 command changes and 148 same-time changes. Conclusion: this is the best-performing KROM tested so far, but its poorer conditioning and the clock defect require follow-up before promotion to a default.

Viewer artifact: `viewer_outputs/rpe-two-eigenfunction-krom-431-433-klqr-controlled.html`.

## Off-equilibrium large-target diagnosis

Target:

```text
r_ref = 1.5 + 0.2 sin(4 tau - pi/2) + 0.3 sin(3 tau - pi/2)
u_ref = d r_ref/dtau
range: r=1.0000..1.95066, u=-1.65768..1.65768
initial target = initial plant = (1,0)
```

The feedback-only implementation returns zero pressure at exact target matching. Nonlinear inversion shows that this target instead requires `p_ff` with RMS `90.578 kPa` and peak `140.912 kPa`, within the `200 kPa` limit.

```text
case                         r RMSE       u RMSE       sat.  max jump
feedback, R=1e-16            1.0968e-2    1.2399e-2    3419  400000 Pa
feedforward, R=1e-16         4.7333e-4    1.4268e-3    3365  400000 Pa
feedforward, R=1e-15         1.7940e-6    9.5070e-6       0     467.94 Pa
```

The final pressure RMS is `90.577 kPa`, essentially the inverse-dynamics requirement. Its SHA-256 is `d7d9249661ed69b3461c9882298b62dd1d4ed561ac50702e47b728344ae5a889`.

Conclusion: the eigenfunction pair is adequate for this in-domain target. Missing reference feedforward causes the offset; excessive feedback gain causes chatter. Production fix: add an explicit reference-dynamics/feedforward contract, increase/normalize `R` after feedforward is present, use exact discrete-time design and an integer-tick clock, and penalize command increments. The current feedforward implementation was experiment-local; production code was not changed.

Viewer artifact: `viewer_outputs/rpe-two-krom-offset1p5-large-feedforward-R1e-15.html`.

## Required measurements for each rerun

- archive SHA-256 and complete metadata;
- radius/velocity RMSE and maximum errors;
- pressure mean, RMS, peak, total variation, and saturation count;
- largest command jumps at update boundaries;
- continuous frozen poles and exact sampled poles along the trajectory;
- `||K||`, `||B_psi||`, and controllability conditioning;
- controller update count, Riccati solve count, and wall time;
- Viewer plot with target, solution, and a separate dimensionless pressure panel.

## MATLAB `soft` parity correction

The off-equilibrium feedforward conclusion above applies to the tested Python
controller only. A subsequent source audit found no automatically applied
pressure feedforward in MATLAB KLQR. `CTRL.u` computes only `u=-K*K_err`.
MATLAB `compute("soft")` calls `ctrl.u(x,t)` from every internal `ode15s` RHS
evaluation, and that call recomputes current `B(x)`, lifted `B_psi(x)`, and the
continuous LQR/CARE gain. It has no independent controller clock or held
command. The ordinary Python controller recomputes `B_psi` and CARE only on its
sample clock while reevaluating physical `B(x)` at RK stages with the held
command. This policy mismatch is now a primary parity hypothesis.

Additional parity facts:

- MATLAB Bubble uses `B(r)=-[0;Eu/r]`; Python currently uses the opposite
  pressure sign. Consistently flipping controller and plant input signs leaves
  ideal feedback state dynamics unchanged, but it changes physical command
  semantics and feedforward/limit interpretation.
- The MATLAB examples reference workspace `kROM{1}` after observable-based
  extraction; two-mode extraction lines are commented out. The exact KROM
  membership and normalization are not yet reproduced by the Python
  two-function KROM.
- MATLAB Bubble tuning loaded from `KLQR_ctrl_tuning.mat` is `Q=5000 I`,
  `R=1e-4` in its dimensionless-input convention.

Feedback-only diagnostic archive:

```text
rpe-continuous-stage-fixed-r2-feedback-only-Q1-R1e-8-tau10__20260826T194822.218528Z.npz
SHA-256: 1bb49733d2679b3c790015468cb01b95714c1c598420606a5a8982fe202533b8
```

This run used the MATLAB scanner's fixed target `[2,0]`, continuous four-stage
RK4 controller recomputation, `Q=I`, `R=1e-8`, `delta_tau=0.005`, and no
feedforward. It completed in 175.24 seconds with 8,000 controller evaluations,
but the selected two-function KROM left radius near equilibrium (radius RMSE
`1.000039`, final radius `0.9999817`, pressure RMS `16.973 Pa`). Harder MATLAB
waveform trials left the valid radius/domain or were prohibitively slow; they
produced no accepted archive.

## Deferred startup-spike study

For strict parity, retain the large initial pressure spikes caused by starting
the plant at equilibrium against a distant/nonzero-velocity target. Do not
silently clip or ramp the parity baseline. After parity is established, compare
smooth target engagement, bumpless transfer, magnitude/slew constraints,
control-increment penalties, and preview control. Preserve post-transient
tracking metrics and archive every mitigated case separately so spike reduction
cannot hide degraded control performance.

Current parity policy is explicit: `control_limit=None`, continuous RHS
state-dependence, no command ramp, no slew limiter, and no startup target blend.
The initial spike is therefore an accepted diagnostic feature, not a production
recommendation. A later mitigation is acceptable only if it reduces peak
pressure and pressure slew while preserving the established post-engagement
tracking error, phase, and stability margins.

## Pressure-input convention migration

The default RPE now uses `external-acoustic-pressure/v2`; positive pressure
compresses the bubble. Archives without the field reconstruct as
`legacy-positive-acceleration/v1`, and no existing archive was modified.
Analytic KROM reuse across this correction requires an explicit copied
control-system rebinding with recorded reason and unchanged spectral lineage.

The corrected-sign `0.05 tau` continuous/sampled diagnostics exactly retained
the earlier state trajectories and pressure magnitudes while reversing command
sign. This confirms the expected LQR invariance when both `B_psi` and physical
`B` are flipped consistently. It also confirms that the two-function KROM's
initial compression is not caused by the old sign convention.

## KROM membership and training-domain audit

Grouped polynomial observability selection across 1,001 nonredundant modes
selected full indices `[433,266,767,1983,272]`. Mode 433 was first, confirming
the prior single-mode choice. The five-mode span reached only radius/velocity
`R^2=[0.3454437181,0.8154874705]` before independent improvement was exhausted.

The training radius is only `[0.23287463,2.3334]`; the hard target begins at
`r=4`. The target is therefore explicit polynomial extrapolation. Controllers
now default to an error outside the archived training domain. Reproduction of
the failing parity case requires `target_domain_policy="ignore"`, which is
stored in controller provenance.

The five-mode KROM reduced initial/peak pressure from `66.158/187.761 kPa` to
`15.049/37.900 kPa` over `0.05 tau`, but still compressed the bubble. Do not
tune this extrapolating case as though it validated KROM tracking. Expand the
training domain and build a validated state-observation map first.

The first-class held-out observation-map gate now confirms this failure:
radius/velocity validation `R^2=[0.3442443,0.8149293]`, normalized RMSE
`[0.1698349,0.0275825]`, condition number `1.3244e16`; no map was promoted.
The all-target dry-run training plan expands only the uncovered upper-radius
axis to `7.493379` (10% margin), retains current velocity coverage, proposes 72
new ICs, and requires screening before any full campaign execution.

## Offset-1.5 clean-pathway follow-up

The critical unchanged-target failure of the current and adjacent minimal
KROMs is retained in
`../case-studies/RPE_OFFSET1P5_MINIMAL_KROM_FAILURE_CASE_STUDY_2026-08-27.md`,
with exact numerical diagnostics in
`RPE_OFFSET1P5_EIGENFUNCTION_FAILURE_DIAGNOSTICS_2026-08-27.json`. Use that
record for future eigenfunction-selection and forced-reference experiments;
do not reconstruct the case from conversational summaries.

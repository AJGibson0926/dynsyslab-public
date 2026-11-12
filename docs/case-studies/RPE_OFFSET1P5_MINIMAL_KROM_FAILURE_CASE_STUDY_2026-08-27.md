# RPE offset-1.5 minimal-KROM KLQR failure case study

**Date:** 2026-08-27  
**Status:** retained critical negative result  
**Scope:** unchanged-target, feedback-only comparison of two minimal KROMs  
**Machine record:**
`../notes/RPE_OFFSET1P5_EIGENFUNCTION_FAILURE_DIAGNOSTICS_2026-08-27.json`

## 1. Executive conclusion

This case is a genuine representation-dependent control failure, not simply a
"harder target." The same clean-isolation controller pathway tracks seven
other fixed or multisine targets well, including much larger excursions, but
both tested low-frequency single-eigenfunction KROMs fail badly on this
offset-1.5 waveform.

The evidence supports three interacting causes:

1. The plant and reference begin at the same state, `(r,u)=(1,0)`, so the pure
   feedback law commands exactly zero at `tau=0`. However, the reference has
   acceleration `d2r_ref/dtau2=5.9` and immediately requires
   `p_control/p0=-0.648876564829`. State coincidence therefore hides a
   nonzero acceleration/forcing mismatch from the control law.
2. The selected quartic eigenfunction map becomes nearly singular along this
   target orbit. It maps a large final physical-state error (`0.740083`) to a
   very small lifted error (`0.00674165`). KLQR can consequently regard the
   trajectory as close in its selected coordinates while the physical radius
   remains badly wrong.
3. The controller produces enough total pressure amplitude, but its waveform
   has the wrong phase/shape. Actual versus inverse-dynamics pressure
   correlation is only `0.419448`, with RMSE `0.937783 p0`.

The adjacent complex pair does not repair the defect. It has similar target-
path folding and slightly worse tracking. The target remains entirely inside
the archived axis-aligned training box, `B_psi` remains nonzero, the lifted
pair remains controllable, CARE succeeds, and the command RMS is not too
small. Those explanations are therefore rejected for this experiment.

This is an evidence-backed diagnosis, not a proof that only these mechanisms
can matter. The historical specialized-KROM result demonstrates that this
target can be controlled accurately, but it used a different representation,
tuning/scaling, sampled policy, and pressure limiting and is not an isolated
one-variable comparison.

## 2. Reproduction contract

The target and all controller/solver settings must remain unchanged when this
case is reproduced. In particular, do not add feedforward, reshape or ramp the
target, cap the command, or substitute a sampled controller.

| Item | Retained value |
|---|---|
| Plant coordinates | dimensionless RPE natural-period coordinates |
| Plant initial state | `(r,u)=(1,0)`; controlled-system equilibrium default |
| Target | `rpe_offset_1p5_prior_feedback_only` |
| Target horizon | `tau=0..10` |
| Published spacing | `dtau=0.001` (10,001 samples) |
| Integrator | SciPy BDF |
| Tolerances | `rtol=atol=1e-9` |
| Control evaluation | continuous RHS; controller called by the adaptive solver |
| State cost | `Q=diag(5000,5000)` in two real lifted coordinates |
| Input cost | `R=[[1e-4]]` |
| Law | `p_control/p0 = -K(x)[Psi(x)-Psi(x_ref)]` |
| `B_psi(x)` | recomputed at every controller evaluation |
| CARE and `K(x)` | recomputed at every controller evaluation |
| Feedforward | none |
| Command cap | none |
| Startup ramp | none |
| Target-domain policy | explicit `ignore`; target is nevertheless in the training box |
| Pressure convention | `external-acoustic-pressure/v2`; positive external pressure compresses |

The archived solution publishes dimensional time and state, but integration,
target evaluation, KROM evaluation, controller state, and controller time all
use the same dimensionless coordinates. The control history is tagged
`coordinate_system=dimensionless`, `units=p_control/p0`, and
`scale_to_dimensional=100000 Pa`.

## 3. Plant and nondimensionalization

The run uses the exact `air_water_20c` RPE preset:

| Quantity | Value |
|---|---:|
| equilibrium radius, `R0` | `1e-5 m` |
| polytropic exponent, `kappa` | `4/3` exactly |
| ambient pressure, `p0` | `100000 Pa` |
| vapor pressure, `pV` | `2330 Pa` |
| liquid density, `rhoL` | `998 kg/m^3` |
| surface tension, `sigma` | `0.0725 N/m` |
| dynamic viscosity, `muL` | `0.001 Pa s` |
| equilibrium gas pressure, `pG0` | `112170 Pa` |

Here `pG0=p0-pV+2 sigma/R0`. There are no `CL` or `CG` parameters.

The state and time scaling are

```text
r   = R/R0
tau = t/T0 = f0 t
u   = dr/dtau = (T0/R0) dR/dt

f0 = sqrt(3 kappa pG0 - 2 sigma/R0)
     / (2 pi R0 sqrt(rhoL))
T0 = 1/f0
```

For the exact `kappa=4/3` run, `f0=331963.06683132035 Hz`,
`T0=3.0123833037972497e-6 s`, and the RPE pressure scale is
`P0=10997.907878457558 Pa`. The archived dimensionless RPE groups are:

| Group | Value |
|---|---:|
| `pG0/P0` | `10.199212544752802` |
| `pV/P0` | `0.21185847578919523` |
| `p0/P0` | `9.092638445888207` |
| `2 sigma/(R0 P0)` | `1.31843257465379` |
| `4 muL/(T0 P0)` | `0.12073680576341679` |

The control coordinate is independently and explicitly defined as
`u*=p_control/p0`. Consequently its dimensionless input column is
`B(r)=[0,-(p0/P0)/r]^T` for the external-pressure convention. This distinction
between the RPE equation's pressure scale `P0` and the controller's declared
pressure ratio `p_control/p0` is intentional and is recorded in the archive.

## 4. Exact control target

The retained target is

```text
r_ref(tau) = 1.5 - 0.2 cos(4 tau) - 0.3 cos(3 tau)
u_ref(tau) =       0.8 sin(4 tau) + 0.9 sin(3 tau)
```

It is represented by
`RPEHarmonicTarget(offset=1.5, amplitudes=(0.2,0.3),
angular_frequencies=(4,3), phases=(-pi/2,-pi/2))`.

At `tau=0`, the plant state and target state are both `(1,0)`, but

```text
d2r_ref/dtau2 = 0.2(4^2) + 0.3(3^2) = 5.9.
```

Direct inversion of the dimensionless RPE input channel gives the pressure
required to realize that acceleration:

```text
u_required*(0) = p_control(0)/p0 = -0.6488765648289958.
```

The retained pure-feedback law instead gives exactly zero because
`Psi(x)-Psi(x_ref)=0`. This is an important distinction from the successful
small equilibrium-centered multisine: that comparison target also begins at
`(1,0)`, but its initial acceleration and required pressure are both
approximately zero.

Over `tau=0..10`, the failed target occupies approximately
`r=[1.0,1.95066]`, `u=[-1.65768,1.65768]`.

## 5. Training and full H-DMD/SINDy lineage

Both tested minimal KROMs come from the same immutable pipeline:

1. One free, dimensionless RPE trajectory, exact IC
   `(2.3333865986017348,0)`, integrated by BDF over `tau=0..25` at published
   `dtau=0.001`, giving 25,001 samples.
2. Radius observable, a `20000 x 2000` window-snapshot H-DMD pair, LAPACK
   `gesvd`, absolute tolerance `1e-10`, full numerical rank, projected modes,
   and descending terminal-amplitude energy sorting. No rank cap was used.
3. Continuous generator conversion `lambda=log(mu)/dtau` before SINDy.
4. Fourth-order polynomial dictionary, no constant, no state normalization,
   coefficient threshold `1e-6`, one sparse-threshold pass, system-RHS
   generator data, conjugate-fit reuse, and coefficient-first storage.

The training-domain envelope is

```text
r in [0.2327207269088546, 2.333386598601735]
u in [-27.09322012592197, 26.94802281811782].
```

The failed target is therefore inside this axis-aligned envelope. This does
not prove dense training coverage of every point on the target orbit, but it
does rule out the gross out-of-box extrapolation that affected earlier large-
offset tests.

| Artifact | SHA-256 |
|---|---|
| `archives/solutions/rayleigh-plesset/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz` | `e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba` |
| `archives/derived/koopman-hankel-dmd/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz` | `c2fabe83b37b94ace1abb5b7d4fdd11cad14dee97c476d2fbb00e5e63eb0bf4b` |
| `archives/derived/koopman-with-sindy/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-with-sindy__20260827T183645.394109Z__c4dde209.npz` | `6e3725beac327fa749eeb0b0250cefa816d158cbbbc2d6fdea0d1b6bfa979119` |

## 6. KROMs compared

Each KROM contains one independent complex quartic polynomial eigenfunction
and its conjugate dual, represented to KLQR as two real coordinates. Neither
contains a state observation/reconstruction map. The coefficient vectors are
stored at full precision in the immutable archives and must be loaded from
there rather than transcribed from rounded report text.

### 6.1 Current clean-isolation minimal KROM

```text
full-spectrum indices: 1/2
mu:     0.9997777985801627 +/- 0.004817601141183756 i
lambda: -0.21061644578341268 +/- 4.818634561556464 i per tau
energy-order retained position: 999 of 1001
```

Archive:
`archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz`

SHA-256:
`8c3f221f563a82271e2995c8cc7e627b37f28b54f9feb93c59ae28050ae4ef86`

At equilibrium, its split-real matrices and gain are

```text
A = [[-0.210616445783, -4.818634561556],
     [ 4.818634561556, -0.210616445783]]

B_psi = [[-0.000858731531],
         [ 0.020922343973]]

K = [[-6825.894006335, 6996.487557274]]
```

The equilibrium controllability rank is 2, controllability condition number
is `4.82825`, gain norm is `9774.64`, and closed-loop poles are approximately
`-4.83348` and `-147.83228`. CARE therefore succeeds and gives stable local
lifted closed-loop poles.

### 6.2 Alternative adjacent minimal KROM

```text
full-spectrum indices: 3/4
lambda: -0.1257962639907903 +/- 5.5962746503729495 i per tau
energy-order retained position: 1000 of 1001
```

Archive:
`archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T193558.860380Z__fbc2035f.npz`

SHA-256:
`6ed276dd97e843d977d3272bdf1ee7aa6e8147e1d0c163f42313a8c28b3699d6`

At equilibrium:

```text
A = [[-0.125796263991, -5.596274650373],
     [ 5.596274650373, -0.125796263991]]

B_psi = [[-0.000704634035],
         [ 0.022141943760]]

K = [[-6889.348608242, 7083.920553783]]
```

Its controllability rank is 2, condition number is `5.59920`, gain norm is
`9881.55`, and closed-loop poles are `-5.61203` and `-156.34580`.

### 6.3 Selection-lineage correction

The current pair was previously described as the first energy-ranked complex
pair. That description was wrong. It is the first positive-frequency pair by
original full-spectrum index, but position 999 in retained energy order. The
alternative is position 1000. The numerical archives were not wrong or
changed; only the selection description was corrected. This matters because
neither experiment tests an energy-leading eigenfunction.

## 7. Controlled results

| Case | r RMSE | u RMSE | r RMSE, second half | u RMSE, second half | peak abs command | command RMS |
|---|---:|---:|---:|---:|---:|---:|
| Current pair 1/2 | `0.333425` | `0.603216` | `0.303291` | `0.557932` | `6.20816` | `1.05114` |
| Alternative pair 3/4 | `0.338447` | `0.622949` | `0.311481` | `0.580473` | `7.02241` | `1.09174` |
| Small equilibrium multisine, current pair | `0.000724051` | `0.00237950` | `0.000672826` | `0.00217474` | `0.707291` | `0.257554` |

The alternative pair is slightly worse in every tracking metric. It is not a
repair, although it is useful evidence that the failure is not unique to one
isolated numerical coefficient vector.

Current-pair solution archive:
`archives/solutions/rayleigh-plesset/rpe-prior-control-clean-rerun-offset-1p5-feedback-only__20260827T191637.359624Z.npz`

SHA-256:
`8666cf0f172b974d3c1b8f51d6b015c68590e48d97f60b40a109f1678c6552bf`

Alternative-pair solution archive:
`archives/solutions/rayleigh-plesset/rpe-offset1p5-alternative-minimal-krom-index3__20260827T193559.632192Z.npz`

SHA-256:
`ec17516e85b2d9328c5a48dd9ffa00e1910d89ec3d61839c2cfe4f272753c27c`

Viewer comparison:
`viewer_outputs/rpe-offset1p5-current-vs-alternative-minimal-krom.html`.

## 8. Detailed diagnosis

### 8.1 Hidden initial forcing mismatch

The feedback error is zero at `tau=0`, hence both KROM controllers initially
command numerical zero (`0` for the current pair and `-4.6e-15` for the
alternative). The nonlinear plant inverse requires `-0.6488766 p0` at that
same instant. A controller using only instantaneous state error has no preview
of the target acceleration. This creates an immediate phase deficit even
though the starting point appears ideal in a phase portrait.

This mechanism alone does not explain the persistent large error: feedback
could recover after an initial lag if the chosen coordinates represented the
subsequent physical error faithfully. The remaining diagnostics show why it
does not.

### 8.2 Target-path folding and error compression

For the current eigenfunction, the real two-coordinate map
`(r,u)->(Re Psi,Im Psi)` was differentiated at all 10,001 target samples.

| Diagnostic | Current pair 1/2 | Alternative pair 3/4 | Successful small target, current pair |
|---|---:|---:|---:|
| minimum `abs(det nabla Psi)` | `4.51309e-8` | `2.93441e-8` | `6.67524e-5` |
| median `cond(nabla Psi)` | `4.63617` | `4.55451` | `23.0382` |
| maximum `cond(nabla Psi)` | `17896.6` | `6953.91` | `52.0210` |
| time of maximum condition | `tau=5.736` | `tau=6.886` | `tau=8.980` |

For the current pair, minimum determinant occurs at `tau=4.314`, near
`(r,u)=(1.22546845,-0.46962052)`. Maximum condition occurs at `tau=5.736`,
near `(1.63708244,-1.54980448)`. For the alternative, minimum determinant is
at `tau=7.932`, near `(1.24008302,-0.62983444)`.

The consequence is visible at final time for the current pair:

```text
target state          = (1.5871121774, -0.2931379333)
controlled state      = (0.8635307540, -0.1377270212)
physical error norm   = 0.7400828520
lifted error norm     = 0.0067416485
```

The alternative ends with physical error norm `0.729190` but lifted error norm
only `0.00620875`. Thus both minimal maps can alias substantially different
physical states. With no validated observation map, the KROM has no contract
requiring small lifted error to imply small radius/velocity error.

### 8.3 Incorrect pressure phase and waveform

For every target sample, the nonlinear RPE inverse pressure was calculated
from the reference acceleration and compared with the archived command.

| Diagnostic | Current pair 1/2 | Alternative pair 3/4 | Successful small target |
|---|---:|---:|---:|
| physical inverse-pressure RMS | `0.900913` | `0.900913` | `0.260652` |
| archived command RMS | `1.05114` | `1.09174` | `0.257554` |
| actual/required correlation | `0.419448` | `0.395818` | `0.999994` |
| actual/required RMSE | `0.937783` | `0.979970` | `0.00322684` |

The failed controllers do not lack gross command magnitude: their RMS exceeds
the target's inverse-pressure RMS. The problem is that the pressure is applied
with the wrong phase and shape. In contrast, the successful small-target
command is almost perfectly correlated with the physical inverse pressure.

The reference derivative was also projected onto `B_psi` by least squares.
For the current pair, inferred lifted input and physical inverse input correlate
at only `0.827489` and differ by RMS `0.745039`. The median relative residual
left outside the one-input `B_psi` direction is `0.359251`, reaching nearly
one. The alternative is worse: correlation `0.811762`, RMSE `0.822281`, and
median residual `0.476609`. The selected lifted dynamics/input direction do
not encode this forced reference path consistently enough for accurate
feedback tracking.

### 8.4 What did not fail

- **Not gross training-box extrapolation:** the entire target envelope lies
  inside both KROMs' archived axis-aligned training domain.
- **Not disappearance of `B_psi`:** current-pair norm stays between
  `0.0209060` and `0.121328`; alternative stays between `0.0221311` and
  `0.105685`.
- **Not loss of lifted controllability:** the two-state controllability matrix
  remains full rank. Its condition is effectively constant at `4.82825` for
  the current pair and `5.59920` for the alternative.
- **Not CARE failure or unstable computed poles:** CARE succeeds and both
  equilibrium closed-loop pole pairs are strictly negative real values.
- **Not insufficient total pressure RMS:** actual command RMS is larger than
  the reference inverse-pressure RMS.
- **Not an ordinary startup-spike case:** the initial command is zero, not a
  large spike. The problem is initially missing required forcing and later
  coordinate aliasing/waveform error.

## 9. Why apparently harder targets can succeed

Target difficulty is not ordered by amplitude alone.

- A fixed `(2,0)` target begins with nonzero state error, immediately exposes
  a feedback direction, and needs no sustained periodic forcing after capture.
- Remote/high-velocity targets begin with a large state error and the uncapped
  controller produces large startup commands. Their physical acceptability is
  a separate question, but they do not hide the initial forcing demand.
- The small equilibrium-centered multisine matches initial state, velocity,
  acceleration, and required pressure. Its orbit also avoids the severe
  Jacobian folding seen here.
- The offset-1.5 target is unusual in the variables seen by pure state
  feedback: it matches initial state and velocity while demanding a large
  immediate acceleration, then passes through regions where both tested
  minimal eigenfunction maps become nearly singular.

This explains why the target can look modest in physical phase space yet be a
much harder task for these particular lifted coordinates.

## 10. Historical comparator and interpretation boundary

An older specialized-KROM experiment reported radius RMSE `0.01096827` and
velocity RMSE `0.01239878` on this target. It is important evidence that the
target itself is not intrinsically uncontrollable. It is not evidence that the
current implementation should match it without further isolation: that case
used another KROM, different scaling/tuning, a sampled control policy, and
pressure limiting.

The correct next comparison changes exactly one factor at a time. Do not use
the older result to justify copying its full pathway wholesale or silently
altering the retained clean-isolation baseline.

## 11. Implications for eigenfunction selection

Eigenvalue frequency, damping, and H-DMD terminal-amplitude energy are not
sufficient control-selection criteria. Candidate eigenfunctions or small
subsets should be screened on the intended reference family for:

1. local injectivity of the split-real map, including determinant/singular
   values and condition number of `nabla Psi` along target paths;
2. separation of physically distinct states relevant to the controller;
3. `B_psi` magnitude and controllability conditioning;
4. forced-reference consistency,
   `dPsi_ref/dtau - A Psi_ref - B_psi u_required`;
5. agreement between lifted least-squares input and nonlinear inverse input;
6. held-out state observation/reconstruction accuracy;
7. robustness across nearby target phases, amplitudes, and offsets.

Thresholds should be derived from successful and failed validation cases, not
chosen arbitrarily. The present comparison already supplies useful anchors:
maximum target-path Jacobian condition near `52` succeeds, while conditions in
the thousands coincide with failure. Those values are evidence for designing
a gate, not yet a universal acceptance threshold.

## 12. Recommended staged follow-up

1. Keep this exact target, plant, solver, controller law, `Q`, and `R` fixed.
2. Screen every independent full-spectrum eigenfunction offline on the target
   path. This requires no controlled integration and is much cheaper than
   blindly solving every candidate.
3. Rank candidates by combined physical-state separation, Jacobian singular
   values, forced-reference residual, input agreement, and `B_psi` health.
4. Re-run only the best individual candidates under this exact contract.
5. Build two- and then small multi-eigenfunction KROMs from complementary
   candidates. Require a held-out state-observation map before interpreting
   small lifted error as successful physical tracking.
6. Compare the best new KROM with the historical specialized KROM using a
   one-factor-at-a-time matrix: representation first, then tuning, evaluation
   cadence, and pressure policy.
7. Only after representation adequacy is established, study optional
   target-preserving engagement or reference-preview mechanisms separately.
   Do not fold them into this baseline.

This sequence preserves the requested law form and isolates whether a better
observable set can recover accurate feedback-only tracking. It also avoids
using arbitrary command clipping to disguise a coordinate failure.

## 13. Reproduction commands

From the repository root:

```text
python scripts/rerun_rpe_prior_control_cases.py
python scripts/try_rpe_alternative_minimal_krom.py
python scripts/diagnose_rpe_offset1p5_eigenfunction_failure.py
```

The first two commands run integrations and publish new immutable archives;
they do not overwrite the artifacts listed above. The third command is a
read-only diagnostic over the retained archives and rewrites only the
machine-readable diagnostic JSON.

For exact reproduction, verify artifact SHA-256 values before comparison. The
JSON is the authoritative compact numeric record; this report supplies the
interpretation and experimental boundaries.

## 14. Related retained reports

- `RPE_PRIOR_CONTROL_CASES_CLEAN_RERUN_2026-08-27.md`: the eight-target
  comparison that exposed this lone remaining failure.
- `RPE_MINIMAL_KROM_KLQR_STABILITY_CASE_STUDY_2026-08-26.md`: earlier KLQR
  stability, tuning, KROM-selection, and historical specialized-target work.
- `../notes/RPE_MINIMAL_KROM_KLQR_STABILITY_REPRO_NOTES_2026-08-26.md`:
  compact reproduction notes for that earlier study.
- `../RPE_CLEAN_ISOLATION_KLQR_SOP.md`: current clean-isolation baseline.
- `RPE_CLEAN_ISOLATION_INITIAL_TRANSIENT_DIAGNOSIS_2026-08-27.md`: separate
  large-startup-transient diagnosis; do not conflate it with this zero-initial-
  command failure.


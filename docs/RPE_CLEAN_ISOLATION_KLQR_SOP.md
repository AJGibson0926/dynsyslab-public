# RPE clean-isolation H-DMD/SINDy/KLQR standard operating procedure

**Policy name:** `rpe-clean-isolation-klqr`  
**Post-processing policy:** `rpe-clean-isolation-hdmd-sindy`  
**Status:** retained reproducibility baseline  
**Purpose:** reproduce the validated large-target, feedback-only RPE/KLQR pathway without making it the default for unrelated systems.

**Current decision:** this is the normal/default workflow for this RPE control
study. Its large initial transient is deliberately retained. Viewer displays
the full command and automatically supplies a smaller-scale detail inset; the
transient is not clipped, ramped, filtered, or removed from the archive.

The smoother production sequence is a deferred, explicit opt-in add-on. It is
documented separately in
`case-studies/RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE_2026-08-27.md` and must
not be selected implicitly by this SOP or its reproduction command.
That report distinguishes the already tested reference-shaped case from the
still-untested, target-preserving controller-modulation path; the former is not
evidence that the latter works.

This SOP is immutable in intent. Experimental changes must receive distinct
policy names, archive labels, and Viewer outputs. Never overwrite any artifact
listed here.

## 1. Validated upstream inputs

| Stage | Artifact | SHA-256 |
|---|---|---|
| Dimensionless BDF training trajectory | `archives/solutions/rayleigh-plesset/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz` | `e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba` |
| Full energy-sorted H-DMD spectrum | `archives/derived/koopman-hankel-dmd/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz` | `c2fabe83b37b94ace1abb5b7d4fdd11cad14dee97c476d2fbb00e5e63eb0bf4b` |

The upstream model uses exact `K=4/3`, the exact seeded dimensionless IC
`(2.3333865986017348,0)`, SciPy BDF with `rtol=atol=1e-9`, `tau=0..25`,
and uniform published spacing `dtau=.001`.

H-DMD uses radius observable 0, a `20000 x 2000` window-snapshot Hankel pair,
`gesvd`, absolute tolerance `1e-10`, full rank, projected modes, and descending
terminal-amplitude energy ordering. It does not apply a rank cap.

## 2. Polynomial eigenfunction fit

Use the named immutable policy in `postprocessing/policies.py`:

- dimensionless integration coordinates;
- the complete energy-sorted H-DMD spectrum;
- continuous generators computed as `lambda=log(mu)/dtau`;
- fourth-order polynomial dictionary;
- no constant feature;
- no state normalization;
- coefficient threshold `1e-6`;
- one sparse-threshold pass;
- system-RHS generator data;
- conjugate-fit reuse;
- guarded Gram-Hermitian eigensolve;
- coefficient-first immutable archive storage.

Retained SINDy archive:

`archives/derived/koopman-with-sindy/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-with-sindy__20260827T183645.394109Z__c4dde209.npz`

SHA-256: `6e3725beac327fa749eeb0b0250cefa816d158cbbbc2d6fdea0d1b6bfa979119`

## 3. Minimal KROM selection

Select the first positive-frequency complex eigenfunction by original full-
spectrum index and include its conjugate dual. In the retained full spectrum
this is indices 1 and 2:

```text
mu     = 0.9997777985801627 +/- 0.004817601141183756 i
lambda = -0.21061644578341268 +/- 4.818634561556464 i per tau
```

Selection audit correction (2026-08-27): this pair is at retained terminal-
amplitude-energy-order position `999` of `1001`, not first in energy order. The
original extraction code and archived eigenvalues selected full-spectrum index
`1` correctly, but the policy label incorrectly described that operation as
energy-first. The immutable archive is retained; this correction changes only
the lineage description and future reproduction metadata.

Retained KROM:

`archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz`

SHA-256: `8c3f221f563a82271e2995c8cc7e627b37f28b54f9feb93c59ae28050ae4ef86`

## 4. Controller and solve contract

The authoritative settings live in `control/policies.py`.

```text
plant/controller coordinates: dimensionless
input: u* = p_control/p0
initial state: (r,u)=(1,0)
target: rpe_multisine_offset_4_mixed
Q: 5000 I2
R: 1e-4
law: u* = -K(x) [Psi(x)-Psi(x_ref)]
Bpsi(x): recomputed at every controller evaluation
CARE/K(x): recomputed at every controller evaluation
evaluation: continuous_rhs
solver: SciPy BDF, rtol=atol=1e-9
horizon: tau=0..10
published grid: dtau=.001
feedforward: none
control cap: none
startup ramp: none
target-domain policy: explicit ignore for this diagnostic extrapolation
```

The no-ramp, out-of-domain target is retained only to preserve the isolation
baseline. For the present research program it remains the established default
so comparisons stay close to the validated pathway. It is not a claim that the
startup pressure is physically suitable for production.

## 5. Required output and verification

Every SOP run must archive solution, target, accepted control history, raw
control history, and exact controller-call history. Control metadata must say:

```text
coordinate_system = dimensionless
units = p_control/p0
scale_to_dimensional = 100000 Pa
evaluation_policy = continuous_rhs
```

The standard Viewer output contains a square phase portrait, separate radius
and velocity time series with the target overlay, and a separate dimensionless
control panel. When the initial transient dwarfs later commands, Viewer retains
the full y-range and adds an automatic inset showing the smaller-scale detail.
This is display-only behavior and does not modify or conceal archived samples.

Retained baseline solution:

`archives/solutions/rayleigh-plesset/rpe-clean-isolation-dimensionless-klqr-offset4__20260827T183646.749498Z.npz`

SHA-256: `76c36044ae81d9c42f66af71d1a6df567c1fccb688d7b8c6c28974a56ce5e771`

Viewer: `viewer_outputs/rpe-clean-isolation-klqr-latest.html`

Summary: `viewer_outputs/rpe-clean-isolation-klqr-latest.json`

## 6. Reproduction route

The case-specific script has been retired under the suite's script-lifecycle
policy. Reproduce the recorded configuration through the Solver and
Post-processor UI or their shared API, using the exact archive paths,
parameters, controller metadata, and hashes above. New outputs remain immutable
downstream archives; the verified trajectory and H-DMD records are inputs and
must not be duplicated or rewritten.

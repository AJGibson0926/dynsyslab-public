# Linear-systems analysis of echo K-ROMs (2026-11-15)

Demonstrates the [linear-systems toolkit](../../../../docs/architecture/LINEAR_SYSTEMS_TOOLKIT.md)
on three existing viscous-RPE echo K-ROMs. The sources are local archives of the
author's checkout, not bundled sample data; edit `sources` to rerun elsewhere.

```powershell
dynsyslab-cli run examples\native\case-studies\rpe-echo-linear-systems-20261115\echo-klqe-krom.json --archive-root archives
```

Results are in `archives/derived/linear-systems/` (and `archives/derived/klqe/`
for section 4); open the linear-systems results in Viewer under
*Control & estimation → Linear systems* (views: Hankel singular values,
Gramians, Frequency response, Poles and zeros).

## 1. Echo KLQE K-ROM (`echo-klqe-krom.json`)

Source: the 5-mode `krom_extract` (observable-energy rule) of the joint
radius/velocity/echo H-DMD on the late unforced baseline, used by the
2026-11-14 echo→state KLQE. Measurement: far-field echo. Balancing input:
process noise on every lifted coordinate.

- The split realization has 10 real states for 5 complex modes: every pole
  appears twice. The structure report gives minimal order **5**, observable
  dimension 5, detectable, with PBH in agreement.
- Hankel singular values: 0.0983, 0.0973, 0.0405, 0.0401, 0.0147, then ≤ 1.4e-9.
  The echo sees five dimensions.
- Order-4 balanced truncation: actual ‖G − G_r‖∞ = 0.02931, bound
  [σ₅, 2Σσ_tail] = [0.0147, 0.0293]. Full H∞ = 0.196 at ω = 6.283 (the
  fundamental), H2 = 0.056.
- The reduced LQE gain (W = I, V = 1e-4), lifted into the full model, gives a
  stable observer.

## 2. Joint H-DMD parent, top 40 modes (`echo-hdmd-top40.json`)

Source: the parent 2000-mode joint H-DMD (modes 0–39 by terminal amplitude
energy), with echo measurement and process-noise balancing.

- 78 real states. With `stability_margin` 1e-6, **5 modes are kept outside
  balancing**: an equilibrium mode λ ≈ 1.9e-8 and a **growing** pair
  λ = +0.131 ± 6.037i, split twice. This matches the earlier finding that growing
  H-DMD modes on the full decay break KLQE.
- PBH: not observable and **not detectable** (smallest observability margin
  9e-15, from the doubled copies), so the growing pair cannot be corrected from
  the echo. The Krylov staircase reports 78 observable states because the copies
  are only nearly equal; `pbh_agrees: false` flags this.
- Stable-part HSVs decay slowly (0.103, 0.102, 0.044, 0.044, 0.027, …): the
  stable part has no sharp numerical order. Order 9 (4 stable + 5 retained):
  actual error 0.055, bound [0.027, 0.47]. Stable-part H∞ = 0.207 at ω = 6.283.

## 3. Echo-jet KLQR K-ROM (`echo-jet-klqr-krom.json`)

Source: the latest (e, e′) echo-jet `input_map` K-ROM with a constant `B_ψ`
(`noisy-d3-q8-rts-bpsi`). Outputs: the lifted split; LQR with Q = 5000 I,
R = 1e-4 (continuous).

- 4 real states for one conjugate pair; minimal order 2, controllable dimension
  2, stabilizable; HSVs 10.40, 10.25, 0, 0.
- λ = −0.0825 ± 5.781i; H∞ = 20.6 at ω = 5.78.
- LQR closed-loop poles −17040 and −5.78 on the controllable pair. The duplicated
  (uncontrollable) copy stays at −0.0825 ± 5.781i. The reduced (order-2) gain
  lifted back to the full model gives the same closed loop.

## 4. Reduced-order KLQE on the held-out record (`echo-klqe-reduced-order-{4,3}.json`)

The offline `klqe` stage of the 2026-11-14 echo→state run (same K-ROM as model 1,
held-out record `trajectory-2`, window [30, 80], W = I, V = 1e-4), now with
`balanced_truncation`. The observer still estimates all 10 lifted coordinates;
only its gain is designed at reduced order.

| Gain | Bound [σ_{r+1}, 2Σσ_tail] | State RMSE (r, ṙ), whole window | RMSE, late half |
|---|---|---|---|
| full order (2026-11-14, `…075717Z`) | — | 0.0096, 0.0414 | 0.00079, 0.00030 |
| order 4 (`…031828Z`) | [0.0147, 0.0293] | 0.0079, 0.0424 | 0.00027, 0.00029 |
| order 3 (`…031829Z`) | [0.0401, 0.1094] | 0.0080, 0.0440 | 0.00025, 0.00031 |

Truth RMS is 1.00 (radius) and 0.197 (velocity). The slowest observer pole
stays at Re λ = −0.0604 in all three; it belongs to the duplicated copies the
echo cannot see. The order-4 and order-3 gains estimate the state about as well
as the full gain (radius slightly better late in the window).

## 5. Echo-jet output feedback (`echo-jet-output-feedback.json`)

Model 3 with both LQR (Q = 5000 I, R = 1e-4) and LQE (W = I, V = 1e-4) weights
and `reduction_order` 2: the coordinated LQE→LQR (LQG) design (`…031830Z`).

- Full 4-state compensator: the closed-loop poles are the LQR poles (−17040,
  −5.78) together with the LQE poles (−100 ± 5.78i). The uncontrollable and
  unobservable copy stays at −0.0825 ± 5.78i (twice). The closed loop is stable.
- The order-2 compensator, designed on the balanced reduced model and closed
  around the full model, places the same poles (−17040, −5.78, −100 ± 5.78i). Only
  one copy of the duplicate pair is left, because the compensator has two
  states instead of four.

## Data roles

Sections 1–3 and 5 are model analyses (no trajectories read). Section 4 filters
the echo of a held-out unforced record (OFF-e); its radius and velocity (OFF-r)
are used only to score the estimate. Inherited roles: models 1
and 2 have eigenvalues and eigenfunctions from H-DMD on unforced radius,
velocity and echo records (OFF-r, OFF-e). Model 3 is an echo-eigenfunction K-ROM
whose `B_ψ` was regressed from forced records (see `CONTROL_DATA_ROLES.md`).

Superseded exploratory runs from the same session: `…005423Z` (model 1 at the
default numerically minimal order 9), and `…005456Z` and `…005718Z` (model 2 at
order 8, and before the PBH fields were added).

# Receiver arrays: retarded echoes, state readouts and ring errors

Several acoustic receivers of one bubble, synchronized on one clock, a readout of the bubble's
state from their pressures, and errors scored oscillation by oscillation. Ported from the
receiver-reconstruction side-car (`tools/experiments/receiver_reconstruction`) into the core and
three post-processing methods; a sensor design (which receivers, which pressure terms, which
readout) is an experiment over them.

| Piece | Core | Method |
|---|---|---|
| Multi-receiver retarded synthesis | `core/ReceiverArray.hpp` `rpeReceiverArray` | `receiver_signals` |
| Linearized delay inverse | `core/ReceiverArray.hpp` `linearDelayInverse` | `state_readout` `mode: linear` |
| k-NN readout with neighbour provenance | `core/NeighborReadout.hpp` `NeighborReadout` | `state_readout` `mode: neighbors` |
| Ring (max-to-max) errors | `core/NeighborReadout.hpp` `oscillationMaxima`, `intervalErrors` | `ring_errors` |

Test: `dynsyslab_receiver_pipeline_tests` (retardation of the echo between receivers, a held-out
ringdown recovered by two receivers, the linear inverse on a small ringdown, per-ring errors).

## Receiver signals

A receiver at range D records at time t the pressure emitted at t − D/c (Lighthill 1952; Ffowcs
Williams and Hawkings 1969). On the nearest receiver's emission clock τ, receiver i holds the
pressure emitted (D_i − D_near)/(c T) earlier, T seconds per model time. Each term of the RPE
receiver catalog (`core/AcousticReceiver.hpp`: monopole, spherical near-field correction, dipole
induction and radiation) is evaluated at the two bracketing source samples and interpolated
linearly; the record starts once every path has history.

```json
{"name":"receiver_signals","settings":{"receivers":[[0.005,0,0],[0,0.006,0]],"sound_speed":1480,
  "stride":40,"window":[0,100],"terms":["monopole_pressure"],"range_corrected":true,"contrast":false}}
```

`values` are the features a readout sees: the sum of the selected terms per receiver, times the
range when `range_corrected` (so the far-field monopole is receiver independent). With `contrast`
(two receivers) a third feature C = D₂p₂(t) − D₁p₁(t − Δ) uses the near receiver's pressure at the
far receiver's emission: under the spherical model the common monopole cancels. Every term per
receiver (`<term>`, `near_at_<term>`), the source states (`source_states`, labels for scoring
only), `ranges_m`, `lag_tau` and `ring_bounds` (radius maxima: radial velocity turning from positive
to nonpositive, a resting release counted) are kept. Only free (unforced) records are accepted:
the catalog uses the intrinsic zero-input acceleration.

## State readout

```json
{"name":"state_readout","settings":{"training":["<features.npz>"],"mode":"neighbors","neighbors":8,"scales":"shared"}}
{"name":"state_readout","settings":{"training":["<features.npz>"],"mode":"linear","equilibrium":[1,0]}}
```

- **neighbors**: k-nearest-neighbour regression (Cover and Hart 1967) in asinh(feature/scale)
  coordinates (linear near zero, logarithmic in the tails; Burbidge, Magee and Robb 1988), weights
  ∝ 1/d² (Shepard 1968), exact k-d tree search (Bentley 1975) with ties broken by row index.
  `"shared"` scales every feature by the training RMS of the first. Each estimate keeps its
  training rows (`neighbor_rows`), weights and nearest distance. Training labels are used only to
  fit; the query's own states are carried untouched for scoring, and a query that is also a
  training record is refused.
- **affine**: x = W(y − ȳ) + x̄ by least squares on the training rows (optional `ridge`; Hoerl and
  Kennard 1970); the weights and bias are kept (`readout_weights`, `readout_bias`).
- **linear**: near an equilibrium x₀ with A = Df(x₀) and feature gradients H_i (the receivers'
  selected terms, times their ranges), receiver i observes H_i exp(−A δ_i)(x − x₀); the stacked
  rows O give x − x₀ = O⁺(y − y₀) (an observability matrix of delayed outputs; Kalman 1960). It is
  refused when O does not observe every state and for features with the contrast column.

## Ring errors

```json
{"name":"ring_errors","settings":{"center":[1,0],"windows":{"early":[0,5],"late":[76,100]},"role":"validation"}}
```

Per half-open interval between consecutive radius maxima: sample count, a partial flag (the first
interval can be cut by the propagation delay), error and reference sums of squares and the relative
RMS percentage per component, the reference deviation measured from `center`. The pooled
whole-record and window relative RMS percentages are under the `role` key, so a `select` stage can
choose a sensor design by them.

## Sensor design

The comparison the side-car's `study.py` ran (cases monopole, near field, induction, radiation and
their sums; one receiver, two receivers, two plus contrast, the linear inverse) is an experiment:
a `sweep` over `terms`, `receivers`, `contrast` and the readout `mode`, each point running
`receiver_signals` on the training and held-out records, `state_readout` and `ring_errors`, and a
`select` stage choosing the design with the smallest held-out error
(`examples/native/experiments/receiver-sensor-design/design.json`).

## References

- M. J. Lighthill, On sound generated aerodynamically I, Proc. R. Soc. A 211:564–587, 1952.
- J. E. Ffowcs Williams and D. L. Hawkings, Sound generation by turbulence and surfaces in
  arbitrary motion, Phil. Trans. R. Soc. A 264:321–342, 1969.
- T. Cover and P. Hart, Nearest neighbor pattern classification, IEEE Trans. Inf. Theory 13:21–27, 1967.
- D. Shepard, A two-dimensional interpolation function for irregularly-spaced data, Proc. 23rd ACM
  National Conference, 1968.
- J. B. Burbidge, L. Magee and A. L. Robb, Alternative transformations to handle extreme values of
  the dependent variable, J. Am. Stat. Assoc. 83:123–127, 1988.
- J. L. Bentley, Multidimensional binary search trees used for associative searching, Commun. ACM
  18:509–517, 1975.
- R. E. Kalman, Contributions to the theory of optimal control, Bol. Soc. Mat. Mexicana 5:102–119, 1960.

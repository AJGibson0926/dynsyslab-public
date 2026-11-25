# Echo-only input-output delay control

`io_delay_model` and `io_delay_lqr` control a system from one scalar measurement and the
controller's own past commands. There is no SVD, latent state, eigenfunction or estimator.

## Model: `io_delay_model`

Post-processing method. The sources are forced training records. It fits, by ridge least squares
on a decimated clock (ARX; Ljung, *System Identification*, 2nd ed., 1999),

$$
y_{k+1}=\sum_{i=0}^{n_y-1}a_i\,y_{k-i}+\sum_{j=0}^{n_u-1}b_j\,u_{k-j}+c .
$$

```json
{"name":"io_delay_model","settings":{"observable":"far_field_acoustic_echo","decimation":5,
 "ny":30,"nu":30,"relative_ridge":1e-10,"fit_offset":true,"validation_sources":["${s14.archive}"]}}
```

- Regression rows never cross a record boundary.
- `u[k]` is held over `[k, k+1)`. A decimated interval whose fine-clock input changed is flagged,
  and every row that uses it is rejected (`rows_rejected_for_hold`).
- The ridge penalty is `relative_ridge * rows * |theta|^2` on regressors scaled to unit deviation.
  The offset is not penalized. The solve uses sequential QR, not normal equations.
- The archive (`derived/io-delay-model`) holds `a`, `b`, `c`, the clock, the companion realization
  `companion_A/B/C/d` for $z_k=(y_k,\dots,y_{k-n_y+1},u_{k-1},\dots,u_{k-n_u+1})$, the open-loop
  poles, training and held-out one-step and free-run errors, and source provenance.

## Controller: `io_delay_lqr`

A discrete LQ tracker (Anderson and Moore, *Optimal Control: Linear Quadratic Methods*, 1990) for
the cost $\sum_k Q(y_k-r_k)^2+Ru_k^2$:

$$
u_k=-Kz_k+f_k,\qquad f_k=S^{-1}B^{\top}(v_{k+1}-Pd),\qquad
v_k=C^{\top}Qr_k+(A-BK)^{\top}(v_{k+1}-Pd),
$$

with $S=R+B^{\top}PB$. The steady $P$ and $K$ come from the suite's DARE (`designLQR`). The
delay-state idea follows Korda and Mezic, *Automatica* 93 (2018).

```json
{"kind":"io_delay_lqr","settings":{"model":"${inputs.model}","Q":1,"R":1e-4,"final_time":300,
 "warmup_input":0,"control_limit":null},
 "target":{"space":"observable","kind":"from_state","state_target":{"kind":"waveform","...":"..."}}}
```

- The solver must use `control_evaluation: sampled`, with `control_sample_time` equal to the model
  clock.
- The state target is lifted to `y_ref` once, offline, through the named observable. Online, the
  law reads only the scalar measurement and its past applied commands. The state argument only
  simulates the sensor.
- `warmup_input` is held until the history window fills. `control_limit` is applied before the
  command enters the history, so the model sees the input the plant received.
- `final_time` sets the known-target horizon. `preview_tail` (default 5) extends the recursion so
  the end of the run is not treated as a terminal condition.
- `model_binding_evidence` records the model path and SHA-256, the orders, clock, weights,
  closed-loop spectral radius and Riccati residual.

The model is linear and time invariant. It carries no guarantee on a nonlinear plant away from the
forced-record operating range.

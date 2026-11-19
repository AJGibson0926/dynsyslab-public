# Integral action and PID in the controls framework

Two general components remove steady offsets that a model-based state feedback leaves when
its model is imperfect (an unmodelled input, a biased map, a readout error).

## Components (core)

| Component | Header | What it does |
|---|---|---|
| `PIDCompensator` | `core/PIDCompensator.hpp` | PID on any error vector (p channels → m outputs), sampled at monotonic times: u = Kp e + Ki ∫e + Kd ė. Options: first-order derivative filter, per-channel integrator limit, output limits, conditional-integration anti-windup. |
| `PIDController` | `core/ClassicalControl.hpp` | State-feedback PID (the `pid` controller kind) on e = x_target − x, built on `PIDCompensator`. |
| `solveContinuousLQI` | `core/LQR.hpp` | LQR with integral action: augments the state with q, q′ = C_i x − y_ref, penalty blkdiag(Q, Q_i). Control u = −K_x(x − x_ref) − K_q q. |
| `KoopmanLQRController::setIntegralAction` | `core/LQR.hpp` | LQI on a lifted output y = y₀ + C z (a KROM readout). The gain is designed on the augmented KROM. |
| `KoopmanLQRController::setPIDCorrection` | `core/LQR.hpp` | Adds PID(y_target − y(z)) to the KLQR input. |

Integral states make a controller stateful: it needs a sampled, monotonic clock and one
instance per trajectory. The factory enforces both.

## Settings

`pid` controller kind (state feedback):

| Key | Meaning |
|---|---|
| `Kp`, `Ki`, `Kd` | 1 × state gains (existing) |
| `derivative_filter_time` | τ of the derivative filter (0 = raw backward difference) |
| `integrator_limit` | |∫e| ≤ limit, scalar or per state |
| `control_limit` | the actuator cap; with it, PID output limits and anti-windup follow the cap |
| `conditional_integration` | anti-windup on/off when a cap is set (default true) |

KLQR adapters (`controller.settings`, both full-state and KLQE-estimated feedback):

```json
"integral_action": {"reconstruction": "state", "components": [0], "penalty": 100,
                    "integrator_limit": [1.0], "target_components": [0]},
"pid_correction":  {"reconstruction": "state", "components": [0], "Kp": 0, "Ki": 2, "Kd": 0,
                    "derivative_filter_time": 0.05, "integrator_limit": [1.0],
                    "output_limits": 0.5, "conditional_integration": true}
```

The controlled output is rows `components` of the named KROM reconstruction (with its
intercept), applied to the lifted state the KLQR acts on: ψ(x) under full-state feedback, the
KLQE estimate ẑ under estimated feedback (then y is the estimated output, e.g. r̂, and the
loop still uses only the estimator's measurement online). It is compared with components
`target_components` of the dynamical-state target. Integral action cannot be combined with
`balanced_truncation`. Archives record the settings in `klqr_augmentations`.

## Tests

`dynsyslab_classical_control_tests`: integrator limit, derivative filter, conditional
integration and recovery, LQI removing a constant-disturbance offset that LQR keeps.
`dynsyslab_koopman_control_tests`: KLQR integral action and PID correction remove an
unmodelled input offset on a scalar KROM.

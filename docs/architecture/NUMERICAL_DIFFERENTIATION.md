# Noise-robust numerical differentiation and measurement noise

An echo-only lab controller has to form time derivatives (ė, and ψ̇ for eigenfunction and B_ψ fits) from a sampled, noisy signal. Plain finite differences amplify white noise roughly in proportion to frequency, by about 1/(√2·Δt) per unit of noise, so even 0.1% noise swamps ė at Δt = 0.001 τ. DynSysLab now has one library for this (`native/include/dynsyslab/core/NumericalDifferentiation.hpp`), used offline and online, plus a matching noise generator.

## Methods

Every smoothing method takes one physical knob, `cutoff` (cycles per τ). Content above `cutoff` is suppressed; at `cutoff` the smoothing gain is about one half. Explicit parameters (`window`, `sigma`, `lambda`, `process_noise`, …) override it.

| Method | Family | Idea | Main parameter |
|---|---|---|---|
| `central` | finite difference | (y₊ − y₋)/2Δt; the baseline | — |
| `savitzky_golay` | local polynomial | least-squares polynomial on a sliding window; derivative of the fit | `order`, `cutoff` or `window` |
| `gaussian` | smoothing kernel | convolution with the derivative of a Gaussian | `cutoff` or `sigma` |
| `holoborodko` | smoothing kernel | Holoborodko's smooth noise-robust differentiator: exact on cubics, maximally flat, guaranteed high-frequency suppression | `points`, `cutoff` or `stride` |
| `weak` | weak form / modulating function | integrate by parts against a compact bump φ = (1−s²)^p: y′ ≈ −Σφ′y / Σφ. The derivative falls on φ, never on the data | `exponent`, `cutoff` or `window` |
| `butterworth` | filter, then difference | zero-phase (forward-backward) Butterworth low-pass, then central difference | `order`, `cutoff` |
| `whittaker` | penalized least squares | smoothing spline: min ‖z−y‖² + λ‖D^d z‖² (banded solve) | `penalty_order`, `cutoff` or `lambda` |
| `kalman` | state estimation | Rauch–Tung–Striebel smoother of a d-times-integrated white-noise model (d = 3: constant acceleration) | `order`, `cutoff` or `process_noise`; `measurement_noise` (estimated if 0) |
| `spectral` | Fourier | FFT derivative of the even extension with a raised-cosine low-pass | `cutoff` |
| `tv_denoise` | total variation | Condat's exact 1-D TV denoising, then central difference | `weight` (default 3 × noise std) |
| `tikhonov_integral` | integral formulation | find u = y′ with min ‖A u − y‖² + α‖Du‖², A = cumulative integral (matrix-free CG). Data are integrated, never differentiated | `cutoff` or `weight` |
| `tvr` | integral formulation + TV | Chartrand's total-variation regularized differentiation (lagged diffusivity) | `cutoff` or `weight`, `iterations` |
| `backward`, `backward2` | causal finite difference | 1st / 2nd-order backward differences | — |
| `causal_savitzky_golay` | causal local polynomial | polynomial fit of the last window, derivative at the newest sample (the algebraic / Jacobi differentiator family) | `order`, `cutoff` or `window` |
| `causal_butterworth` | causal filter | causal low-pass, then backward difference | `order`, `cutoff` |
| `causal_kalman` | causal state estimation | forward Kalman filter of the integrated white-noise model | `order`, `cutoff`, `measurement_noise` |
| `levant` | sliding mode | Levant's 2nd-order robust exact differentiator | `lipschitz` (bound on \|y‴\|) or `cutoff` |

Some of these are close relatives:
- Whittaker, the RTS smoother and the smoothing spline are the same estimator for matching parameters.
- `gaussian` and `weak` differ only in the bump.
- `causal_savitzky_golay` is the polynomial/algebraic differentiator.

## Noise models

| `kind` | Model | Parameters |
|---|---|---|
| `white` | Gaussian, independent samples | `level` (RMS relative to the signal RMS) |
| `uniform` | bounded uniform | `level` |
| `ar1` | red / Ornstein–Uhlenbeck | `level`, `correlation_time` |
| `pink` | 1/f (Kellet filter) | `level` |
| `random_walk` | drift | `level` (RMS accumulated per τ) |
| `quantization` | ADC rounding | `bits` over the signal range, or an absolute step |
| `impulse` | spikes / dropouts | `probability` per sample, `level` |
| `sum` | composition, applied in order | `parts` |

`relative: false` makes `level` absolute. `seed` makes every realization reproducible.

## Where it plugs into the pipeline

| Place | Setting | Effect |
|---|---|---|
| `coordinate_transform` | `{"kind": "measured_derivatives", "observable": "far_field_acoustic_echo", "noise": {...}, "derivative": {...}, "values": "smoothed"\|"measured", "reference_derivatives": [...]}` | Measured coordinates (e, ė) from the noisy echo. Writes `coordinate_derivatives` (ė, ë) for SINDy `gamma_method: measured` and keeps `control_values` of forced records. `reference_derivatives` only scores the result. |
| coordinate SINDy | `"input_space": "coordinates"`, `"gamma_method": "measured"`, or `finite_difference` with `"derivative": {...}` | ψ̇ targets from robust derivatives; works with `sample_stride` / `time_window` |
| `input_map` | measured-coordinate K-ROMs: `fit: controlled_residual`, `"derivative": {...}` | B_ψ from forced measured records. ψ̇ by the chosen differentiator, with the input put through the same smoother |
| `control_target` | `"coordinate_observable": "echo_jet2"` | offline target lift for measured-coordinate K-ROMs |
| KLQE estimator | `"measurement_noise": {..., "scale": [...]}`, `"measurement_derivatives": [0]`, `"measurement_derivative": {"method": "causal_kalman", ...}` | noisy sensors and online causal ė |
| `derivative_benchmark` (post-processing) | `observable`, `reference_derivative`, `noises`, `methods` (each with a `cutoffs` grid), `seeds`, `trace_window` | sweeps methods × noise × cutoff on a record and stores errors and best traces |

## Results

See `docs/reports/NOISE_ROBUST_DIFFERENTIATION_2026-11-15.md`.

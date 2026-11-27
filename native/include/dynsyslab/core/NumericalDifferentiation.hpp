#pragma once
// Noise-robust numerical differentiation of sampled measurements, and measurement-noise models.
//
// Offline (non-causal, uniform samples): differentiate(y, dt, spec) returns the smoothed signal and
// its time derivative. Online (causal): CausalDifferentiator::step(y) returns the current derivative
// estimate from the samples seen so far.
//
// Every smoothing method can be tuned by one physical knob, `cutoff` (cycles per time unit): the
// frequency above which the method suppresses content. Method-specific parameters (window, sigma,
// lambda, q) override it when given.
//
// Offline methods
//   central             2nd-order central difference (no smoothing; the baseline)
//   central_high_order  central difference of accuracy order points-1 (points: odd stencil width,
//                       default 9 = 8th order), one-sided stencils of the same width at the ends;
//                       weights by Fornberg 1988, Math. Comp. 51:699-706
//   savitzky_golay      local least-squares polynomial (window, order), derivative of the fit
//   gaussian            convolution with the derivative of a Gaussian (sigma)
//   holoborodko         smooth noise-robust differentiator (N points, stride): exact on cubics,
//                       maximally flat low-pass response, guaranteed high-frequency suppression
//   weak                weak-form / modulating-function derivative: integrate by parts against a
//                       compact bump phi(s) = (1 - s^2)^p, y' ~ -sum phi'(t_j - t) y_j / sum phi
//   butterworth         zero-phase (forward-backward) Butterworth low-pass, then central difference
//   whittaker           penalized least squares (smoothing spline): min |z-y|^2 + lambda |D^d z|^2
//   kalman              Rauch-Tung-Striebel smoother of a d-times-integrated white-noise model
//                       (d = 2: constant velocity, 3: constant acceleration, 4: constant jerk)
//   spectral            FFT derivative of the even extension with a raised-cosine low-pass
//   tv_denoise          1-D total-variation denoising (Condat's exact algorithm), then central difference
//   tikhonov_integral   integral formulation: find u = y' with  min |A u - (y - y0)|^2 + alpha |D u|^2,
//                       A = cumulative integration (matrix-free conjugate gradients)
//   tvr                 Chartrand's total-variation regularized differentiation: same integral
//                       formulation with alpha |D u|_1 (lagged diffusivity)
// Any offline method can be preceded by Hampel despiking (despike = half-window in samples), which
// removes impulsive outliers that linear smoothers only spread out.
// Causal methods (also usable offline, to measure their lag)
//   backward, backward2           1st / 2nd order backward differences
//   causal_savitzky_golay         polynomial fit of the last `window` samples, derivative at the newest
//   causal_butterworth            causal Butterworth low-pass, then backward difference
//   causal_kalman                 forward Kalman filter of the integrated white-noise model
//   levant                        Levant's 2nd-order robust exact (sliding-mode) differentiator
#include <cstdint>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct DerivativeSpec {
    std::string method{"central"};
    double cutoff{0};            // cycles per time unit; 0 = method default / explicit parameters
    std::size_t window{0};       // samples (savitzky_golay, causal_savitzky_golay, weak); odd for centered
    std::size_t order{3};        // polynomial order (savitzky_golay), Butterworth order, Kalman model order d
    double sigma{0};             // gaussian, in time units
    std::size_t points{0};       // holoborodko N (odd >= 5)
    std::size_t stride{0};       // holoborodko sample stride
    double lambda{0};            // whittaker
    std::size_t penaltyOrder{3}; // whittaker difference order d
    double processNoise{0};      // kalman q (continuous); 0 = from cutoff
    double measurementNoise{0};  // kalman r (per-sample variance); 0 = estimated from the data
    double weight{0};            // tv_denoise / tvr / tikhonov_integral regularization; 0 = from cutoff / noise
    double exponent{4};          // weak bump exponent p
    double lipschitz{0};         // levant: bound on |y'''|; 0 = from the data
    std::size_t iterations{0};   // tvr outer iterations (0 = 20); CG iterations come from the size
    std::size_t despike{0};      // offline only: Hampel filter half-window (samples) applied first; 0 = off
    double despikeThreshold{3};  // replace samples farther than this many robust std (MAD) from the window median
};

// Weights w_j with f^(m)(0) ~ sum_j w_j f(x_j) on distinct nodes x_j (Fornberg 1988); for nodes on an
// integer grid of step h divide by h^m.
[[nodiscard]] std::vector<double> differenceWeights(const std::vector<double>& nodes, std::size_t derivativeOrder);
// dy/dt of a uniformly sampled signal by the central_high_order stencils.
[[nodiscard]] std::vector<double> highOrderCentralDifference(const std::vector<double>& y, double dt, std::size_t points = 9);

// Hampel despiking (median / MAD outlier replacement) used by DerivativeSpec::despike.
[[nodiscard]] std::vector<double> despike(const std::vector<double>& y, std::size_t halfWindow, double threshold = 3);

struct DerivativeResult {
    std::vector<double> values;      // smoothed signal (the input itself for pure differentiators)
    std::vector<double> derivative;  // estimated dy/dt
};

[[nodiscard]] bool isCausalDerivativeMethod(const std::string& method);
[[nodiscard]] std::vector<std::string> derivativeMethods();

// Offline differentiation of a uniformly sampled signal.
[[nodiscard]] DerivativeResult differentiate(const std::vector<double>& y, double dt, const DerivativeSpec& spec);

// Noise standard deviation estimated from the data (MAD of 3rd differences, scaled for white noise).
[[nodiscard]] double estimateWhiteNoiseStd(const std::vector<double>& y);

// Online differentiator: one sample at a time on a fixed clock.
class CausalDifferentiator {
public:
    CausalDifferentiator() = default;
    CausalDifferentiator(const DerivativeSpec& spec, double dt, double signalScale = 1.0);
    double step(double y);              // returns the derivative estimate after this sample
    [[nodiscard]] double value() const { return value_; }
    [[nodiscard]] bool configured() const { return !method_.empty(); }
private:
    std::string method_;
    double dt_{1}, value_{0}, derivative_{0};
    std::size_t count_{0};
    std::vector<double> history_;       // newest last
    std::vector<double> coefficients_;  // causal SG weights (oldest first)
    std::vector<double> b_, a_, xs_, ys_;  // causal Butterworth (direct form, as cascaded biquads)
    std::vector<double> x_, P_, F_, Q_;    // Kalman state, covariance, transition, process covariance
    double r_{0};
    std::size_t d_{0};
    double L_{0};                       // levant
    std::vector<double> z_;
};

// Measurement noise. `level` is the noise RMS relative to `scale` (the signal RMS unless given).
struct NoiseSpec {
    std::string kind{"none"};  // none | white | uniform | ar1 | pink | random_walk | quantization | impulse | sum
    double level{0};
    bool relative{true};        // level relative to the signal RMS (else an absolute RMS / step)
    double correlationTime{0};  // ar1: correlation time (time units)
    double probability{0};      // impulse: probability per sample
    double bits{0};             // quantization: bits over the signal range (or absolute step if !relative and level>0)
    std::uint64_t seed{20261115};
    std::vector<NoiseSpec> parts;  // kind "sum"
};

// Adds noise to a whole record. `scale` overrides the signal RMS for relative levels (0 = RMS of y).
[[nodiscard]] std::vector<double> addMeasurementNoise(const std::vector<double>& y, double dt,
                                                      const NoiseSpec& spec, double scale = 0);

// The same noise, one sample at a time (needs an explicit scale and, for quantization, a range).
class StreamingNoise {
public:
    StreamingNoise() = default;
    StreamingNoise(const NoiseSpec& spec, double dt, double scale, double range = 0);
    double apply(double y);
    [[nodiscard]] bool active() const { return !parts_.empty(); }
private:
    struct Part { NoiseSpec spec; double sigma{0}, a{0}, state{0}, step{0}; double pink[7]{}; };
    std::vector<Part> parts_;
    double dt_{1};
    std::uint64_t rng_{0x9E3779B97F4A7C15ull};
    double normal();
    double uniform();
};

}  // namespace dynsyslab::core

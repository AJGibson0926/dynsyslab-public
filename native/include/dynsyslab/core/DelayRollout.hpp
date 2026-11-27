#pragma once

// Rollout of a delay-coordinate model on its own history (method of steps).
//
// The model gives the current derivative from a window of delayed states,
//   x'(t) = F(x(t), x(t - tau_1), ..., x(t - tau_K)),  0 = tau_0 < tau_1 < ... < tau_K,
// and the rollout advances x with classical RK4 at a fixed step h. Delayed values come from
// the supplied uniformly sampled prehistory and, once the trajectory passes it, from the
// rollout's own output; lags that are not multiples of h (and the half-step RK stages) use
// linear interpolation between neighbouring samples. Each lag must be at least one step, so
// every delayed value is already known: the method of steps (Bellen & Zennaro, Numerical
// Methods for Delay Differential Equations, Oxford 2003, Ch. 3; Driver, Ordinary and Delay
// Differential Equations, Springer 1977). Ported from the tools/experiments
// delay_rollout_kernel (mc_delay_rk4_v1).

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

// window: (K + 1) x dimension row-major [present, tau_1 delayed, ..., tau_K delayed].
using DelayRightHandSide = std::function<void(const double* window, double* derivative)>;

struct DelayRolloutOptions {
    double step{};                       // h
    std::vector<double> lags;            // tau_0 = 0, tau_1, ..., strictly increasing
    std::size_t outputSamples{};         // rows written, the first is the last history row
    std::vector<double> magnitudeLimits; // per-component |x| guard (empty: none)
};

enum class DelayRolloutStatus { complete, magnitudeGuard, nonFinite };

struct DelayRollout {
    std::vector<double> values;          // samples x dimension
    std::size_t samples{};
    DelayRolloutStatus status{DelayRolloutStatus::complete};
    std::string reason;
};

// history: historyRows x dimension, uniformly sampled at `step`, ending at the initial state.
// Invalid input throws std::invalid_argument; a guard or non-finite value stops the rollout
// early with the samples computed so far.
[[nodiscard]] DelayRollout rolloutDelayModel(const DelayRightHandSide& rhs, std::size_t dimension,
                                             const std::vector<double>& history,
                                             const DelayRolloutOptions& options);

}  // namespace dynsyslab::core

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace dynsyslab::core {

// Practical PID options. The defaults reproduce a textbook PID: rectangular integral of
// the current error, raw backward-difference derivative, no limits.
struct PIDOptions {
    double derivativeFilterTime{0.0};       // first-order filter time constant on de/dt (0 = raw)
    std::vector<double> integratorLimit;    // per error channel |integral| <= limit (empty = none)
    double outputMinimum{-std::numeric_limits<double>::infinity()};
    double outputMaximum{std::numeric_limits<double>::infinity()};
    // Anti-windup: hold the integral on a step whose output would saturate and whose
    // integral increment pushes further into saturation (conditional integration).
    bool conditionalIntegration{false};
};

// General PID on an error vector e (p channels) producing m outputs:
// u = Kp e + Ki int(e) + Kd de/dt, gains m x p row major. Sampled: call update() at
// monotonic times. Usable on any signal: a state error, a measured output, an estimate.
class PIDCompensator {
public:
    PIDCompensator(std::vector<double> proportional, std::vector<double> integral,
        std::vector<double> derivative, std::size_t errorDimension,
        std::size_t outputDimension, PIDOptions options = {});
    [[nodiscard]] std::vector<double> update(double time, const std::vector<double>& error);
    void reset() noexcept;
    [[nodiscard]] const std::vector<double>& integral() const noexcept { return accumulated_; }
    [[nodiscard]] std::size_t errorDimension() const noexcept { return errorDimension_; }
    [[nodiscard]] std::size_t outputDimension() const noexcept { return outputDimension_; }

private:
    std::vector<double> combine(const std::vector<double>& error,
        const std::vector<double>& integral, const std::vector<double>& rate) const;
    std::vector<double> proportional_,integralGain_,derivative_;
    std::vector<double> accumulated_,previousError_,filteredRate_;
    std::size_t errorDimension_{},outputDimension_{};
    PIDOptions options_;
    double previousTime_{};
    bool initialized_{};
};

}  // namespace dynsyslab::core

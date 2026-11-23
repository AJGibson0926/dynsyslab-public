#pragma once
#include <dynsyslab/core/Observable.hpp>
namespace dynsyslab::core {
// Simulated sensor: y(x,u) = h(x) + gain * u * scale(x).
// This is separate from the observer's measurement model. Empty gain leaves h
// unchanged. Noise, derivatives and nonlinear measurement transforms follow it.
// Standard input/output state-space observation model; see Astrom & Murray,
// Feedback Systems (2008). No observer approximation is implied by this sensor.
class InputAffineSensor {
public:
    InputAffineSensor(std::vector<Observable> base, std::vector<double> gain = {},
        std::optional<Observable> scale = {}, std::size_t scaleComponent = 0);
    std::vector<double> evaluate(const double* state, double appliedInput) const;
    std::size_t outputDimension() const noexcept { return width_; }
private:
    std::vector<Observable> base_;
    std::vector<double> gain_;
    std::optional<Observable> scale_;
    std::size_t component_{}, width_{};
};
}

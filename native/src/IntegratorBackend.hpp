#pragma once

// Internal stepper interface of the native ODE integrator. One driver (Integrator.cpp)
// owns the output clock, the sampled-control clock, guards, checkpoints and failure
// policy; a Stepper only advances a state across one interval [start, end].

#include <dynsyslab/core/Integrator.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace dynsyslab::core::detail {

// Evaluates the controlled right-hand side of one trajectory and checks accepted states.
struct RhsContext {
    ODESystem& system;
    const ControlPolicy& control;
    const std::vector<double>& heldInput;     // sampled mode: the held input vector
    std::size_t trajectory{};
    IntegrationStatistics& statistics;
    const ExecutionContext* execution{};
    const IntegrationOptions& options;
    std::size_t trajectorySteps{};             // accepted + rejected, this trajectory
    std::vector<double> input;                 // scratch, inputDimension

    void rhs(double time, const double* state, double* derivative);
    // Validity domain and magnitude guards; throws std::runtime_error with the reason.
    void accept(double time, const double* state) const;
    // Counts one internal step against the per-trajectory budget.
    void countStep(bool accepted);
};

class Stepper {
public:
    virtual ~Stepper() = default;
    // Start of a trajectory and after every input discontinuity (sampled-control tick).
    virtual void restart(RhsContext& context, double time, const std::vector<double>& state) = 0;
    // Advance state from start to exactly end.
    virtual void advance(RhsContext& context, double start, double end,
                         std::vector<double>& state) = 0;
};

[[nodiscard]] std::unique_ptr<Stepper> makeExplicitStepper(const IntegrationOptions& options,
                                                           std::size_t dimension);
[[nodiscard]] std::unique_ptr<Stepper> makeRadauStepper(const IntegrationOptions& options,
                                                        std::size_t dimension);
// CVODE: bdf, adams and the lsoda switching driver. Throws when SUNDIALS is not built.
[[nodiscard]] std::unique_ptr<Stepper> makeCvodeStepper(const IntegrationOptions& options,
                                                        std::size_t dimension);

// Mixed absolute/relative tolerance scale used by every adaptive stepper.
[[nodiscard]] inline double toleranceScale(const IntegrationOptions& options, double a, double b) {
    const double magnitude = a < 0.0 ? -a : a, other = b < 0.0 ? -b : b;
    return options.absoluteTolerance
        + options.relativeTolerance * (magnitude > other ? magnitude : other);
}

}  // namespace dynsyslab::core::detail

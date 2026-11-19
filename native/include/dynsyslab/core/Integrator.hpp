#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class IntegrationMethod { euler, rk4, dormandPrince45, bdf };

struct IntegrationOptions {
    double initialTime{0.0};
    double finalTime{1.0};
    double outputStep{0.001};
    double initialStep{0.0};
    double maximumStep{0.0};
    double relativeTolerance{1.0e-6};
    double absoluteTolerance{1.0e-9};
    std::size_t maximumInternalSteps{10000000};
    IntegrationMethod method{IntegrationMethod::dormandPrince45};
    // On an integration failure of one trajectory (left the validity domain, step
    // underflow, a controller error) keep the samples computed so far, fill the rest with
    // NaN, record the reason and continue with the next trajectory. Cancellation still
    // propagates. Default: the failure aborts the whole solve.
    bool truncateOnFailure{false};
};

struct IntegrationProgress {
    double fraction{};
    std::size_t trajectory{};
    std::size_t trajectoryCount{};
    std::size_t acceptedSteps{};
    std::size_t rejectedSteps{};
};

struct IntegrationStatistics {
    std::size_t acceptedSteps{};
    std::size_t rejectedSteps{};
    std::size_t rightHandSideCalls{};
    std::size_t controllerCalls{};
};

using ScalarControlLaw = std::function<double(std::size_t trajectory, double time,
                                               const double* state)>;
enum class ControlEvaluationMode { continuousRhs, sampled };

struct ScalarControlPolicy {
    ScalarControlLaw value;
    ControlEvaluationMode evaluation{ControlEvaluationMode::continuousRhs};
    double sampleTime{};

    explicit operator bool() const noexcept { return static_cast<bool>(value); }
};

using ProgressCallback = std::function<void(const IntegrationProgress&)>;

struct NativeSolution {
    std::vector<double> time;
    std::vector<double> values;
    std::vector<double> controlValues;
    std::size_t trajectoryCount{};
    std::size_t sampleCount{};
    std::size_t stateDimension{};
    IntegrationStatistics statistics;
    std::string coordinateSystem;
    // Per trajectory: samples actually computed (sampleCount unless truncated on failure)
    // and the failure that ended it (empty when complete).
    std::vector<std::size_t> completedSamples;
    std::vector<std::string> terminationReasons;

    [[nodiscard]] const double* state(std::size_t trajectory, std::size_t sample) const;
    [[nodiscard]] double control(std::size_t trajectory, std::size_t sample) const;
    [[nodiscard]] TimeSeries trajectorySeries(std::size_t trajectory) const;
    [[nodiscard]] TimeSeries controlSeries(std::size_t trajectory) const;
};

[[nodiscard]] NativeSolution integrate(ODESystem& system,
                                       const std::vector<double>& initialStates,
                                       std::size_t trajectoryCount,
                                       const IntegrationOptions& options,
                                       ScalarControlLaw controller = {},
                                       ProgressCallback progress = {},
                                       const ExecutionContext* execution = nullptr);

[[nodiscard]] NativeSolution integrateControlled(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t trajectoryCount, const IntegrationOptions& options,
    ScalarControlPolicy controller, ProgressCallback progress = {},
    const ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::core

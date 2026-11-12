#pragma once

#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

enum class FixedMethod { euler, rk4 };

struct FixedIntegrationRequest {
    std::string model;
    std::vector<double> parameters;
    std::vector<double> initialStates;
    std::size_t batch{1};
    double initialTime{0.0};
    double finalTime{1.0};
    double stepSize{0.001};
    FixedMethod method{FixedMethod::rk4};
};

struct TrajectoryBatch {
    std::vector<double> time;
    std::vector<double> values;
    std::size_t batch{};
    std::size_t samples{};
    std::size_t stateDimension{};

    [[nodiscard]] const double* state(std::size_t trajectory, std::size_t sample) const;
    [[nodiscard]] TimeSeries trajectorySeries(std::size_t trajectory) const;
};

[[nodiscard]] std::size_t modelStateDimension(const std::string& model);
[[nodiscard]] std::size_t modelParameterCount(const std::string& model);
[[nodiscard]] std::vector<std::string> fixedStepModels();
[[nodiscard]] TrajectoryBatch integrateFixed(const FixedIntegrationRequest& request);

}  // namespace dynsyslab::core

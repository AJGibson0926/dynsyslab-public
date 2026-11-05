#pragma once

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsyslab::workflow::detail {

struct TrainingPreflightResult {
    bool enabled{};
    dynsyslab::core::json::Value report;
};

using TrainingPreflightProgress = std::function<void(
    double fraction, std::size_t completed, std::size_t planned,
    const std::string& message)>;

void validateTrainingPreflightPlan(
    const dynsyslab::core::json::Value& job,
    std::size_t stateDimension,
    const std::string& coordinateSystem);

[[nodiscard]] TrainingPreflightResult executeTrainingPreflight(
    const dynsyslab::core::json::Value& job,
    dynsyslab::core::ODESystem& system,
    const std::vector<double>& initialStates,
    std::size_t initialConditionCount,
    TrainingPreflightProgress progress,
    const dynsyslab::core::ExecutionContext* execution);

}  // namespace dynsyslab::workflow::detail

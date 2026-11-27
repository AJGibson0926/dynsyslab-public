#pragma once
// Classical in-loop estimators (lqe, luenberger, static_map); see ClassicalEstimator.cpp.
#include <dynsyslab/workflow/EstimatorFactory.hpp>

#include <string>

namespace dynsyslab::workflow::detail {

[[nodiscard]] bool isClassicalEstimator(const std::string& kind);
[[nodiscard]] EstimatorSetup buildClassicalEstimator(const json::Value& spec, double sampleTime,
                                                     const core::ODESystem& system, std::size_t trajectoryCount);

}  // namespace dynsyslab::workflow::detail

#pragma once

#include <string>

namespace dsl::compute {

// Applies job-wide numerical execution and campaign settings without changing
// each System's coordinate/time, controller clock, or PDE-decomposition contract.
[[nodiscard]] std::string applySolverExecutionContract(
    const std::string& jobJson, const std::string& solverJson,
    const std::string& campaignJson);

}  // namespace dsl::compute

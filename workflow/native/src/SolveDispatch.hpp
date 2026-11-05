#pragma once

#include <dynsyslab/workflow/Solve.hpp>

namespace dynsyslab::workflow::detail {

[[nodiscard]] SolveResult executeSingleSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution);

}  // namespace dynsyslab::workflow::detail

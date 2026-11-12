#pragma once

#include <dynsyslab/workflow/Run.hpp>

namespace dynsyslab::workflow::detail {

[[nodiscard]] RunResult executeSolveSet(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution);

}  // namespace dynsyslab::workflow::detail

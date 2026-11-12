#pragma once

#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include <filesystem>

namespace dynsyslab::workflow {

[[nodiscard]] bool isLatticeSolve(const json::Value& job);
[[nodiscard]] SolveResult executeLatticeSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include <filesystem>

namespace dynsyslab::workflow::detail {

// True when system.name is a delay, field or lattice model of a loaded extension.
[[nodiscard]] bool isExtensionModelSolve(const json::Value& job);

[[nodiscard]] SolveResult executeExtensionModelSolve(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath, const dynsyslab::core::ExecutionContext* execution);

}  // namespace dynsyslab::workflow::detail

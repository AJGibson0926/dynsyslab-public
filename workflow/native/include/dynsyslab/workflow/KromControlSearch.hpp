#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <filesystem>

namespace dynsyslab::workflow {

[[nodiscard]] json::Value planKromControlSearch(const json::Value& job);

[[nodiscard]] RunResult executeKromControlSearch(
    const json::Value& job,const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath={},
    const dynsyslab::core::ExecutionContext* execution=nullptr);

}  // namespace dynsyslab::workflow

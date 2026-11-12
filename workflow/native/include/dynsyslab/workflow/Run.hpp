#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct RunResult {
    std::string kind;
    std::wstring summary;
    std::vector<std::filesystem::path> archives;
    std::filesystem::path manifest;
};

void validateJob(const json::Value& job);

[[nodiscard]] RunResult executeJob(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

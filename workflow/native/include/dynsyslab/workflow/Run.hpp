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
    // Method that produced each archive, when the job kind records it
    // (postprocess stages); empty otherwise. Used by pipeline references.
    std::vector<std::string> methods;
};

void validateJob(const json::Value& job);

[[nodiscard]] RunResult executeJob(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

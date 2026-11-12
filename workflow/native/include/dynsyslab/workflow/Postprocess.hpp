#pragma once

#include <dynsyslab/workflow/Json.hpp>
#include <dynsyslab/core/Execution.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct PostprocessResult {
    std::vector<std::filesystem::path> archives;
    std::vector<std::string> methods;
    [[nodiscard]] std::wstring summary() const;
};

[[nodiscard]] PostprocessResult executePostprocess(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

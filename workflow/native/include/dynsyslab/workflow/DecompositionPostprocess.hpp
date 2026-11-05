#pragma once

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow {

[[nodiscard]] std::filesystem::path runDecompositionPostprocess(
    const std::string& method, const std::filesystem::path& source,
    const json::Value& settings, const std::filesystem::path& archiveRoot,
    const std::string& label);

}  // namespace dynsyslab::workflow

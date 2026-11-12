#pragma once
#include <dynsyslab/workflow/Json.hpp>
#include <filesystem>
#include <vector>
namespace dynsyslab::workflow {
std::filesystem::path runHistoryBasis(
    const std::vector<std::filesystem::path>& sources,const json::Value& settings,
    const std::filesystem::path& archiveRoot,const std::string& label);
std::filesystem::path runHistoryProjection(
    const std::filesystem::path& source,const json::Value& settings,
    const std::filesystem::path& archiveRoot,const std::string& label);
}

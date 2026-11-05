#pragma once

#include <dynsyslab/workflow/Solve.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct CampaignPoint {
    std::vector<std::string> coordinates;
    SolveResult solve;
};

struct CampaignResult {
    std::string id;
    std::string label;
    std::filesystem::path manifest;
    std::vector<CampaignPoint> points;

    [[nodiscard]] std::wstring summary() const;
};

[[nodiscard]] bool hasEnabledCampaign(const json::Value& job);
[[nodiscard]] std::size_t campaignPointCount(const json::Value& job);
[[nodiscard]] CampaignResult executeCampaign(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

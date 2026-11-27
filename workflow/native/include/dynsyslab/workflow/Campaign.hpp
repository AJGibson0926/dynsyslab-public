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

// Campaign job keys beyond Cartesian `axes` of explicit `values` (docs/architecture/JOB_QUEUE.md):
//   axes[i] = {path, start, stop, points} (linspace) or a `driver.<field>` path (waveform driver);
//   profile low|medium|high|ultra_high (default points and systems_per_job);
//   systems_per_job (queue array chunk size), point_range [first, count] (one chunk);
//   execution_profile fixed|adaptive_local|adaptive_portable (queue.retry_limit);
//   checkpoint_profile fast|recoverable; time.{transient,retained}_cycles, samples_per_cycle,
//   cycle_profile lyapunov_screen|lyapunov_verification.
// Every run registers its plan and chunk manifest in the campaign registry.
[[nodiscard]] json::Value campaignProfiles();  // dynsyslab.campaign-profiles/v1
[[nodiscard]] bool hasEnabledCampaign(const json::Value& job);
[[nodiscard]] std::size_t campaignPointCount(const json::Value& job);
[[nodiscard]] CampaignResult executeCampaign(
    const json::Value& job, const std::filesystem::path& archiveRoot,
    const std::filesystem::path& progressPath = {},
    const dynsyslab::core::ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::workflow

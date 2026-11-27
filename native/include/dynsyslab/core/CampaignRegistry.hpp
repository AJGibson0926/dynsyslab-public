#pragma once
// Campaign registry: the original campaign plan and which queue job and manifest cover each
// chunk of its points. It lives beside the job store (job_manager/state.sqlite) so campaign
// status and validation can join plans, jobs and published manifests.

#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct CampaignPlanRecord {
    std::string id;
    std::string label;
    std::string kind;
    std::size_t pointCount{};
    std::string planJson;  // the original campaign job, as submitted
    std::string createdUtc;
};

struct CampaignChunkRecord {
    std::string campaignId;
    std::size_t firstPoint{};
    std::size_t pointCount{};
    std::string jobId;         // empty when the chunk ran outside the queue
    std::string manifestPath;  // empty until the chunk publishes its manifest
    std::string updatedUtc;
};

class CampaignRegistry {
public:
    explicit CampaignRegistry(const std::filesystem::path& database);
    ~CampaignRegistry();
    CampaignRegistry(const CampaignRegistry&) = delete;
    CampaignRegistry& operator=(const CampaignRegistry&) = delete;

    // Records a plan once; a second registration of the same id keeps the original.
    void registerPlan(const CampaignPlanRecord& plan);
    // Inserts or updates one chunk; empty job or manifest fields keep their stored value.
    void recordChunk(const CampaignChunkRecord& chunk);
    [[nodiscard]] std::optional<CampaignPlanRecord> plan(const std::string& id) const;
    [[nodiscard]] std::vector<CampaignPlanRecord> plans(std::size_t limit = 100) const;
    [[nodiscard]] std::vector<CampaignChunkRecord> chunks(const std::string& id) const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace dynsyslab::core

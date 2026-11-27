#pragma once
// Campaign catalog (dynsyslab-cli archive campaigns): a derived JSON view of every campaign an
// archive root holds, rebuilt from the two independent records of a campaign, its manifests in
// <root>/campaigns (dynsyslab.campaign/v3, written by native campaigns) and the campaign_id in the
// metadata of its solution archives. It reads only and stores nothing (the queue's campaign
// registry, CampaignRegistry.hpp, is the stored record), so it is safe while a campaign runs.
// dynsyslab.campaign-catalog/v1: per campaign its id, label, kind, created_utc, point count,
// manifests, the solution archives found by metadata, and the archives a manifest lists that no
// longer exist.
#include <filesystem>
#include <string>

namespace dynsyslab::core {

inline constexpr const char* kCampaignCatalogSchema = "dynsyslab.campaign-catalog/v1";

// The catalog of `archiveRoot`, as JSON text.
[[nodiscard]] std::string campaignCatalogJson(const std::filesystem::path& archiveRoot);

}  // namespace dynsyslab::core

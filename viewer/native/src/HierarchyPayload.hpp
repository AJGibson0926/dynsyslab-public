#pragma once

#include "dynsyslab_viewer/PlotPresentation.hpp"

#include <filesystem>
#include <string>

namespace dsl {

// Builds lazy aggregate plots for a Batch, Cluster, or parameter-preserving
// Campaign directly from canonical hierarchy and solution archives.
[[nodiscard]] Payload buildHierarchyPayload(
    const CatalogItem& item, const std::string& view,
    const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot,
    const std::optional<PhaseSpaceDomain>& phaseDomain = std::nullopt);

std::vector<std::filesystem::path> hierarchyTrajectorySources(const CatalogItem& item,
    const std::filesystem::path& projectRoot, const std::filesystem::path& archiveRoot);

}  // namespace dsl

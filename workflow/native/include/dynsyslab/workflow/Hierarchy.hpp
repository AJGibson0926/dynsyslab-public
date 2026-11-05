#pragma once

#include <dynsyslab/workflow/Json.hpp>

#include <dynsyslab/core/Hierarchy.hpp>

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

struct HierarchyIdentity {
    std::string batchId;
    std::string systemRunId;
    std::string parameterSetKey;
    std::string dataRole;
    std::string campaignId;
    std::string campaignKind;
    std::vector<std::string> clusterIds;
    std::vector<std::string> trajectoryIds;
    std::vector<std::size_t> trajectoryClusters;
    std::filesystem::path manifestPath;
};

[[nodiscard]] HierarchyIdentity prepareHierarchy(
    const json::Value& job, const dynsyslab::core::SystemConfiguration& configuration,
    std::size_t trajectoryCount, const std::filesystem::path& archiveRoot,
    const std::string& stamp);

void publishHierarchy(
    const HierarchyIdentity& identity, const json::Value& job,
    const dynsyslab::core::SystemConfiguration& requestedConfiguration,
    const dynsyslab::core::SystemDescriptor& descriptor,
    const dynsyslab::core::SystemConfiguration& resolvedConfiguration,
    const std::vector<double>& initialStates,
    const std::vector<std::filesystem::path>& solutionArchives, bool controlled,
    bool hasControlTarget, const std::string& createdUtc);

}  // namespace dynsyslab::workflow

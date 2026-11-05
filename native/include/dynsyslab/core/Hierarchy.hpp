#pragma once

#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct TimeSeriesReference {
    std::string role;
    std::string timeArray{"time"};
    std::string valueArray;
};

struct TrajectoryRecord {
    std::string id;
    std::string label;
    std::vector<double> initialState;
    std::filesystem::path solutionArchive;
    std::vector<TimeSeriesReference> timeSeries;
    std::vector<std::string> tags;
};

struct ClusterRecord {
    std::string id;
    std::string label;
    std::vector<TrajectoryRecord> trajectories;
};

struct SystemRunRecord {
    std::string id;
    std::string label;
    std::string systemIdentifier;
    std::string systemLabel;
    std::string systemFamily{"ode"};
    std::string preset;
    std::string presetLabel;
    std::string parameterSetKey;
    std::string coordinateSystem;
    std::map<std::string, double> parameterOverrides;
    std::map<std::string, double> resolvedParameters;
    std::map<std::string, std::string> parameterUnits;
    std::map<std::string, double> dimensionlessGroups;
    std::map<std::string, std::string> dimensionlessGroupUnits;
    std::string equationsVersion;
    std::string equationsFingerprint;
    std::vector<std::string> stateNames;
    std::vector<ClusterRecord> clusters;
};

struct BatchRecord {
    std::string id;
    std::string label;
    std::string createdUtc;
    std::string dataRole{"free"};
    std::string campaignId;
    std::string campaignKind;
    std::vector<std::string> tags;
    std::map<std::string, std::string> properties;
    std::vector<SystemRunRecord> systems;
};

[[nodiscard]] std::string stableIdentifier(const std::string& value);
[[nodiscard]] std::string parameterSetKey(const SystemConfiguration& configuration);
[[nodiscard]] std::string hierarchyJson(const BatchRecord& batch);
void writeHierarchyAtomic(const std::filesystem::path& path, const BatchRecord& batch);

}  // namespace dynsyslab::core

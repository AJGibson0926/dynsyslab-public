#pragma once

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct LegacyCsvTrajectory {
    std::size_t systemIndex{};
    std::size_t clusterIndex{};
    std::size_t trajectoryIndex{};
    std::vector<double> time;
    std::vector<double> states;
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t firstSourceLine{};
    std::size_t lastSourceLine{};
};

struct LegacyCsvSystem {
    std::size_t index{};
    std::string dynamicsName;
    double declaredStep{};
    std::size_t declaredTransitions{};
    std::vector<std::string> trailingFields;
    std::vector<double> time;
    std::string forcing{"off"};
    std::string subforcing{"off"};
    std::vector<double> controlTime;
    std::vector<double> controlValues;
    std::vector<LegacyCsvTrajectory> trajectories;
    std::size_t sourceLine{};
};

struct LegacyCsvBatch {
    std::filesystem::path source;
    std::string sourceSha256;
    std::vector<std::string> header;
    std::optional<std::size_t> declaredSystemCount;
    std::vector<LegacyCsvSystem> systems;
    std::vector<std::string> warnings;

    [[nodiscard]] std::size_t trajectoryCount() const;
    [[nodiscard]] json::Value summary() const;
};

[[nodiscard]] LegacyCsvBatch readLegacyDynBatchCsv(
    const std::filesystem::path& path, bool strict = true);

struct LegacyCsvImportOptions {
    std::string label{"legacy-dynbatch-import"};
    bool strict{true};
    double inputScale{1.0};
    std::string coordinateSystem{"unknown"};
    SystemConfiguration mappedSystem;
};

struct LegacyCsvImportResult {
    std::filesystem::path source;
    std::string sourceSha256;
    std::vector<std::filesystem::path> archives;
    std::filesystem::path hierarchyManifest;
    std::filesystem::path importManifest;
    std::size_t systemCount{};
    std::size_t trajectoryCount{};

    [[nodiscard]] json::Value summary() const;
};

[[nodiscard]] LegacyCsvImportResult importLegacyDynBatchCsv(
    const std::filesystem::path& source,
    const std::filesystem::path& archiveRoot,
    const LegacyCsvImportOptions& options = {});

}  // namespace dynsyslab::core

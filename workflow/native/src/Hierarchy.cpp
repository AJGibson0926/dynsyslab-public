#include <dynsyslab/workflow/Hierarchy.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <dynsyslab/core/Text.hpp>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow {
namespace {

std::string safeName(std::string value) {
    for (char& item : value)
        if (!std::isalnum(static_cast<unsigned char>(item)) && item != '-' && item != '_')
            item = '-';
    return value.empty() ? "native-solve" : value;
}

std::size_t clusterCount(const json::Value& job, std::size_t trajectories) {
    const double requested = job.at("initial_conditions").at("clusters").number(1.0);
    return std::min(trajectories, static_cast<std::size_t>(std::max(1.0, requested)));
}

std::string role(const json::Value& job, bool controlled) {
    std::string value = job.at("archive").at("data_role").string(
        job.at("data_role").string());
    if (!value.empty()) return value;
    const std::string campaign = job.at("campaign").at("kind").string();
    if (!campaign.empty()) return "campaign";
    return controlled ? "control_attempt" : "free";
}

void atomicHierarchyJson(const std::filesystem::path& path, const json::Value& value) {
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to update hierarchy manifest");
      stream << json::dump(value) << '\n'; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish merged hierarchy manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

}  // namespace

HierarchyIdentity prepareHierarchy(
    const json::Value& job, const dynsyslab::core::SystemConfiguration& configuration,
    std::size_t trajectoryCount, const std::filesystem::path& archiveRoot,
    const std::string& stamp) {
    if (trajectoryCount == 0) throw std::invalid_argument("hierarchy requires trajectories");
    HierarchyIdentity result;
    const std::string label = job.at("label").string("native-solve");
    result.parameterSetKey = dynsyslab::core::parameterSetKey(configuration);
    result.batchId = job.at("batch").at("id").string(
        dynsyslab::core::stableIdentifier(label + "|" + stamp + "|batch"));
    result.systemRunId = job.at("batch").at("system_run_id").string(
        dynsyslab::core::stableIdentifier(
            result.batchId + "|" + result.parameterSetKey + "|system"));
    const double campaignPoint = job.at("campaign").at("point_index").number(-1.0);
    if (campaignPoint >= 0.0) result.systemRunId = dynsyslab::core::stableIdentifier(
        result.systemRunId + "|point|" + std::to_string(static_cast<std::size_t>(campaignPoint)));
    result.dataRole = role(job, job.at("controller").at("kind").string("none") != "none");
    result.campaignId = job.at("campaign").at("id").string();
    result.campaignKind = job.at("campaign").at("kind").string();
    const std::size_t clusters = clusterCount(job, trajectoryCount);
    result.clusterIds.reserve(clusters);
    for (std::size_t cluster = 0; cluster < clusters; ++cluster)
        result.clusterIds.push_back(dynsyslab::core::stableIdentifier(
            result.systemRunId + "|cluster|" + std::to_string(cluster)));
    result.trajectoryIds.reserve(trajectoryCount);
    result.trajectoryClusters.reserve(trajectoryCount);
    for (std::size_t trajectory = 0; trajectory < trajectoryCount; ++trajectory) {
        const std::size_t cluster = std::min(clusters - 1, trajectory * clusters / trajectoryCount);
        result.trajectoryClusters.push_back(cluster);
        result.trajectoryIds.push_back(dynsyslab::core::stableIdentifier(
            result.clusterIds[cluster] + "|trajectory|" + std::to_string(trajectory)));
    }
    const std::string requestedManifest = job.at("batch").at("manifest_path").string();
    result.manifestPath = requestedManifest.empty()
        ? archiveRoot / L"hierarchies" / dynsyslab::core::wide(configuration.identifier)
            / dynsyslab::core::wide(safeName(label) + "__" + stamp + ".hierarchy.json")
        : std::filesystem::u8path(requestedManifest);
    return result;
}

void publishHierarchy(
    const HierarchyIdentity& identity, const json::Value& job,
    const dynsyslab::core::SystemConfiguration& requestedConfiguration,
    const dynsyslab::core::SystemDescriptor& descriptor,
    const dynsyslab::core::SystemConfiguration& resolvedConfiguration,
    const std::vector<double>& initialStates,
    const std::vector<std::filesystem::path>& solutionArchives, bool controlled,
    bool hasControlTarget, const std::string& createdUtc) {
    const std::size_t count = solutionArchives.size();
    const std::size_t dimension = descriptor.stateDimension;
    if (identity.trajectoryIds.size() != count || initialStates.size() != count * dimension)
        throw std::invalid_argument("hierarchy solve dimensions are inconsistent");
    dynsyslab::core::BatchRecord batch;
    batch.id = identity.batchId;
    batch.label = job.at("batch").at("label").string(
        job.at("label").string("native-solve"));
    batch.createdUtc = createdUtc;
    batch.dataRole = identity.dataRole;
    batch.campaignId = identity.campaignId;
    batch.campaignKind = identity.campaignKind;
    batch.tags = {controlled ? "forced" : "free", identity.dataRole};
    batch.properties.emplace("execution_backend", "native-cpp");
    batch.properties.emplace("python_required", "false");
    if (!identity.campaignId.empty()) batch.properties.emplace("campaign_id", identity.campaignId);
    if (job.at("preflight").at("enabled").boolean(false)) {
        batch.properties.emplace("training_preflight", "passed");
        batch.properties.emplace("training_source_krom",
            job.at("preflight").at("source_krom").string());
    }

    dynsyslab::core::SystemRunRecord run;
    run.id = identity.systemRunId;
    run.label = job.at("ui_system_label").string(batch.label);
    run.systemIdentifier = descriptor.identifier;
    run.systemLabel = descriptor.label;
    run.systemFamily = descriptor.family;
    run.preset = resolvedConfiguration.preset;
    for (const auto& preset : descriptor.presets)
        if (preset.key == resolvedConfiguration.preset) run.presetLabel = preset.description;
    if (run.presetLabel.empty()) run.presetLabel = resolvedConfiguration.preset;
    run.parameterSetKey = identity.parameterSetKey;
    run.coordinateSystem = descriptor.preferredCoordinates;
    run.parameterOverrides = requestedConfiguration.parameters;
    run.resolvedParameters = resolvedConfiguration.parameters;
    run.parameterUnits = resolvedConfiguration.parameterUnits;
    run.dimensionlessGroups = resolvedConfiguration.dimensionlessGroups;
    run.dimensionlessGroupUnits = resolvedConfiguration.dimensionlessGroupUnits;
    run.equationsVersion = resolvedConfiguration.equationsVersion;
    run.equationsFingerprint = dynsyslab::core::equationsFingerprint(resolvedConfiguration);
    run.stateNames = descriptor.stateNames;
    run.clusters.resize(identity.clusterIds.size());
    for (std::size_t cluster = 0; cluster < run.clusters.size(); ++cluster) {
        run.clusters[cluster].id = identity.clusterIds[cluster];
        run.clusters[cluster].label = "cluster-" + std::to_string(cluster);
    }
    for (std::size_t trajectory = 0; trajectory < count; ++trajectory) {
        dynsyslab::core::TrajectoryRecord record;
        record.id = identity.trajectoryIds[trajectory];
        record.label = count == 1 ? batch.label :
            batch.label + "-trajectory-" + std::to_string(trajectory);
        record.initialState.assign(initialStates.begin() + trajectory * dimension,
                                   initialStates.begin() + (trajectory + 1) * dimension);
        record.solutionArchive = solutionArchives[trajectory];
        record.tags = {controlled ? "forced" : "free", identity.dataRole};
        record.timeSeries.push_back({"state", "time", "values"});
        if (controlled) {
            record.timeSeries.push_back({"control", "time", "control_values"});
        }
        if (hasControlTarget) {
            record.timeSeries.push_back({"control_target", "time", "control_target_values"});
        }
        run.clusters[identity.trajectoryClusters[trajectory]].trajectories.push_back(
            std::move(record));
    }
    batch.systems.push_back(std::move(run));
    if (job.at("batch").at("merge").boolean(false)
        && std::filesystem::is_regular_file(identity.manifestPath)) {
        const auto incoming = identity.manifestPath.parent_path()
            / (identity.manifestPath.filename().wstring() + L".incoming");
        dynsyslab::core::writeHierarchyAtomic(incoming, batch);
        auto existing = json::load(identity.manifestPath.wstring());
        const auto next = json::load(incoming.wstring());
        if (existing.at("format").string() != "dynsyslab.hierarchy/v2"
            || existing.at("id").string() != identity.batchId)
            throw std::runtime_error("cannot merge incompatible hierarchy manifests");
        for (const auto& system : next.at("systems").array())
            existing.at("systems").array().push_back(system);
        atomicHierarchyJson(identity.manifestPath, existing);
        std::error_code ignored; std::filesystem::remove(incoming, ignored);
    } else dynsyslab::core::writeHierarchyAtomic(identity.manifestPath, batch);
}

}  // namespace dynsyslab::workflow

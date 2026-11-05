#include "HierarchyCatalog.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <iomanip>
#include <map>
#include <sstream>
#include <system_error>
#include <utility>

namespace dsl {
namespace {

std::string systemLabel(const std::string& identifier) {
    for (const auto& descriptor : dynsyslab::core::availableSystems())
        if (descriptor.identifier == identifier) return descriptor.label;
    if (identifier == "rayleigh-plesset") return "Rayleigh-Plesset";
    return identifier;
}

std::string numberedId(const std::string& parent, const char* level, std::size_t index) {
    std::ostringstream out;
    out << parent << '-' << level << '-' << std::setw(4) << std::setfill('0') << index;
    return out.str();
}

std::string coordinateSystem(const json::Value& system) {
    for (const auto& cluster : system.at("clusters").array())
        for (const auto& trajectory : cluster.at("trajectories").array()) {
            const auto value = trajectory.at("metadata").at("solver_metadata")
                .at("integration_form").string();
            if (!value.empty()) return value;
        }
    return {};
}

std::size_t stateDimension(const json::Value& system) {
    for (const auto& cluster : system.at("clusters").array())
        for (const auto& trajectory : cluster.at("trajectories").array())
            if (!trajectory.at("initial_condition").array().empty())
                return trajectory.at("initial_condition").array().size();
    return 0;
}

std::filesystem::path indexedSourcePath(
    const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot,
    const std::filesystem::path& manifest,
    const std::string& source,
    const std::set<std::string>* indexedSources) {
    if (!indexedSources)
        return resolveHierarchySource(projectRoot, archiveRoot, manifest, source);
    const auto reference = std::filesystem::u8path(source);
    if (reference.is_absolute()) return reference.lexically_normal();
    std::vector<std::filesystem::path> candidates{manifest.parent_path() / reference,
        projectRoot / reference, archiveRoot / reference};
    bool copy = false;
    std::filesystem::path archiveRelative;
    for (const auto& part : reference) {
        const auto name = part.generic_u8string();
        if (!copy && (name == "solutions" || name == "derived")) copy = true;
        if (copy) archiveRelative /= part;
    }
    if (!archiveRelative.empty()) candidates.push_back(archiveRoot / archiveRelative);
    for (const auto& candidate : candidates)
        if (indexedSources->count(normalizedSource(candidate)))
            return candidate.lexically_normal();
    return candidates.front().lexically_normal();
}

bool completed(const json::Value& trajectory, const std::filesystem::path& projectRoot,
               const std::filesystem::path& archiveRoot,
               const std::filesystem::path& manifest,
               const std::set<std::string>* indexedSources) {
    if (trajectory.at("metadata").at("status").string() != "complete") return false;
    const auto source = trajectory.at("solution_archive").string();
    if (source.empty()) return false;
    const auto path = indexedSourcePath(projectRoot, archiveRoot, manifest, source,
                                        indexedSources);
    return indexedSources ? indexedSources->count(normalizedSource(path)) != 0
                          : std::filesystem::is_regular_file(path);
}

std::size_t completedCount(const json::Value& system,
                           const std::filesystem::path& projectRoot,
                           const std::filesystem::path& archiveRoot,
                           const std::filesystem::path& manifest,
                           const std::set<std::string>* indexedSources) {
    std::size_t count = 0;
    for (const auto& cluster : system.at("clusters").array())
        for (const auto& trajectory : cluster.at("trajectories").array())
            if (completed(trajectory, projectRoot, archiveRoot, manifest,
                          indexedSources)) ++count;
    return count;
}

void aggregateViews(CatalogItem& item, std::size_t dimension, bool forced) {
    item.views = {"ensemble_time"};
    if (dimension >= 2) item.views.insert(item.views.begin(), "ensemble_phase");
    if (forced) item.views.push_back("ensemble_control");
}

void appendConfiguration(CatalogItem& item, const json::Value& system) {
    const auto& configuration = system.at("dynamical_system");
    const auto identifier = configuration.at("name").string();
    item.properties["system id"] = identifier;
    item.properties["coordinate system"] = coordinateSystem(system);
    item.properties["parameter set key"] = dynsyslab::core::stableIdentifier(
        identifier + json::dump(configuration.at("preset")) +
        json::dump(configuration.at("ambient_pressure_pa")) +
        json::dump(configuration.at("equilibrium_radius_m")));
    for (const auto& value : configuration.at("preset").object()) {
        const auto text = value.second.string();
        item.properties["parameter · " + value.first] = text.empty()
            ? json::dump(value.second) : text;
    }
    for (const auto& value : configuration.at("driver").object())
        item.properties["control · " + value.first] = json::dump(value.second);
}

}  // namespace

void appendLegacyHierarchyCatalog(Catalog& catalog,
                                  const std::filesystem::path& projectRoot,
                                  const std::filesystem::path& archiveRoot,
                                  std::set<std::string>& ownedSources,
                                  const std::vector<dynsyslab::core::ArchiveIndexEntry>*
                                      indexEntries) {
    const auto hierarchyRoot = archiveRoot / L"hierarchies";
    std::vector<std::pair<std::filesystem::path, json::Value>> manifests;
    std::set<std::string> indexedSources;
    std::map<std::string, const dynsyslab::core::ArchiveIndexEntry*> indexedEntries;
    if (indexEntries) {
        for (const auto& entry : *indexEntries) {
            indexedSources.insert(normalizedSource(entry.path));
            indexedEntries.emplace(normalizedSource(entry.path), &entry);
            if (entry.path.parent_path().lexically_normal() != hierarchyRoot.lexically_normal()
                || entry.path.extension() != L".npz") continue;
            try { manifests.emplace_back(entry.path, json::parse(entry.summaryJson)); }
            catch (...) { /* A malformed summary remains a generic indexed archive. */ }
        }
    } else {
        if (!std::filesystem::is_directory(hierarchyRoot)) return;
        for (const auto& entry : std::filesystem::directory_iterator(hierarchyRoot)) {
            if (!entry.is_regular_file() || entry.path().extension() != L".npz") continue;
            try {
                const auto index = dynsyslab::core::readNpzIndex(entry.path());
                if (!index.metadataJson.empty())
                    manifests.emplace_back(entry.path(), json::parse(index.metadataJson));
            } catch (...) { /* A malformed archive remains available to generic scan. */ }
        }
    }
    const auto* available = indexEntries ? &indexedSources : nullptr;
    for (const auto& manifest : manifests) {
        const auto& entryPath = manifest.first;
        const auto& root = manifest.second;
        if (root.at("format").string() != "dynsyslab.hierarchy/v1") continue;
        ownedSources.insert(normalizedSource(entryPath));

        std::size_t systemCount = 0, clusterCount = 0, trajectoryCount = 0;
        bool forced = false; std::size_t maximumDimension = 0;
        for (const auto& system : root.at("systems").array()) {
            const auto count = completedCount(system, projectRoot, archiveRoot,
                                              entryPath, available);
            if (!count) continue;
            ++systemCount; trajectoryCount += count;
            maximumDimension = std::max(maximumDimension, stateDimension(system));
            forced = forced || !system.at("dynamical_system").at("driver").isNull();
            for (const auto& cluster : system.at("clusters").array())
                if (std::any_of(cluster.at("trajectories").array().begin(),
                    cluster.at("trajectories").array().end(), [&](const auto& trajectory) {
                        return completed(trajectory, projectRoot, archiveRoot, entryPath,
                                         available);
                    })) ++clusterCount;
        }
        if (!trajectoryCount) continue;
        const auto& job = root.at("job");
        const auto jobId = job.at("job_id").string(
            dynsyslab::core::stableIdentifier(entryPath.generic_u8string()));
        const std::string batchId = "legacy-batch-" + jobId;
        const auto campaignId = job.at("tags").at("campaign_id").string();
        const auto campaignKind = job.at("tags").at("campaign_kind").string();
        const std::string role = forced ? "forced" : "free";
        const std::string manifestSource = entryPath.generic_u8string();

        CatalogItem batch;
        batch.id = batchId; batch.kind = "batch"; batch.category = "Archived solves";
        batch.label = job.at("label").string(entryPath.stem().u8string());
        batch.subtitle = "Completed archive batch (legacy hierarchy)";
        batch.source = manifestSource; batch.lineageId = batchId;
        batch.role = role; batch.tags = {role, "legacy hierarchy"};
        if (!campaignId.empty()) batch.tags.push_back("campaign");
        if (!campaignKind.empty()) batch.tags.push_back(campaignKind);
        batch.properties = {{"archive schema", "dynsyslab.hierarchy/v1"},
            {"system count", std::to_string(systemCount)},
            {"cluster count", std::to_string(clusterCount)},
            {"trajectory count", std::to_string(trajectoryCount)},
            {"campaign id", campaignId}, {"campaign kind", campaignKind}};
        aggregateViews(batch, maximumDimension, forced);
        if (!campaignId.empty()) {
            batch.views.insert(batch.views.begin(), "campaign_summary");
            batch.views.insert(batch.views.begin() + 1, "campaign_bifurcation");
        }
        catalog.items.push_back(std::move(batch));

        std::size_t systemIndex = 0;
        for (const auto& system : root.at("systems").array()) {
            const auto systemTrajectories = completedCount(
                system, projectRoot, archiveRoot, entryPath, available);
            const auto currentSystem = systemIndex++;
            if (!systemTrajectories) continue;
            const std::string systemId = numberedId(batchId, "system", currentSystem);
            const auto identifier = system.at("dynamical_system").at("name").string();
            CatalogItem systemItem;
            systemItem.id = systemId; systemItem.kind = "system_run";
            systemItem.category = "Archived solves";
            systemItem.label = system.at("label").string(systemLabel(identifier));
            systemItem.system = systemLabel(identifier);
            systemItem.subtitle = "Completed system run"; systemItem.source = manifestSource;
            systemItem.parentId = batchId; systemItem.lineageId = batchId;
            systemItem.role = role; systemItem.tags = {role, coordinateSystem(system)};
            systemItem.properties["trajectory count"] = std::to_string(systemTrajectories);
            systemItem.properties["legacy system index"] = std::to_string(currentSystem);
            appendConfiguration(systemItem, system);
            aggregateViews(systemItem, stateDimension(system), forced);
            catalog.items.push_back(std::move(systemItem));

            std::size_t clusterIndex = 0;
            for (const auto& cluster : system.at("clusters").array()) {
                const auto currentCluster = clusterIndex++;
                std::vector<const json::Value*> archived;
                for (const auto& trajectory : cluster.at("trajectories").array())
                    if (completed(trajectory, projectRoot, archiveRoot, entryPath,
                                  available))
                        archived.push_back(&trajectory);
                if (archived.empty()) continue;
                const std::string clusterId = numberedId(systemId, "cluster", currentCluster);
                CatalogItem clusterItem;
                clusterItem.id = clusterId; clusterItem.kind = "cluster";
                clusterItem.category = "Archived solves";
                clusterItem.label = cluster.at("label").string("Cluster " +
                    std::to_string(currentCluster + 1));
                clusterItem.system = systemLabel(identifier);
                clusterItem.subtitle = "Archived trajectory cluster";
                clusterItem.source = manifestSource; clusterItem.parentId = systemId;
                clusterItem.lineageId = batchId; clusterItem.role = role;
                clusterItem.tags = {role};
                clusterItem.properties = {{"trajectory count", std::to_string(archived.size())},
                    {"legacy system index", std::to_string(currentSystem)},
                    {"legacy cluster index", std::to_string(currentCluster)}};
                aggregateViews(clusterItem, stateDimension(system), forced);
                catalog.items.push_back(std::move(clusterItem));

                for (std::size_t trajectoryIndex = 0; trajectoryIndex < archived.size(); ++trajectoryIndex) {
                    const auto& trajectory = *archived[trajectoryIndex];
                    const auto source = indexedSourcePath(projectRoot, archiveRoot,
                        entryPath, trajectory.at("solution_archive").string(), available);
                    CatalogItem item;
                    item.id = numberedId(clusterId, "trajectory", trajectoryIndex);
                    item.kind = "trajectory"; item.category = "Archived solves";
                    item.label = trajectory.at("label").string("Trajectory " +
                        std::to_string(trajectoryIndex + 1));
                    item.system = systemLabel(identifier);
                    item.subtitle = "Archived dynamical trajectory";
                    item.source = source.generic_u8string(); item.parentId = clusterId;
                    item.lineageId = batchId; item.role = role; item.tags = {role};
                    item.views = {"data_overview", "time"};
                    if (stateDimension(system) >= 2) item.views.push_back("phase");
                    if (stateDimension(system) >= 3) item.views.push_back("three_d");
                    item.views.push_back("fourier");
                    item.views.push_back("time_frequency");
                    if (forced) {
                        item.views.push_back("control");
                        item.views.push_back("control_fourier");
                    }
                    item.properties["IC"] = json::dump(
                        trajectory.at("initial_condition"));
                    const auto indexed = indexedEntries.find(normalizedSource(source));
                    if (indexed != indexedEntries.end()) {
                        item.bytes = static_cast<double>(indexed->second->sizeBytes);
                        item.modified = static_cast<double>(indexed->second->modifiedTicks);
                    } else {
                        std::error_code fileError;
                        item.bytes = static_cast<double>(std::filesystem::file_size(source, fileError));
                        if (fileError) item.bytes = 0.0;
                        fileError.clear();
                        const auto modified = std::filesystem::last_write_time(source, fileError);
                        item.modified = fileError ? 0.0 : static_cast<double>(
                            modified.time_since_epoch().count());
                    }
                    ownedSources.insert(normalizedSource(source));
                    catalog.items.push_back(std::move(item));
                }
            }
        }
    }
}

}  // namespace dsl

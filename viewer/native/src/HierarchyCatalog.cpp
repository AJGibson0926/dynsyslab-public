#include "HierarchyCatalog.hpp"
#include "ArchiveCatalog.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <map>
#include <sstream>
#include <system_error>

namespace dsl {
namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::vector<std::string> strings(const json::Value& value) {
    std::vector<std::string> result;
    for (const auto& item : value.array()) result.push_back(item.string());
    return result;
}

bool contains(const std::vector<std::string>& values, const std::string& value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

std::string joined(const std::vector<double>& values) {
    std::string result;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) result += ", ";
        result += std::to_string(values[index]);
    }
    return result;
}

std::string measured(double value, const std::string& unit) {
    std::ostringstream out;
    out << std::setprecision(12) << value;
    if (!unit.empty()) out << ' ' << unit;
    return out.str();
}

void appendSystemProperties(CatalogItem& item, const json::Value& system) {
    item.properties.emplace("system id", system.at("system").string());
    item.properties.emplace("system family", system.at("system_family").string("ode"));
    item.properties.emplace("preset", system.at("preset").string());
    const auto presetLabel = system.at("preset_label").string();
    if (!presetLabel.empty()) item.properties.emplace("physics preset", presetLabel);
    item.properties.emplace("parameter set key", system.at("parameter_set_key").string());
    item.properties.emplace("coordinate system", system.at("coordinate_system").string());
    item.properties.emplace("equations version", system.at("equations_version").string());
    item.properties.emplace("equations fingerprint", system.at("equations_fingerprint").string());
    for (const auto& parameter : system.at("resolved_parameters").object())
        item.properties.emplace("parameter · " + parameter.first, measured(
            parameter.second.number(), system.at("parameter_units").at(parameter.first).string()));
    for (const auto& group : system.at("dimensionless_groups").object())
        item.properties.emplace("group · " + group.first, measured(
            group.second.number(), system.at("dimensionless_group_units").at(group.first).string()));
    for (const auto& parameter : system.at("parameter_overrides").object())
        item.properties.emplace("override · " + parameter.first,
            std::to_string(parameter.second.number()));
}

CatalogItem systemItem(const json::Value& system, const json::Value& root,
                       const std::filesystem::path& manifest,
                       const std::string& batchId) {
    CatalogItem item;
    item.id = system.at("id").string(batchId + "-system-" +
        system.at("parameter_set_key").string());
    item.kind = "system_run";
    item.category = "Archived solves";
    item.label = system.at("label").string(system.at("system_label").string());
    item.system = system.at("system_label").string(system.at("system").string());
    item.subtitle = "Completed system run and exact physical parameter set";
    item.source = manifest.generic_u8string();
    item.parentId = batchId;
    item.role = root.at("data_role").string("free");
    item.lineageId = batchId;
    item.tags = {system.at("coordinate_system").string(), system.at("preset").string()};
    for (const auto& tag : strings(root.at("tags")))
        if (!contains(item.tags, tag)) item.tags.push_back(tag);
    item.properties.emplace("cluster count", std::to_string(system.at("clusters").array().size()));
    appendSystemProperties(item, system);
    item.views = {"ensemble_time"};
    if (system.at("state_names").array().size() >= 2)
        item.views.insert(item.views.begin(), "ensemble_phase");
    if (contains(item.tags, "forced")) item.views.push_back("ensemble_control");
    return item;
}

CatalogItem batchItem(const json::Value& root, const std::filesystem::path& manifest) {
    CatalogItem item;
    item.id = root.at("id").string();
    item.kind = "batch";
    item.category = "Archived solves";
    item.label = root.at("label").string("batch");
    item.subtitle = "Completed archive batch";
    item.source = manifest.generic_u8string();
    item.lineageId = item.id;
    item.role = root.at("data_role").string("free");
    item.tags = strings(root.at("tags"));
    if (item.tags.empty()) item.tags = {item.role};
    const std::string campaign = root.at("campaign_kind").string();
    if (!campaign.empty()) item.tags.push_back(campaign);
    item.properties.emplace("data role", item.role);
    item.properties.emplace("created UTC", root.at("created_utc").string());
    item.properties.emplace("campaign id", root.at("campaign_id").string());
    item.properties.emplace("campaign kind", campaign);
    std::size_t clusters = 0, trajectories = 0;
    for (const auto& system : root.at("systems").array()) {
        clusters += system.at("clusters").array().size();
        for (const auto& cluster : system.at("clusters").array())
            trajectories += cluster.at("trajectories").array().size();
    }
    item.properties.emplace("system count", std::to_string(root.at("systems").array().size()));
    item.properties.emplace("cluster count", std::to_string(clusters));
    item.properties.emplace("trajectory count", std::to_string(trajectories));
    for (const auto& property : root.at("properties").object())
        item.properties.emplace(property.first, property.second.string());
    item.views = {"ensemble_time"};
    bool phase = false;
    for (const auto& system : root.at("systems").array())
        phase = phase || system.at("state_names").array().size() >= 2;
    if (phase) item.views.insert(item.views.begin(), "ensemble_phase");
    if (contains(item.tags, "forced")) item.views.push_back("ensemble_control");
    if (!root.at("campaign_id").string().empty()) {
        item.views.insert(item.views.begin(), "campaign_summary");
        item.views.insert(item.views.begin() + 1, "campaign_bifurcation");
    }
    return item;
}

// Campaign type label for the catalog: bifurcation, lyapunov, trajectory ensemble, parameter
// sweep, or generic, from the manifest's kind, label and goal.
std::string campaignType(const json::Value& campaign) {
    std::string text;
    for (const char* key : {"kind", "label", "goal", "campaign_type"}) text += " " + campaign.at(key).string();
    for (auto& character : text) character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    if (text.find("bifurcation") != std::string::npos) return "bifurcation";
    if (text.find("lyapunov") != std::string::npos) return "lyapunov";
    if (text.find("ensemble") != std::string::npos || text.find("unforced") != std::string::npos)
        return "trajectory ensemble";
    if (text.find("sweep") != std::string::npos) return "parameter sweep";
    return "generic";
}

}  // namespace

std::filesystem::path resolveHierarchySource(
    const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot,
    const std::filesystem::path& manifest,
    const std::string& source) {
    const auto reference = std::filesystem::u8path(source);
    if (reference.is_absolute()) return std::filesystem::weakly_canonical(reference);
    std::vector<std::filesystem::path> candidates{
        manifest.parent_path() / reference,
        projectRoot / reference,
        archiveRoot / reference};
    bool copy = false;
    std::filesystem::path archiveRelative;
    for (const auto& part : reference) {
        const auto name = part.generic_u8string();
        if (!copy && (name == "solutions" || name == "derived" || name == "objects")) copy = true;
        if (copy) archiveRelative /= part;
    }
    if (!archiveRelative.empty()) candidates.push_back(archiveRoot / archiveRelative);
    for (const auto& candidate : candidates)
        if (std::filesystem::is_regular_file(candidate))
            return std::filesystem::weakly_canonical(candidate);
    return std::filesystem::weakly_canonical(candidates.front());
}

std::string normalizedSource(const std::filesystem::path& path) {
    std::error_code ignored;
    auto absolute = std::filesystem::absolute(path, ignored);
    if (ignored) absolute = path;
    auto normalized = absolute.lexically_normal().generic_u8string();
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    return normalized;
}

void appendHierarchyCatalog(Catalog& catalog, const std::filesystem::path& projectRoot,
                            const std::filesystem::path& archiveRoot,
                            std::set<std::string>& ownedSources,
                            const std::vector<dynsyslab::core::ArchiveIndexEntry>*
                                indexEntries) {
    std::map<std::string, const dynsyslab::core::ArchiveIndexEntry*> indexedSources;
    if (indexEntries)
        for (const auto& indexed : *indexEntries)
            indexedSources.emplace(normalizedSource(indexed.path), &indexed);
    const auto hierarchyRoot = archiveRoot / L"hierarchies";
    if (!std::filesystem::is_directory(hierarchyRoot)) return;
    std::set<std::string> representedCampaigns;
    std::set<std::string> representedBatches;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(
             hierarchyRoot, std::filesystem::directory_options::skip_permission_denied)) {
        if (!entry.is_regular_file() || entry.path().extension() != L".json") continue;
        json::Value root;
        try { root = json::load(entry.path().wstring()); }
        catch (...) { continue; }
        if (root.at("schema").string() != "dynsyslab.hierarchy/v2") continue;
        ownedSources.insert(normalizedSource(entry.path()));
        const std::string batchId = root.at("id").string();
        if (batchId.empty() || !representedBatches.insert(batchId).second) continue;
        catalog.items.push_back(batchItem(root, entry.path()));
        const std::string campaignId = root.at("campaign_id").string();
        if (!campaignId.empty()) representedCampaigns.insert(campaignId);
        for (const auto& system : root.at("systems").array()) {
            CatalogItem archivedSystem = systemItem(system, root, entry.path(), batchId);
            const std::string systemId = archivedSystem.id;
            catalog.items.push_back(std::move(archivedSystem));
            for (const auto& cluster : system.at("clusters").array()) {
                CatalogItem clusterItem;
                clusterItem.id = cluster.at("id").string();
                clusterItem.kind = "cluster";
                clusterItem.category = "Archived solves";
                clusterItem.label = cluster.at("label").string("cluster");
                clusterItem.system = system.at("system_label").string();
                clusterItem.subtitle = "IC trajectory group";
                clusterItem.source = entry.path().generic_u8string();
                clusterItem.parentId = systemId;
                clusterItem.lineageId = batchId;
                clusterItem.role = root.at("data_role").string();
                clusterItem.tags = strings(root.at("tags"));
                clusterItem.properties.emplace("trajectory count",
                    std::to_string(cluster.at("trajectories").array().size()));
                clusterItem.properties.emplace("parameter set key",
                    system.at("parameter_set_key").string());
                clusterItem.properties.emplace("system run id", systemId);
                clusterItem.views = {"ensemble_time"};
                if (system.at("state_names").array().size() >= 2)
                    clusterItem.views.insert(clusterItem.views.begin(), "ensemble_phase");
                if (contains(clusterItem.tags, "forced"))
                    clusterItem.views.push_back("ensemble_control");
                const std::string clusterId = clusterItem.id;
                catalog.items.push_back(std::move(clusterItem));
                for (const auto& trajectory : cluster.at("trajectories").array()) {
                    CatalogItem item;
                    item.id = trajectory.at("id").string();
                    item.kind = "trajectory";
                    item.category = "Archived solves";
                    item.label = trajectory.at("label").string("trajectory");
                    item.system = system.at("system_label").string();
                    item.subtitle = "Archived dynamical trajectory";
                    const auto source = resolveHierarchySource(projectRoot, archiveRoot, entry.path(),
                        trajectory.at("solution_archive").string());
                    item.source = source.generic_u8string();
                    ownedSources.insert(normalizedSource(source));
                    item.parentId = clusterId;
                    item.lineageId = batchId;
                    item.role = root.at("data_role").string();
                    item.tags = strings(trajectory.at("tags"));
                    if (item.tags.empty()) item.tags = {item.role};
                    item.views = {"data_overview", "time"};
                    const std::string family = system.at("system_family").string("ode");
                    if (family == "pde") item.views.push_back("field");
                    else if (system.at("state_names").array().size() >= 2)
                        item.views.push_back("phase");
                    if (family == "ode" && system.at("state_names").array().size() >= 3)
                        item.views.push_back("three_d");
                    if (family == "ode") {
                        item.views.push_back("fourier");
                        item.views.push_back("time_frequency");
                    }
                    if (contains(item.tags, "forced")) {
                        item.views.push_back("control");
                        item.views.push_back("control_fourier");
                    }
                    item.properties.emplace("IC",
                        joined(numbers(trajectory.at("initial_state"))));
                    item.properties.emplace("parameter set key",
                        system.at("parameter_set_key").string());
                    const auto indexed = indexedSources.find(normalizedSource(source));
                    if (indexed != indexedSources.end()) {
                        item.bytes = static_cast<double>(indexed->second->sizeBytes);
                        item.modified = static_cast<double>(indexed->second->modifiedTicks);
                        enrichArchiveCatalogItem(item, {indexed->second->arrays,
                            indexed->second->summaryJson});
                    } else {
                        std::error_code fileError;
                        item.bytes = static_cast<double>(std::filesystem::file_size(source, fileError));
                        if (fileError) item.bytes = 0.0;
                        fileError.clear();
                        const auto modified = std::filesystem::last_write_time(source, fileError);
                        item.modified = fileError ? 0.0 : static_cast<double>(
                            modified.time_since_epoch().count());
                        enrichArchiveCatalogItem(item, source);
                    }
                    catalog.items.push_back(std::move(item));
                }
            }
        }
    }
    const auto campaignRoot = archiveRoot / L"campaigns";
    if (!std::filesystem::is_directory(campaignRoot)) return;
    for (const auto& entry : std::filesystem::directory_iterator(campaignRoot)) {
        if (!entry.is_regular_file() || entry.path().extension() != L".json") continue;
        try {
            const auto campaign = json::load(entry.path().wstring());
            if (campaign.at("schema").string() == "dynsyslab.campaign/v3" &&
                representedCampaigns.count(campaign.at("id").string())) {
                ownedSources.insert(normalizedSource(entry.path()));
                const std::string campaignId = campaign.at("id").string();
                for (auto& item : catalog.items) {
                    const auto identity = item.properties.find("campaign id");
                    if (item.kind != "batch" || identity == item.properties.end() ||
                        identity->second != campaignId) continue;
                    item.properties["global campaign points"] = std::to_string(
                        campaign.at("points").array().size());
                    item.properties["campaign type"] = campaignType(campaign);
                    if (!campaign.at("bifurcation_tail_fraction").isNull())
                        item.properties["bifurcation tail fraction"] =
                            json::dump(campaign.at("bifurcation_tail_fraction"));
                    for (const auto& axis : campaign.at("axes").array()) {
                        std::string values;
                        for (const auto& value : axis.at("values").array()) {
                            if (!values.empty()) values += ", ";
                            values += json::dump(value);
                        }
                        item.properties["axis · " + axis.at("path").string()] = values;
                    }
                }
            }
        } catch (...) { /* Unrecognized JSON remains available to the generic scan. */ }
    }
}

}  // namespace dsl

// `datapack trace`: follows the lineage recorded in archive metadata and lists what a pack author
// still has to supply. It reports gaps and never guesses a setting.
#include <dynsyslab/workflow/DataPack.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>

#include <map>
#include <set>
#include <variant>

namespace dynsyslab::workflow {
namespace {

void collectParents(const json::Value& metadata, std::set<std::string>& parents) {
    for (const auto* key : {"source_archive", "source_krom", "source_spectrum_archive", "history_transform"}) {
        const auto text = metadata.at(key).string();
        if (!text.empty()) parents.insert(text);
    }
    for (const auto& item : metadata.at("dataset_references").at("references").array()) {
        const auto text = item.at("path").string();
        if (!text.empty()) parents.insert(text);
    }
    const auto krom = metadata.at("controller_configuration").at("settings").at("krom").string();
    if (!krom.empty()) parents.insert(krom);
}

void walk(const std::filesystem::path& path, std::map<std::string, json::Value>& nodes, json::Value::Array& gaps) {
    const auto key = path.generic_u8string();
    if (nodes.count(key)) return;
    nodes[key] = json::Value();
    if (!std::filesystem::is_regular_file(path)) {
        gaps.emplace_back("missing parent archive: " + key);
        return;
    }
    if (path.extension() == L".png") {
        gaps.emplace_back("a PNG records no lineage; trace its plotted archive: " + key);
        return;
    }
    const auto index = core::readNpzIndex(path);
    json::Value metadata;
    try { metadata = json::parse(index.metadataJson); }
    catch (const std::exception&) { gaps.emplace_back("unreadable metadata: " + key); return; }
    const auto method = metadata.at("method").string(metadata.at("schema").string() == "dynsyslab.solution/v1" ? "solve" : "");
    std::set<std::string> parents;
    collectParents(metadata, parents);
    json::Value::Object node{{"method", json::Value(method)}, {"label", metadata.at("label")}};
    if (metadata.at("run_provenance").isNull())
        gaps.emplace_back(key + ": no recorded resolved job or source commit (archive predates run provenance)");
    else node["resolved_job"] = metadata.at("run_provenance").at("resolved_job");
    if (method.empty()) gaps.emplace_back(key + ": producing method not recorded");
    if (method == "solve" && metadata.at("run_provenance").isNull()) {
        node["solve_settings"] = json::Value(json::Value::Object{{"time", metadata.at("time_configuration")},
            {"solver", metadata.at("solver_configuration")}, {"controller", metadata.at("controller_configuration")},
            {"initial_state", metadata.at("initial_state")}, {"system_parameters", metadata.at("system_parameters")}});
    }
    json::Value::Array parentList;
    for (const auto& parent : parents) parentList.emplace_back(parent);
    node["parents"] = json::Value(std::move(parentList));
    nodes[key] = json::Value(std::move(node));
    for (const auto& parent : parents) walk(std::filesystem::u8path(parent), nodes, gaps);
}

}  // namespace

json::Value traceArchiveLineage(const std::filesystem::path& archive) {
    std::map<std::string, json::Value> nodes;
    json::Value::Array gaps;
    walk(archive, nodes, gaps);
    json::Value::Object out;
    for (auto& [key, node] : nodes) out[key] = std::move(node);
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.data-pack-trace/v1")},
        {"root", json::Value(archive.generic_u8string())}, {"node_count", json::Value(double(out.size()))},
        {"nodes", json::Value(std::move(out))}, {"gaps", json::Value(std::move(gaps))}});
}

}  // namespace dynsyslab::workflow

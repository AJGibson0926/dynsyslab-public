#include "dynsyslab_viewer/ArchiveFocus.hpp"

#include "HierarchyCatalog.hpp"

#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

namespace dsl {
namespace {

std::filesystem::path resolved(const std::filesystem::path& projectRoot,
                               const std::string& source) {
    auto path = std::filesystem::u8path(source);
    if (path.is_relative()) path = projectRoot / path;
    return path;
}

}  // namespace

void applyArchiveFocus(Catalog& catalog,
                       const std::filesystem::path& manifest,
                       const std::filesystem::path& projectRoot) {
    namespace json = dynsyslab::core::json;
    if (manifest.empty() || !std::filesystem::exists(manifest)) return;

    json::Value root;
    try { root = json::load(manifest.wstring()); }
    catch (const std::exception& error) {
        throw std::runtime_error("Unable to read Viewer Mask: " + std::string(error.what()));
    }
    if (root.at("schema").string() != "dynsyslab.viewer.archive-focus/v1"
        || !std::holds_alternative<bool>(root.at("active").storage()))
        throw std::runtime_error("Invalid Viewer Mask: " + manifest.generic_u8string());
    if (!root.at("active").boolean()) return;
    if (!std::holds_alternative<json::Value::Array>(root.at("allowed_sources").storage()))
        throw std::runtime_error("Viewer Mask requires an allowed_sources array");

    std::set<std::string> allowedSources;
    for (const auto& value : root.at("allowed_sources").array()) {
        if (!std::holds_alternative<std::string>(value.storage()))
            throw std::runtime_error("Viewer Mask sources must be paths");
        const auto source = value.string();
        if (!source.empty())
            allowedSources.insert(normalizedSource(resolved(projectRoot, source)));
    }

    std::set<std::string> retainedIds;
    std::map<std::string, const CatalogItem*> byId;
    for (const auto& item : catalog.items) {
        byId.emplace(item.id, &item);
        if (allowedSources.count(normalizedSource(resolved(projectRoot, item.source))))
            retainedIds.insert(item.id);
    }
    // An active Mask remains active even when none of its records are present.

    bool changed = true;
    while (changed) {
        changed = false;
        const auto snapshot = retainedIds;
        for (const auto& id : snapshot) {
            const auto found = byId.find(id);
            if (found == byId.end() || found->second->parentId.empty()) continue;
            changed = retainedIds.insert(found->second->parentId).second || changed;
        }
    }

    const std::size_t originalCount = catalog.items.size();
    catalog.items.erase(std::remove_if(catalog.items.begin(), catalog.items.end(),
        [&](const CatalogItem& item) { return !retainedIds.count(item.id); }),
        catalog.items.end());
    for (auto& item : catalog.items) {
        item.relations.erase(std::remove_if(item.relations.begin(), item.relations.end(),
            [&](const Relation& relation) { return !retainedIds.count(relation.target); }),
            item.relations.end());
    }
    catalog.focusLabel = root.at("label").string("Archive focus");
    if (catalog.focusLabel.empty()) catalog.focusLabel = "Archive focus";
    catalog.hiddenItemCount = originalCount - catalog.items.size();
}

}  // namespace dsl

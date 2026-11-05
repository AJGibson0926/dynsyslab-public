#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

namespace dsl {
namespace {

using dynsyslab::core::json::Value;

std::vector<std::string> strings(const Value& value) {
    std::vector<std::string> result;
    for (const auto& item : value.array()) result.push_back(item.string());
    return result;
}

Value stringArray(const std::vector<std::string>& values) {
    Value::Array result;
    for (const auto& value : values) result.emplace_back(value);
    return Value(std::move(result));
}

std::string propertyString(const Value& value) {
    const auto text = value.string();
    if (!text.empty()) return text;
    if (value.isNull()) return {};
    constexpr double sentinel = -9.87654321e307;
    const double number = value.number(sentinel);
    if (number != sentinel) return std::to_string(number);
    return value.boolean(false) ? "yes" : "no";
}

CatalogItem itemFrom(const Value& value) {
    CatalogItem result;
    result.id = value.at("id").string();
    result.kind = value.at("kind").string();
    result.category = value.at("category").string();
    result.label = value.at("label").string();
    result.source = value.at("source").string();
    result.system = value.at("system").string();
    result.subtitle = value.at("subtitle").string();
    result.modified = value.at("modified").number();
    result.bytes = value.at("bytes").number();
    result.views = strings(value.at("views"));
    result.tags = strings(value.at("tags"));
    result.role = value.at("role").string("primary");
    result.parentId = value.at("parent").string();
    result.lineageId = value.at("lineage").string(result.id);
    result.primary = value.at("primary").boolean(true);
    result.variantCount = static_cast<int>(value.at("variants").number());
    for (const auto& entry : value.at("properties").object())
        result.properties.emplace(entry.first, propertyString(entry.second));
    for (const auto& relation : value.at("relations").array())
        result.relations.push_back({relation.at("target").string(),
            relation.at("kind").string(), relation.at("label").string()});
    return result;
}

Value itemValue(const CatalogItem& item) {
    Value::Object value;
    value["id"] = item.id; value["kind"] = item.kind;
    value["category"] = item.category; value["label"] = item.label;
    value["source"] = item.source; value["system"] = item.system;
    value["subtitle"] = item.subtitle; value["modified"] = item.modified;
    value["bytes"] = item.bytes; value["views"] = stringArray(item.views);
    value["tags"] = stringArray(item.tags); value["role"] = item.role;
    value["parent"] = item.parentId; value["lineage"] = item.lineageId;
    value["primary"] = item.primary;
    value["variants"] = static_cast<double>(item.variantCount);
    Value::Object properties;
    for (const auto& entry : item.properties) properties[entry.first] = entry.second;
    value["properties"] = Value(std::move(properties));
    Value::Array relations;
    for (const auto& relation : item.relations) {
        Value::Object encoded;
        encoded["target"] = relation.target; encoded["kind"] = relation.kind;
        encoded["label"] = relation.label;
        relations.emplace_back(std::move(encoded));
    }
    value["relations"] = Value(std::move(relations));
    return Value(std::move(value));
}

std::string systemToken(std::string value) {
    std::string result;
    for (const unsigned char item : value)
        if (std::isalnum(item)) result.push_back(
            static_cast<char>(std::tolower(item)));
    return result;
}

void canonicalize(Catalog& catalog) {
    for (auto& item : catalog.items) {
        const auto previous = item.system;
        item.system = canonicalSystemLabel(previous);
        for (auto& tag : item.tags)
            if (systemToken(tag) == systemToken(previous)) tag = item.system;
        if (item.kind == "system_run" &&
            systemToken(item.label) == systemToken(previous)) item.label = item.system;
    }
}

}  // namespace

Catalog Catalog::fromJson(const std::string& text) {
    const auto root = dynsyslab::core::json::parse(text);
    if (root.at("schema").string() != "dynsyslab.viewer.catalog/v1")
        throw std::runtime_error("Unsupported Viewer catalog schema");
    Catalog result;
    result.projectName = root.at("project").at("name").string();
    result.projectRoot = root.at("project").at("root").string();
    for (const auto& entry : root.at("view_labels").object())
        result.viewLabels.emplace(entry.first, entry.second.string());
    const auto rootPath = std::filesystem::u8path(result.projectRoot);
    for (const auto& item : root.at("items").array()) {
        auto loaded = itemFrom(item);
        auto source = std::filesystem::u8path(loaded.source);
        if (source.is_relative() && !rootPath.empty()) {
            source = rootPath / source;
            loaded.source = source.lexically_normal().generic_u8string();
        }
        result.items.push_back(std::move(loaded));
    }
    canonicalize(result);
    return result;
}

std::string Catalog::toJson() const {
    Value::Object root;
    root["schema"] = "dynsyslab.viewer.catalog/v1";
    Value::Object project;
    project["name"] = projectName; project["root"] = projectRoot;
    root["project"] = Value(std::move(project));
    Value::Object labels;
    for (const auto& entry : viewLabels) labels[entry.first] = entry.second;
    root["view_labels"] = Value(std::move(labels));
    Value::Array encodedItems;
    for (const auto& item : items) encodedItems.push_back(itemValue(item));
    root["items"] = Value(std::move(encodedItems));
    return dynsyslab::core::json::dump(Value(std::move(root)));
}

}  // namespace dsl

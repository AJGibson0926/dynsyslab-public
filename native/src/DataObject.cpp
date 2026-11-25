#include <dynsyslab/core/DataObject.hpp>

namespace dynsyslab::core {

json::Value dataObjectJson(const DataObjectInfo& info) {
    return json::Value::Object{{"schema", kDataObjectSchema}, {"kind", info.kind},
                               {"produced_by", dataOriginJson(info.producedBy)},
                               {"properties", info.properties.isNull() ? json::Value(json::Value::Object{}) : info.properties}};
}

std::optional<DataObjectInfo> dataObjectFromMetadata(const std::string& metadataJson) {
    if (metadataJson.find("\"data_object\"") == std::string::npos) return std::nullopt;
    const auto block = json::parse(metadataJson).at("data_object");
    if (block.at("schema").string() != kDataObjectSchema) return std::nullopt;
    return DataObjectInfo{block.at("kind").string(), dataOriginFromJson(block.at("produced_by")), block.at("properties")};
}

std::string dataObjectKind(const std::string& metadataJson) {
    const auto info = dataObjectFromMetadata(metadataJson);
    return info ? info->kind : std::string();
}

std::filesystem::path objectDirectory(const std::filesystem::path& archiveRoot, const std::string& kind) {
    return archiveRoot / L"objects" / std::filesystem::u8path(kind);
}

}  // namespace dynsyslab::core

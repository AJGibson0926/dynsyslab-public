#include "CliArchiveTools.hpp"
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/OutputPath.hpp>

#include <iostream>
#include <stdexcept>

namespace dynsyslab::cli {
namespace json = dynsyslab::core::json;

void materializeArchive(const std::filesystem::path& recorded, const std::filesystem::path& outputArgument) {
    const auto output = dynsyslab::core::prepareOutputFile(outputArgument, "materialize output");
    auto record = dynsyslab::core::readNpz(recorded);
    auto meta = json::parse(record.metadataJson());
    if (meta.at("recorded_archive").isNull()) throw std::invalid_argument("the archive already stores its arrays");
    meta.at("materialized_from") = json::Value(std::filesystem::absolute(recorded).generic_u8string());
    meta.at("materialized_record") = meta.at("recorded_archive");
    meta.object().erase("recorded_archive");
    record.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(json::dump(meta));
    dynsyslab::core::writeNpzAtomicNew(output, record);
    std::cout << "{\"schema\":\"dynsyslab.materialized-archive/v1\",\"output\":"
              << json::dump(json::Value(output.generic_u8string())) << "}\n";
}

}  // namespace dynsyslab::cli

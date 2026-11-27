#pragma once
// Archives hold data objects, never a method's state. Every archive a processing step writes carries
// one envelope in metadata_json["data_object"]:
//
//   {"schema": "dynsyslab.data-object/v1",
//    "kind": "krom" | "dynamical_model" | "time_series" | ...,
//    "produced_by": {"method": "hdmd", "settings": {...}, "inputs": [{"path", "sha256", "role"}]},
//    "properties": {...}}
//
// The kind says what the object is; produced_by says how it was made (a tag, not an owner). Archives
// are filed by kind under <archive root>/objects/<kind>/.
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <filesystem>
#include <optional>
#include <string>

namespace dynsyslab::core {

inline constexpr const char* kDataObjectSchema = "dynsyslab.data-object/v1";

// Object kinds the suite writes. A kind names the contract of the archive's arrays.
namespace object_kind {
// A Koopman reduced-order model: spectrum and eigenfunctions, optionally input maps B and readouts C;
// one that extends another references it as its parent (KoopmanArchive.hpp).
inline constexpr const char* krom = "krom";
inline constexpr const char* dynamicalModel = "dynamical_model";   // an identified vector field
inline constexpr const char* timeSeries = "time_series";
inline constexpr const char* coordinateMap = "coordinate_map";
inline constexpr const char* matrix = "matrix";
inline constexpr const char* field = "field";
inline constexpr const char* figure = "figure";                    // animations and figure datasets
inline constexpr const char* report = "report";                    // diagnostics and comparisons
}  // namespace object_kind

struct DataObjectInfo {
    std::string kind;
    DataOrigin producedBy;
    json::Value properties;
};

[[nodiscard]] json::Value dataObjectJson(const DataObjectInfo& info);
// The envelope of an archive's metadata text; nullopt when it has none.
[[nodiscard]] std::optional<DataObjectInfo> dataObjectFromMetadata(const std::string& metadataJson);
// The declared kind, or empty.
[[nodiscard]] std::string dataObjectKind(const std::string& metadataJson);
[[nodiscard]] std::filesystem::path objectDirectory(const std::filesystem::path& archiveRoot, const std::string& kind);

}  // namespace dynsyslab::core

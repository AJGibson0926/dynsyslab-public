#pragma once

#include <dynsyslab/core/Json.hpp>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Dataset use is relative to a particular model; it is independent of forcing.
// References contain no numerical arrays. A separate input archive is optional.
struct KoopmanDatasetReference {
    std::string path;
    std::string role;       // training, validation, test
    std::string dynamics{"unknown"}; // unforced, forced, unknown (legacy only)
    std::string inputPath;
};

void mergeKoopmanDatasetReference(std::vector<KoopmanDatasetReference>& references,
                                  KoopmanDatasetReference reference);
[[nodiscard]] json::Value koopmanDatasetReferencesJson(
    const std::vector<KoopmanDatasetReference>& references);
// Reads the versioned contract and recognized legacy provenance keys. Unknown
// forcing remains unknown; an algorithm name is not evidence of free dynamics.
[[nodiscard]] std::vector<KoopmanDatasetReference> koopmanDatasetReferencesFromMetadata(
    const json::Value& metadata);

} // namespace dynsyslab::core

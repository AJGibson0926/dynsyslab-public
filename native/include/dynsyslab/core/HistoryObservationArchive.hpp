#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
namespace dynsyslab::core {
// Embed the same transform contract in a KROM or a standalone derived archive.
void writeFrozenHistoryTransform(NpzRecord& archive,const FrozenHistoryTransform& transform);
FrozenHistoryTransform readFrozenHistoryTransform(const NpzRecord& archive);
FrozenHistoryTransform loadFrozenHistoryTransform(const std::filesystem::path& path);
void saveFrozenHistoryTransform(const std::filesystem::path& path,
    const FrozenHistoryTransform& transform,const std::string& metadataJson = "{}");
}

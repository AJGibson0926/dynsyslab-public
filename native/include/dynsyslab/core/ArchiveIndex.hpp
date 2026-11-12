#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

inline constexpr std::uint32_t kArchiveSummaryVersion = 4;

struct ArchiveIndexEntry {
    std::filesystem::path path;
    std::string relativePath;
    std::uint64_t modifiedTicks{};
    std::uint64_t sizeBytes{};
    std::string sha256;
    std::vector<std::string> arrays;
    std::uint32_t summaryVersion{};
    std::string summaryJson;
};

struct ArchiveIndexUpdate {
    std::size_t discovered{};
    std::size_t unchanged{};
    std::size_t indexed{};
    std::size_t removed{};
    std::size_t failed{};
};

// The Archiver-owned metadata index lives beside the archives, never in a UI
// cache. File size and modification time are the cheap change detector; a
// SHA-256 digest is computed only when a record is new or changed.
[[nodiscard]] std::filesystem::path archiveIndexPath(
    const std::filesystem::path& archiveRoot);
[[nodiscard]] std::string summarizeArchiveMetadata(
    const std::string& metadataJson);
[[nodiscard]] std::vector<ArchiveIndexEntry> loadArchiveIndex(
    const std::filesystem::path& archiveRoot);
// Consumer snapshots are derived views of the Archiver-owned index.  They are
// returned only while their schema and index generation are current.
[[nodiscard]] std::optional<std::string> loadArchiveCatalogSnapshot(
    const std::filesystem::path& archiveRoot, const std::string& consumer,
    std::uint32_t schemaVersion);
[[nodiscard]] bool storeArchiveCatalogSnapshot(
    const std::filesystem::path& archiveRoot, const std::string& consumer,
    std::uint32_t schemaVersion, std::uint64_t sourceGeneration,
    const std::string& payloadJson);
[[nodiscard]] std::uint64_t archiveIndexGeneration(
    const std::filesystem::path& archiveRoot);
ArchiveIndexUpdate indexArchiveFiles(
    const std::filesystem::path& archiveRoot,
    const std::vector<std::filesystem::path>& paths);
ArchiveIndexUpdate reconcileArchiveIndex(
    const std::filesystem::path& archiveRoot);

}  // namespace dynsyslab::core

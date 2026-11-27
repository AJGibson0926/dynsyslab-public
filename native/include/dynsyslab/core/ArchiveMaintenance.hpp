#pragma once
// Lossless storage maintenance of archives (dynsyslab-cli archive recompress | share-time-grids |
// edit-metadata). Scientific content is never changed: recompression and grid sharing rewrite a
// file only in how it stores its arrays, and verify that readNpz returns identical arrays before
// the rewrite is published; a metadata edit writes a new archive and leaves its source as it was.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ArchiveAudit.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ArchiveRewrite {
    std::filesystem::path path;
    std::uintmax_t bytesBefore{};
    std::uintmax_t bytesAfter{};
    bool rewritten{};
    std::string error;  // non-empty: the file was left unchanged
};

struct ArchiveMaintenanceResult {
    std::string schema;  // dynsyslab.archive-recompress/v1 | dynsyslab.time-grid-sharing/v1
    std::filesystem::path root;
    bool dryRun{true};
    std::string compression;
    std::vector<ArchiveRewrite> files;
    std::vector<std::filesystem::path> grids;  // shared grids written (or that would be)
};

// Re-encodes each archive with `compression`, in place and atomically, keeping every stored array
// and metadata_json byte for byte. Without `apply`, reports the files only.
[[nodiscard]] ArchiveMaintenanceResult recompressArchives(
    const std::vector<std::filesystem::path>& archives, ArchiveCompression compression, bool apply);

// Stores each final time mesh shared by two or more of the archives once, as
// <root>/time_grids/<sha256>.npz, and rewrites those archives without their `time` array
// (TimeGrid.hpp). Raw adaptive meshes (`raw_time`) stay with their solve. Archives that already
// reference a grid, are recorded archives, or have no `time` are left alone.
[[nodiscard]] ArchiveMaintenanceResult shareTimeGrids(
    const std::filesystem::path& root, const std::vector<std::filesystem::path>& archives, bool apply);

[[nodiscard]] std::string archiveMaintenanceJson(const ArchiveMaintenanceResult& result);

// Writes `output`, a new archive with the arrays of `source` and its metadata changed by
// `patchJson` (a JSON object: each key is set, a null value removes the key; RFC 7386 at the top
// level). The edit is recorded in metadata "metadata_edits" (source path and SHA-256, keys, UTC
// time) and "metadata_revision" is incremented. Keys that define how the arrays are read
// (recorded_archive, time_grid_*, run_provenance, metadata_edits) cannot be patched.
// Fails when `output` exists or is `source`. Returns the output path.
std::filesystem::path editArchiveMetadata(const std::filesystem::path& source,
                                          const std::filesystem::path& output,
                                          const std::string& patchJson);

}  // namespace dynsyslab::core

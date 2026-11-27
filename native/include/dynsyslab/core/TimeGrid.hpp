#pragma once
// Shared time grids. Solutions that sample the same final time mesh may store it once, in
// <archive root>/time_grids/<sha256>.npz (one `time` array), instead of each storing `time`.
// Such an archive has no `time` (or `raw_time`) array; its metadata names the grid:
//   "time_grid_ref":     "time_grids/<sha256>.npz"   (relative, generic separators)
//   "time_grid_root":    absolute archive root the reference is relative to
//   "time_grid_shared":  true
//   "time_grid_samples": sample count, checked on read (absent in older archives)
// ("raw_time_grid_ref" does the same for `raw_time`). readNpz / readNpzIndex restore the array,
// so every reader sees an ordinary archive. The reference resolves against time_grid_root, then
// archive_root, then each folder above the archive, so a relocated archive tree still reads.
// Writing shared grids: core/ArchiveMaintenance.hpp shareTimeGrids.
#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::core {

inline constexpr const char* kTimeGridSchema = "dynsyslab.time-grid/v1";
inline constexpr const char* kTimeGridFolder = "time_grids";

// The grid file a reference names (key "time_grid_ref" or "raw_time_grid_ref"); empty when
// metadata has no such reference. Throws when the reference is unsafe (absolute or with "..")
// or names no existing file.
[[nodiscard]] std::filesystem::path resolveTimeGrid(const std::filesystem::path& archive,
                                                    const std::string& metadataJson,
                                                    const std::string& referenceKey);
// Adds `time` / `raw_time` from referenced grids to a record read from `archive`.
void attachTimeGrids(const std::filesystem::path& archive, NpzRecord& record);
// Lists `time` / `raw_time` in an index whose metadata references a grid.
void addTimeGridNames(NpzIndex& index);

}  // namespace dynsyslab::core

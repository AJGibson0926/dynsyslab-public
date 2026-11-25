#pragma once
// dynsyslab-cli materialize: a recorded archive (core/Archive.hpp recorded archives, e.g. a
// coordinate change stored as its record) written out with its recomputed arrays, for tools that
// cannot replay it. Prints a dynsyslab.materialized-archive/v1 JSON line.
#include <filesystem>

namespace dynsyslab::cli {

void materializeArchive(const std::filesystem::path& recorded, const std::filesystem::path& output);

}  // namespace dynsyslab::cli

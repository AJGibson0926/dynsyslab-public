#pragma once
#include <dynsyslab/core/EdmdControl.hpp>
#include <filesystem>

namespace dynsyslab::core {
// Read the canonical immutable EDMDc artifact without loading its training data
// or constructing a physical system. Requires the archive-I/O build option.
[[nodiscard]] EdmdControlResult loadEdmdControlModel(const std::filesystem::path& path);
}

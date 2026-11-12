#pragma once

#include <filesystem>
#include <string>

namespace dynsyslab::core {

// Lower-case SHA-256 of the exact file bytes. Kept in the archive layer so
// imported evidence has the same durable identity regardless of UI or adapter.
[[nodiscard]] std::string sha256Text(const std::string& text);
[[nodiscard]] std::string sha256File(const std::filesystem::path& path);

}  // namespace dynsyslab::core

#pragma once

#include <cstddef>
#include <filesystem>
#include <string>

namespace dsl {

// Persists the small amount of state needed to diagnose installed-app startup.
// This never scans or modifies the scientific archive.
void writeStartupDiagnostics(
    const std::string& application,
    const std::filesystem::path& cacheRoot,
    const std::filesystem::path& projectRoot,
    const std::filesystem::path& archiveRoot,
    const std::filesystem::path& catalogPath,
    std::size_t catalogItems,
    std::size_t navigationNodes,
    const std::filesystem::path& archiveFocus = {},
    const std::string& focusLabel = {},
    std::size_t hiddenItems = 0) noexcept;

}  // namespace dsl

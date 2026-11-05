#pragma once

#include <filesystem>
#include <string>

namespace dsl {

struct RuntimePaths {
    std::filesystem::path executable;
    std::filesystem::path applicationRoot;
    std::filesystem::path projectRoot;
    std::filesystem::path archiveRoot;
    std::filesystem::path cacheRoot;
    // Persistent application preferences do not follow --cache overrides.
    std::filesystem::path preferencesRoot;
};

// Deterministic resolver used by the native applications and their tests.
// userDataRoot is the shared DynSysLab directory, not its parent.
[[nodiscard]] RuntimePaths resolveRuntimePaths(
    const std::wstring& applicationName,
    const std::filesystem::path& executable,
    const std::filesystem::path& userDataRoot,
    const std::filesystem::path& projectOverride = {},
    const std::filesystem::path& archiveOverride = {},
    const std::filesystem::path& cacheOverride = {});

[[nodiscard]] RuntimePaths discoverRuntimePaths(
    const std::wstring& applicationName,
    const std::filesystem::path& projectOverride = {},
    const std::filesystem::path& archiveOverride = {},
    const std::filesystem::path& cacheOverride = {});

}  // namespace dsl

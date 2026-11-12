#include "dynsyslab_viewer/RuntimePaths.hpp"

#include <dynsyslab/core/Json.hpp>

#include <shlobj.h>
#include <windows.h>

#include <stdexcept>
#include <utility>
#include <vector>

namespace dsl {
namespace {

std::filesystem::path executablePath() {
    std::vector<wchar_t> buffer(32768);
    const DWORD count = GetModuleFileNameW(nullptr, buffer.data(),
                                           static_cast<DWORD>(buffer.size()));
    if (!count || count >= buffer.size())
        throw std::runtime_error("Unable to resolve the DynSysLab executable path");
    return std::filesystem::path(std::wstring(buffer.data(), count));
}

std::filesystem::path localAppData() {
    PWSTR value = nullptr;
    const HRESULT status = SHGetKnownFolderPath(FOLDERID_LocalAppData, KF_FLAG_DEFAULT,
                                                 nullptr, &value);
    if (FAILED(status) || !value)
        throw std::runtime_error("Unable to resolve the local application-data directory");
    std::filesystem::path result(value);
    CoTaskMemFree(value);
    return result;
}

bool sourceRoot(const std::filesystem::path& candidate) {
    return std::filesystem::is_directory(candidate / L"native") &&
           std::filesystem::is_directory(candidate / L"archives") &&
           std::filesystem::is_directory(candidate / L"docs");
}

struct SavedPaths {
    std::filesystem::path projectRoot;
    std::filesystem::path archiveRoot;
};

SavedPaths loadSavedPaths(const std::filesystem::path& userRoot,
                          const std::filesystem::path& executable) {
    namespace json = dynsyslab::core::json;
    const std::pair<std::filesystem::path, std::string> candidates[]{
        {executable.parent_path().parent_path() / L"resources" / L"runtime.json",
         "dynsyslab.runtime-paths/v1"},
        {userRoot / L"runtime.json", "dynsyslab.runtime-paths/v1"},
        {userRoot / L"Compute" / L"launch.json", "dynsyslab.compute-launch/v1"}};
    for (const auto& candidate : candidates) {
        if (!std::filesystem::is_regular_file(candidate.first)) continue;
        try {
            const auto root = json::load(candidate.first.wstring());
            if (root.at("schema").string() != candidate.second) continue;
            SavedPaths saved{
                std::filesystem::u8path(root.at("project_root").string()),
                std::filesystem::u8path(root.at("archive_root").string())};
            if (!saved.projectRoot.empty()
                && !std::filesystem::is_directory(saved.projectRoot))
                saved.projectRoot.clear();
            if (!saved.archiveRoot.empty()
                && !std::filesystem::is_directory(saved.archiveRoot))
                saved.archiveRoot.clear();
            if (!saved.projectRoot.empty() || !saved.archiveRoot.empty()) return saved;
        } catch (...) {
            // A stale or malformed preference must never prevent the apps from
            // falling back to their self-contained user-data location.
        }
    }
    return {};
}

std::filesystem::path discoverProject(const std::filesystem::path& executable) {
    auto candidate = executable.parent_path();
    while (!candidate.empty()) {
        if (sourceRoot(candidate)) return candidate;
        if (candidate == candidate.parent_path()) break;
        candidate = candidate.parent_path();
    }
    const auto installed = executable.parent_path() / L"resources";
    if (std::filesystem::is_directory(installed)) return installed;
    const auto distribution = executable.parent_path().parent_path() / L"resources";
    return std::filesystem::is_directory(distribution) ? distribution : executable.parent_path();
}

}  // namespace

RuntimePaths resolveRuntimePaths(const std::wstring& applicationName,
                                 const std::filesystem::path& executable,
                                 const std::filesystem::path& userRoot,
                                 const std::filesystem::path& projectOverride,
                                 const std::filesystem::path& archiveOverride,
                                 const std::filesystem::path& cacheOverride) {
    RuntimePaths result;
    result.executable = executable;
    result.applicationRoot = executable.parent_path();
    const auto discoveredProject = discoverProject(executable);
    const auto saved = sourceRoot(discoveredProject)
        ? SavedPaths{} : loadSavedPaths(userRoot, executable);
    result.projectRoot = !projectOverride.empty() ? projectOverride
        : (!saved.projectRoot.empty() ? saved.projectRoot : discoveredProject);
    result.archiveRoot = !archiveOverride.empty() ? archiveOverride
        : (!projectOverride.empty() && sourceRoot(result.projectRoot)
            ? result.projectRoot / L"archives"
            : (!saved.archiveRoot.empty() ? saved.archiveRoot
            : (sourceRoot(result.projectRoot) ? result.projectRoot / L"archives"
                                              : userRoot / L"archives")));
    result.preferencesRoot = userRoot / applicationName;
    result.cacheRoot = cacheOverride.empty() ? result.preferencesRoot : cacheOverride;
    std::filesystem::create_directories(result.archiveRoot);
    std::filesystem::create_directories(result.cacheRoot);
    return result;
}

RuntimePaths discoverRuntimePaths(const std::wstring& applicationName,
                                  const std::filesystem::path& projectOverride,
                                  const std::filesystem::path& archiveOverride,
                                  const std::filesystem::path& cacheOverride) {
    const auto userRoot = localAppData() / L"DynSysLab";
    return resolveRuntimePaths(applicationName, executablePath(), userRoot,
        projectOverride, archiveOverride, cacheOverride);
}

}  // namespace dsl

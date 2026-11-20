#include "CliRuntimePaths.hpp"

#include <dynsyslab/core/Json.hpp>

#include <cstdlib>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::cli {
namespace {

bool sourceCheckout(const std::filesystem::path& candidate) {
    std::error_code ignored;
    return std::filesystem::is_directory(candidate / L"native", ignored)
        && std::filesystem::is_directory(candidate / L"archives", ignored)
        && std::filesystem::is_directory(candidate / L"docs", ignored);
}

std::filesystem::path savedArchiveRoot(const std::filesystem::path& file) {
    std::error_code ignored;
    if (!std::filesystem::is_regular_file(file, ignored)) return {};
    try {
        const auto root = core::json::load(file.wstring());
        if (root.at("schema").string() != "dynsyslab.runtime-paths/v1") return {};
        const auto path = std::filesystem::u8path(root.at("archive_root").string());
        return !path.empty() && std::filesystem::is_directory(path, ignored) ? path : std::filesystem::path{};
    } catch (...) {
        return {};  // a stale or malformed preference never blocks the CLI
    }
}

}  // namespace

std::filesystem::path cliExecutablePath() {
#ifdef _WIN32
    std::vector<wchar_t> buffer(32768);
    const DWORD count = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    if (count && count < buffer.size()) return std::filesystem::path(std::wstring(buffer.data(), count));
#endif
    return std::filesystem::current_path() / L"dynsyslab-cli";
}

std::filesystem::path sourceCheckoutContaining(const std::filesystem::path& executable) {
    for (auto candidate = executable.parent_path(); !candidate.empty(); candidate = candidate.parent_path()) {
        if (sourceCheckout(candidate)) return candidate;
        if (candidate == candidate.parent_path()) break;
    }
    return {};
}

std::filesystem::path userDataRoot() {
#ifdef _WIN32
    std::vector<wchar_t> local(32768);
    const DWORD size = GetEnvironmentVariableW(L"LOCALAPPDATA", local.data(), static_cast<DWORD>(local.size()));
    if (size && size < local.size()) return std::filesystem::path(std::wstring(local.data(), size)) / L"DynSysLab";
#else
    if (const char* home = std::getenv("HOME"); home && *home)
        return std::filesystem::path(home) / ".local" / "share" / "DynSysLab";
#endif
    return std::filesystem::current_path() / L"DynSysLab";
}

std::filesystem::path installedProjectRoot() {
    const auto file = cliExecutablePath().parent_path().parent_path() / L"resources" / L"runtime.json";
    std::error_code ignored;
    if (!std::filesystem::is_regular_file(file, ignored)) return {};
    try {
        const auto root = core::json::load(file.wstring());
        if (root.at("schema").string() != "dynsyslab.runtime-paths/v1") return {};
        return std::filesystem::u8path(root.at("project_root").string());
    } catch (...) { return {}; }
}

ArchiveRootChoice resolveArchiveRoot(const std::filesystem::path& executable,
                                     const std::filesystem::path& userData) {
    // An installation's own runtime.json wins, even when it is installed inside a checkout.
    if (auto saved = savedArchiveRoot(executable.parent_path().parent_path() / L"resources" / L"runtime.json");
        !saved.empty()) return {saved, "install_runtime"};
    if (const auto checkout = sourceCheckoutContaining(executable); !checkout.empty())
        return {checkout / L"archives", "source_checkout"};
    if (auto saved = savedArchiveRoot(userData / L"runtime.json"); !saved.empty())
        return {saved, "user_runtime"};
    return {userData / L"archives", "user_default"};
}

ArchiveRootChoice defaultArchiveRoot() { return resolveArchiveRoot(cliExecutablePath(), userDataRoot()); }

}  // namespace dynsyslab::cli

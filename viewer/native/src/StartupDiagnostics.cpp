#include "dynsyslab_viewer/StartupDiagnostics.hpp"

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Text.hpp>

#include <windows.h>

#include <fstream>
#include <vector>

namespace dsl {
namespace {

std::filesystem::path executablePath() {
    std::vector<wchar_t> value(32768);
    const DWORD size = GetModuleFileNameW(nullptr, value.data(),
                                          static_cast<DWORD>(value.size()));
    return size && size < value.size()
        ? std::filesystem::path(std::wstring(value.data(), size))
        : std::filesystem::path{};
}

std::string pathText(const std::filesystem::path& path) {
    return path.generic_u8string();
}

}  // namespace

void writeStartupDiagnostics(const std::string& application,
                             const std::filesystem::path& cacheRoot,
                             const std::filesystem::path& projectRoot,
                             const std::filesystem::path& archiveRoot,
                             const std::filesystem::path& catalogPath,
                             std::size_t catalogItems,
                             std::size_t navigationNodes,
                             const std::filesystem::path& archiveFocus,
                             const std::string& focusLabel,
                             std::size_t hiddenItems) noexcept {
    try {
        namespace json = dynsyslab::core::json;
        std::error_code ignored;
        const auto workingDirectory = std::filesystem::current_path(ignored);
        ignored.clear();
        const bool projectExists = std::filesystem::is_directory(projectRoot, ignored);
        ignored.clear();
        const bool archiveExists = std::filesystem::is_directory(archiveRoot, ignored);
        ignored.clear();
        const bool catalogExists = !catalogPath.empty() &&
            std::filesystem::is_regular_file(catalogPath, ignored);
        const DWORD processId = GetCurrentProcessId();
        json::Value::Object root;
        root["schema"] = "dynsyslab.native-startup/v1";
        root["application"] = application;
        root["process_id"] = static_cast<double>(processId);
        root["executable"] = pathText(executablePath());
        root["working_directory"] = pathText(workingDirectory);
        root["command_line"] = dynsyslab::core::narrow(GetCommandLineW());
        root["project_root"] = pathText(projectRoot);
        root["project_exists"] = projectExists;
        root["archive_root"] = pathText(archiveRoot);
        root["archive_exists"] = archiveExists;
        root["catalog_path"] = pathText(catalogPath);
        root["catalog_exists"] = catalogExists;
        root["catalog_items"] = static_cast<double>(catalogItems);
        root["navigation_nodes"] = static_cast<double>(navigationNodes);
        if (!archiveFocus.empty()) {
            root["archive_focus_path"] = pathText(archiveFocus);
            root["archive_focus_applied"] = !focusLabel.empty();
            root["archive_focus_label"] = focusLabel;
            root["archive_focus_hidden_items"] = static_cast<double>(hiddenItems);
        }
        std::filesystem::create_directories(cacheRoot, ignored);
        const auto destination = cacheRoot / L"startup-diagnostics.json";
        const auto temporary = cacheRoot /
            (L"startup-diagnostics-" + std::to_wstring(processId) + L".tmp");
        const std::string document = json::dump(json::Value(std::move(root)));
        {
            std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
            stream << document;
        }
        ignored.clear();
        std::filesystem::remove(destination, ignored);
        ignored.clear();
        std::filesystem::rename(temporary, destination, ignored);
    } catch (...) {
        // Diagnostics must never prevent the UI from starting.
    }
}

}  // namespace dsl

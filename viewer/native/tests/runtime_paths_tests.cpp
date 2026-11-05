#include "dynsyslab_viewer/RuntimePaths.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Json.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void writeConfiguration(const std::filesystem::path& path,
                        const std::filesystem::path& project,
                        const std::filesystem::path& archives) {
    namespace json = dynsyslab::core::json;
    json::Value::Object root;
    root["schema"] = "dynsyslab.runtime-paths/v1";
    root["project_root"] = project.generic_u8string();
    root["archive_root"] = archives.generic_u8string();
    std::filesystem::create_directories(path.parent_path());
    std::ofstream stream(path, std::ios::binary);
    stream << json::dump(json::Value(std::move(root)));
}

void writeCatalog(const std::filesystem::path& path,
                  const std::filesystem::path& project) {
    namespace json = dynsyslab::core::json;
    json::Value::Object item;
    item["id"] = "relative-source";
    item["kind"] = "trajectory";
    item["category"] = "Training data";
    item["label"] = "Relative source";
    item["source"] = "archives/example.npz";
    json::Value::Array items;
    items.emplace_back(std::move(item));
    json::Value::Object projectObject;
    projectObject["name"] = "DynSysLab";
    projectObject["root"] = project.generic_u8string();
    json::Value::Object root;
    root["schema"] = "dynsyslab.viewer.catalog/v1";
    root["project"] = json::Value(std::move(projectObject));
    root["items"] = json::Value(std::move(items));
    std::ofstream stream(path, std::ios::binary);
    stream << json::dump(json::Value(std::move(root)));
}

}  // namespace

int main() {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto root = std::filesystem::temp_directory_path()
        / ("dynsyslab-runtime-paths-" + std::to_string(stamp));
    try {
        const auto project = root / "project";
        const auto archives = project / "archives";
        const auto userData = root / "user" / "DynSysLab";
        const auto installedExecutable = root / "install" / "bin" / "DynSysLab.exe";
        std::filesystem::create_directories(project / "native");
        std::filesystem::create_directories(project / "docs");
        std::filesystem::create_directories(archives);
        std::filesystem::create_directories(installedExecutable.parent_path());
        writeConfiguration(userData / "runtime.json", project, archives);

        const auto installationProject = root / "installation-project";
        const auto installationArchives = installationProject / "archives";
        std::filesystem::create_directories(installationProject);
        std::filesystem::create_directories(installationArchives);
        writeConfiguration(root / "install" / "resources" / "runtime.json",
            installationProject, installationArchives);

        const auto installed = dsl::resolveRuntimePaths(
            L"Viewer", installedExecutable, userData);
        require(installed.projectRoot == installationProject,
            "installed app did not prefer its colocated project binding");
        require(installed.archiveRoot == installationArchives,
            "installed app did not prefer its colocated archive binding");

        const auto fallbackExecutable = root / "fallback-install" / "bin" /
            "DynSysLab.exe";
        std::filesystem::create_directories(fallbackExecutable.parent_path());
        const auto fallback = dsl::resolveRuntimePaths(
            L"Viewer", fallbackExecutable, userData);
        require(fallback.projectRoot == project,
            "installed app did not restore the per-user project fallback");
        require(fallback.archiveRoot == archives,
            "installed app did not restore the per-user archive fallback");

        const auto sourceExecutable = project / "build" / "bin" / "DynSysLab.exe";
        const auto source = dsl::resolveRuntimePaths(L"Compute", sourceExecutable, userData);
        require(source.projectRoot == project,
            "source build did not discover its source root");
        require(source.archiveRoot == archives,
            "source build did not prefer its colocated archive root");

        const auto explicitProject = root / "explicit-project";
        const auto explicitArchives = root / "explicit-archives";
        const auto explicitCache = root / "explicit-cache";
        std::filesystem::create_directories(explicitProject);
        const auto explicitPaths = dsl::resolveRuntimePaths(L"Viewer", installedExecutable,
            userData, explicitProject, explicitArchives, explicitCache);
        require(explicitPaths.projectRoot == explicitProject,
            "explicit project root did not win over saved configuration");
        require(explicitPaths.archiveRoot == explicitArchives,
            "explicit archive root did not win over saved configuration");
        require(explicitPaths.cacheRoot == explicitCache,
            "explicit cache root did not win over saved configuration");
        require(explicitPaths.preferencesRoot == installed.preferencesRoot
            && explicitPaths.preferencesRoot == userData / "Viewer",
            "Viewer preferences moved with a cache or project override");

        // Compute can use a temporary cache without moving any explicit Viewer
        // preferences. The launcher's sibling cache is not preferences.
        const auto computePaths = dsl::resolveRuntimePaths(L"Compute", installedExecutable,
            userData, project, archives, root / "temporary" / "Compute");
        const auto launchedViewer = dsl::resolveRuntimePaths(L"Viewer", installedExecutable,
            userData, project, archives, computePaths.cacheRoot.parent_path() / "Viewer");
        require(launchedViewer.cacheRoot != installed.cacheRoot
            && launchedViewer.preferencesRoot == installed.preferencesRoot,
            "Compute's Viewer launch changed the preferences location");

        const auto catalogPath = root / "catalog.json";
        writeCatalog(catalogPath, project);
        const auto catalog = dsl::Catalog::load(catalogPath.wstring());
        require(catalog.items.size() == 1,
            "catalog fixture did not load its single item");
        require(std::filesystem::u8path(catalog.items.front().source) ==
                (project / "archives" / "example.npz").lexically_normal(),
            "catalog source still depends on the launcher working directory");
    } catch (const std::exception& error) {
        std::error_code ignored;
        std::filesystem::remove_all(root, ignored);
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::error_code ignored;
    std::filesystem::remove_all(root, ignored);
    return 0;
}

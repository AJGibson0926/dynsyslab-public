// dynsyslab_cli default archive root: DYNSYSLAB_ARCHIVE_DIR, then the installer's runtime.json,
// source checkout, user runtime.json, per-user default; never the working directory.
#include "../src/CliRuntimePaths.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace fs = std::filesystem;
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }

void runtimeJson(const fs::path& file, const fs::path& archives, const char* schema = "dynsyslab.runtime-paths/v1") {
    fs::create_directories(file.parent_path());
    std::ofstream(file) << "{\"schema\":\"" << schema << "\",\"project_root\":\"\",\"archive_root\":\""
                        << archives.generic_u8string() << "\"}";
}
}  // namespace

int main() {
    const auto root = fs::temp_directory_path() / ("dynsyslab-cli-paths-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        const auto install = root / "install", user = root / "user", exe = install / "bin" / "dynsyslab-cli.exe";
        fs::create_directories(exe.parent_path());
        using dynsyslab::cli::resolveArchiveRoot;

        auto choice = resolveArchiveRoot(exe, user);
        require(choice.source == "user_default" && choice.path == user / "archives", "no configuration: per-user default");

        const auto userArchives = root / "user-archives";
        fs::create_directories(userArchives);
        runtimeJson(user / "runtime.json", userArchives);
        choice = resolveArchiveRoot(exe, user);
        require(choice.source == "user_runtime" && fs::equivalent(choice.path, userArchives), "user runtime.json");

        const auto installArchives = root / "Documents" / "DynSysLab" / "archives";
        fs::create_directories(installArchives);
        runtimeJson(install / "resources" / "runtime.json", installArchives);
        choice = resolveArchiveRoot(exe, user);
        require(choice.source == "install_runtime" && fs::equivalent(choice.path, installArchives), "install runtime.json wins");

        runtimeJson(install / "resources" / "runtime.json", root / "missing");  // stale folder: skipped
        choice = resolveArchiveRoot(exe, user);
        require(choice.source == "user_runtime", "missing archive folder falls through");
        runtimeJson(user / "runtime.json", userArchives, "other/v1");          // wrong schema: skipped
        choice = resolveArchiveRoot(exe, user);
        require(choice.source == "user_default", "wrong schema falls through");

        const auto checkout = root / "checkout";
        for (const char* dir : {"native", "archives", "docs", "build/bin"}) fs::create_directories(checkout / dir);
        choice = resolveArchiveRoot(checkout / "build" / "bin" / "dynsyslab-cli.exe", user);
        require(choice.source == "source_checkout" && choice.path == checkout / "archives", "source checkout archives");

        const auto nested = checkout / "build" / "install-test" / "app";  // an install inside a checkout
        fs::create_directories(nested / "bin");
        runtimeJson(nested / "resources" / "runtime.json", installArchives);
        choice = resolveArchiveRoot(nested / "bin" / "dynsyslab-cli.exe", user);
        require(choice.source == "install_runtime" && fs::equivalent(choice.path, installArchives),
                "an installation's runtime.json wins inside a checkout");

        const auto chosen = root / "from-environment";  // DYNSYSLAB_ARCHIVE_DIR wins over runtime.json
        choice = resolveArchiveRoot(nested / "bin" / "dynsyslab-cli.exe", user, chosen);
        require(choice.source == "environment" && choice.path == chosen.lexically_normal(),
                "DYNSYSLAB_ARCHIVE_DIR wins over every configured default");

        fs::remove_all(root);
        std::cout << "dynsyslab_cli_runtime_paths_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; fs::remove_all(root, ignored);
        std::cerr << "dynsyslab_cli_runtime_paths_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}

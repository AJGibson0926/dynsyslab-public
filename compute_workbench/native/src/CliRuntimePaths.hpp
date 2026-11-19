#pragma once
// Default archive root of dynsyslab_cli when --archive-root is not given. Same order as the
// desktop apps (viewer/native/src/RuntimePaths.cpp), never the current working directory:
//   1. source checkout containing the executable   -> <checkout>/archives
//   2. <install>/resources/runtime.json             -> its archive_root (written by install.ps1)
//   3. %LOCALAPPDATA%/DynSysLab/runtime.json        -> its archive_root
//   4. %LOCALAPPDATA%/DynSysLab/archives
// A runtime.json is used only with schema dynsyslab.runtime-paths/v1 and an existing folder.
#include <filesystem>
#include <string>

namespace dynsyslab::cli {

struct ArchiveRootChoice {
    std::filesystem::path path;
    std::string source;  // "source_checkout", "install_runtime", "user_runtime", "user_default"
};

[[nodiscard]] ArchiveRootChoice resolveArchiveRoot(const std::filesystem::path& executable,
                                                   const std::filesystem::path& userDataRoot);
[[nodiscard]] ArchiveRootChoice defaultArchiveRoot();

// The running executable, and the source checkout containing it (empty for an install).
[[nodiscard]] std::filesystem::path cliExecutablePath();
[[nodiscard]] std::filesystem::path sourceCheckoutContaining(const std::filesystem::path& executable);

}  // namespace dynsyslab::cli

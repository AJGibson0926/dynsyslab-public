// Viewer command line: smoke/background checks and backend roots, resolved against the
// discovered runtime paths.
#include "dynsyslab_viewer/Workbench.hpp"
#include "dynsyslab_viewer/RuntimePaths.hpp"

#include <shellapi.h>

#include <filesystem>
#include <stdexcept>
#include <string>

namespace dsl {

WorkbenchArguments parseArguments() {
    int count = 0;
    wchar_t** values = CommandLineToArgvW(GetCommandLineW(), &count);
    WorkbenchArguments result;
    for (int index = 1; index < count; ++index) {
        const std::wstring key = values[index];
        if (key == L"--smoke-test") { result.smokeTest = true; continue; }
        if (key == L"--background-check") { result.backgroundCheck = true; continue; }
        if (key == L"--background-krom-check") { result.backgroundCheck = result.backgroundKromCheck = true; continue; }
        if (key == L"--background-workspace-check") {
            result.backgroundCheck = result.backgroundWorkspaceCheck = true; continue;
        }
        if (index + 1 >= count) break;
        const std::filesystem::path value = values[++index];
        if (key == L"--project-root") result.backend.projectRoot = value;
        else if (key == L"--archive-root") result.backend.archiveRoot = value;
        else if (key == L"--catalog") result.backend.catalog = value;
        else if (key == L"--cache") result.backend.cache = value;
        else if (key == L"--archive-focus") result.backend.archiveFocus = value;
        else if (key == L"--source") result.initialSource = value;
        else if (key == L"--background-refresh-count") {
            result.backgroundRefreshCount = std::stoul(value.wstring());
            if (result.backgroundRefreshCount > 20)
                throw std::runtime_error("Background refresh count must be at most 20");
        }
        else --index;
    }
    LocalFree(values);
    const auto paths = discoverRuntimePaths(L"Viewer", result.backend.projectRoot,
        result.backend.archiveRoot, result.backend.cache);
    result.backend.projectRoot = paths.projectRoot;
    result.backend.archiveRoot = paths.archiveRoot;
    result.backend.cache = paths.cacheRoot;
    return result;
}

}  // namespace dsl

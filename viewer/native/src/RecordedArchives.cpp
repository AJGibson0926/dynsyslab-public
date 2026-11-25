// Recorded archives (core/Archive.hpp) are replayed by the workflow library, which registers its
// readers when it loads. Viewer does not link it, so it loads it once from its own folder; without
// it, opening a recorded archive reports how to materialize one.
#include "RecordedArchives.hpp"

#include <filesystem>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dsl {

void loadRecordedArchiveReaders() noexcept {
#ifdef _WIN32
    static const bool loaded = [] {
        std::vector<wchar_t> filename(32768);
        const auto length = GetModuleFileNameW(nullptr, filename.data(), DWORD(filename.size()));
        if (!length || length >= filename.size()) return false;
        const auto dll = std::filesystem::path(std::wstring(filename.data(), length)).parent_path() / L"dynsyslab_workflow.dll";
        return LoadLibraryExW(dll.c_str(), nullptr, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS) != nullptr;
    }();
    (void)loaded;
#endif
}

}  // namespace dsl

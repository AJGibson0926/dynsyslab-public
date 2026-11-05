#include "PhaseContext.hpp"
#include <stdexcept>
#include <vector>
#include <windows.h>
namespace dsl {
namespace {
// The independently packaged workflow module is optional in a standalone Viewer.
// Resolve its C ABI from the executable directory, never a working-directory DLL.
struct Callbacks {
    int (*cancellation)(void *);
    void (*progress)(const char *, const char *, double, size_t, size_t, void *);
    void *userData;
    size_t maximumThreads;
    size_t maximumMemoryBytes;
};
struct Bridge {
    HMODULE module{};
    using Query = int (*)(const char *, const char *, const Callbacks *, void **);
    using Read = int (*)(const void *, char *, size_t, size_t *);
    using Error = int (*)(char *, size_t, size_t *);
    using Destroy = void (*)(void *);
    Query query{};
    Read read{};
    Error error{};
    Destroy destroy{};
    Bridge() {
        std::vector<wchar_t> filename(32768);
        const auto length = GetModuleFileNameW(nullptr, filename.data(), DWORD(filename.size()));
        if (!length || length >= filename.size())
            throw std::runtime_error("Cannot resolve Viewer installation");
        const auto dll =
            std::filesystem::path(std::wstring(filename.data(), length)).parent_path() /
            L"dynsyslab_workflow_c_api.dll";
        module =
            LoadLibraryExW(dll.c_str(), nullptr,
                           LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
        if (!module)
            throw std::runtime_error(
                "Phase-space analysis requires the installed DynSysLab solver module");
        query = reinterpret_cast<Query>(GetProcAddress(module, "dsl_workflow_phase_space_json"));
        read = reinterpret_cast<Read>(GetProcAddress(module, "dsl_workflow_result_json"));
        error = reinterpret_cast<Error>(GetProcAddress(module, "dsl_workflow_last_error"));
        destroy = reinterpret_cast<Destroy>(GetProcAddress(module, "dsl_workflow_result_destroy"));
        if (!query || !read || !error || !destroy) {
            FreeLibrary(module);
            module = nullptr;
            throw std::runtime_error("Installed solver does not support the phase-space API");
        }
    }
    ~Bridge() {
        if (module)
            FreeLibrary(module);
    }
};
} // namespace
json::Value callPhaseSidecar(const json::Value &request, const std::filesystem::path &archives,
                             const std::atomic_bool *cancelled) {
    static Bridge bridge;
    Callbacks callbacks{};
    callbacks.userData = const_cast<std::atomic_bool *>(cancelled);
    callbacks.cancellation = [](void *state) {
        return state && static_cast<std::atomic_bool *>(state)->load() ? 1 : 0;
    };
    void *raw = nullptr;
    const auto status =
        bridge.query(json::dump(request).c_str(), archives.u8string().c_str(), &callbacks, &raw);
    struct Result {
        void *handle;
        Bridge::Destroy destroy;
        ~Result() {
            if (handle)
                destroy(handle);
        }
    } result{raw, bridge.destroy};
    size_t size = 0;
    if (status != 0) {
        bridge.error(nullptr, 0, &size);
        std::vector<char> buffer(size ? size : 1);
        bridge.error(buffer.data(), buffer.size(), &size);
        throw std::runtime_error(buffer.data());
    }
    bridge.read(raw, nullptr, 0, &size);
    if (!size || size > 256 * 1024 * 1024)
        throw std::runtime_error("Invalid phase-space response size");
    std::vector<char> buffer(size);
    if (bridge.read(raw, buffer.data(), buffer.size(), &size) != 0)
        throw std::runtime_error("Cannot read phase-space response");
    return json::parse(buffer.data());
}
} // namespace dsl

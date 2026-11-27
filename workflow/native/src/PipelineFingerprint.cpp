#include "PipelineFingerprint.hpp"

#include <dynsyslab/core/BuildInfo.hpp>
#include <dynsyslab/core/FileDigest.hpp>

#include <algorithm>
#include <filesystem>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace dynsyslab::workflow::detail {
namespace {

std::filesystem::path thisModule() {
#ifdef _WIN32
    HMODULE module = nullptr;
    if (GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           reinterpret_cast<LPCWSTR>(&thisModule), &module)) {
        std::wstring path(32768, L'\0');
        const DWORD length = GetModuleFileNameW(module, path.data(), DWORD(path.size()));
        if (length > 0 && length < path.size()) { path.resize(length); return path; }
    }
#endif
    return {};
}

std::string compute() {
    const auto module = thisModule();
    if (module.empty()) return "build:" + core::sha256Text(core::buildIdentityJson());  // no binary to read
    std::vector<std::filesystem::path> files{module};
    std::error_code error;
    for (const auto& entry : std::filesystem::directory_iterator(module.parent_path(), error)) {
        const auto name = entry.path().filename().u8string();
        if (entry.path() != module && name.rfind("dynsyslab", 0) == 0 && entry.path().extension() == L".dll")
            files.push_back(entry.path());
    }
    std::sort(files.begin() + 1, files.end());
    core::Sha256Stream digest;
    for (const auto& file : files) {
        digest.update(file.filename().u8string());
        digest.update(core::sha256File(file));
    }
    return digest.finish();
}

}  // namespace

const std::string& pipelineFingerprint() {
    static const std::string value = compute();
    return value;
}

}  // namespace dynsyslab::workflow::detail

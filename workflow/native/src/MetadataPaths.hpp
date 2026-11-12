#pragma once
#include <filesystem>
#include <string>

namespace dynsyslab::workflow::detail {
// Metadata uses nested content hashes and can exceed Win32's legacy MAX_PATH.
// Keep ordinary paths in public metadata; use this form only for filesystem I/O.
inline std::filesystem::path metadataIoPath(const std::filesystem::path &path) {
#ifdef _WIN32
    auto value = std::filesystem::absolute(path).lexically_normal().make_preferred().wstring();
    if (value.rfind(LR"(\\?\)", 0) == 0)
        return value;
    if (value.rfind(LR"(\\)", 0) == 0)
        return std::wstring(LR"(\\?\UNC\)") + value.substr(2);
    return std::wstring(LR"(\\?\)") + value;
#else
    return path;
#endif
}
} // namespace dynsyslab::workflow::detail

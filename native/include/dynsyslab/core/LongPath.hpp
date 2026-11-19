#pragma once
// Windows extended-length form (\\?\C:\..., \\?\UNC\server\...) of a path, so file APIs accept
// archive and manifest paths past MAX_PATH (260 characters). Other platforms: unchanged.
#include <filesystem>
#include <string>

namespace dynsyslab::core {

[[nodiscard]] inline std::filesystem::path extendedLengthPath(const std::filesystem::path& source) {
#ifdef _WIN32
    auto value = std::filesystem::absolute(source).lexically_normal().make_preferred().wstring();
    if (value.rfind(LR"(\\?\)", 0) == 0) return value;
    if (value.rfind(LR"(\\)", 0) == 0) return std::wstring(LR"(\\?\UNC\)") + value.substr(2);
    return std::wstring(LR"(\\?\)") + value;
#else
    return source;
#endif
}

}  // namespace dynsyslab::core

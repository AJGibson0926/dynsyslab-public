#include <dynsyslab/core/OutputPath.hpp>

#include <dynsyslab/core/LongPath.hpp>

#include <utility>

namespace dynsyslab::core {
namespace {

std::string describe(const std::string& option, const std::filesystem::path& path, const std::string& what,
                     const std::string& reason) {
    return option + ": " + what + " \"" + path.u8string() + "\": " + reason;
}

void createDirectory(const std::filesystem::path& directory, const std::filesystem::path& reported,
                     const std::string& option, const char* what) {
    std::error_code error;
    const auto io = extendedLengthPath(directory);
    if (std::filesystem::is_directory(io, error)) return;
    if (std::filesystem::exists(io, error))
        throw OutputPathError(option, reported, describe(option, reported, what, "a file is in the way"));
    std::filesystem::create_directories(io, error);
    if (error && !std::filesystem::is_directory(io))
        throw OutputPathError(option, reported, describe(option, reported, what, error.message()));
}

}  // namespace

OutputPathError::OutputPathError(std::string option, std::filesystem::path path, const std::string& reason)
    : std::runtime_error(reason), option_(std::move(option)), path_(std::move(path)) {}

std::filesystem::path resolveOutputPath(const std::filesystem::path& path) {
    if (path.empty()) return path;
    std::error_code error;
    const auto absolute = std::filesystem::absolute(path, error);
    return (error ? path : absolute).lexically_normal();
}

std::filesystem::path prepareOutputFile(const std::filesystem::path& path, const std::string& option) {
    if (path.empty()) throw OutputPathError(option, path, option + ": the output path is empty");
    const auto resolved = resolveOutputPath(path);
    if (!resolved.has_filename())
        throw OutputPathError(option, resolved, describe(option, resolved, "expected a file path, got", "a folder"));
    const auto parent = resolved.parent_path();
    if (!parent.empty()) createDirectory(parent, resolved, option, "cannot create parent directory of");
    return resolved;
}

std::filesystem::path prepareOutputDirectory(const std::filesystem::path& path, const std::string& option) {
    if (path.empty()) throw OutputPathError(option, path, option + ": the output folder is empty");
    const auto resolved = resolveOutputPath(path);
    createDirectory(resolved, resolved, option, "cannot create directory");
    return resolved;
}

}  // namespace dynsyslab::core

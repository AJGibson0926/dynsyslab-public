#pragma once
// User-supplied output paths (progress files, manifests, figures, archives). Relative paths
// resolve against the current working directory; a bare file name ("progress.jsonl") is valid.
// Parent folders are created only when missing. Failures name the option and the resolved path.
#include <filesystem>
#include <stdexcept>
#include <string>

namespace dynsyslab::core {

class OutputPathError : public std::runtime_error {
public:
    OutputPathError(std::string option, std::filesystem::path path, const std::string& reason);
    [[nodiscard]] const std::string& option() const noexcept { return option_; }
    [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }

private:
    std::string option_;
    std::filesystem::path path_;
};

// Absolute, normalized form of `path` (empty stays empty).
[[nodiscard]] std::filesystem::path resolveOutputPath(const std::filesystem::path& path);

// Resolves an output file and makes sure its parent folder exists. Returns the resolved path.
// Throws OutputPathError("<option>: cannot create parent directory of \"<path>\": <reason>").
std::filesystem::path prepareOutputFile(const std::filesystem::path& path, const std::string& option = "output");

// Resolves an output folder and creates it when missing.
std::filesystem::path prepareOutputDirectory(const std::filesystem::path& path, const std::string& option = "output");

}  // namespace dynsyslab::core

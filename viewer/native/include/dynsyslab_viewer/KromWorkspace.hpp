#pragma once
#include "dynsyslab_viewer/Model.hpp"
#include <optional>
namespace dsl {
struct KromWorkspaceOptions {
    // Local retained-mode indices. Unspecified selects the first mode; an empty list selects none.
    std::optional<std::vector<std::size_t>> modes;
    bool allModes = false;
    std::string complexDisplay = "real-imaginary";
};
Payload buildKromWorkspace(const std::filesystem::path& source,
    const KromWorkspaceOptions& options = {});
}

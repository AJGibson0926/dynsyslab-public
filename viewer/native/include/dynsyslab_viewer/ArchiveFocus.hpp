#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <filesystem>

namespace dsl {

// Applies a Viewer-only, fail-open allowlist. Scientific archives are never
// moved or rewritten; hierarchy ancestors required to display allowed leaves
// are retained automatically and links to masked objects are removed.
void applyArchiveFocus(Catalog& catalog,
                       const std::filesystem::path& manifest,
                       const std::filesystem::path& projectRoot);

}  // namespace dsl

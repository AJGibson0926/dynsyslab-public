#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <filesystem>
#include <optional>
#include <string>

namespace dsl::compute {

struct KoopmanPipelineSettings {
    std::string hdmd;
    std::string sindy;
};

// Edits the standard scalar H-DMD -> supplied-spectrum SINDy policy. Unknown
// advanced keys are preserved in both returned JSON objects.
std::optional<KoopmanPipelineSettings> showKoopmanPipelineDialog(
    HWND owner, const std::filesystem::path& source,
    const std::string& currentHdmd, const std::string& currentSindy,
    const dsl::theme::Fonts& fonts);

}  // namespace dsl::compute

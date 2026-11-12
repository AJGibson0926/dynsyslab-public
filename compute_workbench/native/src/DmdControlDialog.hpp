#pragma once
#include "dynsyslab_viewer/Theme.hpp"
#include <optional>
#include <string>

namespace dsl::compute {
// Presentation only: edits the public postprocess method=settings contract.
std::optional<std::string> showDmdControlDialog(HWND owner,
    const std::string& current,const dsl::theme::Fonts& fonts,bool extended = false);
}

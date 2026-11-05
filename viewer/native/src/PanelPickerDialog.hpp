#pragma once

#include "dynsyslab_viewer/PanelCatalog.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <vector>

namespace dsl {

std::optional<std::size_t> showPanelPickerDialog(
    HWND owner, const std::vector<PanelChoice>& choices,
    const theme::Fonts& fonts);

}  // namespace dsl

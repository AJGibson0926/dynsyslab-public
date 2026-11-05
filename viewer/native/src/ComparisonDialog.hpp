#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

#include <optional>
#include <string>

namespace dsl {

// Selects a second archive object that supports the currently visible view.
std::optional<int> showComparisonDialog(HWND owner, const Catalog& catalog,
                                        int currentIndex, const std::string& view,
                                        const theme::Fonts& fonts);

}  // namespace dsl

#pragma once

#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>

namespace dsl {

// Edits the visible page's axes without mutating archived scientific data.
bool showAxisRangeDialog(HWND owner, PlotCanvas& canvas, const theme::Fonts& fonts);

}  // namespace dsl

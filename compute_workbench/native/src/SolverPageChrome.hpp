#pragma once

#include <windows.h>

namespace dsl::compute {

void drawSolverFolderTab(const DRAWITEMSTRUCT& item, bool active, HFONT font);
void drawSolverEditorFrame(HDC dc, const RECT& bounds);

}  // namespace dsl::compute

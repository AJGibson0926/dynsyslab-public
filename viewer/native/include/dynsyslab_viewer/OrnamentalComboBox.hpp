#pragma once

#include <windows.h>

namespace dsl {

// Applies DynSysLab's owner-painted closed-face treatment while retaining the
// native ComboBox window, popup, keyboard navigation, focus, and accessibility.
void styleOrnamentalComboBox(HWND comboBox);

}  // namespace dsl

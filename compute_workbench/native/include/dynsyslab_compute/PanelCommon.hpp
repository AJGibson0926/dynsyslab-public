#pragma once

#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>
#include <string>
#include <vector>

namespace dsl::compute {

HWND makeLabel(HWND parent, const wchar_t* text, HFONT font);
HWND makeCenteredLabel(HWND parent, const wchar_t* text, HFONT font);
HWND makeEdit(HWND parent, int id, const wchar_t* text, HFONT font, bool multiline = false);
// Read-mostly numeric/value field: paints as compact scientific text and only
// becomes an edit control while the user is actively changing it.
HWND makeInlineValue(HWND parent, int id, const wchar_t* text, HFONT font);
// Scalar variant used beside mathematical labels. Its resting state is
// borderless and uses the same scientific-number notation as the parameter book.
HWND makeScientificInlineValue(HWND parent, int id, const wchar_t* text, HFONT font);
// Compact TeX-like label supporting Greek commands and braced subscripts.
HWND makeMathLabel(HWND parent, const wchar_t* latex, HFONT body, HFONT caption);
HWND makeCombo(HWND parent, int id, HFONT font);
HWND makeButton(HWND parent, int id, const wchar_t* text, HFONT font);
HWND makeCheck(HWND parent, int id, const wchar_t* text, HFONT font, bool checked = false);
HWND makeProgress(HWND parent, HFONT font);
std::wstring windowText(HWND window);
std::string utf8Text(HWND window);
std::string jsonEscape(const std::string& value);
std::string jsonObject(const std::string& value, const char* fallback = "{}");
std::string selectedText(HWND combo);
void setComboItems(HWND combo, const std::vector<std::string>& values, int selection = 0);
bool checked(HWND window);
int scalePixels(HWND window, int logicalPixels);
void moveFieldControl(HWND control, int x, int y, int width, int visibleHeight,
                      int dropDownRows = 8);

}  // namespace dsl::compute

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Model.hpp"

#include <algorithm>
#include <commctrl.h>
#include <sstream>
#include <windowsx.h>

namespace dsl::compute {

namespace {
HWND create(HWND parent, const wchar_t* type, const wchar_t* text, DWORD style,
            DWORD extended, int id, HFONT font) {
    HWND result = CreateWindowExW(extended, type, text, WS_CHILD | WS_VISIBLE | style,
        0, 0, 0, 0, parent, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr), nullptr);
    dsl::theme::styleControl(result, font); return result;
}
}

HWND makeLabel(HWND parent, const wchar_t* text, HFONT font) {
    return create(parent, L"STATIC", text, SS_LEFT, 0, 0, font);
}
HWND makeCenteredLabel(HWND parent, const wchar_t* text, HFONT font) {
    return create(parent, L"STATIC", text, SS_CENTER, 0, 0, font);
}
HWND makeEdit(HWND parent, int id, const wchar_t* text, HFONT font, bool multiline) {
    DWORD style = ES_AUTOHSCROLL | WS_TABSTOP;
    if (multiline) style = ES_MULTILINE | ES_AUTOVSCROLL | ES_WANTRETURN |
        WS_VSCROLL | WS_TABSTOP;
    return create(parent, L"EDIT", text, style, WS_EX_CLIENTEDGE, id, font);
}
HWND makeCombo(HWND parent, int id, HFONT font) {
    return create(parent, L"COMBOBOX", L"", CBS_DROPDOWNLIST | WS_VSCROLL | WS_TABSTOP,
                  0, id, font);
}
HWND makeButton(HWND parent, int id, const wchar_t* text, HFONT font) {
    return create(parent, L"BUTTON", text, BS_OWNERDRAW | WS_TABSTOP, 0, id, font);
}
HWND makeCheck(HWND parent, int id, const wchar_t* text, HFONT font, bool selected) {
    HWND result = create(parent, L"BUTTON", text, BS_AUTOCHECKBOX | WS_TABSTOP, 0, id, font);
    Button_SetCheck(result, selected ? BST_CHECKED : BST_UNCHECKED); return result;
}
HWND makeProgress(HWND parent, HFONT font) {
    HWND result = create(parent, PROGRESS_CLASSW, L"", PBS_SMOOTH | PBS_MARQUEE,
                         0, 0, font);
    SendMessageW(result, PBM_SETRANGE32, 0, 1000);
    return result;
}
std::wstring windowText(HWND window) {
    const int count = GetWindowTextLengthW(window);
    std::wstring result(static_cast<std::size_t>(count + 1), L'\0');
    GetWindowTextW(window, result.data(), count + 1); result.resize(static_cast<std::size_t>(count)); return result;
}
std::string utf8Text(HWND window) { return dsl::narrow(windowText(window)); }
std::string jsonEscape(const std::string& value) {
    std::ostringstream result;
    for (const unsigned char character : value) {
        switch (character) {
        case '\\': result << "\\\\"; break; case '"': result << "\\\""; break;
        case '\n': result << "\\n"; break; case '\r': result << "\\r"; break;
        case '\t': result << "\\t"; break;
        default: if (character >= 0x20) result << static_cast<char>(character); break;
        }
    }
    return result.str();
}
std::string jsonObject(const std::string& value, const char* fallback) {
    const auto first = value.find_first_not_of(" \t\r\n");
    return first == std::string::npos ? fallback : value;
}
std::string selectedText(HWND combo) {
    const LRESULT selected = SendMessageW(combo, CB_GETCURSEL, 0, 0);
    if (selected == CB_ERR) return {};
    const LRESULT count = SendMessageW(combo, CB_GETLBTEXTLEN, selected, 0);
    std::wstring value(static_cast<std::size_t>(count + 1), L'\0');
    SendMessageW(combo, CB_GETLBTEXT, selected, reinterpret_cast<LPARAM>(value.data()));
    value.resize(static_cast<std::size_t>(count)); return dsl::narrow(value);
}
void setComboItems(HWND combo, const std::vector<std::string>& values, int selection) {
    SendMessageW(combo, CB_RESETCONTENT, 0, 0);
    for (const auto& value : values) SendMessageW(combo, CB_ADDSTRING, 0,
        reinterpret_cast<LPARAM>(dsl::wide(value).c_str()));
    if (!values.empty()) SendMessageW(combo, CB_SETCURSEL,
        static_cast<WPARAM>(std::clamp(selection, 0, static_cast<int>(values.size()) - 1)), 0);
}
bool checked(HWND window) { return Button_GetCheck(window) == BST_CHECKED; }
int scalePixels(HWND window, int logicalPixels) {
    return dsl::theme::scalePixels(window, logicalPixels);
}
void moveFieldControl(HWND control, int x, int y, int width, int visibleHeight,
                      int dropDownRows) {
    wchar_t className[32]{};
    GetClassNameW(control, className, static_cast<int>(std::size(className)));
    int windowHeight = visibleHeight;
    if (lstrcmpiW(className, L"ComboBox") == 0) {
        const LRESULT countResult = SendMessageW(control, CB_GETCOUNT, 0, 0);
        const int count = countResult == CB_ERR ? 1 : static_cast<int>(countResult);
        const int rows = std::clamp(count, 1, std::max(1, dropDownRows));
        const LRESULT itemHeightResult = SendMessageW(control, CB_GETITEMHEIGHT, 0, 0);
        const int itemHeight = itemHeightResult == CB_ERR
            ? scalePixels(control, 18) : static_cast<int>(itemHeightResult);
        windowHeight += rows * itemHeight + scalePixels(control, 4);
    }
    MoveWindow(control, x, y, width, windowHeight, TRUE);
}

}  // namespace dsl::compute

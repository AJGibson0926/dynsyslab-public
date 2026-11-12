#include "dynsyslab_viewer/Theme.hpp"

#include <commctrl.h>

namespace dsl::theme {
namespace {
constexpr UINT_PTR kToggleSubclass = 1;

void paintToggle(HWND window, HDC target) {
    RECT area{}; GetClientRect(window, &area);
    const bool enabled = IsWindowEnabled(window) != FALSE;
    const bool selected = SendMessageW(window, BM_GETCHECK, 0, 0) == BST_CHECKED;
    const bool pressed = (SendMessageW(window, BM_GETSTATE, 0, 0) & BST_PUSHED) != 0;
    const auto fillColor = !enabled ? controlBackground
        : pressed ? accent : selected ? selection : panelRaised;
    const auto fill = CreateSolidBrush(fillColor);
    FillRect(target, &area, fill); DeleteObject(fill);
    const auto outline = CreateSolidBrush(selected && enabled ? accent : border);
    FrameRect(target, &area, outline); DeleteObject(outline);
    SetBkMode(target, TRANSPARENT);
    SetTextColor(target, !enabled ? disabled
        : pressed ? panelRaised : selected ? accent : text);
    const auto oldFont = SelectObject(target,
        reinterpret_cast<HFONT>(SendMessageW(window, WM_GETFONT, 0, 0)));
    const auto length = GetWindowTextLengthW(window);
    std::wstring label(static_cast<std::size_t>(length + 1), L'\0');
    GetWindowTextW(window, label.data(), length + 1);
    DrawTextW(target, label.c_str(), length, &area,
        DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    if (GetFocus() == window && !(SendMessageW(window, WM_QUERYUISTATE, 0, 0) & UISF_HIDEFOCUS)) {
        InflateRect(&area, -3, -3);
        const auto focusPen = CreatePen(PS_SOLID, 2, focus);
        const auto oldPen = SelectObject(target, focusPen);
        const auto oldBrush = SelectObject(target, GetStockObject(NULL_BRUSH));
        Rectangle(target, area.left, area.top, area.right, area.bottom);
        SelectObject(target, oldBrush); SelectObject(target, oldPen);
        DeleteObject(focusPen);
    }
    SelectObject(target, oldFont);
}

LRESULT CALLBACK toggleProcedure(HWND window, UINT message, WPARAM wparam,
                                 LPARAM lparam, UINT_PTR, DWORD_PTR) {
    if (message == WM_PAINT) {
        PAINTSTRUCT state{}; const auto target = BeginPaint(window, &state);
        paintToggle(window, target); EndPaint(window, &state); return 0;
    }
    if (message == WM_PRINTCLIENT) { paintToggle(window, reinterpret_cast<HDC>(wparam)); return 0; }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_NCDESTROY) RemoveWindowSubclass(window, toggleProcedure, kToggleSubclass);
    const auto result = DefSubclassProc(window, message, wparam, lparam);
    if (message == BM_SETCHECK || message == BM_SETSTATE || message == WM_SETFOCUS
        || message == WM_KILLFOCUS || message == WM_ENABLE || message == WM_UPDATEUISTATE)
        InvalidateRect(window, nullptr, FALSE);
    return result;
}
}

void styleToggleButton(HWND control, HFONT font) {
    styleControl(control, font);
    SetWindowSubclass(control, toggleProcedure, kToggleSubclass, 0);
}

}  // namespace dsl::theme

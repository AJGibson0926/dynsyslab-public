#include "dynsyslab_viewer/OrnamentalComboBox.hpp"

#include "dynsyslab_viewer/GdiBackBuffer.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <algorithm>
#include <commctrl.h>
#include <iterator>
#include <string>
#include <uxtheme.h>

namespace dsl {
namespace {

constexpr UINT_PTR kSubclassId = 0xd51c;

struct ComboState {
    GdiBackBuffer backBuffer;
    bool hovered{};
};

int scaled(HWND window, int value) {
    const UINT dpi = GetDpiForWindow(window);
    return MulDiv(value, dpi ? static_cast<int>(dpi) : 96, 96);
}

void fill(HDC dc, const RECT& area, COLORREF color) {
    HBRUSH brush = CreateSolidBrush(color);
    FillRect(dc, &area, brush); DeleteObject(brush);
}

void frame(HDC dc, const RECT& area, COLORREF color) {
    HBRUSH brush = CreateSolidBrush(color);
    FrameRect(dc, &area, brush); DeleteObject(brush);
}

std::wstring selectedText(HWND combo) {
    const LRESULT selected = SendMessageW(combo, CB_GETCURSEL, 0, 0);
    if (selected == CB_ERR) return {};
    const LRESULT length = SendMessageW(combo, CB_GETLBTEXTLEN, selected, 0);
    if (length == CB_ERR) return {};
    std::wstring result(static_cast<std::size_t>(length + 1), L'\0');
    SendMessageW(combo, CB_GETLBTEXT, selected,
                 reinterpret_cast<LPARAM>(result.data()));
    result.resize(static_cast<std::size_t>(length));
    return result;
}

void diamond(HDC dc, int centerX, int centerY, int radius) {
    POINT points[]{{centerX, centerY - radius}, {centerX + radius, centerY},
        {centerX, centerY + radius}, {centerX - radius, centerY},
        {centerX, centerY - radius}};
    Polyline(dc, points, static_cast<int>(std::size(points)));
}

void paintFiligree(HWND combo, ComboState& state, HDC dc) {
    RECT area{}; GetClientRect(combo, &area);
    const bool enabled = IsWindowEnabled(combo) != FALSE;
    const bool focused = GetFocus() == combo;
    const bool dropped = SendMessageW(combo, CB_GETDROPPEDSTATE, 0, 0) != FALSE;
    const bool active = enabled && (focused || dropped || state.hovered);
    fill(dc, area, enabled ? theme::panelRaised : theme::controlBackground);

    const int buttonWidth = std::min(static_cast<int>(area.right), scaled(combo, 29));
    RECT button{area.right - buttonWidth, area.top, area.right, area.bottom};
    fill(dc, button, active ? theme::selection : theme::controlBackground);
    HPEN accent = CreatePen(PS_SOLID, 1,
        active ? theme::focus : theme::borderStrong);
    HGDIOBJ previous = SelectObject(dc, accent);
    MoveToEx(dc, button.left, button.top + 1, nullptr);
    LineTo(dc, button.left, button.bottom - 1);

    const int centerX = (button.left + button.right) / 2;
    const int centerY = (button.top + button.bottom) / 2;
    const int radius = scaled(combo, 3);
    POINT arrow[]{{centerX - radius, centerY - 1},
        {centerX, centerY + radius - 1},
        {centerX + radius, centerY - 1}};
    Polyline(dc, arrow, static_cast<int>(std::size(arrow)));
    diamond(dc, centerX, button.top + scaled(combo, 4), scaled(combo, 1));
    diamond(dc, centerX, button.bottom - scaled(combo, 4), scaled(combo, 1));
    const int curl = scaled(combo, 4);
    POINT flourish[]{{button.left + 3, centerY},
        {button.left + curl, centerY - curl},
        {centerX - curl, centerY - curl}, {centerX - radius, centerY - 1}};
    PolyBezier(dc, flourish, static_cast<DWORD>(std::size(flourish)));
    SelectObject(dc, previous); DeleteObject(accent);

    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, enabled ? theme::text : theme::disabled);
    if (HFONT font = reinterpret_cast<HFONT>(SendMessageW(combo, WM_GETFONT, 0, 0)))
        SelectObject(dc, font);
    RECT text{area.left + scaled(combo, 8), area.top,
        button.left - scaled(combo, 6), area.bottom};
    const std::wstring value = selectedText(combo);
    DrawTextW(dc, value.c_str(), -1, &text,
        DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_END_ELLIPSIS | DT_NOPREFIX);
    frame(dc, area, active ? theme::blue : theme::border);
    if (focused) {
        RECT inner = area; InflateRect(&inner, -scaled(combo, 2), -scaled(combo, 2));
        frame(dc, inner, theme::focus);
    }
}

LRESULT CALLBACK comboProcedure(HWND combo, UINT message, WPARAM wparam,
                                LPARAM lparam, UINT_PTR, DWORD_PTR reference) {
    auto* state = reinterpret_cast<ComboState*>(reference);
    if (message == WM_PAINT && state) {
        PAINTSTRUCT paint{}; HDC target = BeginPaint(combo, &paint);
        RECT area{}; GetClientRect(combo, &area);
        HDC frameDc = state->backBuffer.target(target,
            static_cast<int>(area.right), static_cast<int>(area.bottom));
        paintFiligree(combo, *state, frameDc);
        if (frameDc != target) state->backBuffer.present(target, paint.rcPaint);
        EndPaint(combo, &paint); return 0;
    }
    if (message == WM_PRINTCLIENT && state) {
        paintFiligree(combo, *state, reinterpret_cast<HDC>(wparam)); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_MOUSEMOVE && state && !state->hovered) {
        state->hovered = true;
        TRACKMOUSEEVENT tracking{sizeof(tracking), TME_LEAVE, combo, 0};
        TrackMouseEvent(&tracking); InvalidateRect(combo, nullptr, FALSE);
    }
    if (message == WM_MOUSELEAVE && state) {
        state->hovered = false; InvalidateRect(combo, nullptr, FALSE); return 0;
    }
    if (message == WM_NCDESTROY) {
        const LRESULT result = DefSubclassProc(combo, message, wparam, lparam);
        delete state; return result;
    }
    const LRESULT result = DefSubclassProc(combo, message, wparam, lparam);
    if (message == WM_SETFOCUS || message == WM_KILLFOCUS
        || message == WM_ENABLE || message == WM_SETFONT
        || message == CB_SETCURSEL || message == CB_RESETCONTENT
        || message == CB_ADDSTRING || message == CB_SHOWDROPDOWN)
        InvalidateRect(combo, nullptr, FALSE);
    return result;
}

}  // namespace

void styleOrnamentalComboBox(HWND comboBox) {
    if (!comboBox) return;
    const LONG_PTR style = GetWindowLongPtrW(comboBox, GWL_STYLE);
    if ((style & 0x3) != CBS_DROPDOWNLIST) return;
    DWORD_PTR existing = 0;
    if (GetWindowSubclass(comboBox, comboProcedure, kSubclassId, &existing)) return;

    SetWindowTheme(comboBox, L"", L"");
    const LONG_PTR extended = GetWindowLongPtrW(comboBox, GWL_EXSTYLE);
    SetWindowLongPtrW(comboBox, GWL_EXSTYLE, extended & ~WS_EX_CLIENTEDGE);
    SetWindowPos(comboBox, nullptr, 0, 0, 0, 0,
        SWP_NOMOVE | SWP_NOSIZE | SWP_NOZORDER | SWP_NOACTIVATE
        | SWP_FRAMECHANGED | SWP_NOREDRAW);
    auto* state = new ComboState;
    if (!SetWindowSubclass(comboBox, comboProcedure, kSubclassId,
                           reinterpret_cast<DWORD_PTR>(state))) {
        delete state; return;
    }
    COMBOBOXINFO info{sizeof(info)};
    if (GetComboBoxInfo(comboBox, &info) && info.hwndList)
        SetWindowTheme(info.hwndList, L"Explorer", nullptr);
    InvalidateRect(comboBox, nullptr, FALSE);
}

}  // namespace dsl

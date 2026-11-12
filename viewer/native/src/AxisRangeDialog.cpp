#include "AxisRangeDialog.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cwchar>
#include <cwctype>
#include <string>

namespace dsl {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabAxisRangeDialog";
constexpr int kPanelId = 501;
constexpr int kFirstEditId = 510;
constexpr int kAutoId = 520;
constexpr int kApplyId = 521;
constexpr int kCloseId = 522;

struct DialogState {
    HWND owner = nullptr;
    HWND window = nullptr;
    HWND panel = nullptr;
    std::array<HWND, 4> edits{};
    PlotCanvas* canvas = nullptr;
    const theme::Fonts* fonts = nullptr;
    std::vector<AxisPanelState> panels;
    bool changed = false;
};

HWND makeControl(HWND parent, const wchar_t* type, const wchar_t* text,
                 DWORD style, int id, int x, int y, int width, int height,
                 HFONT font) {
    HWND result = CreateWindowExW(type == std::wstring(L"EDIT") ? WS_EX_CLIENTEDGE : 0,
        type, text, WS_CHILD | WS_VISIBLE | style, x, y, width, height, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(result, font);
    return result;
}

void putNumber(HWND edit, double value) {
    wchar_t buffer[64]{};
    swprintf_s(buffer, L"%.12g", value);
    SetWindowTextW(edit, buffer);
}

void loadPanel(DialogState& state) {
    const auto selected = static_cast<std::size_t>(std::max<LRESULT>(
        0, SendMessageW(state.panel, CB_GETCURSEL, 0, 0)));
    if (selected >= state.panels.size()) return;
    const auto& panel = state.panels[selected];
    const AxisLimits limits = panel.manual.value_or(panel.automatic);
    putNumber(state.edits[0], limits.xMin); putNumber(state.edits[1], limits.xMax);
    putNumber(state.edits[2], limits.yMin); putNumber(state.edits[3], limits.yMax);
}

bool readNumber(HWND edit, double& value) {
    wchar_t buffer[96]{};
    GetWindowTextW(edit, buffer, 96);
    wchar_t* end = nullptr;
    value = std::wcstod(buffer, &end);
    while (end && std::iswspace(*end)) ++end;
    return end && end != buffer && *end == L'\0' && std::isfinite(value);
}

void apply(DialogState& state) {
    const auto selected = static_cast<std::size_t>(std::max<LRESULT>(
        0, SendMessageW(state.panel, CB_GETCURSEL, 0, 0)));
    AxisLimits limits;
    if (selected >= state.panels.size() || !readNumber(state.edits[0], limits.xMin)
        || !readNumber(state.edits[1], limits.xMax)
        || !readNumber(state.edits[2], limits.yMin)
        || !readNumber(state.edits[3], limits.yMax)
        || limits.xMin >= limits.xMax || limits.yMin >= limits.yMax) {
        MessageBoxW(state.window,
            L"Enter four finite limits with minimum strictly below maximum.",
            L"Invalid axis range", MB_OK | MB_ICONWARNING);
        return;
    }
    if (!state.canvas->setAxisLimits(selected, limits)) {
        MessageBoxW(state.window,L"Logarithmic Y limits must be positive. Use Symmetric log for ranges including zero or negative values.",
            L"Invalid logarithmic range",MB_OK | MB_ICONWARNING); return;
    }
    state.panels[selected].manual = limits;
    state.changed = true;
}

void useAutomatic(DialogState& state) {
    const auto selected = static_cast<std::size_t>(std::max<LRESULT>(
        0, SendMessageW(state.panel, CB_GETCURSEL, 0, 0)));
    if (selected >= state.panels.size()) return;
    state.canvas->setAxisLimits(selected, std::nullopt);
    state.panels[selected].manual.reset();
    loadPanel(state); state.changed = true;
}

void createChildren(DialogState& state) {
    const auto body = state.fonts->body;
    makeControl(state.window, L"STATIC", L"Panel", SS_LEFT, 0, 20, 17, 80, 20, body);
    state.panel = makeControl(state.window, L"COMBOBOX", L"",
        CBS_DROPDOWNLIST | WS_VSCROLL, kPanelId, 105, 13, 295, 190, body);
    for (const auto& panel : state.panels)
        SendMessageW(state.panel, CB_ADDSTRING, 0,
            reinterpret_cast<LPARAM>(panel.title.c_str()));
    SendMessageW(state.panel, CB_SETCURSEL, 0, 0);
    constexpr std::array<const wchar_t*, 4> labels{L"X minimum", L"X maximum",
                                                   L"Y minimum", L"Y maximum"};
    for (int index = 0; index < 4; ++index) {
        const int row = index / 2, column = index % 2;
        const int x = 20 + column * 194, y = 58 + row * 64;
        makeControl(state.window, L"STATIC", labels[index], SS_LEFT, 0,
                    x, y, 90, 20, body);
        state.edits[static_cast<std::size_t>(index)] = makeControl(state.window,
            L"EDIT", L"", ES_AUTOHSCROLL, kFirstEditId + index,
            x, y + 23, 180, 25, state.fonts->mono);
    }
    makeControl(state.window, L"BUTTON", L"Use auto", BS_PUSHBUTTON,
                kAutoId, 110, 197, 90, 29, body);
    makeControl(state.window, L"BUTTON", L"Apply", BS_DEFPUSHBUTTON,
                kApplyId, 207, 197, 90, 29, body);
    makeControl(state.window, L"BUTTON", L"Close", BS_PUSHBUTTON,
                kCloseId, 304, 197, 90, 29, body);
    loadPanel(state);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<DialogState*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        state->window = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_CREATE) { createChildren(*state); return 0; }
    if (message == WM_COMMAND) {
        const int id = LOWORD(wparam);
        if (id == kPanelId && HIWORD(wparam) == CBN_SELCHANGE) { loadPanel(*state); return 0; }
        if (id == kAutoId) { useAutomatic(*state); return 0; }
        if (id == kApplyId || id == IDOK) { apply(*state); return 0; }
        if (id == kCloseId || id == IDCANCEL) { DestroyWindow(window); return 0; }
    }
    if (message == WM_CLOSE) { DestroyWindow(window); return 0; }
    if (message == WM_CTLCOLORSTATIC || message == WM_CTLCOLOREDIT) {
        static HBRUSH background = CreateSolidBrush(theme::panel);
        HDC dc = reinterpret_cast<HDC>(wparam); SetBkColor(dc, theme::panel);
        SetTextColor(dc, theme::text); return reinterpret_cast<LRESULT>(background);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerDialogClass() {
    WNDCLASSEXW type{sizeof(type)};
    type.lpfnWndProc = procedure;
    type.hInstance = GetModuleHandleW(nullptr);
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(theme::panel);
    type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

bool showAxisRangeDialog(HWND owner, PlotCanvas& canvas, const theme::Fonts& fonts) {
    DialogState state;
    state.owner = owner; state.canvas = &canvas; state.fonts = &fonts;
    state.panels = canvas.currentAxisPanels();
    if (state.panels.empty() || !registerDialogClass()) return false;
    RECT ownerArea{}; GetWindowRect(owner, &ownerArea);
    const int width = 430, height = 285;
    const int x = ownerArea.left + (ownerArea.right - ownerArea.left - width) / 2;
    const int y = ownerArea.top + (ownerArea.bottom - ownerArea.top - height) / 2;
    state.window = CreateWindowExW(WS_EX_DLGMODALFRAME, kClassName,
        L"Plot axis ranges", WS_CAPTION | WS_SYSMENU | WS_POPUP,
        x, y, width, height, owner, nullptr, GetModuleHandleW(nullptr), &state);
    if (!state.window) return false;
    theme::configureTitleBar(state.window);
    EnableWindow(owner, FALSE); ShowWindow(state.window, SW_SHOW); UpdateWindow(state.window);
    MSG event{};
    while (IsWindow(state.window) && GetMessageW(&event, nullptr, 0, 0) > 0) {
        if (!IsDialogMessageW(state.window, &event)) {
            TranslateMessage(&event); DispatchMessageW(&event);
        }
    }
    EnableWindow(owner, TRUE); SetForegroundWindow(owner);
    return state.changed;
}

}  // namespace dsl

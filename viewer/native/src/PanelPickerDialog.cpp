#include "PanelPickerDialog.hpp"

#include <algorithm>

namespace dsl {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabPanelPickerDialog";
constexpr int kFamilyId = 701;
constexpr int kPanelId = 702;
constexpr int kAddId = 703;
constexpr int kCancelId = 704;

struct DialogState {
    HWND owner{};
    HWND window{};
    HWND family{};
    HWND panel{};
    HWND detail{};
    const theme::Fonts* fonts{};
    const std::vector<PanelChoice>* choices{};
    std::vector<std::string> families;
    std::vector<std::size_t> filtered;
    std::optional<std::size_t> selection;
};

int scaled(HWND window, int value) {
    return theme::scalePixels(window, value);
}

HWND control(DialogState& state, const wchar_t* type, const wchar_t* text,
             DWORD style, int id, int x, int y, int width, int height,
             HFONT font) {
    HWND result = CreateWindowExW(0, type, text, WS_CHILD | WS_VISIBLE | style,
        scaled(state.window, x), scaled(state.window, y),
        scaled(state.window, width), scaled(state.window, height), state.window,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(result, font);
    return result;
}

void updateDetail(DialogState& state) {
    const LRESULT selected = SendMessageW(state.panel, CB_GETCURSEL, 0, 0);
    if (selected == CB_ERR || static_cast<std::size_t>(selected) >= state.filtered.size()) {
        SetWindowTextW(state.detail, L""); return;
    }
    const auto& choice = state.choices->at(
        state.filtered[static_cast<std::size_t>(selected)]);
    SetWindowTextW(state.detail, wide(choice.detail).c_str());
}

void populatePanels(DialogState& state) {
    SendMessageW(state.panel, CB_RESETCONTENT, 0, 0);
    state.filtered.clear();
    const LRESULT selected = SendMessageW(state.family, CB_GETCURSEL, 0, 0);
    if (selected == CB_ERR || static_cast<std::size_t>(selected) >= state.families.size())
        return;
    const auto& family = state.families[static_cast<std::size_t>(selected)];
    for (std::size_t index = 0; index < state.choices->size(); ++index) {
        if (state.choices->at(index).family != family) continue;
        state.filtered.push_back(index);
        const auto label = wide(state.choices->at(index).label);
        SendMessageW(state.panel, CB_ADDSTRING, 0,
            reinterpret_cast<LPARAM>(label.c_str()));
    }
    SendMessageW(state.panel, CB_SETCURSEL, 0, 0);
    updateDetail(state);
}

void createChildren(DialogState& state) {
    control(state, L"STATIC", L"PLOT FAMILY", SS_LEFT, 0,
        20, 18, 150, 20, state.fonts->caption);
    state.family = control(state, L"COMBOBOX", L"",
        CBS_DROPDOWNLIST | WS_VSCROLL, kFamilyId,
        20, 40, 430, 230, state.fonts->body);
    control(state, L"STATIC", L"INFORMATION TO DISPLAY", SS_LEFT, 0,
        20, 82, 220, 20, state.fonts->caption);
    state.panel = control(state, L"COMBOBOX", L"",
        CBS_DROPDOWNLIST | WS_VSCROLL, kPanelId,
        20, 104, 430, 250, state.fonts->body);
    state.detail = control(state, L"STATIC", L"", SS_LEFT | SS_NOPREFIX, 0,
        20, 150, 430, 46, state.fonts->body);
    control(state, L"STATIC",
        L"The new panel is session-local. Drag its header to reorder it; "
        L"drag its lower-right corner to resize it.",
        SS_LEFT | SS_NOPREFIX, 0, 20, 202, 430, 42, state.fonts->caption);
    control(state, L"BUTTON", L"Add panel", BS_DEFPUSHBUTTON,
        kAddId, 260, 258, 90, 30, state.fonts->body);
    control(state, L"BUTTON", L"Cancel", BS_PUSHBUTTON,
        kCancelId, 360, 258, 90, 30, state.fonts->body);

    for (const auto& choice : *state.choices)
        if (std::find(state.families.begin(), state.families.end(), choice.family)
            == state.families.end()) state.families.push_back(choice.family);
    for (const auto& family : state.families) {
        const auto label = wide(family);
        SendMessageW(state.family, CB_ADDSTRING, 0,
            reinterpret_cast<LPARAM>(label.c_str()));
    }
    SendMessageW(state.family, CB_SETCURSEL, 0, 0);
    populatePanels(state);
}

void accept(DialogState& state) {
    const LRESULT selected = SendMessageW(state.panel, CB_GETCURSEL, 0, 0);
    if (selected != CB_ERR && static_cast<std::size_t>(selected) < state.filtered.size())
        state.selection = state.filtered[static_cast<std::size_t>(selected)];
    DestroyWindow(state.window);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<DialogState*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<DialogState*>(
            reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        state->window = window;
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_CREATE) { createChildren(*state); return 0; }
    if (message == WM_COMMAND) {
        const int id = LOWORD(wparam);
        if (id == kFamilyId && HIWORD(wparam) == CBN_SELCHANGE) {
            populatePanels(*state); return 0;
        }
        if (id == kPanelId && HIWORD(wparam) == CBN_SELCHANGE) {
            updateDetail(*state); return 0;
        }
        if (id == kAddId || id == IDOK) { accept(*state); return 0; }
        if (id == kCancelId || id == IDCANCEL) { DestroyWindow(window); return 0; }
    }
    if (message == WM_CLOSE) { DestroyWindow(window); return 0; }
    if (message == WM_CTLCOLORSTATIC) {
        static HBRUSH background = CreateSolidBrush(theme::panel);
        HDC dc = reinterpret_cast<HDC>(wparam);
        SetBkColor(dc, theme::panel); SetTextColor(dc, theme::text);
        return reinterpret_cast<LRESULT>(background);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerDialogClass() {
    WNDCLASSEXW type{sizeof(type)};
    type.lpfnWndProc = procedure; type.hInstance = GetModuleHandleW(nullptr);
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(theme::panel);
    type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

std::optional<std::size_t> showPanelPickerDialog(
    HWND owner, const std::vector<PanelChoice>& choices,
    const theme::Fonts& fonts) {
    if (choices.empty() || !registerDialogClass()) return std::nullopt;
    DialogState state;
    state.owner = owner; state.choices = &choices; state.fonts = &fonts;
    RECT ownerArea{}; GetWindowRect(owner, &ownerArea);
    const int width = theme::scalePixels(owner, 490);
    const int height = theme::scalePixels(owner, 350);
    const int x = ownerArea.left + (ownerArea.right - ownerArea.left - width) / 2;
    const int y = ownerArea.top + (ownerArea.bottom - ownerArea.top - height) / 2;
    state.window = CreateWindowExW(WS_EX_DLGMODALFRAME | WS_EX_CONTROLPARENT,
        kClassName, L"Add Viewer panel", WS_CAPTION | WS_SYSMENU | WS_POPUP,
        x, y, width, height, owner, nullptr, GetModuleHandleW(nullptr), &state);
    if (!state.window) return std::nullopt;
    theme::configureTitleBar(state.window);
    EnableWindow(owner, FALSE); ShowWindow(state.window, SW_SHOW); UpdateWindow(state.window);
    MSG event{};
    while (IsWindow(state.window) && GetMessageW(&event, nullptr, 0, 0) > 0) {
        if (!IsDialogMessageW(state.window, &event)) {
            TranslateMessage(&event); DispatchMessageW(&event);
        }
    }
    EnableWindow(owner, TRUE); SetForegroundWindow(owner);
    return state.selection;
}

}  // namespace dsl

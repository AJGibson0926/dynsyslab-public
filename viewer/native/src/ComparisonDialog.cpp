#include "ComparisonDialog.hpp"

#include <algorithm>
#include <vector>

namespace dsl {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabComparisonDialog";
constexpr int kListId = 601;
constexpr int kCompareId = 602;
constexpr int kCancelId = 603;

struct DialogState {
    HWND owner = nullptr;
    HWND window = nullptr;
    HWND list = nullptr;
    const theme::Fonts* fonts = nullptr;
    const Catalog* catalog = nullptr;
    std::vector<int> candidates;
    std::optional<int> selection;
};

HWND control(HWND parent, const wchar_t* type, const wchar_t* text,
             DWORD style, int id, int x, int y, int width, int height, HFONT font) {
    HWND result = CreateWindowExW(type == std::wstring(L"LISTBOX") ? WS_EX_CLIENTEDGE : 0,
        type, text, WS_CHILD | WS_VISIBLE | style, x, y, width, height, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    theme::styleControl(result, font); return result;
}

std::wstring itemLabel(const CatalogItem& item) {
    std::wstring result = wide(item.system.empty() ? item.category : item.system);
    if (!result.empty()) result += L"  /  ";
    return result + wide(item.label);
}

void createChildren(DialogState& state) {
    control(state.window, L"STATIC",
        L"Choose another archive object supporting this same view:", SS_LEFT,
        0, 18, 15, 455, 22, state.fonts->body);
    state.list = control(state.window, L"LISTBOX", L"",
        LBS_NOTIFY | WS_VSCROLL, kListId, 18, 43, 455, 230, state.fonts->body);
    for (const int index : state.candidates) {
        const auto label = itemLabel(state.catalog->items[static_cast<std::size_t>(index)]);
        SendMessageW(state.list, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(label.c_str()));
    }
    SendMessageW(state.list, LB_SETCURSEL, 0, 0);
    control(state.window, L"BUTTON", L"Compare", BS_DEFPUSHBUTTON,
        kCompareId, 282, 286, 90, 29, state.fonts->body);
    control(state.window, L"BUTTON", L"Cancel", BS_PUSHBUTTON,
        kCancelId, 383, 286, 90, 29, state.fonts->body);
}

void accept(DialogState& state) {
    const LRESULT selected = SendMessageW(state.list, LB_GETCURSEL, 0, 0);
    if (selected != LB_ERR && static_cast<std::size_t>(selected) < state.candidates.size())
        state.selection = state.candidates[static_cast<std::size_t>(selected)];
    DestroyWindow(state.window);
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
        if (id == kCompareId || (id == kListId && HIWORD(wparam) == LBN_DBLCLK)) {
            accept(*state); return 0;
        }
        if (id == kCancelId || id == IDCANCEL) { DestroyWindow(window); return 0; }
    }
    if (message == WM_CLOSE) { DestroyWindow(window); return 0; }
    if (message == WM_CTLCOLORSTATIC) {
        static HBRUSH background = CreateSolidBrush(theme::panel);
        HDC dc = reinterpret_cast<HDC>(wparam); SetBkColor(dc, theme::panel);
        SetTextColor(dc, theme::text); return reinterpret_cast<LRESULT>(background);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerDialogClass() {
    WNDCLASSEXW type{sizeof(type)};
    type.lpfnWndProc = procedure; type.hInstance = GetModuleHandleW(nullptr);
    type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(theme::panel); type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

std::optional<int> showComparisonDialog(HWND owner, const Catalog& catalog,
                                        int currentIndex, const std::string& view,
                                        const theme::Fonts& fonts) {
    DialogState state; state.owner = owner; state.fonts = &fonts; state.catalog = &catalog;
    for (std::size_t index = 0; index < catalog.items.size(); ++index) {
        const auto& item = catalog.items[index];
        if (static_cast<int>(index) != currentIndex
            && std::find(item.views.begin(), item.views.end(), view) != item.views.end())
            state.candidates.push_back(static_cast<int>(index));
    }
    if (state.candidates.empty() || !registerDialogClass()) return std::nullopt;
    RECT ownerArea{}; GetWindowRect(owner, &ownerArea);
    const int width = 510, height = 365;
    state.window = CreateWindowExW(WS_EX_DLGMODALFRAME, kClassName,
        L"Compare archive objects", WS_CAPTION | WS_SYSMENU | WS_POPUP,
        ownerArea.left + (ownerArea.right - ownerArea.left - width) / 2,
        ownerArea.top + (ownerArea.bottom - ownerArea.top - height) / 2,
        width, height, owner, nullptr, GetModuleHandleW(nullptr), &state);
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

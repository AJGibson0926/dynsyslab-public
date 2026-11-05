#include "SystemParameterDialog.hpp"

#include "dynsyslab_viewer/Json.hpp"
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>

#include <algorithm>
#include <cmath>
#include <cwchar>
#include <cwctype>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabSystemParameterDialog";
constexpr int kParameterListId = 701, kValueId = 702, kSetId = 703;
constexpr int kResetId = 704, kResetAllId = 705, kDoneId = 706, kCancelId = 707;

struct DialogState {
    HWND window = nullptr, list = nullptr, value = nullptr;
    const dsl::theme::Fonts* fonts = nullptr;
    std::map<std::string, double> defaults, overrides, groups;
    std::map<std::string, std::string> units, groupUnits;
    std::vector<std::string> keys;
    std::optional<std::string> output;
    std::wstring context;
};

HWND control(HWND parent, const wchar_t* type, const wchar_t* text, DWORD style,
             int id, int x, int y, int width, int height, HFONT font) {
    const bool bordered = type == std::wstring(L"EDIT") || type == std::wstring(L"LISTBOX");
    HWND result = CreateWindowExW(bordered ? WS_EX_CLIENTEDGE : 0, type, text,
        WS_CHILD | WS_VISIBLE | style, x, y, width, height, parent,
        reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    dsl::theme::styleControl(result, font); return result;
}

std::wstring number(double value) {
    std::wostringstream text; text << std::setprecision(12) << value; return text.str();
}

std::size_t selection(const DialogState& state) {
    const auto value = SendMessageW(state.list, LB_GETCURSEL, 0, 0);
    return value == LB_ERR ? state.keys.size() : static_cast<std::size_t>(value);
}

void loadValue(DialogState& state) {
    const auto index = selection(state); if (index >= state.keys.size()) return;
    const auto& key = state.keys[index];
    const auto override = state.overrides.find(key);
    const double value = override == state.overrides.end() ? state.defaults.at(key) : override->second;
    SetWindowTextW(state.value, number(value).c_str());
}

void refresh(DialogState& state, std::size_t selected = 0) {
    SendMessageW(state.list, LB_RESETCONTENT, 0, 0);
    for (const auto& key : state.keys) {
        const auto override = state.overrides.find(key);
        const double value = override == state.overrides.end() ? state.defaults.at(key) : override->second;
        std::wstring label = override == state.overrides.end() ? L"  " : L"● ";
        label += dynsyslab::core::wide(key) + L" = " + number(value);
        const auto unit = state.units.find(key);
        if (unit != state.units.end() && !unit->second.empty())
            label += L"  " + dynsyslab::core::wide(unit->second);
        SendMessageW(state.list, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(label.c_str()));
    }
    if (!state.keys.empty()) SendMessageW(state.list, LB_SETCURSEL,
        static_cast<WPARAM>(std::min(selected, state.keys.size() - 1)), 0);
    loadValue(state);
}

bool readValue(HWND edit, double& result) {
    wchar_t buffer[96]{}; GetWindowTextW(edit, buffer, 96);
    wchar_t* end = nullptr; result = std::wcstod(buffer, &end);
    while (end && std::iswspace(*end)) ++end;
    return end && end != buffer && *end == L'\0' && std::isfinite(result);
}

void setOverride(DialogState& state) {
    const auto index = selection(state); double value = 0.0;
    if (index >= state.keys.size() || !readValue(state.value, value)) {
        MessageBoxW(state.window, L"Enter one finite numeric parameter value.",
                    L"Invalid parameter", MB_OK | MB_ICONWARNING); return;
    }
    state.overrides[state.keys[index]] = value; refresh(state, index);
}

std::string encoded(const std::map<std::string, double>& values) {
    dsl::json::Value::Object object;
    for (const auto& entry : values) object[entry.first] = entry.second;
    return dsl::json::dump(dsl::json::Value(std::move(object)));
}

void createChildren(DialogState& state) {
    control(state.window, L"STATIC", state.context.c_str(), SS_LEFT, 0,
        18, 13, 590, 22, state.fonts->heading);
    control(state.window, L"STATIC",
        L"Resolved parameters  (● = explicit override)", SS_LEFT, 0,
        18, 43, 390, 20, state.fonts->body);
    state.list = control(state.window, L"LISTBOX", L"",
        LBS_NOTIFY | WS_VSCROLL, kParameterListId, 18, 65, 590, 220, state.fonts->mono);
    control(state.window, L"STATIC", L"Selected value", SS_LEFT, 0,
        18, 296, 110, 20, state.fonts->body);
    state.value = control(state.window, L"EDIT", L"", ES_AUTOHSCROLL,
        kValueId, 130, 291, 190, 27, state.fonts->mono);
    control(state.window, L"BUTTON", L"Set override", BS_PUSHBUTTON,
        kSetId, 330, 291, 95, 28, state.fonts->body);
    control(state.window, L"BUTTON", L"Reset selected", BS_PUSHBUTTON,
        kResetId, 433, 291, 105, 28, state.fonts->body);
    control(state.window, L"BUTTON", L"Reset all", BS_PUSHBUTTON,
        kResetAllId, 546, 291, 62, 28, state.fonts->body);
    std::wstring groups = L"Dimensionless groups:";
    if (state.groups.empty()) groups += L" none declared";
    for (const auto& group : state.groups) {
        groups += L"\r\n" + dynsyslab::core::wide(group.first) + L" = " + number(group.second);
        const auto unit = state.groupUnits.find(group.first);
        if (unit != state.groupUnits.end() && !unit->second.empty())
            groups += L" " + dynsyslab::core::wide(unit->second);
    }
    control(state.window, L"EDIT", groups.c_str(),
        ES_LEFT | ES_MULTILINE | ES_READONLY | WS_VSCROLL, 0,
        18, 334, 590, 80, state.fonts->caption);
    control(state.window, L"BUTTON", L"Keep overrides", BS_DEFPUSHBUTTON,
        kDoneId, 380, 426, 112, 29, state.fonts->body);
    control(state.window, L"BUTTON", L"Cancel", BS_PUSHBUTTON,
        kCancelId, 501, 426, 107, 29, state.fonts->body);
    refresh(state);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<DialogState*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<DialogState*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        state->window = window; SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_CREATE) { createChildren(*state); return 0; }
    if (message == WM_COMMAND) {
        const int id = LOWORD(wparam);
        if (id == kParameterListId && HIWORD(wparam) == LBN_SELCHANGE) { loadValue(*state); return 0; }
        if (id == kSetId) { setOverride(*state); return 0; }
        if (id == kResetId) { const auto index=selection(*state); if(index<state->keys.size()){
            state->overrides.erase(state->keys[index]);refresh(*state,index);} return 0; }
        if (id == kResetAllId) { state->overrides.clear(); refresh(*state); return 0; }
        if (id == kDoneId || id == IDOK) { state->output=encoded(state->overrides);DestroyWindow(window);return 0; }
        if (id == kCancelId || id == IDCANCEL) { DestroyWindow(window); return 0; }
    }
    if (message == WM_CLOSE) { DestroyWindow(window); return 0; }
    if (message == WM_CTLCOLORSTATIC || message == WM_CTLCOLOREDIT || message == WM_CTLCOLORLISTBOX) {
        static HBRUSH background=CreateSolidBrush(dsl::theme::panel);HDC dc=reinterpret_cast<HDC>(wparam);
        SetBkColor(dc,dsl::theme::panel);SetTextColor(dc,dsl::theme::text);return reinterpret_cast<LRESULT>(background);
    }
    return DefWindowProcW(window,message,wparam,lparam);
}

bool registerDialogClass() {
    WNDCLASSEXW type{sizeof(type)};type.lpfnWndProc=procedure;type.hInstance=GetModuleHandleW(nullptr);
    type.hCursor=LoadCursorW(nullptr,IDC_ARROW);type.hbrBackground=CreateSolidBrush(dsl::theme::panel);
    type.lpszClassName=kClassName;return RegisterClassExW(&type)||GetLastError()==ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

std::optional<std::string> showSystemParameterDialog(
    HWND owner, const std::string& system, const std::string& preset,
    const std::string& currentOverrides, const dsl::theme::Fonts& fonts) {
    DialogState state;state.fonts=&fonts;
    try {
        dynsyslab::core::SystemConfiguration requested;requested.identifier=system;requested.preset=preset;
        const auto resolved=dynsyslab::workflow::resolveSystemDefinition(requested).configuration;
        state.defaults=resolved.parameters;state.units=resolved.parameterUnits;
        state.groups=resolved.dimensionlessGroups;state.groupUnits=resolved.dimensionlessGroupUnits;
        const auto existing=dsl::json::parse(currentOverrides.empty()?"{}":currentOverrides);
        for(const auto& item:existing.object()) state.overrides[item.first]=item.second.number();
        for(const auto& item:state.overrides) if(!state.defaults.count(item.first)) state.defaults[item.first]=item.second;
        for(const auto& item:state.defaults) state.keys.push_back(item.first);
        state.context=dynsyslab::core::wide(resolved.identifier)
            +(preset.empty()?L" / model default":L" / "+dynsyslab::core::wide(preset));
    } catch(const std::exception& error) {
        MessageBoxW(owner,dynsyslab::core::wide(error.what()).c_str(),
                    L"Unable to resolve parameters",MB_OK|MB_ICONERROR);
        return std::nullopt;
    }
    if(!registerDialogClass()) return std::nullopt;
    RECT area{};GetWindowRect(owner,&area);const int width=645,height=515;
    state.window=CreateWindowExW(WS_EX_DLGMODALFRAME,kClassName,L"System parameters",
        WS_CAPTION|WS_SYSMENU|WS_POPUP,area.left+(area.right-area.left-width)/2,
        area.top+(area.bottom-area.top-height)/2,width,height,owner,nullptr,GetModuleHandleW(nullptr),&state);
    if(!state.window) return std::nullopt;
    dsl::theme::configureTitleBar(state.window);EnableWindow(owner,FALSE);
    ShowWindow(state.window,SW_SHOW);UpdateWindow(state.window);MSG event{};
    while(IsWindow(state.window)&&GetMessageW(&event,nullptr,0,0)>0){
        if(!IsDialogMessageW(state.window,&event)){TranslateMessage(&event);DispatchMessageW(&event);}}
    EnableWindow(owner,TRUE);SetForegroundWindow(owner);return state.output;
}

}  // namespace dsl::compute

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_compute/ScientificParameterFormatting.hpp"

#include <algorithm>
#include <commctrl.h>
#include <cwchar>
#include <cwctype>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <uxtheme.h>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabInlineValue";
constexpr wchar_t kMathLabelClassName[] = L"DynSysLabMathLabel";

struct State {
    std::wstring value;
    HFONT font = nullptr;
    HWND editor = nullptr;
    int id = 0;
    bool scientific = false;
};

struct MathLabelState {
    std::wstring latex;
    HFONT math = nullptr;
    HFONT caption = nullptr;
};

void closeEditor(State& state) {
    HWND editor = state.editor; state.editor = nullptr;
    if (editor) DestroyWindow(editor);
}

void notifyChange(HWND window, const State& state) {
    HWND parent = GetParent(window);
    if (parent && state.id)
        SendMessageW(parent, WM_COMMAND, MAKEWPARAM(state.id, EN_CHANGE),
                     reinterpret_cast<LPARAM>(window));
}

void commit(HWND window, State& state) {
    if (!state.editor) return;
    const int length = GetWindowTextLengthW(state.editor);
    std::wstring value(static_cast<std::size_t>(length + 1), L'\0');
    GetWindowTextW(state.editor, value.data(), length + 1);
    value.resize(static_cast<std::size_t>(length));
    const auto first = value.find_first_not_of(L" \t\r\n");
    const auto last = value.find_last_not_of(L" \t\r\n");
    if (first == std::wstring::npos) value.clear();
    else value = value.substr(first, last - first + 1);
    if (state.scientific) {
        double numeric = 0.0;
        if (!parseScientificParameterInput(value, numeric)) {
            MessageBeep(MB_ICONWARNING); SetFocus(state.editor);
            SendMessageW(state.editor, EM_SETSEL, 0, -1); return;
        }
        std::wostringstream canonical;
        canonical << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << std::defaultfloat << numeric;
        value = canonical.str();
    }
    const bool changed = value != state.value;
    state.value = std::move(value); closeEditor(state);
    InvalidateRect(window, nullptr, FALSE);
    if (changed) notifyChange(window, state);
}

LRESULT CALLBACK editorProcedure(HWND editor, UINT message, WPARAM wparam,
                                 LPARAM lparam, UINT_PTR, DWORD_PTR reference) {
    HWND valueControl = reinterpret_cast<HWND>(reference);
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(valueControl, GWLP_USERDATA));
    if (state && message == WM_KEYDOWN && wparam == VK_RETURN) {
        commit(valueControl, *state); return 0;
    }
    if (state && message == WM_KEYDOWN && wparam == VK_ESCAPE) {
        closeEditor(*state); SetFocus(valueControl);
        InvalidateRect(valueControl, nullptr, FALSE); return 0;
    }
    if (state && message == WM_KILLFOCUS
        && reinterpret_cast<HWND>(wparam) != valueControl)
        PostMessageW(valueControl, WM_APP + 72, 0, 0);
    return DefSubclassProc(editor, message, wparam, lparam);
}

void beginEdit(HWND window, State& state) {
    if (state.editor) return;
    RECT area{}; GetClientRect(window, &area);
    state.editor = CreateWindowExW(WS_EX_CLIENTEDGE, L"EDIT", state.value.c_str(),
        WS_CHILD | WS_VISIBLE | ES_AUTOHSCROLL
            | (state.scientific ? ES_RIGHT : 0), 0, 0, area.right, area.bottom,
        window, nullptr, GetModuleHandleW(nullptr), nullptr);
    SendMessageW(state.editor, WM_SETFONT, reinterpret_cast<WPARAM>(state.font), TRUE);
    SetWindowTheme(state.editor, L"Explorer", nullptr);
    SetWindowSubclass(state.editor, editorProcedure, 1,
                      reinterpret_cast<DWORD_PTR>(window));
    SetFocus(state.editor); SendMessageW(state.editor, EM_SETSEL, 0, -1);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<State*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_SETTEXT) {
        state->value = lparam ? reinterpret_cast<const wchar_t*>(lparam) : L"";
        if (state->editor) SetWindowTextW(state->editor, state->value.c_str());
        InvalidateRect(window, nullptr, FALSE); return TRUE;
    }
    if (message == WM_GETTEXTLENGTH) return static_cast<LRESULT>(state->value.size());
    if (message == WM_GETTEXT) {
        const std::size_t capacity = static_cast<std::size_t>(wparam);
        if (!capacity || !lparam) return 0;
        wchar_t* output = reinterpret_cast<wchar_t*>(lparam);
        const std::size_t copied = std::min(capacity - 1, state->value.size());
        std::wmemcpy(output, state->value.data(), copied); output[copied] = L'\0';
        return static_cast<LRESULT>(copied);
    }
    if (message == WM_LBUTTONUP || (message == WM_KEYDOWN && (wparam == VK_RETURN || wparam == VK_SPACE))) {
        beginEdit(window, *state); return 0;
    }
    if (message == WM_SETFOCUS) { InvalidateRect(window, nullptr, FALSE); return 0; }
    if (message == WM_KILLFOCUS) { InvalidateRect(window, nullptr, FALSE); return 0; }
    if (message == WM_APP + 72) {
        if (state->editor && GetFocus() != state->editor) commit(window, *state);
        return 0;
    }
    if (message == WM_SIZE && state->editor) {
        MoveWindow(state->editor, 0, 0, LOWORD(lparam), HIWORD(lparam), TRUE); return 0;
    }
    if (message == WM_SETCURSOR) { SetCursor(LoadCursorW(nullptr, IDC_IBEAM)); return TRUE; }
    if (message == WM_PAINT) {
        PAINTSTRUCT paint{}; HDC dc = BeginPaint(window, &paint);
        RECT area{}; GetClientRect(window, &area);
        HBRUSH background = CreateSolidBrush(state->scientific
            ? dsl::theme::background : dsl::theme::panelRaised);
        FillRect(dc, &area, background); DeleteObject(background);
        if (!state->scientific) {
            HBRUSH border = CreateSolidBrush(GetFocus() == window
                ? dsl::theme::accent : dsl::theme::border);
            FrameRect(dc, &area, border); DeleteObject(border);
        }
        SetBkMode(dc, TRANSPARENT); SetTextColor(dc, dsl::theme::text);
        SelectObject(dc, state->font);
        RECT text{area.left + (state->scientific ? 2 : 8), area.top,
                  area.right - (state->scientific ? 2 : 7), area.bottom};
        std::wstring rendered = state->value;
        if (state->scientific) {
            double numeric = 0.0;
            if (parseScientificParameterInput(state->value, numeric))
                rendered = scientificParameterNumber(numeric, false);
        }
        const wchar_t* display = rendered.empty() ? L"—" : rendered.c_str();
        DrawTextW(dc, display, -1, &text,
                  (state->scientific ? DT_RIGHT : DT_LEFT) | DT_SINGLELINE
                  | DT_VCENTER | DT_END_ELLIPSIS | DT_NOPREFIX);
        EndPaint(window, &paint); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_NCDESTROY) {
        closeEditor(*state); delete state; SetWindowLongPtrW(window, GWLP_USERDATA, 0);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerInlineValueClass() {
    WNDCLASSEXW type{sizeof(type)}; type.lpfnWndProc = procedure;
    type.hInstance = GetModuleHandleW(nullptr); type.hCursor = LoadCursorW(nullptr, IDC_IBEAM);
    type.hbrBackground = CreateSolidBrush(dsl::theme::panelRaised); type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

std::wstring greek(std::wstring value) {
    const std::pair<const wchar_t*, const wchar_t*> replacements[]{
        {L"\\Delta", L"\u0394"}, {L"\\delta", L"\u03b4"},
        {L"\\tau", L"\u03c4"}, {L"\\omega", L"\u03c9"}};
    for (const auto& replacement : replacements) {
        std::size_t position = 0;
        while ((position = value.find(replacement.first, position))
               != std::wstring::npos) {
            value.replace(position, std::wcslen(replacement.first), replacement.second);
            position += std::wcslen(replacement.second);
        }
    }
    return value;
}

void splitLatex(const std::wstring& source, std::wstring& base,
                std::wstring& subscript) {
    const auto marker = source.find(L"_{");
    if (marker == std::wstring::npos || source.empty() || source.back() != L'}') {
        base = greek(source); return;
    }
    base = greek(source.substr(0, marker));
    subscript = greek(source.substr(marker + 2, source.size() - marker - 3));
}

LRESULT CALLBACK mathLabelProcedure(HWND window, UINT message, WPARAM wparam,
                                    LPARAM lparam) {
    auto* state = reinterpret_cast<MathLabelState*>(
        GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<MathLabelState*>(
            reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_SETTEXT) {
        state->latex = lparam ? reinterpret_cast<const wchar_t*>(lparam) : L"";
        InvalidateRect(window, nullptr, FALSE); return TRUE;
    }
    if (message == WM_GETTEXTLENGTH)
        return static_cast<LRESULT>(state->latex.size());
    if (message == WM_GETTEXT) {
        const std::size_t capacity = static_cast<std::size_t>(wparam);
        if (!capacity || !lparam) return 0;
        wchar_t* output = reinterpret_cast<wchar_t*>(lparam);
        const std::size_t copied = std::min(capacity - 1, state->latex.size());
        std::wmemcpy(output, state->latex.data(), copied); output[copied] = L'\0';
        return static_cast<LRESULT>(copied);
    }
    if (message == WM_PAINT) {
        PAINTSTRUCT paint{}; HDC dc = BeginPaint(window, &paint);
        RECT area{}; GetClientRect(window, &area);
        HBRUSH background = CreateSolidBrush(dsl::theme::background);
        FillRect(dc, &area, background); DeleteObject(background);
        SetBkMode(dc, TRANSPARENT); SetTextColor(dc, dsl::theme::accent);
        std::wstring base, subscript; splitLatex(state->latex, base, subscript);
        SelectObject(dc, state->math);
        SIZE baseSize{}; GetTextExtentPoint32W(dc, base.c_str(),
            static_cast<int>(base.size()), &baseSize);
        SIZE subscriptSize{};
        if (!subscript.empty()) {
            SelectObject(dc, state->caption);
            GetTextExtentPoint32W(dc, subscript.c_str(),
                static_cast<int>(subscript.size()), &subscriptSize);
        }
        const int left = std::max(0,
            (static_cast<int>(area.right) - static_cast<int>(baseSize.cx)
                - static_cast<int>(subscriptSize.cx)) / 2);
        SelectObject(dc, state->math);
        RECT main{left, area.top, area.right, area.bottom};
        DrawTextW(dc, base.c_str(), -1, &main,
            DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
        if (!subscript.empty()) {
            SelectObject(dc, state->caption);
            RECT sub{left + baseSize.cx, area.top + 6, area.right, area.bottom};
            DrawTextW(dc, subscript.c_str(), -1, &sub,
                DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
        }
        EndPaint(window, &paint); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_NCDESTROY) {
        if (state->math) DeleteObject(state->math);
        delete state; SetWindowLongPtrW(window, GWLP_USERDATA, 0);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerMathLabelClass() {
    WNDCLASSEXW type{sizeof(type)}; type.lpfnWndProc = mathLabelProcedure;
    type.hInstance = GetModuleHandleW(nullptr); type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(dsl::theme::background);
    type.lpszClassName = kMathLabelClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

HWND makeInlineValueImpl(HWND parent, int id, const wchar_t* text, HFONT font,
                         bool scientific) {
    if (!registerInlineValueClass()) return nullptr;
    auto* state = new State; state->value = text ? text : L"";
    state->font = font; state->id = id; state->scientific = scientific;
    HWND result = CreateWindowExW(0, kClassName, L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP,
        0, 0, 0, 0, parent, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
        GetModuleHandleW(nullptr), state);
    if (!result) delete state;
    return result;
}

}  // namespace

HWND makeInlineValue(HWND parent, int id, const wchar_t* text, HFONT font) {
    return makeInlineValueImpl(parent, id, text, font, false);
}

HWND makeScientificInlineValue(HWND parent, int id, const wchar_t* text,
                               HFONT font) {
    return makeInlineValueImpl(parent, id, text, font, true);
}

HWND makeMathLabel(HWND parent, const wchar_t* latex, HFONT body,
                   HFONT caption) {
    if (!registerMathLabelClass()) return nullptr;
    auto* state = new MathLabelState; state->latex = latex ? latex : L"";
    state->caption = caption;
    LOGFONTW specification{}; GetObjectW(body, sizeof(specification), &specification);
    specification.lfItalic = TRUE;
    lstrcpynW(specification.lfFaceName, L"Cambria Math", LF_FACESIZE);
    state->math = CreateFontIndirectW(&specification);
    HWND result = CreateWindowExW(0, kMathLabelClassName, L"",
        WS_CHILD | WS_VISIBLE, 0, 0, 0, 0, parent, nullptr,
        GetModuleHandleW(nullptr), state);
    if (!result) { if (state->math) DeleteObject(state->math); delete state; }
    return result;
}

}  // namespace dsl::compute

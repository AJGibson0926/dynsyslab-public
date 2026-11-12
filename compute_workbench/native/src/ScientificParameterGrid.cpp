#include "dynsyslab_compute/ScientificParameterGrid.hpp"

#include "dynsyslab_compute/ScientificParameterFormatting.hpp"
#include "dynsyslab_viewer/GdiBackBuffer.hpp"
#include "dynsyslab_viewer/OrnamentalScrollBar.hpp"

#include <algorithm>
#include <array>
#include <commctrl.h>
#include <utility>
#include <vector>
#include <uxtheme.h>
#include <windowsx.h>

namespace dsl::compute {
namespace {

constexpr wchar_t kClassName[] = L"DynSysLabScientificParameterGrid";
constexpr std::size_t kNone = static_cast<std::size_t>(-1);

struct Item {
    std::string key;
    std::wstring latex;
    std::wstring value;
    std::wstring unit;
    double numericValue{};
};

struct Section {
    std::wstring title;
    std::vector<Item> items;
    std::wstring message;
    bool editable{};
    bool showUnits{};
};

struct State {
    std::array<Section, 3> sections;
    HFONT body{};
    HFONT caption{};
    HFONT math{};
    ScientificParameterEditHandler onEdit;
    HWND editor{};
    std::size_t editingSection{kNone};
    std::size_t editingItem{kNone};
    bool longPrecision{};
    dsl::OrnamentalScrollBar scroll;
    dsl::GdiBackBuffer backBuffer;
};

struct Cell {
    std::size_t item{};
    RECT area{};
    int valueLeft{};
    int valueRight{};
    int unitLeft{};
};

struct SectionLayout {
    RECT title{};
    int bottom{};
    std::vector<Cell> cells;
};

struct Layout {
    std::array<SectionLayout, 3> sections;
    int contentHeight{};
};

int scaled(HWND window, int value) {
    const UINT dpi = GetDpiForWindow(window);
    return MulDiv(value, dpi ? static_cast<int>(dpi) : 96, 96);
}

std::size_t sectionIndex(ScientificParameterSection section) {
    return static_cast<std::size_t>(section);
}

std::wstring greek(std::wstring value) {
    const std::pair<const wchar_t*, const wchar_t*> replacements[]{
        {L"\\alpha", L"\u03b1"}, {L"\\beta", L"\u03b2"}, {L"\\gamma", L"\u03b3"},
        {L"\\delta", L"\u03b4"}, {L"\\epsilon", L"\u03b5"}, {L"\\kappa", L"\u03ba"},
        {L"\\lambda", L"\u03bb"}, {L"\\mu", L"\u03bc"}, {L"\\nu", L"\u03bd"},
        {L"\\rho", L"\u03c1"}, {L"\\sigma", L"\u03c3"}, {L"\\tau", L"\u03c4"},
        {L"\\omega", L"\u03c9"}, {L"\\chi", L"\u03c7"}, {L"\\Pi", L"\u03a0"}};
    for (const auto& replacement : replacements) {
        std::size_t position = 0;
        while ((position = value.find(replacement.first, position)) != std::wstring::npos) {
            value.replace(position, std::wcslen(replacement.first), replacement.second);
            position += std::wcslen(replacement.second);
        }
    }
    return value;
}

void splitLatex(const std::wstring& source, std::wstring& base,
                std::wstring& subscript) {
    const auto marker = source.find(L"_{");
    if (marker == std::wstring::npos || source.back() != L'}') {
        base = greek(source); return;
    }
    base = greek(source.substr(0, marker));
    subscript = greek(source.substr(marker + 2, source.size() - marker - 3));
}

int textWidth(HDC dc, HFONT font, const std::wstring& text) {
    HGDIOBJ previous = SelectObject(dc, font);
    SIZE size{};
    GetTextExtentPoint32W(dc, text.c_str(), static_cast<int>(text.size()), &size);
    SelectObject(dc, previous);
    return size.cx;
}

int mathWidth(HDC dc, const Item& item, const State& state) {
    std::wstring base, subscript;
    splitLatex(item.latex, base, subscript);
    return textWidth(dc, state.math, base)
        + (subscript.empty() ? 0 : textWidth(dc, state.caption, subscript));
}

Layout geometry(HWND window, HDC dc, const State& state) {
    RECT client{}; GetClientRect(window, &client);
    const int padding = scaled(window, 8), displayGap = scaled(window, 18);
    const int itemPadding = scaled(window, 6), beforeEquals = scaled(window, 3);
    const int afterEquals = scaled(window, 5), unitGap = scaled(window, 7);
    const int titleH = scaled(window, 24), rowH = scaled(window, 24);
    const int sectionGap = scaled(window, 8);
    const int clientWidth = static_cast<int>(client.right)
        - dsl::OrnamentalScrollBar::preferredWidth(window);
    const int usableWidth = std::max(1, clientWidth - 2 * padding);
    const int equalsWidth = textWidth(dc, state.body, L"=");
    Layout result; int top = scaled(window, 3);
    for (std::size_t sectionIndexValue = 0; sectionIndexValue < state.sections.size();
         ++sectionIndexValue) {
        const auto& section = state.sections[sectionIndexValue];
        auto& layout = result.sections[sectionIndexValue];
        layout.title = {padding, top, clientWidth - padding, top + titleH};
        top += titleH;
        if (!section.message.empty()) {
            top += scaled(window, 34);
        } else if (!section.items.empty()) {
            int left = padding, rowTop = top;
            for (std::size_t item = 0; item < section.items.size(); ++item) {
                const auto& value = section.items[item];
                const int symbolWidth = mathWidth(dc, value, state);
                const int valueWidth = textWidth(dc, state.body, value.value);
                const int unitWidth = section.showUnits
                    ? textWidth(dc, state.caption, value.unit) : 0;
                const int unitSlot = unitWidth > 0 ? unitGap + unitWidth : 0;
                const int intrinsicWidth = 2 * itemPadding + symbolWidth
                    + beforeEquals + equalsWidth + afterEquals + valueWidth + unitSlot;
                const int displayWidth = std::min(usableWidth, intrinsicWidth);
                if (left > padding && left + displayWidth > clientWidth - padding) {
                    left = padding; rowTop += rowH;
                }
                const int valueLeft = left + itemPadding + symbolWidth
                    + beforeEquals + equalsWidth + afterEquals;
                const int valueRight = std::min(valueLeft + valueWidth,
                    left + displayWidth - itemPadding - unitSlot);
                layout.cells.push_back({item,
                    {left, rowTop, left + displayWidth, rowTop + rowH},
                    valueLeft, valueRight, valueRight + unitGap});
                left += displayWidth + displayGap;
            }
            top = rowTop + rowH;
        }
        layout.bottom = top;
        top += sectionGap;
    }
    result.contentHeight = top;
    return result;
}

void cancelEdit(State& state) {
    HWND editor = state.editor; state.editor = nullptr;
    state.editingSection = kNone; state.editingItem = kNone;
    if (editor) DestroyWindow(editor);
}

void updateScroll(HWND window, State& state) {
    RECT area{}; GetClientRect(window, &area);
    HDC dc = GetDC(window); const auto layout = geometry(window, dc, state);
    ReleaseDC(window, dc);
    state.scroll.configure(layout.contentHeight, static_cast<int>(area.bottom));
}

bool parseNumber(HWND editor, double& value) {
    wchar_t buffer[128]{}; GetWindowTextW(editor, buffer, 128);
    return parseScientificParameterInput(buffer, value);
}

void commitEdit(HWND window, State& state) {
    if (!state.editor || state.editingSection >= state.sections.size()) return;
    auto& items = state.sections[state.editingSection].items;
    if (state.editingItem >= items.size()) return;
    double value = 0.0;
    if (!parseNumber(state.editor, value)) {
        MessageBeep(MB_ICONWARNING); SetFocus(state.editor);
        SendMessageW(state.editor, EM_SETSEL, 0, -1); return;
    }
    const std::string key = items[state.editingItem].key;
    const double previous = items[state.editingItem].numericValue;
    const auto callback = state.onEdit;
    cancelEdit(state); InvalidateRect(window, nullptr, FALSE);
    if (callback && value != previous) callback(key, value);
}

LRESULT CALLBACK editorProcedure(HWND editor, UINT message, WPARAM wparam,
                                 LPARAM lparam, UINT_PTR, DWORD_PTR reference) {
    HWND grid = reinterpret_cast<HWND>(reference);
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(grid, GWLP_USERDATA));
    if (state && message == WM_KEYDOWN && wparam == VK_RETURN) {
        commitEdit(grid, *state); return 0;
    }
    if (state && message == WM_KEYDOWN && wparam == VK_ESCAPE) {
        cancelEdit(*state); SetFocus(grid); InvalidateRect(grid, nullptr, FALSE); return 0;
    }
    if (state && message == WM_KILLFOCUS && reinterpret_cast<HWND>(wparam) != grid)
        PostMessageW(grid, WM_APP + 71, 0, 0);
    return DefSubclassProc(editor, message, wparam, lparam);
}

void beginEdit(HWND window, State& state, int x, int y) {
    if (!state.onEdit) return;
    if (state.editor) commitEdit(window, state);
    HDC dc = GetDC(window); const auto layout = geometry(window, dc, state);
    const int contentY = y + state.scroll.position();
    for (std::size_t sectionIndexValue = 0; sectionIndexValue < state.sections.size();
         ++sectionIndexValue) {
        const auto& section = state.sections[sectionIndexValue];
        if (!section.editable || !section.message.empty()) continue;
        for (const auto& cell : layout.sections[sectionIndexValue].cells) {
            if (x < cell.valueLeft || x >= cell.valueRight
                || contentY < cell.area.top || contentY >= cell.area.bottom) continue;
            const int renderedWidth = textWidth(dc, state.body,
                section.items[cell.item].value);
            if (x < cell.valueRight - renderedWidth) continue;
            ReleaseDC(window, dc);
            state.editingSection = sectionIndexValue; state.editingItem = cell.item;
            state.editor = CreateWindowExW(WS_EX_CLIENTEDGE, L"EDIT",
                section.items[cell.item].value.c_str(),
                WS_CHILD | WS_VISIBLE | ES_AUTOHSCROLL | ES_RIGHT,
                cell.valueLeft, cell.area.top - state.scroll.position() + 1,
                cell.valueRight - cell.valueLeft, cell.area.bottom - cell.area.top - 2,
                window, nullptr, GetModuleHandleW(nullptr), nullptr);
            SendMessageW(state.editor, WM_SETFONT,
                reinterpret_cast<WPARAM>(state.body), TRUE);
        SetWindowTheme(state.editor, L"Explorer", nullptr);
            SetWindowSubclass(state.editor, editorProcedure, 1,
                reinterpret_cast<DWORD_PTR>(window));
            SetFocus(state.editor); SendMessageW(state.editor, EM_SETSEL, 0, -1);
            return;
        }
    }
    ReleaseDC(window, dc);
}

int drawMath(HWND window, HDC dc, const Item& item, const State& state, RECT cell) {
    std::wstring base, subscript; splitLatex(item.latex, base, subscript);
    SelectObject(dc, state.math); SetTextColor(dc, dsl::theme::accent);
    RECT main{cell.left, cell.top + scaled(window, 1), cell.right, cell.bottom};
    DrawTextW(dc, base.c_str(), -1, &main,
        DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
    SIZE size{}; GetTextExtentPoint32W(dc, base.c_str(), static_cast<int>(base.size()), &size);
    if (subscript.empty()) return size.cx;
    SelectObject(dc, state.caption);
    RECT sub{cell.left + size.cx, cell.top + scaled(window, 8), cell.right, cell.bottom};
    DrawTextW(dc, subscript.c_str(), -1, &sub,
        DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
    SIZE subscriptSize{};
    GetTextExtentPoint32W(dc, subscript.c_str(), static_cast<int>(subscript.size()),
                         &subscriptSize);
    return size.cx + subscriptSize.cx;
}

void paint(HWND window, State& state, HDC dc) {
    RECT area{}; GetClientRect(window, &area);
    HBRUSH background = CreateSolidBrush(dsl::theme::panelRaised);
    FillRect(dc, &area, background); DeleteObject(background);
    SetBkMode(dc, TRANSPARENT); const auto layout = geometry(window, dc, state);
    const int offset = state.scroll.position();
    for (std::size_t sectionIndexValue = 0; sectionIndexValue < state.sections.size();
         ++sectionIndexValue) {
        const auto& section = state.sections[sectionIndexValue];
        const auto& sectionLayout = layout.sections[sectionIndexValue];
        RECT title = sectionLayout.title; OffsetRect(&title, 0, -offset);
        SelectObject(dc, state.caption); SetTextColor(dc, dsl::theme::muted);
        DrawTextW(dc, section.title.c_str(), -1, &title,
            DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
        if (sectionIndexValue > 0) {
            HPEN divider = CreatePen(PS_SOLID, 1, dsl::theme::border);
            HGDIOBJ old = SelectObject(dc, divider);
            MoveToEx(dc, title.left, title.top - scaled(window, 3), nullptr);
            LineTo(dc, title.right, title.top - scaled(window, 3));
            SelectObject(dc, old); DeleteObject(divider);
        }
        if (!section.message.empty()) {
            SelectObject(dc, state.body); SetTextColor(dc, dsl::theme::warning);
            RECT message{title.left, title.bottom,
                title.right, sectionLayout.bottom - offset};
            DrawTextW(dc, section.message.c_str(), -1, &message,
                DT_LEFT | DT_WORDBREAK | DT_NOPREFIX);
            continue;
        }
        for (const auto& cell : sectionLayout.cells) {
            RECT line = cell.area; OffsetRect(&line, 0, -offset);
            if (line.bottom <= 0 || line.top >= area.bottom) continue;
            RECT symbol{line.left + scaled(window, 6), line.top,
                cell.valueLeft - scaled(window, 5), line.bottom};
            const int symbolWidth = drawMath(window, dc, section.items[cell.item],
                state, symbol);
            SelectObject(dc, state.body); SetTextColor(dc, dsl::theme::muted);
            RECT equals{symbol.left + symbolWidth + scaled(window, 3), line.top,
                cell.valueLeft, line.bottom};
            DrawTextW(dc, L"=", 1, &equals,
                DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_NOPREFIX);
            if (state.editor && state.editingSection == sectionIndexValue
                && state.editingItem == cell.item) continue;
            SelectObject(dc, state.body); SetTextColor(dc, dsl::theme::text);
            RECT value{cell.valueLeft, line.top, cell.valueRight, line.bottom};
            DrawTextW(dc, section.items[cell.item].value.c_str(), -1, &value,
                DT_RIGHT | DT_SINGLELINE | DT_VCENTER | DT_END_ELLIPSIS | DT_NOPREFIX);
            if (section.showUnits && !section.items[cell.item].unit.empty()) {
                SelectObject(dc, state.caption); SetTextColor(dc, dsl::theme::muted);
                RECT unit{cell.unitLeft, line.top,
                    line.right - scaled(window, 6), line.bottom};
                DrawTextW(dc, section.items[cell.item].unit.c_str(), -1, &unit,
                    DT_LEFT | DT_SINGLELINE | DT_VCENTER | DT_END_ELLIPSIS | DT_NOPREFIX);
            }
        }
    }
    state.scroll.paint(window, dc, area);
    HBRUSH border = CreateSolidBrush(dsl::theme::border);
    FrameRect(dc, &area, border); DeleteObject(border);
}

LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(window, GWLP_USERDATA));
    if (message == WM_NCCREATE) {
        state = static_cast<State*>(reinterpret_cast<CREATESTRUCTW*>(lparam)->lpCreateParams);
        SetWindowLongPtrW(window, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(state));
    }
    if (!state) return DefWindowProcW(window, message, wparam, lparam);
    if (message == WM_SIZE) { cancelEdit(*state); updateScroll(window, *state);
        InvalidateRect(window, nullptr, FALSE); return 0; }
    if (message == WM_LBUTTONDOWN) {
        RECT area{}; GetClientRect(window, &area);
        if (state->scroll.buttonDown(window, area,
                GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam))) {
            cancelEdit(*state); InvalidateRect(window, nullptr, FALSE); return 0;
        }
    }
    if (message == WM_LBUTTONUP) {
        RECT area{}; GetClientRect(window, &area);
        if (state->scroll.buttonUp(window, area,
                GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam))) {
            InvalidateRect(window, nullptr, FALSE); return 0;
        }
        beginEdit(window, *state, GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam));
        return 0;
    }
    if (message == WM_MOUSEMOVE) {
        RECT area{}; GetClientRect(window, &area);
        if (state->scroll.mouseMove(window, area,
                GET_X_LPARAM(lparam), GET_Y_LPARAM(lparam))) {
            cancelEdit(*state); InvalidateRect(window, nullptr, FALSE); return 0;
        }
    }
    if (message == WM_APP + 71) { if (state->editor && GetFocus() != state->editor)
        commitEdit(window, *state); return 0; }
    if (message == WM_MOUSEWHEEL) {
        const int pixels = -GET_WHEEL_DELTA_WPARAM(wparam) / WHEEL_DELTA
            * scaled(window, 58);
        if (state->scroll.scrollBy(pixels)) {
            cancelEdit(*state); InvalidateRect(window, nullptr, FALSE);
        }
        return 0;
    }
    if (message == WM_CAPTURECHANGED) state->scroll.cancelDrag(window);
    if (message == WM_PAINT) {
        PAINTSTRUCT paintState{};
        HDC target = BeginPaint(window, &paintState);
        RECT area{}; GetClientRect(window, &area);
        HDC frame = state->backBuffer.target(target,
            static_cast<int>(area.right), static_cast<int>(area.bottom));
        paint(window, *state, frame);
        if (frame != target) state->backBuffer.present(target, paintState.rcPaint);
        EndPaint(window, &paintState); return 0;
    }
    if (message == WM_PRINTCLIENT) {
        paint(window, *state, reinterpret_cast<HDC>(wparam)); return 0;
    }
    if (message == WM_ERASEBKGND) return 1;
    if (message == WM_NCDESTROY) {
        if (state->math) DeleteObject(state->math); delete state;
        SetWindowLongPtrW(window, GWLP_USERDATA, 0);
    }
    return DefWindowProcW(window, message, wparam, lparam);
}

bool registerClass() {
    WNDCLASSEXW type{sizeof(type)}; type.lpfnWndProc = procedure;
    type.hInstance = GetModuleHandleW(nullptr); type.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    type.hbrBackground = CreateSolidBrush(dsl::theme::panelRaised);
    type.lpszClassName = kClassName;
    return RegisterClassExW(&type) || GetLastError() == ERROR_CLASS_ALREADY_EXISTS;
}

}  // namespace

HWND makeScientificParameterGrid(HWND parent, const dsl::theme::Fonts& fonts,
                                 ScientificParameterEditHandler onEdit) {
    if (!registerClass()) return nullptr;
    auto* state = new State;
    state->sections[0] = {L"PHYSICAL / MODEL PARAMETERS  \u2022  click a value to edit",
        {}, {}, true, false};
    state->sections[1] = {L"DERIVED PHYSICAL SCALES", {}, {}, false, false};
    state->sections[2] = {L"DIMENSIONLESS GROUPS", {}, {}, false, false};
    state->body = fonts.body; state->caption = fonts.caption;
    state->onEdit = std::move(onEdit);
    LOGFONTW specification{}; GetObjectW(fonts.body, sizeof(specification), &specification);
    specification.lfItalic = TRUE;
    lstrcpynW(specification.lfFaceName, L"Cambria Math", LF_FACESIZE);
    state->math = CreateFontIndirectW(&specification);
    HWND result = CreateWindowExW(0, kClassName, L"",
        WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0, 0, 0,
        parent, nullptr, GetModuleHandleW(nullptr), state);
    if (!result) { if (state->math) DeleteObject(state->math); delete state; }
    return result;
}

void setScientificParameterGridSection(
    HWND grid, ScientificParameterSection section,
    const std::map<std::string, double>& values,
    const std::map<std::string, std::string>& units) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(grid, GWLP_USERDATA));
    if (!state) return;
    cancelEdit(*state); auto& target = state->sections[sectionIndex(section)];
    target.items.clear(); target.message.clear(); target.showUnits = false;
    for (const auto& entry : values) {
        const auto unit = units.find(entry.first);
        const std::string sourceUnit = unit == units.end() ? std::string{} : unit->second;
        if (!sourceUnit.empty() && sourceUnit != "dimensionless") target.showUnits = true;
        target.items.push_back({entry.first, scientificParameterLatex(entry.first),
            scientificParameterNumber(entry.second, state->longPrecision),
            scientificParameterUnit(sourceUnit), entry.second});
    }
    updateScroll(grid, *state); InvalidateRect(grid, nullptr, FALSE);
}

void setScientificParameterGridSectionMessage(
    HWND grid, ScientificParameterSection section, const std::wstring& message) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(grid, GWLP_USERDATA));
    if (!state) return;
    cancelEdit(*state); auto& target = state->sections[sectionIndex(section)];
    target.items.clear(); target.message = message;
    updateScroll(grid, *state); InvalidateRect(grid, nullptr, FALSE);
}

void setScientificParameterGridLongPrecision(HWND grid, bool enabled) {
    auto* state = reinterpret_cast<State*>(GetWindowLongPtrW(grid, GWLP_USERDATA));
    if (!state || state->longPrecision == enabled) return;
    cancelEdit(*state); state->longPrecision = enabled;
    for (auto& section : state->sections)
        for (auto& item : section.items)
            item.value = scientificParameterNumber(item.numericValue, enabled);
    updateScroll(grid, *state); InvalidateRect(grid, nullptr, FALSE);
}

}  // namespace dsl::compute

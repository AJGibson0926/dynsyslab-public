#include "dynsyslab_viewer/FolderTab.hpp"
#include "dynsyslab_viewer/Theme.hpp"
#include <string>
namespace dsl {
namespace {
std::wstring windowText(HWND window) {
    std::wstring result(GetWindowTextLengthW(window) + 1, L'\0');
    const auto size = GetWindowTextW(window, result.data(), static_cast<int>(result.size()));
    result.resize(size); return result;
}
void selectAndDelete(HDC dc, HGDIOBJ object, HGDIOBJ prior) {
    SelectObject(dc, prior); DeleteObject(object);
}
}
void drawFolderTab(const DRAWITEMSTRUCT& item, bool active, HFONT font) {
    RECT bounds = item.rcItem;
    HBRUSH background = CreateSolidBrush(dsl::theme::background);
    FillRect(item.hDC, &bounds, background);
    DeleteObject(background);

    const int left = bounds.left;
    const int top = bounds.top;
    const int right = bounds.right - 1;
    const int bottom = bounds.bottom - 1;
    const int diagonal = 12;
    const int radius = 8;
    const COLORREF fillColor = active ? dsl::theme::panelRaised : dsl::theme::panel;
    const COLORREF borderColor = active ? dsl::theme::accent : dsl::theme::border;

    BeginPath(item.hDC);
    MoveToEx(item.hDC, left, bottom, nullptr);
    LineTo(item.hDC, left, top + radius);
    POINT shoulder[]{{left, top + 4}, {left + 3, top + 1},
                     {left + radius, top + 1}};
    PolyBezierTo(item.hDC, shoulder, 3);
    LineTo(item.hDC, right - diagonal - 4, top + 1);
    LineTo(item.hDC, right - diagonal, top + 3);
    LineTo(item.hDC, right, bottom);
    LineTo(item.hDC, left, bottom);
    CloseFigure(item.hDC);
    EndPath(item.hDC);
    HBRUSH fill = CreateSolidBrush(fillColor);
    HPEN noPen = CreatePen(PS_NULL, 0, fillColor);
    HGDIOBJ oldBrush = SelectObject(item.hDC, fill);
    HGDIOBJ oldPen = SelectObject(item.hDC, noPen);
    FillPath(item.hDC);
    SelectObject(item.hDC, oldBrush);
    SelectObject(item.hDC, oldPen);
    DeleteObject(fill);
    DeleteObject(noPen);

    HPEN border = CreatePen(PS_SOLID, 1, borderColor);
    oldPen = SelectObject(item.hDC, border);
    MoveToEx(item.hDC, left, bottom, nullptr);
    LineTo(item.hDC, left, top + radius);
    PolyBezierTo(item.hDC, shoulder, 3);
    LineTo(item.hDC, right - diagonal - 4, top + 1);
    LineTo(item.hDC, right - diagonal, top + 3);
    LineTo(item.hDC, right, bottom);
    if (!active) LineTo(item.hDC, left, bottom);
    selectAndDelete(item.hDC, border, oldPen);

    SetBkMode(item.hDC, TRANSPARENT);
    SetTextColor(item.hDC, dsl::theme::text);
    SelectObject(item.hDC, font);
    RECT text = bounds;
    text.left += radius;
    text.right -= diagonal;
    DrawTextW(item.hDC, windowText(item.hwndItem).c_str(), -1, &text,
              DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
}

}

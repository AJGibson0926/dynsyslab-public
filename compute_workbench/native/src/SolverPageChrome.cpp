#include "SolverPageChrome.hpp"
#include "dynsyslab_viewer/FolderTab.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "dynsyslab_viewer/Theme.hpp"

namespace dsl::compute {
namespace {

void selectAndDelete(HDC dc, HGDIOBJ object, HGDIOBJ prior) {
    SelectObject(dc, prior);
    DeleteObject(object);
}

}  // namespace

void drawSolverFolderTab(const DRAWITEMSTRUCT& item, bool active, HFONT font) {
    dsl::drawFolderTab(item, active, font);
}

void drawSolverEditorFrame(HDC dc, const RECT& requestedBounds) {
    RECT bounds = requestedBounds;
    HBRUSH fill = CreateSolidBrush(dsl::theme::background);
    FillRect(dc, &bounds, fill);
    DeleteObject(fill);

    HPEN border = CreatePen(PS_SOLID, 1, dsl::theme::border);
    HGDIOBJ oldPen = SelectObject(dc, border);
    HGDIOBJ oldBrush = SelectObject(dc, GetStockObject(NULL_BRUSH));
    Rectangle(dc, bounds.left, bounds.top, bounds.right, bounds.bottom);
    SelectObject(dc, oldBrush);
    selectAndDelete(dc, border, oldPen);
}

}  // namespace dsl::compute

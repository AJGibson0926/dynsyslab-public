#include "dynsyslab_compute/Workbench.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "SolverPageChrome.hpp"

namespace dsl::compute {
namespace {

constexpr int kRefreshId = 5003;
constexpr int kSolverTabId = 5004;
constexpr int kPostTabId = 5005;
constexpr int kHelpId = 5007;
constexpr int kStageSystemId = 5017;
constexpr int kStageInitialConditionsId = 5018;
constexpr int kStageControlId = 5019;
constexpr int kStageIntegratorId = 5020;
constexpr int kStageReviewId = 5021;

bool isStage(int id) {
    return id >= kStageSystemId && id <= kStageReviewId;
}

SolverStage stageFor(int id) {
    if (id == kStageInitialConditionsId) return SolverStage::initialConditions;
    if (id == kStageControlId) return SolverStage::control;
    if (id == kStageIntegratorId) return SolverStage::integrator;
    if (id == kStageReviewId) return SolverStage::review;
    return SolverStage::system;
}

}  // namespace

void Workbench::drawButton(const DRAWITEMSTRUCT& item) {
    if (item.CtlID == kSolverSystemPageId || item.CtlID == kSolverControlPageId) {
        drawSolverFolderTab(item,
            solver_.tabButtonActive(static_cast<int>(item.CtlID)), fonts_.body);
        return;
    }
    const bool pressed = (item.itemState & ODS_SELECTED) != 0;
    const bool enabled = (item.itemState & ODS_DISABLED) == 0;
    const bool active = (item.CtlID == kSolverTabId && solverActive_)
        || (item.CtlID == kPostTabId && !solverActive_)
        || (isStage(static_cast<int>(item.CtlID))
            && stageFor(static_cast<int>(item.CtlID)) == activeSolverStage_)
        || solver_.tabButtonActive(static_cast<int>(item.CtlID));
    const bool primary = item.CtlID == kSolverRunId || item.CtlID == kPostRunId;
    const COLORREF fill = item.CtlID == kSolverAddControlId
        ? dsl::theme::background
        : !enabled ? dsl::theme::controlBackground
        : pressed ? dsl::theme::accent
        : primary ? dsl::theme::blue
        : active ? dsl::theme::selection : dsl::theme::panelRaised;
    HBRUSH brush = CreateSolidBrush(fill);
    FillRect(item.hDC, &item.rcItem, brush);
    DeleteObject(brush);
    if (item.CtlID != kSolverAddControlId) {
        HBRUSH border = CreateSolidBrush(
            primary ? dsl::theme::blue : active ? dsl::theme::accent : dsl::theme::border);
        FrameRect(item.hDC, &item.rcItem, border);
        DeleteObject(border);
    }
    SetBkMode(item.hDC, TRANSPARENT);
    SetTextColor(item.hDC, !enabled ? dsl::theme::disabled
        : item.CtlID == kSolverAddControlId ? dsl::theme::accentStrong
        : (pressed || primary) ? dsl::theme::panelRaised
        : active ? dsl::theme::accentStrong : dsl::theme::text);
    SelectObject(item.hDC, fonts_.body);
    std::wstring label = windowText(item.hwndItem);
    if (item.CtlID == kRefreshId) label = L"\u21BB";
    else if (item.CtlID == kHelpId) label = L"?";
    RECT area = item.rcItem;
    DrawTextW(item.hDC, label.c_str(), -1, &area,
        (item.CtlID == kSolverAddControlId ? DT_LEFT : DT_CENTER)
        | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS | DT_NOPREFIX);
    if (item.itemState & ODS_FOCUS) {
        RECT focusBounds = item.rcItem;
        InflateRect(&focusBounds, -scalePixels(window_, 3),
            -scalePixels(window_, 3));
        HPEN focusPen = CreatePen(PS_SOLID, scalePixels(window_, 2),
            dsl::theme::focus);
        HGDIOBJ oldPen = SelectObject(item.hDC, focusPen);
        HGDIOBJ oldBrush = SelectObject(item.hDC, GetStockObject(NULL_BRUSH));
        Rectangle(item.hDC, focusBounds.left, focusBounds.top,
            focusBounds.right, focusBounds.bottom);
        SelectObject(item.hDC, oldBrush);
        SelectObject(item.hDC, oldPen);
        DeleteObject(focusPen);
    }
}

}  // namespace dsl::compute

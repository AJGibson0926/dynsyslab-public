#include "dynsyslab_compute/SolverPanel.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"
#include "SolverPageChrome.hpp"

#include <algorithm>

namespace dsl::compute {
namespace {

void visible(HWND window, bool show) {
    // A hidden root makes IsWindowVisible false even for locally visible
    // controls. Honor the requested page state before the root is shown too.
    if (((GetWindowLongPtrW(window,GWL_STYLE)&WS_VISIBLE)!=0) != show)
        ShowWindow(window, show ? SW_SHOW : SW_HIDE);
}

}  // namespace

bool SolverPanel::hasControl() const {
    return controller_ && selectedText(controller_) != "none";
}

bool SolverPanel::tabButtonActive(int id) const {
    if (!hasControl()) return false;
    return (id == kSolverSystemPageId && !controlPageActive_)
        || (id == kSolverControlPageId && controlPageActive_);
}

void SolverPanel::setPanelPage(bool control) {
    controlPageActive_ = control && hasControl();
    updatePanelPageVisibility();
    if (layoutWidth_ > 0)
        layoutEditor(layoutX_, layoutY_, layoutWidth_, layoutHeight_);
    RECT upper{layoutX_, layoutY_, layoutX_ + layoutWidth_,
        layoutY_ + layoutHeight_};
    RedrawWindow(parent_, &upper, nullptr,
        RDW_INVALIDATE | RDW_ALLCHILDREN | RDW_UPDATENOW);
}

void SolverPanel::updatePanelPageVisibility() {
    for (HWND item : controls_) visible(item, false);
    if (!editorVisible_) return;
    const bool controlled = hasControl();
    visible(addControl_, !controlled);
    visible(systemPage_, controlled); visible(controlPage_, controlled);
    if (controlled && controlPageActive_) {
        for (HWND item : {controllerSection_, controllerSeparator_, labels_[13],
             controller_}) visible(item, true);
        updateControllerControls();
        return;
    }
    for (HWND item : {systemSection_, systemSeparator_, labels_[1],
         system_, labels_[2], preset_, parameterSummary_, longParameters_,
         labels_[4], coordinates_, labels_[5], start_,
         labels_[6], stop_, labels_[7], step_, timeLabelDash_,
         timeRangeColon_, timeValueDash_, timeStepColon_}) visible(item, true);
    const bool pde = systemIndex() < registry_->systems.size()
        && registry_->systems[systemIndex()].family == "pde";
    for (HWND item : {pdePartition_, pdePartitionShapeLabel_, pdePartitionShape_,
         pdeGhostWidthLabel_, pdeGhostWidth_, pdeLocalWorkersLabel_,
         pdeLocalWorkers_}) visible(item, pde);
}

void SolverPanel::showEditor(bool visible) {
    editorVisible_ = visible;
    updatePanelPageVisibility();
    if (layoutWidth_ > 0) {
        RECT upper{layoutX_, layoutY_, layoutX_ + layoutWidth_,
            layoutY_ + layoutHeight_};
        RedrawWindow(parent_, &upper, nullptr,
            RDW_INVALIDATE | RDW_ALLCHILDREN);
    }
}

void SolverPanel::showExecution(bool visible, bool enabled) {
    executionVisible_ = visible;
    executionEnabled_ = enabled;
    execution_.show(visible, enabled);
}

void SolverPanel::layoutEditor(int x, int y, int width, int height) {
    layoutX_ = x; layoutY_ = y; layoutWidth_ = width; layoutHeight_ = height;
    const int gap = scalePixels(parent_, 8);
    const int inner = scalePixels(parent_, 10);
    const int labelH = scalePixels(parent_, 15);
    const int headerH = scalePixels(parent_, 19);
    const int rowH = scalePixels(parent_, 27);
    const int tabH = scalePixels(parent_, 29);
    const int tabW = scalePixels(parent_, 104);
    auto inlineField = [&](int label, HWND window, int px, int py, int pw,
                           int captionW) {
        MoveWindow(labels_[static_cast<std::size_t>(label)], px,
            py + (rowH - labelH) / 2, captionW, labelH, TRUE);
        moveFieldControl(window, px + captionW, py, pw - captionW, rowH);
    };

    const bool controlled = hasControl();
    if (controlled) {
        MoveWindow(systemPage_, x, y, tabW, tabH, TRUE);
        MoveWindow(controlPage_, x + tabW + scalePixels(parent_, 3), y,
            tabW, tabH, TRUE);
    } else {
        MoveWindow(addControl_, x + inner, y, scalePixels(parent_, 104), tabH, TRUE);
    }
    const int frameY = y + tabH - 1;
    editorBounds_ = RECT{x, frameY, x + width, y + height};
    const int bodyX = x + inner;
    const int bodyY = frameY + inner;
    const int bodyW = std::max(1, width - 2 * inner);
    const int bodyH = std::max(1, y + height - inner - bodyY);

    if (controlled && controlPageActive_) {
        layoutControllerArea(bodyX, bodyY, bodyW, bodyH);
        return;
    }

    MoveWindow(systemSection_, bodyX, bodyY, bodyW, headerH, TRUE);
    const int systemRow = bodyY + headerH + scalePixels(parent_, 3);
    const int systemW = (bodyW - gap) / 2;
    const int presetW = bodyW - systemW - gap;
    inlineField(1, system_, bodyX, systemRow, systemW,
        scalePixels(parent_, 70));
    inlineField(2, preset_, bodyX + systemW + gap, systemRow, presetW,
        scalePixels(parent_, 58));
    const int gridY = systemRow + rowH + gap;
    const int longH = scalePixels(parent_, 18);
    const int longW = scalePixels(parent_, 52);
    const int timeGap = scalePixels(parent_, 7);
    const int timeBlockH = rowH;
    const bool pde = systemIndex() < registry_->systems.size()
        && registry_->systems[systemIndex()].family == "pde";
    const int pdeH = pde ? rowH + timeGap : 0;
    const int gridH = std::max(scalePixels(parent_, 34),
        bodyY + bodyH - gridY - longH - timeGap - timeBlockH
            - pdeH - scalePixels(parent_, 2));
    MoveWindow(parameterSummary_, bodyX, gridY, bodyW, gridH, TRUE);
    MoveWindow(longParameters_, bodyX + bodyW - longW, gridY + gridH,
        longW, longH, TRUE);
    int timeRow = gridY + gridH + longH + timeGap;
    if (pde) {
        const int toggleW = bodyW * 24 / 100;
        MoveWindow(pdePartition_, bodyX, timeRow, toggleW, rowH, TRUE);
        const int remainingX = bodyX + toggleW + gap;
        const int remainingW = bodyW - toggleW - gap;
        const int third = (remainingW - 2 * gap) / 3;
        auto pdeField = [&](HWND label, HWND field, int px, int pw,
                            int captionW) {
            MoveWindow(label, px, timeRow + (rowH - labelH) / 2,
                captionW, labelH, TRUE);
            moveFieldControl(field, px + captionW, timeRow,
                pw - captionW, rowH);
        };
        pdeField(pdePartitionShapeLabel_, pdePartitionShape_, remainingX,
            third, scalePixels(parent_, 78));
        pdeField(pdeGhostWidthLabel_, pdeGhostWidth_, remainingX + third + gap,
            third, scalePixels(parent_, 62));
        pdeField(pdeLocalWorkersLabel_, pdeLocalWorkers_,
            remainingX + 2 * (third + gap), remainingW - 2 * (third + gap),
            scalePixels(parent_, 48));
        timeRow += rowH + timeGap;
    }
    const int timeGroupW = std::min(scalePixels(parent_, 404),
        std::max(scalePixels(parent_, 330), bodyW * 64 / 100));
    const int timeGroupX = bodyX + bodyW - timeGroupW;
    const int coordinateW = std::max(scalePixels(parent_, 150),
        timeGroupX - bodyX - gap);
    inlineField(4, coordinates_, bodyX, timeRow, coordinateW,
        scalePixels(parent_, 70));

    const int symbolW = scalePixels(parent_, 25);
    const int finalSymbolW = scalePixels(parent_, 31);
    const int punctuationW = scalePixels(parent_, 13);
    const int valueW = scalePixels(parent_, 72);
    const int expressionGap = scalePixels(parent_, 18);
    const int deltaW = scalePixels(parent_, 37);
    int fieldX = timeGroupX;
    MoveWindow(labels_[5], fieldX, timeRow, symbolW, rowH, TRUE);
    fieldX += symbolW;
    MoveWindow(timeLabelDash_, fieldX, timeRow, punctuationW, rowH, TRUE);
    fieldX += punctuationW;
    MoveWindow(labels_[6], fieldX, timeRow, finalSymbolW, rowH, TRUE);
    fieldX += finalSymbolW;
    MoveWindow(timeRangeColon_, fieldX, timeRow, punctuationW, rowH, TRUE);
    fieldX += punctuationW;
    MoveWindow(start_, fieldX, timeRow, valueW, rowH, TRUE);
    fieldX += valueW;
    MoveWindow(timeValueDash_, fieldX, timeRow, punctuationW, rowH, TRUE);
    fieldX += punctuationW;
    MoveWindow(stop_, fieldX, timeRow, valueW, rowH, TRUE);
    fieldX += valueW + expressionGap;
    MoveWindow(labels_[7], fieldX, timeRow, deltaW, rowH, TRUE);
    fieldX += deltaW;
    MoveWindow(timeStepColon_, fieldX, timeRow, punctuationW, rowH, TRUE);
    fieldX += punctuationW;
    MoveWindow(step_, fieldX, timeRow,
        std::max(1, bodyX + bodyW - fieldX), rowH, TRUE);
    MoveWindow(systemSeparator_, bodyX, bodyY + bodyH, bodyW, 1, TRUE);
}

void SolverPanel::paintEditorFrame(HDC dc) const {
    if (editorVisible_) drawSolverEditorFrame(dc, editorBounds_);
}

void SolverPanel::layoutExecution(int x, int y, int width, int height) {
    executionX_ = x; executionY_ = y;
    executionWidth_ = width; executionHeight_ = height;
    execution_.layout(x, y, width, height);
}

void SolverPanel::activateStage(SolverStage stage) {
    switch (stage) {
    case SolverStage::system:
        setPanelPage(false);
        if (system_) SetFocus(system_);
        break;
    case SolverStage::initialConditions:
        setPanelPage(false);
        editInitialConditions();
        break;
    case SolverStage::control:
        if (!hasControl()) {
            SendMessageW(controller_, CB_SETCURSEL, 1, 0);
            controlPageActive_ = true;
            changeControllerKind();
        } else {
            setPanelPage(true);
        }
        if (controller_) SetFocus(controller_);
        break;
    case SolverStage::integrator:
        execution_.focus();
        break;
    case SolverStage::review:
        if (plan_) SetFocus(plan_);
        break;
    }
}

}  // namespace dsl::compute

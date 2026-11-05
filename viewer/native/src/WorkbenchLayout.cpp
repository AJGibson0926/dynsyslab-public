#include "dynsyslab_viewer/Workbench.hpp"

#include <algorithm>
#include <vector>

namespace dsl {
namespace {

struct Placement {
    HWND window{};
    int x{}, y{}, width{}, height{};
    UINT visibility{};
};

void applyLayout(HWND parent, const std::vector<Placement>& placements) {
    constexpr UINT common = SWP_NOZORDER | SWP_NOACTIVATE | SWP_NOOWNERZORDER |
                            SWP_NOREDRAW;
    HDWP transaction = BeginDeferWindowPos(static_cast<int>(placements.size()));
    for (const auto& item : placements) {
        if (!transaction) break;
        transaction = DeferWindowPos(transaction, item.window, nullptr,
            item.x, item.y, std::max(1, item.width), std::max(1, item.height),
            common | item.visibility);
    }
    const bool committed = transaction && EndDeferWindowPos(transaction);
    if (!committed) {
        for (const auto& item : placements)
            SetWindowPos(item.window, nullptr, item.x, item.y,
                std::max(1, item.width), std::max(1, item.height),
                (common & ~SWP_NOREDRAW) | item.visibility);
    }
    RedrawWindow(parent, nullptr, nullptr,
        RDW_INVALIDATE | RDW_ALLCHILDREN | RDW_UPDATENOW);
}

}  // namespace

void Workbench::layout() {
    if (!window_) return;
    RECT area{};
    GetClientRect(window_, &area);
    const int width = area.right;
    const int height = area.bottom;
    const auto ui = theme::metrics(window_);
    const int margin = ui.margin;
    const int headerHeight = ui.headerHeight;
    const int toolbarHeight = ui.toolbarHeight;
    const int plotToolbarHeight = toolbarHeight * plotDragToolbar_.rowCount() + ui.compactGap;
    const int searchHeight = ui.toolbarHeight;
    const int sortHeight = ui.toolbarHeight;
    const int tabHeight = theme::scalePixels(window_, 34);
    const int viewTop = headerHeight + tabHeight + ui.compactGap;
    const int searchGap = ui.compactGap;
    const int bottomMargin = ui.margin;
    const int topControlY = theme::scalePixels(window_, 8);
    const int topControlHeight = ui.toolbarHeight;
    constexpr int splitterWidth = 8;
    HorizontalPaneConstraints constraints;
    constraints.totalWidth = width;
    constraints.margin = margin;
    constraints.desiredNavigationWidth = navigationWidth_;
    constraints.desiredInspectorWidth = inspectorWidth_;
    constraints.minimumNavigationWidth = 210;
    constraints.maximumNavigationWidth = 500;
    constraints.minimumCenterWidth = 420;
    constraints.minimumInspectorWidth = 250;
    constraints.maximumInspectorWidth = 520;
    constraints.inspectorThreshold = 1180;
    const auto panes = resolveHorizontalPaneLayout(constraints);
    const int leftColumnWidth = panes.navigationWidth +
        (panes.inspectorVisible ? constraints.gap + panes.inspectorWidth : 0);
    const int centerX = margin + leftColumnWidth + constraints.gap;
    const int centerWidth = std::max(0, width - centerX - margin);

    const int viewWidth = std::max(theme::scalePixels(window_, 100),
        std::min(theme::scalePixels(window_, 250),
            centerWidth - theme::scalePixels(window_, 252)));
    const int contentHeight = height - headerHeight - bottomMargin;
    const int inspectorGap = ui.compactGap;
    const int inspectorHeight = panes.inspectorVisible
        ? std::clamp(contentHeight * 28 / 100,
            theme::scalePixels(window_, 190), theme::scalePixels(window_, 340)) : 0;
    const int treeHeight = contentHeight - searchHeight - searchGap
        - sortHeight - searchGap
        - (panes.inspectorVisible ? inspectorGap + inspectorHeight : 0);
    const int searchY = headerHeight;
    const int sortY = searchY + searchHeight + searchGap;
    const int treeY = sortY + sortHeight + searchGap;
    const int inspectorY = treeY + treeHeight + inspectorGap;
    const int refreshWidth = theme::scalePixels(window_, 40);
    const int sortWidth = theme::scalePixels(window_, 40);
    const int sortControlGap = theme::scalePixels(window_, 6);
    const int helpWidth = theme::scalePixels(window_, 110);
    const int saveWidth = theme::scalePixels(window_, 86);
    const int copyWidth = theme::scalePixels(window_, 88);
    const int axesWidth = theme::scalePixels(window_, 64);
    const int addPanelWidth = theme::scalePixels(window_, 78);
    const int compareWidth = theme::scalePixels(window_, 76);
    const int computeWidth = theme::scalePixels(window_, 104);
    const int actionGap = theme::scalePixels(window_, 8);
    int actionRight = width - margin;
    auto actionX = [&](int itemWidth) {
        actionRight -= itemWidth;
        const int result = actionRight;
        actionRight -= actionGap;
        return result;
    };
    const int helpX = actionX(helpWidth);
    const int saveX = actionX(saveWidth);
    const int copyX = actionX(copyWidth);
    const int axesX = actionX(axesWidth);
    const int addPanelX = actionX(addPanelWidth);
    const int compareX = actionX(compareWidth);
    const int computeX = actionX(computeWidth);
    const int pageWidth = theme::scalePixels(window_, 170);
    const int pageRightReserve = theme::scalePixels(window_, 244);
    const int previousWidth = theme::scalePixels(window_, 30);
    const int nextWidth = theme::scalePixels(window_, 30);
    const int previousOffset = theme::scalePixels(window_, 68);
    const int nextOffset = theme::scalePixels(window_, 33);
    const bool krom = selected_ >= 0 && catalog_.items[static_cast<std::size_t>(selected_)].kind == "krom";
    const UINT pageVisibility = plots_.verticallyScrollable() ? SWP_HIDEWINDOW : SWP_SHOWWINDOW;
    const UINT inspectorVisibility = panes.inspectorVisible ? SWP_SHOWWINDOW : SWP_HIDEWINDOW;
    const int leftSplitterCenter = margin + leftColumnWidth + constraints.gap / 2;
    applyLayout(window_, {
        {help_, helpX, topControlY, helpWidth, topControlHeight},
        {save_, saveX, topControlY, saveWidth, topControlHeight},
        {copy_, copyX, topControlY, copyWidth, topControlHeight},
        {axes_, axesX, topControlY, axesWidth, topControlHeight},
        {addPanel_, addPanelX, topControlY, addPanelWidth, topControlHeight},
        {compare_, compareX, topControlY, compareWidth, topControlHeight},
        {compute_, computeX, topControlY, computeWidth, topControlHeight},
        {sortDate_, margin, sortY, sortWidth, sortHeight},
        {sortTrajectories_, margin + sortWidth + sortControlGap, sortY,
            sortWidth, sortHeight},
        {sortSize_, margin + (sortWidth + sortControlGap) * 2, sortY,
            sortWidth, sortHeight},
        {navigation_.handle(), margin, treeY, leftColumnWidth, treeHeight},
        {search_, margin, searchY,
            leftColumnWidth - refreshWidth - searchGap, searchHeight},
        {refresh_, margin + leftColumnWidth - refreshWidth, searchY,
            refreshWidth, searchHeight},
        {leftSplitter_.handle(), leftSplitterCenter - splitterWidth / 2,
            headerHeight, splitterWidth,
            height - headerHeight - bottomMargin},
        {phaseContext_, centerX + centerWidth - theme::scalePixels(window_, 120), headerHeight,
            theme::scalePixels(window_, 120), tabHeight, static_cast<UINT>(krom ? SWP_HIDEWINDOW : SWP_SHOWWINDOW)},
        {mainTab_, centerX, headerHeight, theme::scalePixels(window_, 118), tabHeight},
        {observablesTab_, centerX + theme::scalePixels(window_, 125), headerHeight,
            theme::scalePixels(window_, 166), tabHeight, static_cast<UINT>(krom ? SWP_HIDEWINDOW : SWP_SHOWWINDOW)},
        {selectionContext_, centerX + theme::scalePixels(window_, 304), headerHeight,
            std::max(1, centerWidth - theme::scalePixels(window_, 436)), tabHeight,
            static_cast<UINT>(krom ? SWP_HIDEWINDOW : SWP_SHOWWINDOW)},
        {viewSelect_, centerX, viewTop, viewWidth, toolbarHeight * 8},
        {pageSelect_, centerX + centerWidth - pageRightReserve, viewTop, pageWidth,
            toolbarHeight * 7, pageVisibility},
        {previous_, centerX + centerWidth - previousOffset, viewTop, previousWidth,
            toolbarHeight - theme::scalePixels(window_, 4), pageVisibility},
        {next_, centerX + centerWidth - nextOffset, viewTop, nextWidth,
            toolbarHeight - theme::scalePixels(window_, 4), pageVisibility},
        {plotDragToolbar_.handle(), centerX, viewTop + toolbarHeight,
            centerWidth, plotToolbarHeight},
        {plots_.handle(), centerX, viewTop + toolbarHeight + plotToolbarHeight,
            centerWidth, height - viewTop - toolbarHeight - plotToolbarHeight - bottomMargin},
        {inspector_.handle(), margin, inspectorY, leftColumnWidth,
            inspectorHeight, inspectorVisibility}
    });
}

}  // namespace dsl

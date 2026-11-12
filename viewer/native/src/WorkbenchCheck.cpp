#include "dynsyslab_viewer/Workbench.hpp"
#include <commctrl.h>
#include <gdiplus.h>
#include <algorithm>
#include <fstream>
#include <iterator>

namespace dsl {
namespace {

void renderHiddenControl(HWND window) {
    RECT area{}; GetClientRect(window, &area);
    Gdiplus::Bitmap bitmap(std::max(1L, area.right), std::max(1L, area.bottom), PixelFormat32bppARGB);
    Gdiplus::Graphics graphics(&bitmap);
    HDC dc = graphics.GetHDC();
    SendMessageW(window, WM_PRINTCLIENT, reinterpret_cast<WPARAM>(dc), PRF_CLIENT);
    graphics.ReleaseHDC(dc);
}

LRESULT CALLBACK checkTreePaint(HWND window, UINT message, WPARAM wparam,
                               LPARAM lparam, UINT_PTR, DWORD_PTR inspector) {
    // Reproduce the synchronous sibling paint seen in the refresh crash dump,
    // even though this workbench never becomes visible.
    if (message == TVM_DELETEITEM)
        renderHiddenControl(reinterpret_cast<HWND>(inspector));
    return DefSubclassProc(window, message, wparam, lparam);
}

}  // namespace

void Workbench::installBackgroundPaintCheck() {
    SetWindowSubclass(navigation_.handle(), checkTreePaint, 1,
        reinterpret_cast<DWORD_PTR>(inspector_.handle()));
}

void Workbench::finishBackgroundCheck(int status,const std::wstring& error) {
    if (status == 0 && arguments_.backgroundKromCheck && !advanceKromCheck()) return;
    if (status == 0 && arguments_.backgroundWorkspaceCheck && checkedWorkspaceStage_ < 2) {
        if (checkedWorkspaceStage_ == 0) {
            RECT toolbar{}, plot{}, tab{}, views{};
            GetWindowRect(plotDragToolbar_.handle(),&toolbar);GetWindowRect(plots_.handle(),&plot);
            GetWindowRect(mainTab_,&tab);GetWindowRect(viewSelect_,&views);
            if (toolbar.bottom > plot.top || tab.bottom > views.top) {
                finishBackgroundCheck(1,L"Workspace tabs or toolbar overlap plot content"); return;
            }
            if (!plots_.hasPayload() || mainSources_.empty() ||
                !supportsObservableWorkspace(mainSources_.front())) {
                finishBackgroundCheck(1,L"No observable workspace for the selected Main trajectory"); return;
            }
            const auto axes = plots_.currentAxisPanels();
            if (!axes.empty()) plots_.setAxisLimits(0,axes.front().automatic);
            checkedMainWorkspace_ = plots_.captureState();
            saveBackgroundWorkspaceSnapshot(L"workspace-main.png");
            checkedWorkspaceStage_ = 1;
            SendMessageW(window_,WM_COMMAND,115,0);
            return;
        }
        const auto observableState = plots_.captureState();
        if (!observablesActive_ || !observableState.payload || observableState.payload->figures.empty()
            || std::any_of(observableState.payload->figures.begin(), observableState.payload->figures.end(),
                [](const auto& figure) { return figure.kind != "observable_time" && figure.kind != "observable_heatmap"; })) {
            finishBackgroundCheck(1,L"Observables tab did not load its time-series and phase-space panels"); return;
        }
        if (!plots_.verticallyScrollable() || plots_.pageCount() != 1
            || (GetWindowLongPtrW(pageSelect_,GWL_STYLE)&WS_VISIBLE)
            || (GetWindowLongPtrW(previous_,GWL_STYLE)&WS_VISIBLE)
            || (GetWindowLongPtrW(next_,GWL_STYLE)&WS_VISIBLE)) {
            finishBackgroundCheck(1,L"Observables must scroll without paging controls"); return;
        }
        checkedObservableFigures_ = observableState.payload->figures.size();
        saveBackgroundWorkspaceSnapshot(L"workspace-observables.png");
        plots_.setScrollPosition(plots_.scrollMaximum());
        saveBackgroundWorkspaceSnapshot(L"workspace-observables-bottom.png");
        const auto scroll = plots_.scrollPosition();
        const auto visible = plots_.visibleFigureIndices();
        if (std::find(visible.begin(),visible.end(),checkedObservableFigures_-1) == visible.end()) {
            finishBackgroundCheck(1,L"Final observable cannot be reached by scrolling"); return;
        }
        renderHiddenControl(plots_.handle());
        SendMessageW(window_,WM_COMMAND,114,0);
        const auto restored = plots_.captureState();
        if (observablesActive_ || restored.payload != checkedMainWorkspace_.payload
            || restored.axisLimits.size() != checkedMainWorkspace_.axisLimits.size()
            || restored.page != checkedMainWorkspace_.page) {
            finishBackgroundCheck(1,L"Main workspace changed after switching tabs"); return;
        }
        SendMessageW(window_,WM_COMMAND,115,0);
        if (taskActive_ || plots_.captureState().payload != observableState.payload || plots_.scrollPosition() != scroll) {
            finishBackgroundCheck(1,L"Observables did not preserve their workspace"); return;
        }
        SendMessageW(window_,WM_COMMAND,114,0);
        checkedWorkspaceStage_ = 2;
    }
    if (status == 0) {
        renderHiddenControl(inspector_.handle());
        if (plots_.hasPayload()) renderHiddenControl(plots_.handle());
        if (checkedRefreshes_ < arguments_.backgroundRefreshCount) {
            ++checkedRefreshes_;
            refreshCatalog();
            return;
        }
    }
    json::Value::Object result{{"schema","dynsyslab.viewer.workbench-check/v1"},
        {"process_id",static_cast<double>(GetCurrentProcessId())},
        {"status",static_cast<double>(status)},{"error",narrow(error)},
        {"catalog_items",static_cast<double>(catalog_.items.size())},
        {"navigation_nodes",static_cast<double>(navigation_.nodeCount())},
        {"selected_id",selected_ >= 0 ? catalog_.items[static_cast<std::size_t>(selected_)].id : ""},
        {"selected_label",selected_ >= 0 ? catalog_.items[static_cast<std::size_t>(selected_)].label : ""},
        {"training_data_selected",isTrainingDataItem(catalog_,selected_)},
        {"focus_label",catalog_.focusLabel},{"hidden_items",static_cast<double>(catalog_.hiddenItemCount)},
        {"focus_path",arguments_.backend.archiveFocus.generic_u8string()},
        {"refreshes",static_cast<double>(checkedRefreshes_)},
        {"workspace_tabs_checked",checkedWorkspaceStage_ == 2},
        {"observable_figures",static_cast<double>(checkedObservableFigures_)},
        {"krom_stages_checked",static_cast<double>(checkedKromStage_)},
        {"visible",IsWindowVisible(window_)!=FALSE},{"has_plot",plots_.hasPayload()}};
    std::filesystem::create_directories(arguments_.backend.cache);
    RECT sortBounds{}, treeBounds{};
    GetWindowRect(sortDate_, &sortBounds);
    GetWindowRect(navigation_.handle(), &treeBounds);
    wchar_t dateName[64]{}, trajectoryName[64]{}, sizeName[64]{};
    GetWindowTextW(sortDate_, dateName, static_cast<int>(std::size(dateName)));
    GetWindowTextW(sortTrajectories_, trajectoryName,
        static_cast<int>(std::size(trajectoryName)));
    GetWindowTextW(sortSize_, sizeName, static_cast<int>(std::size(sizeName)));
    const bool sortControlsChecked = std::wstring(dateName) == L"Sort by date" &&
        std::wstring(trajectoryName) == L"Sort by trajectory count" &&
        std::wstring(sizeName) == L"Sort by total size" &&
        treeSort_ == CatalogTreeSort::date &&
        SendMessageW(sortTooltip_, TTM_GETTOOLCOUNT, 0, 0) >= 3 &&
        sortBounds.bottom <= treeBounds.top;
    result["tree_sort_controls_checked"] = sortControlsChecked;
    result["tree_sort_icon_lengths"] = static_cast<double>(
        3);
    result["tree_sort_tooltip_count"] = static_cast<double>(
        SendMessageW(sortTooltip_, TTM_GETTOOLCOUNT, 0, 0));
    result["tree_sort_date_selected"] = treeSort_ == CatalogTreeSort::date;
    result["tree_sort_above_tree"] = sortBounds.bottom <= treeBounds.top;
    std::ofstream output(arguments_.backend.cache / L"workbench-check.json",std::ios::binary);
    output << json::dump(result);
    PostQuitMessage(status);
}
}

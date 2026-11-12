#include "dynsyslab_viewer/Workbench.hpp"
#include "dynsyslab_viewer/ViewerSession.hpp"
#include <fstream>
#include <algorithm>
namespace dsl {
bool Workbench::advanceKromCheck() {
    if (checkedKromStage_ >= 5) return true;
    const auto state = plots_.captureState();
    auto fail = [&](const wchar_t* message) { finishBackgroundCheck(1,message); return false; };
    if (!state.payload || !plots_.isKrom()) return fail(L"Selected KROM did not load the split workspace");
    const auto& payload = *state.payload;
    const auto modes = payload.figures.front().traces.size()-1;
    if (checkedKromStage_ == 0) {
        RECT toolbar{},plot{}; GetWindowRect(plotDragToolbar_.handle(),&toolbar);GetWindowRect(plots_.handle(),&plot);
        if (toolbar.bottom > plot.top) return fail(L"KROM radio toolbar overlaps plot content");
        if ((GetWindowLongPtrW(observablesTab_,GWL_STYLE)&WS_VISIBLE)
            || (GetWindowLongPtrW(pageSelect_,GWL_STYLE)&WS_VISIBLE)
            || (GetWindowLongPtrW(previous_,GWL_STYLE)&WS_VISIBLE)
            || (GetWindowLongPtrW(next_,GWL_STYLE)&WS_VISIBLE)) return fail(L"KROM exposes Observables or paging controls");
        wchar_t tab[32]{};GetWindowTextW(mainTab_,tab,32);
        if (std::wstring(tab) != L"KROM") return fail(L"KROM workspace tab has the wrong label");
        for (int id : {5,6,7,10,11,12,13}) if ((GetWindowLongPtrW(GetDlgItem(plotDragToolbar_.handle(),id),GWL_STYLE)&BS_TYPEMASK) != BS_AUTORADIOBUTTON)
            return fail(L"Complex-display buttons are not round radio controls");
        checkedKromStage_ = 1;
        SendMessageW(GetDlgItem(plotDragToolbar_.handle(),8),BM_CLICK,0,0);
        return false;
    }
    if (checkedKromStage_ <= 3) {
        const auto expected = checkedKromStage_ == 3 ? "domain-coloring" : checkedKromStage_ == 2 ? "magnitude-phase" : "real-imaginary";
        if (payload.complexDisplay != expected || payload.selectedModes.size() != modes)
            return fail(L"KROM radio change or All selection did not rebuild the workspace");
        const auto components = checkedKromStage_ == 3 ? 1u : 2u;
        if (payload.figures.size() != 1 + modes*components*2)
            return fail(L"KROM time/map component count is wrong for the selected representation");
        const auto stem = L"krom-"+std::to_wstring(checkedKromStage_);
        plots_.setScrollPosition(0);
        saveBackgroundWorkspaceSnapshot((stem+L".png").c_str());
        if (plots_.pageCount() != 1 || !plots_.verticallyScrollable()) return fail(L"KROM workspace is paginated");
        plots_.setScrollPosition(plots_.scrollMaximum());
        saveBackgroundWorkspaceSnapshot((stem+L"-bottom.png").c_str());
        const auto visible = plots_.visibleFigureIndices();
        if (std::find(visible.begin(),visible.end(),payload.figures.size()-1) == visible.end())
            return fail(L"Final KROM eigenfunction cannot be reached by scrolling");
        plots_.setScrollPosition(0);
        std::ofstream output(arguments_.backend.cache/(stem+L".json"),std::ios::binary);
        output<<json::dump(inspectPlot(payload));
        ++checkedKromStage_;
        SendMessageW(GetDlgItem(plotDragToolbar_.handle(),checkedKromStage_ == 2 ? 6 : checkedKromStage_ == 3 ? 7 : 9),BM_CLICK,0,0);
        return false;
    }
    if (payload.figures.size() != 1 || !payload.selectedModes.empty()) return fail(L"None should leave the spectrum visible");
    checkedKromStage_ = 5; return true;
}
}

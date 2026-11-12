#include "dynsyslab_viewer/Workbench.hpp"
#include <algorithm>
namespace dsl {
void Workbench::resetWorkspace() {
    observablesActive_ = false;
    mainWorkspace_ = {}; observablesWorkspace_ = {}; mainSources_.clear(); observablesDomain_.reset();
    updateWorkspaceTabs();
}
void Workbench::updateWorkspaceTabs() {
    const bool krom = selected_ >= 0 && catalog_.items[static_cast<std::size_t>(selected_)].kind == "krom";
    SetWindowTextW(mainTab_,krom ? L"KROM" : L"Main");
    ShowWindow(observablesTab_,krom ? SW_HIDE : SW_SHOWNA);
    const bool available = !mainSources_.empty() && std::all_of(mainSources_.begin(), mainSources_.end(),
        [](const auto& item) { return supportsObservableWorkspace(item); });
    SetWindowTextW(observablesTab_, observablesOpened_ ? L"Observables" : L"+ observables");
    EnableWindow(mainTab_, !taskActive_);
    EnableWindow(observablesTab_, !taskActive_ && available);
    EnableWindow(viewSelect_, !taskActive_ && !observablesActive_ && !availableViews_.empty());
    InvalidateRect(mainTab_, nullptr, FALSE); InvalidateRect(observablesTab_, nullptr, FALSE);
}
void Workbench::activateWorkspace(bool observables) {
    if (taskActive_ || observables == observablesActive_ || plots_.isKrom()) return;
    if (observables) {
        if (mainSources_.empty() || !std::all_of(mainSources_.begin(), mainSources_.end(),
            [](const auto& item) { return supportsObservableWorkspace(item); })) return;
        const auto domain = plots_.phaseSpaceDomain();
        if (domain != observablesDomain_) observablesWorkspace_ = {};
        observablesDomain_ = domain;
        mainWorkspace_ = plots_.captureState();
    } else observablesWorkspace_ = plots_.captureState();
    observablesActive_ = observables;
    if (!observables) plots_.restoreState(mainWorkspace_);
    else if (observablesWorkspace_.payload) plots_.restoreState(observablesWorkspace_);
    else {
        observablesOpened_ = true;
        plots_.setMessage(L"Evaluating system observables", L"Using the same archived trajectories as Main…");
        const auto sources = mainSources_; const auto requestedIndex = selected_;
        const auto domain = observablesDomain_;
        taskActive_ = true;
        task_ = std::async(std::launch::async, [this, sources, requestedIndex, domain]() {
            AsyncResult result; result.itemIndex = requestedIndex;
            result.observables = true; result.view = "observable_workspace";
            try { result.payload = std::make_shared<Payload>(backend_.loadObservables(sources, domain)); }
            catch (const std::exception& error) { result.error = wide(error.what()); }
            return result;
        });
    }
    updateExportButtons(); updatePageSelector(); layout();
}
}

#include "dynsyslab_viewer/Workbench.hpp"
#include "PanelPickerDialog.hpp"
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <chrono>
#include <algorithm>
namespace dsl {
void Workbench::requestView(const std::string& view) {
    if (selected_ < 0) return;
    if (taskActive_) { pendingView_ = view; return; }
    pendingView_.clear(); pendingKrom_.reset();
    resetWorkspace();
    const CatalogItem item = catalog_.items[static_cast<std::size_t>(selected_)];
    activeView_ = view;
    if (view == "open") {
        plots_.setMessage(L"No native scientific view is available",
            L"Viewer does not hand archive objects to Windows file associations.");
        updateExportButtons();
        updatePageSelector();
        return;
    }
    plots_.setMessage(L"Loading " + wide(catalog_.viewLabels.count(view) ?
        catalog_.viewLabels.at(view) : view), L"Preparing a compact plot payload…");
    updateExportButtons();
    updatePageSelector();
    InvalidateRect(window_, nullptr, FALSE);
    taskActive_ = true;
    updateExportButtons();
    const int requestedIndex = selected_;
    task_ = std::async(std::launch::async, [this, item, view, requestedIndex]() {
        AsyncResult result; result.itemIndex = requestedIndex; result.view = view; result.sources = {item};
        try { result.payload = std::make_shared<Payload>(backend_.loadPayload(item, view)); }
        catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
}

void Workbench::requestKrom(KromWorkspaceOptions options) {
    if (selected_ < 0 || catalog_.items[selected_].kind != "krom") return;
    if (taskActive_) { pendingKrom_ = std::move(options); return; }
    pendingKrom_.reset();
    const auto item = catalog_.items[selected_]; const auto index = selected_;
    taskActive_ = true; updateExportButtons();
    task_ = std::async(std::launch::async,[this,item,index,options]() {
        AsyncResult result; result.itemIndex = index; result.sources = {item};
        try { result.payload = std::make_shared<Payload>(buildKromWorkspace(backend_.resolveSource(item),options)); }
        catch (const std::exception& e) { result.error = wide(e.what()); }
        return result;
    });
}

void Workbench::refreshCatalog() {
    if (taskActive_) return;
    plots_.setMessage(L"Refreshing archive catalogue",
        L"Reconciling the Archiver index with files on disk…");
    taskActive_ = true;
    task_ = std::async(std::launch::async, [this]() {
        AsyncResult result; result.kind = AsyncResult::Kind::Catalog;
        try { result.catalog = std::make_unique<Catalog>(backend_.refreshCatalog()); }
        catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
}

void Workbench::watchArchive() {
    if (taskActive_) return;
    std::uint64_t generation = 0;
    try { generation = dynsyslab::core::archiveIndexGeneration(arguments_.backend.archiveRoot); }
    catch (const std::exception&) { return; }  // no index yet: nothing to watch
    if (!archiveGeneration_ || generation == archiveGeneration_) { archiveGeneration_ = generation; return; }
    archiveGeneration_ = generation;
    if (selected_ >= 0 && static_cast<std::size_t>(selected_) < catalog_.items.size())
        arguments_.initialSource = backend_.resolveSource(catalog_.items[static_cast<std::size_t>(selected_)]);
    taskActive_ = true;
    task_ = std::async(std::launch::async, [this]() {
        AsyncResult result; result.kind = AsyncResult::Kind::Catalog;
        try { result.catalog = std::make_unique<Catalog>(backend_.loadCatalog()); }
        catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
}

void Workbench::pollTask() {
    if (!taskActive_ || task_.wait_for(std::chrono::milliseconds(0)) != std::future_status::ready) return;
    AsyncResult result = task_.get(); taskActive_ = false;
    if(result.phaseContext) phaseCancellation_.reset();
    if (!pendingView_.empty()) { const auto view = pendingView_; requestView(view); return; }
    if (pendingKrom_) { const auto options = *pendingKrom_; requestKrom(options); return; }
    if (!result.error.empty()) {
        plots_.setMessage(L"Unable to prepare this view", result.error);
        updateExportButtons();
        updatePageSelector();
        if (arguments_.backgroundCheck) finishBackgroundCheck(1,result.error);
        return;
    }
    if (result.kind == AsyncResult::Kind::Catalog) {
        if (!result.catalog) {
            plots_.setMessage(L"Unable to load archive catalogue",
                L"The Archiver index did not return a catalogue.");
            return;
        }
        applyCatalog(std::move(*result.catalog));
        if (arguments_.backgroundCheck && !taskActive_) finishBackgroundCheck(0);
        return;
    }
    if (result.kind == AsyncResult::Kind::PanelChoices) {
        if (result.itemIndex != selected_) { updateExportButtons(); return; }
        if (result.panelChoices.empty()) {
            MessageBoxW(window_, L"This archive object has no additional native plot panels.",
                L"Add Viewer panel", MB_OK | MB_ICONINFORMATION);
        } else {
            const auto selected = showPanelPickerDialog(window_, result.panelChoices, fonts_);
            if (selected && *selected < result.panelChoices.size())
                plots_.appendPanel(std::move(result.panelChoices[*selected].figure));
        }
        updateExportButtons(); updatePageSelector();
        return;
    }
    if (result.itemIndex != selected_) { updateExportButtons(); return; }
    if (!result.observables && !result.phaseContext) {
        mainSources_ = std::move(result.sources);
        observablesWorkspace_ = {};
    }
    auto previous = plots_.captureState();
    const bool kromUpdate = plots_.isKrom() && result.payload && result.payload->layout == "krom"
        && previous.payload->source == result.payload->source;
    plots_.setPayload(result.payload);
    if (result.phaseContext) {
        previous.payload = result.payload; plots_.restoreState(std::move(previous));
    } else if (kromUpdate) {
        auto next = plots_.captureState(); next.scrollY = previous.scrollY; next.page = std::min(previous.page,plots_.pageCount()-1);
        if (previous.axisLimits.count(0)) next.axisLimits[0] = previous.axisLimits.at(0);
        plots_.restoreState(std::move(next));
    }
    layout();
    updateExportButtons();
    updatePageSelector();
    if (arguments_.backgroundCheck) finishBackgroundCheck(0);
}

}

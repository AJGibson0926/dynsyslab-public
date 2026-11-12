#include "dynsyslab_viewer/Workbench.hpp"
#include "dynsyslab_viewer/ComparisonPayload.hpp"
#include "dynsyslab_viewer/PlotExport.hpp"
#include "dynsyslab_viewer/Process.hpp"
#include "ComparisonDialog.hpp"

#include <algorithm>
#include <filesystem>

namespace dsl {

void Workbench::copyPlots() {
    (void)copyPlotToClipboard(window_, plots_.handle());
}

void Workbench::savePlots() {
    std::wstring stem = L"dynsyslab-plot";
    if (selected_ >= 0) stem = wide(catalog_.items[static_cast<std::size_t>(selected_)].label);
    if (observablesActive_) stem += L"-observables";
    else if (!activeView_.empty()) stem += L"-" + wide(activeView_);
    (void)savePlotToDisk(window_, plots_.handle(), stem);
}

void Workbench::openPostprocessor() {
    if (selected_<0) return;
    wchar_t executableBuffer[32768]{};
    const DWORD count=GetModuleFileNameW(nullptr,executableBuffer,32768);
    const auto executable=std::filesystem::path(
        std::wstring(executableBuffer,count)).parent_path()/L"DynSysLab.exe";
    const auto source=backend_.resolveSource(
        catalog_.items[static_cast<std::size_t>(selected_)]);
    const std::vector<std::wstring> command{executable.wstring(),L"--project-root",
        arguments_.backend.projectRoot.wstring(),L"--archive-root",
        arguments_.backend.archiveRoot.wstring(),L"--cache",
        (arguments_.backend.cache.parent_path()/L"Compute").wstring(),L"--source",source.wstring(),
        L"--postprocessor"};
    const auto result=launchDetached(command,arguments_.backend.projectRoot.wstring());
    if(result.exitCode) plots_.setMessage(L"Unable to open Post-processor",result.error);
}

void Workbench::startAddPanel() {
    if (selected_ < 0 || taskActive_) return;
    const int requestedIndex = selected_;
    const CatalogItem item = catalog_.items[static_cast<std::size_t>(selected_)];
    taskActive_ = true; updateExportButtons();
    task_ = std::async(std::launch::async, [this, item, requestedIndex]() {
        AsyncResult result; result.kind = AsyncResult::Kind::PanelChoices;
        result.itemIndex = requestedIndex;
        try { result.panelChoices = backend_.loadPanelChoices(item); }
        catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
}

void Workbench::startComparison() {
    if (selected_ < 0 || activeView_.empty() || taskActive_ || observablesActive_) return;
    const auto selected = showComparisonDialog(
        window_, catalog_, selected_, activeView_, fonts_);
    if (!selected) return;
    const int firstIndex = selected_;
    const CatalogItem first = catalog_.items[static_cast<std::size_t>(firstIndex)];
    const CatalogItem second = catalog_.items[static_cast<std::size_t>(*selected)];
    const std::string view = activeView_;
    plots_.setMessage(L"Building synchronized comparison",
        wide(first.label) + L"  ↔  " + wide(second.label));
    updateExportButtons(); updatePageSelector();
    taskActive_ = true;
    task_ = std::async(std::launch::async, [this, first, second, firstIndex, view]() {
        AsyncResult result; result.itemIndex = firstIndex; result.view = view + " comparison"; result.sources = {first, second};
        try {
            const auto left = backend_.loadPayload(first, view);
            const auto right = backend_.loadPayload(second, view);
            result.payload = std::make_shared<Payload>(
                comparePayload(left, first.label, right, second.label));
        } catch (const std::exception& error) { result.error = wide(error.what()); }
        return result;
    });
}

}  // namespace dsl

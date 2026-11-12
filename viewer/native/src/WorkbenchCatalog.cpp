#include "dynsyslab_viewer/Workbench.hpp"
#include "dynsyslab_viewer/StartupDiagnostics.hpp"

namespace dsl {

void Workbench::applyCatalog(Catalog catalog) {
    // Tree deletion can synchronously paint sibling controls. Clear borrowed
    // item pointers before releasing the old catalog, and ignore its selection.
    updatingCatalog_ = true;
    inspector_.setItem(nullptr, nullptr);
    selected_ = -1;
    resetWorkspace();
    catalog_ = std::move(catalog);
    const auto title = catalog_.focusLabel.empty() ? L"DynSysLab Viewer"
        : L"DynSysLab Viewer — Mask: " + wide(catalog_.focusLabel);
    SetWindowTextW(window_, title.c_str());
    navigation_.setCatalog(&catalog_);
    updatingCatalog_ = false;
    writeStartupDiagnostics("Viewer", arguments_.backend.cache,
        arguments_.backend.projectRoot, arguments_.backend.archiveRoot,
        arguments_.backend.catalog, catalog_.items.size(), navigation_.nodeCount(),
        arguments_.backend.archiveFocus, catalog_.focusLabel, catalog_.hiddenItemCount);
    // Record the actual tree before plot preparation can fail or dispatch paints.
    initialSelection();
    if (catalog_.items.empty()) {
        plots_.setMessage(catalog_.focusLabel.empty() ? L"No archived objects"
            : L"No archived objects match the active Mask", wide(catalog_.focusLabel));
        updateViewSelector();
        updateExportButtons();
        updatePageSelector();
    }
}

}  // namespace dsl

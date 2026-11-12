#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/ArchiveFocus.hpp"
#include "dynsyslab_viewer/NativePayload.hpp"
#include "HierarchyPayload.hpp"
#include "CatalogPresentation.hpp"
#include "dynsyslab_viewer/KromWorkspace.hpp"
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <stdexcept>

namespace dsl {

Backend::Backend(BackendConfig config) : config_(std::move(config)) {
    std::filesystem::create_directories(config_.cache);
}

Catalog Backend::loadCatalog() const {
    const auto presented = [&](Catalog catalog) {
        applyArchiveFocus(catalog, config_.archiveFocus, config_.projectRoot);
        applyCatalogPresentation(catalog);
        return catalog;
    };
    auto indexed = Catalog::indexed(config_.projectRoot, config_.archiveRoot);
    if (!indexed.items.empty()) return presented(std::move(indexed));
    if (!config_.catalog.empty() && std::filesystem::is_regular_file(config_.catalog)) {
        Catalog catalog;
        try {
            catalog = Catalog::load(config_.catalog.wstring());
        }
        catch (const std::exception&) { catalog = Catalog::scan(config_.projectRoot, config_.archiveRoot); }
        return presented(std::move(catalog));
    }
    auto catalog = Catalog::scan(config_.projectRoot, config_.archiveRoot);
    return presented(std::move(catalog));
}

Catalog Backend::refreshCatalog() const {
    (void)dynsyslab::core::reconcileArchiveIndex(config_.archiveRoot);
    return loadCatalog();
}

std::filesystem::path Backend::resolveSource(const CatalogItem& item) const {
    std::filesystem::path source = wide(item.source);
    if (source.is_relative()) source = config_.projectRoot / source;
    return std::filesystem::weakly_canonical(source);
}

std::filesystem::path Backend::payloadPath(const CatalogItem& item,
                                            const std::string& view) const {
    std::string safe = "presentation-v9-" + item.id + "-" + view;
    for (char& value : safe) if (!std::isalnum(static_cast<unsigned char>(value))) value = '-';
    return config_.cache / wide(safe + ".json");
}

Payload Backend::loadPayload(const CatalogItem& item, const std::string& view) const {
    const auto finish = [&](Payload payload) { return applyPhaseSpace(std::move(payload), {item}); };
    const std::filesystem::path source = resolveSource(item);
    if (item.kind == "krom" && (view == "koopman_spectrum" || view == "krom_dashboard"))
        return buildKromWorkspace(source);
    const std::filesystem::path output = payloadPath(item, view);
    std::error_code ignored;
    if (std::filesystem::is_regular_file(output, ignored) &&
        std::filesystem::last_write_time(output, ignored) >=
        std::filesystem::last_write_time(source, ignored)) {
        try { return finish(Payload::load(output.wstring())); }
        catch (const std::exception&) { /* Regenerate an invalid cache entry. */ }
    }
    if (item.kind == "batch" || item.kind == "system_run" || item.kind == "cluster")
        return finish(buildHierarchyPayload(item, view, config_.projectRoot, config_.archiveRoot));
    if (source.extension() == L".json") {
        if (item.kind == "forcing_campaign")
            return finish(buildHierarchyPayload(item, view, config_.projectRoot, config_.archiveRoot));
        return finish(Payload::load(source.wstring()));
    }
    if (source.extension() == L".npz") return finish(buildNativePayload(source, view));
    throw std::runtime_error("This file type does not provide a native plot payload");
}

std::vector<PanelChoice> Backend::loadPanelChoices(const CatalogItem& item) const {
    std::vector<std::pair<std::string, Payload>> payloads;
    const auto source = resolveSource(item);
    if (source.extension() == L".npz") {
        const auto archive = dynsyslab::core::readNpz(source);
        dynsyslab::core::json::Value metadata;
        if (!archive.metadataJson().empty())
            metadata = dynsyslab::core::json::parse(archive.metadataJson());
        for (const auto& view : panelSourceViews(item)) {
            try { payloads.emplace_back(view,
                buildNativePayload(source, view, archive, metadata)); }
            catch (const std::exception&) { /* A single unavailable panel is non-fatal. */ }
        }
        return composePanelChoices(std::move(payloads));
    }
    for (const auto& view : panelSourceViews(item)) {
        try { payloads.emplace_back(view, loadPayload(item, view)); }
        catch (const std::exception&) { /* A single unavailable panel is non-fatal. */ }
    }
    return composePanelChoices(std::move(payloads));
}

}  // namespace dsl

#pragma once

#include "dynsyslab_viewer/Model.hpp"
#include "dynsyslab_viewer/PanelCatalog.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"

#include <filesystem>
#include <atomic>
#include <string>

namespace dsl {

struct BackendConfig {
    std::filesystem::path projectRoot;
    std::filesystem::path archiveRoot;
    std::filesystem::path catalog;
    std::filesystem::path cache;
    std::filesystem::path archiveFocus;
};

bool supportsObservableWorkspace(const CatalogItem& item);

class Backend {
public:
    explicit Backend(BackendConfig config);

    const BackendConfig& config() const { return config_; }
    Catalog loadCatalog() const;
    Catalog refreshCatalog() const;
    Payload loadPayload(const CatalogItem& item, const std::string& view) const;
    // Evaluates every system observable on the same selected Main trajectories.
    Payload loadObservables(const std::vector<CatalogItem>& items,
        std::optional<PhaseSpaceDomain> phaseDomain = std::nullopt) const;
    // Explicit compute or read-only cache lookup, for the selected trajectories' native systems.
    json::Value queryPhaseSpace(const CatalogItem& item, const std::string& mode = "lookup",
        std::optional<PhaseSpaceDomain> domain = std::nullopt, const json::Value& options = {},
        const std::atomic_bool* cancelled = nullptr) const;
    Payload applyPhaseSpace(Payload payload, const std::vector<CatalogItem>& items,
        const std::string& mode = "lookup", std::optional<PhaseSpaceDomain> domain = std::nullopt,
        const std::atomic_bool* cancelled = nullptr) const;
    std::vector<PanelChoice> loadPanelChoices(const CatalogItem& item) const;
    std::filesystem::path resolveSource(const CatalogItem& item) const;

private:
    BackendConfig config_;
    std::filesystem::path payloadPath(const CatalogItem& item,
                                      const std::string& view) const;
};

}  // namespace dsl

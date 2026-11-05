#pragma once
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/KromWorkspace.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"

namespace dsl {
// Background navigation uses the same full catalog, optional explicit Mask,
// and payload services as Workbench. It never controls a desktop application.
class ViewerSession {
public:
    explicit ViewerSession(BackendConfig config);
    const Catalog& catalog() const { return catalog_; }
    json::Value tree(const std::string& search = {},
                     const std::optional<std::string>& parent = std::nullopt) const;
    const CatalogItem& select(const std::string& id);
    const CatalogItem& selectSource(const std::filesystem::path& source);
    const CatalogItem& selected() const;
    Payload plot(const std::string& view = {}, const std::string& coordinates = {}) const;
    Payload workspace(const std::string& workspace = "main", const std::string& view = {},
        const std::string& coordinates = {}, const std::string& comparisonId = {},
        const KromWorkspaceOptions& krom = {}) const;
    json::Value phaseSpace(const std::string& mode = "lookup", const json::Value& options = {}) const;
    std::vector<PanelChoice> panels() const;
private:
    Backend backend_;
    Catalog catalog_;
    int selection_ = -1;
};

// A compact inspection includes labels, sample counts/ranges, reference points,
// and available coordinate transforms; includeSamples adds the full trace arrays.
json::Value inspectPlot(const Payload& payload, bool includeSamples = false);

struct PlotRenderOptions {
    int width = 1100, height = 760, page = 0;
    std::optional<std::size_t> figure;
    bool toolbar = true;
    int scrollY = 0;
    int benchmarkScrollSteps = 0;
};
// Uses PlotCanvas's exact drawing path on an invisible, non-activating surface.
// Returns the rendered page's axis bounds. No dialogs, clipboard, or screen capture.
json::Value renderPlotPng(const Payload& payload, const std::filesystem::path& output,
                          const PlotRenderOptions& options = {});
}

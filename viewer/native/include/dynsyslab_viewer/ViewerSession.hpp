#pragma once
#include "dynsyslab_viewer/Backend.hpp"
#include "dynsyslab_viewer/KromWorkspace.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"

namespace dsl {
// One headless view request. `view` is a route (as `panels` reports it) or a panel
// label; `panel` keeps only that panel of the route. Comparisons overlay further runs
// (catalog ids or archive paths), optionally relabelled and recolored in order
// (selected run first). The window is in the archive's stored time units.
struct ViewRequest {
    std::string workspace = "main", view, panel, coordinates;
    std::vector<std::string> compareIds;
    std::vector<std::filesystem::path> compareSources;
    std::vector<std::string> labels, colors;
    TimeWindow window;
    KromWorkspaceOptions krom;
};

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
    Payload workspace(const ViewRequest& request) const;
    // Routes this item offers (catalog views plus composed panel routes).
    std::vector<std::string> routes() const;
    json::Value phaseSpace(const std::string& mode = "lookup", const json::Value& options = {}) const;
    std::vector<PanelChoice> panels() const;
private:
    int indexOfSource(const std::filesystem::path& source) const;
    std::pair<std::string, std::string> resolveView(const CatalogItem& item, const std::string& view,
                                                    const std::string& panel) const;
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
    double fontScale = defaultPlotFontScale;  // 0.5..4: scales every plot font (and the font-derived margins)
    // Fixed axis ranges for every panel on the page. A missing range is fitted to the
    // data (the ordinate to the samples inside the abscissa range).
    std::optional<std::pair<double, double>> xRange, yRange;
};
// Uses PlotCanvas's exact drawing path on an invisible, non-activating surface.
// Returns the rendered page's axis bounds. No dialogs, clipboard, or screen capture.
json::Value renderPlotPng(const Payload& payload, const std::filesystem::path& output,
                          const PlotRenderOptions& options = {});
}

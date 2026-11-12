#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/ComparisonPayload.hpp"
#include <algorithm>
#include <cwctype>
#include <stdexcept>

namespace dsl {
namespace {
std::wstring folded(const std::string& value) {
    auto text = wide(value);
    std::transform(text.begin(), text.end(), text.begin(), [](wchar_t c) { return std::towlower(c); });
    return text;
}
}
ViewerSession::ViewerSession(BackendConfig config)
    : backend_(std::move(config)), catalog_(backend_.loadCatalog()) {}

json::Value ViewerSession::tree(const std::string& search,
                                const std::optional<std::string>& parent) const {
    if (parent && !parent->empty() && catalog_.indexOf(*parent) < 0)
        throw std::runtime_error("Parent is not in the visible data tree: " + *parent);
    Catalog filtered = catalog_;
    const auto query = folded(search);
    filtered.items.erase(std::remove_if(filtered.items.begin(), filtered.items.end(),
        [&](const CatalogItem& item) {
            if (parent && item.parentId != *parent) return true;
            auto text = item.id + " " + item.label + " " + item.system + " " +
                item.category + " " + item.role + " " + item.source;
            for (const auto& tag : item.tags) text += " " + tag;
            return !query.empty() && folded(text).find(query) == std::wstring::npos;
        }), filtered.items.end());
    auto result = json::parse(filtered.toJson());
    result.object()["schema"] = "dynsyslab.viewer.tree/v1";
    result.object()["focus_label"] = catalog_.focusLabel;
    result.object()["hidden_item_count"] = static_cast<double>(catalog_.hiddenItemCount);
    result.object()["visible_item_count"] = static_cast<double>(catalog_.items.size());
    result.object()["returned_item_count"] = static_cast<double>(filtered.items.size());
    return result;
}
const CatalogItem& ViewerSession::select(const std::string& id) {
    const int index = catalog_.indexOf(id);
    if (index < 0) throw std::runtime_error("Object is not in the visible data tree: " + id);
    if (catalog_.items[static_cast<std::size_t>(index)].kind == "presentation_group")
        throw std::runtime_error("Navigation groups do not contain scientific payloads: " + id);
    selection_ = index;
    return selected();
}
const CatalogItem& ViewerSession::selectSource(const std::filesystem::path& source) {
    const auto& config = backend_.config();
    const auto path = std::filesystem::weakly_canonical(
        source.is_relative() ? config.projectRoot / source : source);
    for (std::size_t i = 0; i < catalog_.items.size(); ++i) {
        if (folded(backend_.resolveSource(catalog_.items[i]).generic_u8string())
            == folded(path.generic_u8string())) {
            selection_ = static_cast<int>(i); return selected();
        }
    }
    throw std::runtime_error("Source is not in the visible data tree: " + path.u8string());
}
const CatalogItem& ViewerSession::selected() const {
    if (selection_ < 0) throw std::runtime_error("Select a data-tree item with --id or --source");
    return catalog_.items[static_cast<std::size_t>(selection_)];
}
Payload ViewerSession::plot(const std::string& view, const std::string& coordinates) const {
    const auto& item = selected();
    const auto route = view.empty() ? (item.views.empty() ? "data_overview" : item.views.front()) : view;
    const auto panelViews = panelSourceViews(item);
    if (!item.views.empty() && std::find(item.views.begin(), item.views.end(), route) == item.views.end()
        && std::find(panelViews.begin(),panelViews.end(),route) == panelViews.end())
        throw std::runtime_error("View is not available for this object: " + route);
    auto payload = backend_.loadPayload(item, route);
    if (!coordinates.empty()) {
        if (coordinates != "dimensional" && coordinates != "dimensionless")
            throw std::runtime_error("Coordinates must be dimensional or dimensionless");
        for (const auto& figure : payload.figures)
            if (!canDisplayCoordinates(figure, coordinates))
                throw std::runtime_error("This plot does not supply the requested coordinate scale");
        for (auto& figure : payload.figures) changeDisplayCoordinates(figure, coordinates);
    }
    return payload;
}
Payload ViewerSession::workspace(const std::string& workspace, const std::string& view,
    const std::string& coordinates, const std::string& comparisonId, const KromWorkspaceOptions& krom) const {
    if (workspace != "main" && workspace != "observables")
        throw std::runtime_error("Workspace must be main or observables");
    std::vector<CatalogItem> items{selected()};
    if (!comparisonId.empty()) {
        const auto index = catalog_.indexOf(comparisonId);
        if (index < 0) throw std::runtime_error("Comparison is not in the visible data tree: " + comparisonId);
        items.push_back(catalog_.items[static_cast<std::size_t>(index)]);
    }
    if (selected().kind == "krom" && !comparisonId.empty())
        throw std::runtime_error("Select KROM modes in the eigenvalue plot to compare eigenfunctions");
    auto payload = workspace == "observables" ? backend_.loadObservables(items) :
        selected().kind == "krom" && (view.empty() || view == "koopman_spectrum" || view == "krom_dashboard")
            ? buildKromWorkspace(backend_.resolveSource(selected()),krom) : plot(view);
    if (workspace == "main" && items.size() > 1) {
        const auto route = view.empty() ? items.front().views.front() : view;
        payload = comparePayload(payload, items.front().label,
            backend_.loadPayload(items.back(), route), items.back().label);
    }
    if (!coordinates.empty()) for (auto& figure : payload.figures) {
        if (!canDisplayCoordinates(figure, coordinates))
            throw std::runtime_error("This workspace does not supply the requested coordinate scale");
        changeDisplayCoordinates(figure, coordinates);
    }
    return payload;
}
json::Value ViewerSession::phaseSpace(const std::string& mode,const json::Value& options) const {
    return backend_.queryPhaseSpace(selected(),mode,std::nullopt,options);
}
std::vector<PanelChoice> ViewerSession::panels() const {
    return backend_.loadPanelChoices(selected());
}
}

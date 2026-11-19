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
int ViewerSession::indexOfSource(const std::filesystem::path& source) const {
    const auto& config = backend_.config();
    const auto path = std::filesystem::weakly_canonical(
        source.is_relative() ? config.projectRoot / source : source);
    for (std::size_t i = 0; i < catalog_.items.size(); ++i)
        if (folded(backend_.resolveSource(catalog_.items[i]).generic_u8string())
            == folded(path.generic_u8string())) return static_cast<int>(i);
    throw std::runtime_error("Source is not in the visible data tree: " + path.u8string());
}
const CatalogItem& ViewerSession::selectSource(const std::filesystem::path& source) {
    selection_ = indexOfSource(source);
    return selected();
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
    ViewRequest request;
    request.workspace = workspace; request.view = view; request.coordinates = coordinates;
    if (!comparisonId.empty()) request.compareIds.push_back(comparisonId);
    request.krom = krom;
    return this->workspace(request);
}

std::vector<std::string> ViewerSession::routes() const {
    const auto& item = selected();
    auto result = item.views;
    for (const auto& route : panelSourceViews(item))
        if (std::find(result.begin(), result.end(), route) == result.end()) result.push_back(route);
    return result;
}

namespace {
std::string panelLabel(const Figure& figure) { return figure.title.empty() ? figure.yLabel : figure.title; }

std::string listed(const std::vector<std::string>& values) {
    std::string result;
    for (const auto& value : values) result += (result.empty() ? "" : ", ") + value;
    return result;
}

void keepPanel(Payload& payload, const std::string& panel) {
    if (panel.empty()) return;
    std::vector<Figure> kept;
    std::vector<std::string> available;
    for (auto& figure : payload.figures) {
        available.push_back(panelLabel(figure));
        if (folded(panelLabel(figure)) == folded(panel)) kept.push_back(std::move(figure));
    }
    if (kept.empty())
        throw std::runtime_error("Panel \"" + panel + "\" is not in this view; its panels are: " + listed(available));
    for (std::size_t i = 0; i < kept.size(); ++i) { kept[i].dashboardPage = static_cast<int>(i / 2); kept[i].gridSpan = 2; }
    payload.figures = std::move(kept);
    payload.pageLabels.clear();
}
}  // namespace

std::pair<std::string, std::string> ViewerSession::resolveView(const CatalogItem& item,
    const std::string& view, const std::string& panel) const {
    const auto available = routes();
    if (view.empty() || std::find(available.begin(), available.end(), view) != available.end()
        || item.views.empty()) return {view, panel};
    // A catalog view label ("Phase space", "Frequency spectra") names its route.
    for (const auto& [route, label] : catalog_.viewLabels)
        if (folded(label) == folded(view)
            && std::find(available.begin(), available.end(), route) != available.end())
            return {route, panel};
    // A panel label as `panels` reports it ("r versus time", "Control history").
    std::vector<std::string> matches;
    std::vector<std::string> labels;
    for (const auto& choice : panels()) {
        labels.push_back(choice.label + " [" + choice.route + "]");
        if (folded(choice.label) == folded(view)
            && std::find(matches.begin(), matches.end(), choice.route) == matches.end())
            matches.push_back(choice.route);
    }
    if (matches.size() == 1) return {matches.front(), panel.empty() ? view : panel};
    if (matches.size() > 1)
        throw std::runtime_error("Panel \"" + view + "\" exists in several routes (" + listed(matches)
            + "); use --view <route> --panel \"" + view + "\"");
    throw std::runtime_error("View is not available for this object: " + view + ". Routes: " + listed(available)
        + ". Panels: " + listed(labels));
}

Payload ViewerSession::workspace(const ViewRequest& request) const {
    if (request.workspace != "main" && request.workspace != "observables")
        throw std::runtime_error("Workspace must be main or observables");
    std::vector<CatalogItem> items{selected()};
    for (const auto& id : request.compareIds) {
        const auto index = catalog_.indexOf(id);
        if (index < 0) throw std::runtime_error("Comparison is not in the visible data tree: " + id);
        items.push_back(catalog_.items[static_cast<std::size_t>(index)]);
    }
    for (const auto& source : request.compareSources)
        items.push_back(catalog_.items[static_cast<std::size_t>(indexOfSource(source))]);
    if (selected().kind == "krom" && items.size() > 1)
        throw std::runtime_error("Select KROM modes in the eigenvalue plot to compare eigenfunctions");
    if (!request.labels.empty() && request.labels.size() != items.size())
        throw std::runtime_error("Give one --label per run (" + std::to_string(items.size()) + "), selected run first");
    if (!request.colors.empty() && request.colors.size() != items.size())
        throw std::runtime_error("Give one --color per run (" + std::to_string(items.size()) + "), selected run first");
    const auto convert = [&](Payload& payload) {
        if (request.coordinates.empty()) return;
        if (request.coordinates != "dimensional" && request.coordinates != "dimensionless")
            throw std::runtime_error("Coordinates must be dimensional or dimensionless");
        for (auto& figure : payload.figures) {
            if (!canDisplayCoordinates(figure, request.coordinates))
                throw std::runtime_error("This plot does not supply the requested coordinate scale");
            changeDisplayCoordinates(figure, request.coordinates);
        }
    };
    if (request.workspace == "observables") {
        if (request.window.active()) throw std::runtime_error("Time windows apply to the main workspace");
        auto payload = backend_.loadObservables(items);
        convert(payload);
        return payload;
    }
    if (selected().kind == "krom" && !request.window.active() && request.panel.empty()
        && (request.view.empty() || request.view == "koopman_spectrum" || request.view == "krom_dashboard")) {
        auto payload = buildKromWorkspace(backend_.resolveSource(selected()), request.krom);
        convert(payload);
        return payload;
    }
    auto [route, panel] = resolveView(selected(), request.view, request.panel);
    if (route.empty()) route = selected().views.empty() ? "data_overview" : selected().views.front();
    std::vector<ComparisonEntry> entries;
    for (std::size_t k = 0; k < items.size(); ++k) {
        const auto& item = items[k];
        const auto offered = panelSourceViews(item);
        if (!item.views.empty() && std::find(item.views.begin(), item.views.end(), route) == item.views.end()
            && std::find(offered.begin(), offered.end(), route) == offered.end())
            throw std::runtime_error("View is not available for " + item.label + ": " + route);
        auto payload = backend_.loadPayload(item, route, request.window);
        keepPanel(payload, panel);
        convert(payload);
        entries.push_back({request.labels.empty() ? item.label : request.labels[k], std::move(payload),
                           request.colors.empty() ? std::string{} : request.colors[k]});
    }
    if (entries.size() == 1) return std::move(entries.front().payload);
    return comparePayloads(entries);
}

json::Value ViewerSession::phaseSpace(const std::string& mode,const json::Value& options) const {
    return backend_.queryPhaseSpace(selected(),mode,std::nullopt,options);
}
std::vector<PanelChoice> ViewerSession::panels() const {
    return backend_.loadPanelChoices(selected());
}
}

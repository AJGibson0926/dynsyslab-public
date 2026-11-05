#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"

#include <algorithm>
#include <map>

namespace dsl {
namespace {

std::vector<std::size_t> visibleIndices(const Payload& payload, int page) {
    const bool dashboard = std::any_of(payload.figures.begin(), payload.figures.end(),
        [](const Figure& figure) { return figure.dashboardPage >= 0; });
    std::vector<std::size_t> result;
    if (dashboard) {
        for (std::size_t index = 0; index < payload.figures.size(); ++index)
            if (payload.figures[index].dashboardPage == page) result.push_back(index);
        return result;
    }
    constexpr std::size_t perPage = 4;
    const std::size_t start = static_cast<std::size_t>(page) * perPage;
    for (std::size_t index = start;
         index < std::min(payload.figures.size(), start + perPage); ++index)
        result.push_back(index);
    return result;
}

template <typename Value>
void eraseIndexed(std::map<std::size_t, Value>& values, std::size_t removed) {
    std::map<std::size_t, Value> adjusted;
    for (auto& entry : values) {
        if (entry.first == removed) continue;
        adjusted[entry.first > removed ? entry.first - 1 : entry.first]
            = std::move(entry.second);
    }
    values = std::move(adjusted);
}

}  // namespace

void PlotCanvas::setPage(int page) {
    const int selected = std::clamp(page, 0, pageCount() - 1);
    if (selected == page_) return;
    page_ = selected; hoverProbe_.reset(); InvalidateRect(window_, nullptr, FALSE);
    PostMessageW(GetParent(window_), kPlotPageChangedMessage,
                 static_cast<WPARAM>(page_), 0);
}

void PlotCanvas::nextPage() { setPage(page_ + 1); }
void PlotCanvas::previousPage() { setPage(page_ - 1); }

void PlotCanvas::notifyCompositionChanged() {
    PostMessageW(GetParent(window_), kPlotCompositionChangedMessage, 0, 0);
}

void PlotCanvas::setBackgroundDragMode(PlotBackgroundDragMode mode) {
    if (backgroundDragMode_ == mode) return;
    pointerUp(); backgroundDragMode_ = mode;
}

void PlotCanvas::ensureComposition() {
    if (composed_ || !payload_ || verticallyScrollable()) return;
    const auto indices = visibleIndices(*payload_, page_);
    std::vector<Figure> visible;
    std::map<std::size_t, CurveProbe> visibleProbes;
    std::map<std::size_t, AxisLimits> visibleLimits;
    std::map<std::size_t, Plot3DCameraState> visibleCameras;
    for (std::size_t visualIndex = 0; visualIndex < indices.size(); ++visualIndex) {
        const auto index = indices[visualIndex];
        visible.push_back(std::move(payload_->figures[index]));
        const auto pinned = pinnedProbes_.find(index);
        if (pinned != pinnedProbes_.end()) visibleProbes[visualIndex] = pinned->second;
        const auto limits = axisLimits_.find(index);
        if (limits != axisLimits_.end()) visibleLimits[visualIndex] = limits->second;
        const auto camera = cameras3D_.find(index);
        if (camera != cameras3D_.end()) visibleCameras[visualIndex] = camera->second;
    }
    rasterCache_.reset(); payload_->figures = std::move(visible);
    for (auto& figure : payload_->figures) {
        figure.dashboardPage = -1;
        figure.gridSpan = std::max(1, figure.gridSpan);
        figure.rowSpan = std::max(1, figure.rowSpan);
    }
    payload_->pageLabels.clear(); payload_->columns = 0;
    axisLimits_ = std::move(visibleLimits); pinnedProbes_ = std::move(visibleProbes);
    cameras3D_ = std::move(visibleCameras);
    hoverProbe_.reset(); page_ = 0; composed_ = true;
}

void PlotCanvas::appendPanel(Figure figure) {
    const auto current = displayCoordinates();
    if (!current.empty()) changeDisplayCoordinates(figure, current);
    if (!payload_) payload_ = std::make_shared<Payload>();
    ensureComposition(); figure.dashboardPage = -1;
    figure.gridSpan = std::max(1, figure.gridSpan);
    figure.rowSpan = std::max(1, figure.rowSpan);
    rasterCache_.reset(); payload_->figures.push_back(std::move(figure));
    hoverProbe_.reset(); InvalidateRect(window_, nullptr, FALSE);
    notifyCompositionChanged();
}

void PlotCanvas::removePanel(std::size_t visualIndex) {
    ensureComposition();
    if (!payload_ || visualIndex >= payload_->figures.size()) return;
    rasterCache_.reset();
    payload_->figures.erase(payload_->figures.begin()
        + static_cast<std::ptrdiff_t>(visualIndex));
    eraseIndexed(axisLimits_, visualIndex); eraseIndexed(pinnedProbes_, visualIndex);
    eraseIndexed(cameras3D_, visualIndex);
    if (isKrom()) page_ = std::min(page_,pageCount()-1);
    hoverProbe_.reset();
    if (payload_->figures.empty()) {
        messageTitle_ = L"Plot workspace is empty";
        messageDetail_ = L"Use + Panel to add a scientific view.";
    }
    InvalidateRect(window_, nullptr, FALSE); notifyCompositionChanged();
}

std::optional<std::size_t> PlotCanvas::hitPanel(POINT point) const {
    for (const auto& hit : hitRegions_)
        if (PtInRect(&hit.card, point)) return hit.visualIndex;
    return std::nullopt;
}

}  // namespace dsl

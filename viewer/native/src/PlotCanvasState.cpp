#include "dynsyslab_viewer/PlotCanvas.hpp"
namespace dsl {
std::optional<PhaseSpaceDomain> PlotCanvas::phaseSpaceDomain() const {
    if (!payload_) return std::nullopt;
    for (std::size_t i = 0; i < payload_->figures.size(); ++i) {
        const auto& figure = payload_->figures[i];
        if (figure.kind != "phase") continue;
        const auto manual = axisLimits_.find(i);
        return PhaseSpaceDomain{manual == axisLimits_.end() ? automaticAxisLimits(figure) : manual->second,
            figure.coordinates};
    }
    return std::nullopt;
}

PlotCanvasState PlotCanvas::captureState() const {
    return {payload_, axisLimits_, pinnedProbes_, cameras3D_, page_, composed_,
            scrollBar_.position()};
}
void PlotCanvas::restoreState(PlotCanvasState state) {
    if (GetCapture() == window_) ReleaseCapture();
    rasterCache_.reset();
    payload_ = std::move(state.payload);
    axisLimits_ = std::move(state.axisLimits);
    pinnedProbes_ = std::move(state.pinnedProbes);
    cameras3D_ = std::move(state.cameras3D);
    page_ = state.page; composed_ = state.composed;
    if (!preferredValueScale_.empty()) setValueScaling(preferredValueScale_);
    hitRegions_.clear(); hoverProbe_.reset(); hoveredPanel_.reset();
    active3DFigure_.reset(); hovered3DFigure_.reset();
    hovered3DControl_ = Plot3DControl::none; dragMode_ = DragMode::none;
    if (!preferredCoordinates_.empty() && canDisplayCoordinates(preferredCoordinates_))
        setDisplayCoordinates(preferredCoordinates_);
    setScrollPosition(state.scrollY);
    InvalidateRect(window_, nullptr, FALSE);
}
}

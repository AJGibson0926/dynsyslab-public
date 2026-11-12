#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PanelWorkspace.hpp"
#include "MathTextDrawing.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
std::string PlotCanvas::valueScaling() const {
    return payload_ && !payload_->figures.empty() ? payload_->figures.front().valueScale
        : preferredValueScale_.empty() ? "auto" : preferredValueScale_;
}
void PlotCanvas::setValueScaling(const std::string& requested) {
    if (payload_ && !payload_->figures.empty() && std::all_of(payload_->figures.begin(),payload_->figures.end(),
        [&](const auto& f) { return f.valueScale == requested; })) { preferredValueScale_ = requested; return; }
    rasterCache_.reset();
    if (payload_) {
        payload_ = std::make_shared<Payload>(*payload_);
        setValueScale(*payload_,requested);
        for (std::size_t i = 0; i < payload_->figures.size(); ++i)
            if (supportsValueScale(payload_->figures[i])) axisLimits_.erase(i);
    } else { Payload empty; setValueScale(empty,requested); }
    preferredValueScale_ = requested;
    hoverProbe_.reset(); hitRegions_.clear(); InvalidateRect(window_,nullptr,FALSE);
}
bool PlotCanvas::verticallyScrollable() const {
    return payload_ && (isKrom() || payload_->view == "observable_workspace");
}
void PlotCanvas::updateScrollLayout() {
    GetClientRect(window_, &scrollViewport_);
    if (isKrom()) scrollViewport_.right = static_cast<LONG>(scrollViewport_.right * .66);
    if (!verticallyScrollable()) { scrollBar_.configure(0,scrollViewport_.bottom); scrollMaximum_ = 0; return; }
    const int width = scrollViewport_.right - OrnamentalScrollBar::preferredWidth(window_);
    const auto scale = scientificTextScale(fonts_->body);
    layoutColumns_ = width >= 900 * scale ? 2 : 1;
    std::vector<PanelSpan> spans;
    for (const auto& figure : payload_->figures)
        if (figure.kind != "krom_spectrum") spans.push_back({figure.gridSpan,figure.rowSpan});
    const auto placements = packPanelGrid(spans,layoutColumns_);
    layoutRows_ = packedPanelRows(placements);
    scrollRowHeight_ = std::max(340.f * scale, static_cast<float>(scrollViewport_.bottom-16) / std::max(1,std::min(2,layoutRows_)) - 4);
    const int content = static_cast<int>(16 + layoutRows_ * (scrollRowHeight_+8) - 8);
    scrollMaximum_ = std::max(0,content-static_cast<int>(scrollViewport_.bottom));
    scrollBar_.configure(content,scrollViewport_.bottom);
}
void PlotCanvas::setScrollPosition(int pixels) {
    stopScroll();
    updateScrollLayout();
    if (scrollBar_.scrollBy(std::clamp(pixels,0,scrollMaximum_) - scrollBar_.position())) {
        hoverProbe_.reset(); hitRegions_.clear(); InvalidateRect(window_,nullptr,FALSE);
    }
}
std::vector<std::size_t> PlotCanvas::visibleFigureIndices() const {
    std::vector<std::size_t> result;
    for (const auto& hit : hitRegions_) result.push_back(hit.figureIndex);
    return result;
}
void PlotCanvas::stopScroll() {
    if (window_) KillTimer(window_,scrollTimerId_);
    scrollTick_ = 0; scrollTarget_ = scrollBar_.position();
}
void PlotCanvas::wheelScroll(int delta) {
    KillTimer(window_,cacheTimerId_);
    updateScrollLayout();
    UINT lines = 3; SystemParametersInfoW(SPI_GETWHEELSCROLLLINES,0,&lines,0);
    const double distance = lines == WHEEL_PAGESCROLL ? std::max(1L,scrollViewport_.bottom-40)
        : 40.0*lines;
    if (!scrollTick_) scrollTarget_ = scrollBar_.position();
    scrollTarget_ = std::clamp(scrollTarget_-static_cast<int>(std::lround(delta*distance/WHEEL_DELTA)),0,scrollMaximum_);
    if (scrollTarget_ == scrollBar_.position()) { stopScroll(); return; }
    if (!scrollTick_) scrollTick_ = GetTickCount64()-16;
    advanceScroll();
    if (scrollTick_) SetTimer(window_,scrollTimerId_,16,nullptr);
}
void PlotCanvas::advanceScroll() {
    if (!scrollTick_) return;
    scrollTarget_ = std::clamp(scrollTarget_,0,scrollMaximum_);
    const auto now = GetTickCount64();
    const double elapsed = static_cast<double>(now-scrollTick_); scrollTick_ = now;
    const int remaining = scrollTarget_-scrollBar_.position();
    const int step = std::abs(remaining) <= 1 ? remaining : static_cast<int>(std::copysign(
        std::max(1.0,std::round(std::abs(remaining)*(1-std::exp(-elapsed/32.0)))),remaining));
    if (scrollBar_.scrollBy(step)) {
        hoverProbe_.reset(); hoveredPanel_.reset(); hitRegions_.clear();
        InvalidateRect(window_,&scrollViewport_,FALSE);
    }
    if (scrollBar_.position() == scrollTarget_) stopScroll();
}

}

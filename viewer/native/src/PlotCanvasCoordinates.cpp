#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include <algorithm>
namespace dsl {
bool PlotCanvas::canDisplayCoordinates(const std::string& coordinates) const {
    return hasPayload()&&std::all_of(payload_->figures.begin(),payload_->figures.end(),
        [&](const auto& figure){return dsl::canDisplayCoordinates(figure,coordinates);});
}
std::string PlotCanvas::displayCoordinates() const {
    if(!hasPayload()) return {};
    const auto& first=payload_->figures.front().coordinates;
    return std::all_of(payload_->figures.begin(),payload_->figures.end(),
        [&](const auto& figure){return figure.coordinates==first;})?first:std::string{};
}
void PlotCanvas::setDisplayCoordinates(const std::string& coordinates) {
    if(!canDisplayCoordinates(coordinates)) return;
    preferredCoordinates_=coordinates;
    if (displayCoordinates() == coordinates) return;
    rasterCache_.reset();
    pointerUp();
    payload_=std::make_shared<Payload>(*payload_);
    for(std::size_t i=0;i<payload_->figures.size();++i) {
        auto& figure=payload_->figures[i];
        if(figure.coordinates==coordinates) continue;
        const auto scale=*figure.alternateCoordinates;
        const auto found=axisLimits_.find(i);
        if(found!=axisLimits_.end()) {
            found->second.xMin*=scale.xScale;found->second.xMax*=scale.xScale;
            found->second.yMin*=scale.yScale;found->second.yMax*=scale.yScale;
        }
        changeDisplayCoordinates(figure,coordinates);
    }
    hoverProbe_.reset();pinnedProbes_.clear();hitRegions_.clear();
    InvalidateRect(window_,nullptr,FALSE);
}
}

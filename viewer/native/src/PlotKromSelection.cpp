#include "dynsyslab_viewer/PlotCanvas.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
void PlotCanvas::setComplexDisplay(const std::string& display) {
    if (!isKrom() || payload_->complexDisplay == display) return;
    if (display != "real-imaginary" && display != "magnitude-phase" && display != "domain-coloring") return;
    payload_->complexDisplay = display;
    PostMessageW(GetParent(window_),kPlotKromChangedMessage,0,0);
}
void PlotCanvas::selectKromModes(std::vector<std::size_t> modes) {
    if (!isKrom()) return;
    std::sort(modes.begin(),modes.end()); modes.erase(std::unique(modes.begin(),modes.end()),modes.end());
    if (modes == payload_->selectedModes) return;
    rasterCache_.reset();
    payload_->selectedModes = std::move(modes);
    for (auto& f : payload_->figures) if (f.kind == "krom_spectrum")
        for (auto& t : f.traces) if (t.selectionId >= 0)
            t.selected = std::binary_search(payload_->selectedModes.begin(),payload_->selectedModes.end(),t.selectionId);
    InvalidateRect(window_,nullptr,FALSE);
    PostMessageW(GetParent(window_),kPlotKromChangedMessage,0,0);
}
void PlotCanvas::selectAllKromModes(bool all) {
    std::vector<std::size_t> modes;
    if (all && isKrom()) for (const auto& f : payload_->figures) if (f.kind == "krom_spectrum")
        for (const auto& t : f.traces) if (t.selectionId >= 0) modes.push_back(t.selectionId);
    selectKromModes(std::move(modes));
}
bool PlotCanvas::zoomKromAt(POINT point, double factor) {
    if (!isKrom()) return false;
    for (const auto& hit : hitRegions_) {
        if (payload_->figures[hit.figureIndex].kind != "krom_spectrum" || !PtInRect(&hit.card,point)) continue;
        if (factor == 0) axisLimits_.erase(hit.figureIndex);
        else {
            const auto& t = hit.transform;
            const double x = t.xMin+(point.x-t.left)*(t.xMax-t.xMin)/t.width;
            const double y = t.yMax-(point.y-t.top)*(t.yMax-t.yMin)/t.height;
            if ((t.xMax-t.xMin)*factor < 1e-12 || (t.xMax-t.xMin)*factor > 1e12) return true;
            axisLimits_[hit.figureIndex] = {x+(t.xMin-x)*factor,x+(t.xMax-x)*factor,
                y+(t.yMin-y)*factor,y+(t.yMax-y)*factor};
        }
        InvalidateRect(window_,nullptr,FALSE); return true;
    }
    return false;
}
bool PlotCanvas::selectKromAt(POINT point, bool control) {
    if (!isKrom()) return false;
    for (const auto& hit : hitRegions_) {
        const auto& f = payload_->figures[hit.figureIndex];
        if (f.kind != "krom_spectrum" || !PtInRect(&hit.card,point)) continue;
        const auto& t = hit.transform;
        int selected = -1; double distance = 12;
        for (const auto& trace : f.traces) if (trace.selectionId >= 0 && !trace.x.empty() && !trace.y.empty()) {
            const double px = t.left+(trace.x[0]-t.xMin)*t.width/(t.xMax-t.xMin);
            const double py = t.top+(t.yMax-trace.y[0])*t.height/(t.yMax-t.yMin);
            const double d = std::hypot(px-point.x,py-point.y);
            if (d < distance) { distance = d; selected = trace.selectionId; }
        }
        if (selected >= 0) {
            auto modes = control ? payload_->selectedModes : std::vector<std::size_t>{};
            const auto found = std::find(modes.begin(),modes.end(),selected);
            if (found == modes.end()) modes.push_back(selected); else modes.erase(found);
            selectKromModes(std::move(modes));
        }
        return true;
    }
    return false;
}
}

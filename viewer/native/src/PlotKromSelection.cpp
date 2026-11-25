#include "dynsyslab_viewer/PlotCanvas.hpp"
#include "ComplexHeatmap.hpp"
#include <algorithm>
#include <cmath>
namespace dsl {
void PlotCanvas::setComplexDisplay(const std::string& display) {
    if (!payload_ || payload_->complexDisplay.empty() || payload_->complexDisplay == display) return;
    if (!validComplexDisplay(display)) return;
    payload_->complexDisplay = display;
    if (isKrom()) { PostMessageW(GetParent(window_),kPlotKromChangedMessage,0,0); return; }
    // Other complex heat maps keep their values: re-encode in place.
    for (auto& figure : payload_->figures) encodeComplexHeatmap(figure, display);
    rasterCache_.reset();
    InvalidateRect(window_,nullptr,FALSE);
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

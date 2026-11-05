#include "dynsyslab_viewer/PlotCanvas.hpp"

namespace dsl {
namespace {
constexpr unsigned modifierControl = 1u;
constexpr unsigned modifierAlt = 2u;
constexpr unsigned modifierShift = 4u;
}

bool PlotCanvas::zoom3DAt(POINT point, double factor) {
    if (!payload_) return false;
    for (const auto& hit : hitRegions_) {
        if (!PtInRect(&hit.plot, point) || hit.figureIndex >= payload_->figures.size()
            || payload_->figures[hit.figureIndex].kind != "three_d") continue;
        Plot3DCamera camera(cameras3D_[hit.figureIndex]);
        camera.zoomBy(factor); cameras3D_[hit.figureIndex] = camera.state();
        active3DFigure_ = hit.figureIndex; rasterCache_.reset();
        InvalidateRect(window_, nullptr, FALSE); return true;
    }
    return false;
}

bool PlotCanvas::key3D(WPARAM key, unsigned modifiers) {
    if (!payload_) return false;
    std::optional<std::size_t> index = active3DFigure_;
    if (!index || *index >= payload_->figures.size()
        || payload_->figures[*index].kind != "three_d") {
        for (const auto& hit : hitRegions_)
            if (hit.figureIndex < payload_->figures.size()
                && payload_->figures[hit.figureIndex].kind == "three_d") {
                index = hit.figureIndex; break;
            }
    }
    if (!index) return false;
    Plot3DCamera camera(cameras3D_[*index]);
    bool handled = true;
    const bool control = (modifiers & modifierControl) != 0;
    const bool alt = (modifiers & modifierAlt) != 0;
    const bool fine = (modifiers & modifierShift) != 0;
    const double step = fine ? 7.0 : 22.0;
    if (key == VK_HOME || key == L'0' || key == VK_NUMPAD0) camera.reset();
    else if (key == L'P') camera.toggleProjection();
    else if (key == VK_ADD || key == VK_OEM_PLUS) camera.zoomBy(fine ? 1.06 : 1.18);
    else if (key == VK_SUBTRACT || key == VK_OEM_MINUS) camera.zoomBy(fine ? 1.0 / 1.06 : 1.0 / 1.18);
    else if (alt && (key == VK_UP || key == VK_RIGHT)) camera.zoomBy(fine ? 1.06 : 1.18);
    else if (alt && (key == VK_DOWN || key == VK_LEFT)) camera.zoomBy(fine ? 1.0 / 1.06 : 1.0 / 1.18);
    else if (control && key == VK_LEFT) camera.apply(Plot3DControl::pan_left);
    else if (control && key == VK_RIGHT) camera.apply(Plot3DControl::pan_right);
    else if (control && key == VK_UP) camera.apply(Plot3DControl::pan_up);
    else if (control && key == VK_DOWN) camera.apply(Plot3DControl::pan_down);
    else if (key == VK_LEFT) camera.orbit(-step, 0.0);
    else if (key == VK_RIGHT) camera.orbit(step, 0.0);
    else if (key == VK_UP) camera.orbit(0.0, -step);
    else if (key == VK_DOWN) camera.orbit(0.0, step);
    else handled = false;
    if (!handled) return false;
    cameras3D_[*index] = camera.state(); active3DFigure_ = index;
    rasterCache_.reset(); InvalidateRect(window_, nullptr, FALSE); return true;
}

}  // namespace dsl

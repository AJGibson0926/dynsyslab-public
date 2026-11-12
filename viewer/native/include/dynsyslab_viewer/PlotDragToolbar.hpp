#pragma once

#include "dynsyslab_viewer/PlotCanvas.hpp"

namespace dsl {

// Session-local, mutually exclusive plot background drag tools.
class PlotDragToolbar {
public:
    bool create(HWND parent, PlotCanvas& canvas, const theme::Fonts& fonts);
    HWND handle() const { return window_; }
    void setMode(PlotBackgroundDragMode mode);
    void updateCoordinates();
    int rowCount() const { return canvas_ && canvas_->isKrom() ? 3 : 2; }

private:
    HWND window_ = nullptr;
    HWND arrange_ = nullptr;
    HWND pan_ = nullptr;
    HWND realImaginary_ = nullptr, magnitudePhase_ = nullptr, domainColoring_ = nullptr;
    HWND scales_[4]{};
    HWND allModes_ = nullptr, noModes_ = nullptr;
    HWND dimensional_ = nullptr;
    HWND dimensionless_ = nullptr;
    PlotCanvas* canvas_ = nullptr;
    const theme::Fonts* fonts_ = nullptr;

    void layout();
    void paint(HDC target);
    LRESULT message(UINT message, WPARAM wparam, LPARAM lparam);
    static LRESULT CALLBACK procedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam);
};

}  // namespace dsl

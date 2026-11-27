#pragma once
// One vertical colour bar to the right of the plot for every colour-encoded figure: time-coloured
// lines, scalar heat maps, and domain colouring (hue = phase in [-pi, pi]).
#include "dynsyslab_viewer/PlotInteractionModel.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <gdiplus.h>

#include <optional>
#include <string>

namespace dsl {

struct PlotColorbar {
    enum class Kind { time, scalar, phase } kind = Kind::time;
    double low = 0.0, high = 1.0;  // values at the bottom and top of the bar
    std::string label;             // what the colour encodes (math text)
};

// The figure's colour bar, or none when no trace is colour-encoded.
[[nodiscard]] std::optional<PlotColorbar> plotColorbar(const Figure&);
// Width (device-independent pixels before text scaling) to reserve right of the plot.
[[nodiscard]] float plotColorbarReserve(const Figure&);
void drawPlotColorbar(Gdiplus::Graphics&, const Figure&, const Gdiplus::RectF& plot, const theme::Fonts&);

}  // namespace dsl

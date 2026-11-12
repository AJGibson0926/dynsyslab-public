#pragma once
#include "Model.hpp"

namespace dsl {
struct TimeSeriesCurve {
    std::string key, label, color, lineStyle;
    double lineWidthPoints=0.;
};

struct TimeSeriesOverlay {
    std::vector<TimeSeriesCurve> curves;
    std::vector<std::string> componentLabels;
    std::string timeKey = "time", timeLabel = "Time", title;
    std::optional<double> start, stop;
};

// Read aligned N or N-by-C arrays. A zero bin count retains every source vertex.
Payload timeSeriesOverlay(const std::filesystem::path& source, const TimeSeriesOverlay& request,
                          std::size_t horizontalBins = 0);
}

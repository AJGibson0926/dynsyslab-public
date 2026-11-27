#pragma once
#include "Model.hpp"

namespace dsl {
struct TimeSeriesCurve {
    std::string key, label, color, lineStyle;
    std::vector<std::string> componentColors;  // optional: one colour per component (panel)
    double lineWidthPoints=0.;
    bool rightAxis=false;  // fitted to its own ordinate on the right (TimeSeriesOverlay::stackAxes)
};

struct TimeSeriesOverlay {
    std::vector<TimeSeriesCurve> curves;
    std::vector<std::string> componentLabels;
    std::string timeKey = "time", timeLabel = "Time", title;
    std::optional<double> start, stop;
    // 0 or 1: one full-width panel per row, two per page (default). N >= 2: tile the
    // components N per row, N x N panels per page (e.g. 2 -> a 2x2 tile).
    std::size_t columns = 0;
    // Optional ordinate label shared by every panel; the component label then titles the panel.
    std::string valueLabel;
    // Phase overlays: a curve key "xkey,ykey" plots y against x (aligned N or N-by-C arrays,
    // restricted to the time window); such panels become phase panels with this abscissa label.
    std::string xLabel;
    // Spectra: each (non-phase) curve becomes its one-sided amplitude spectrum over the
    // time window (mean removed, Hann window, the archive Fourier route's convention);
    // frequency = cycles per clock unit x frequencyScale, without the zero bin.
    bool spectrum = false;
    double frequencyScale = 1.0;
    // Spectra: only the bins in this (scaled) frequency range are kept, before display binning, so
    // a zoom on a low band keeps its full resolution (the CLI's --x-range with --spectrum).
    std::optional<std::pair<double, double>> frequencyRange;
    // Curves with lineStyle "stem" are drawn as vertical stems from a baseline to each
    // sample with a marker at the tip (every sample, no binning). The baseline is
    // stemBaseline when given; otherwise 0, or for logBaseline one decade below the
    // smallest positive sample of any stem curve, shared by the figure, so stems stay
    // visible and comparable on a log axis.
    std::optional<double> stemBaseline;
    bool logBaseline = false;
    bool integerTicks = false;  // index abscissa: ticks at integers (Figure::integerXTicks)
    // Right ordinate showing value = (plotted - offset) / scale on every panel. With curves marked
    // rightAxis the scale and offset are fitted instead (only the label is used): the right curves
    // span the left curves' range, or with stackAxes a band of the same height below them, so
    // each side keeps its own scale and the two sets do not overlap.
    std::optional<Figure::RightAxis> rightAxis;
    bool stackAxes = false;
};

// Read aligned N or N-by-C arrays. A zero bin count retains every source vertex.
Payload timeSeriesOverlay(const std::filesystem::path& source, const TimeSeriesOverlay& request,
                          std::size_t horizontalBins = 0);
}

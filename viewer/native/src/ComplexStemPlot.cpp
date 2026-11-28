#include "ComplexStemPlot.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dsl {

Figure complexPlaneStemFigure(const ComplexStemPlot& plot) {
    if (plot.bases.size() != plot.heights.size())
        throw std::invalid_argument("complex stem plot needs one height per base point");
    constexpr double pi = 3.14159265358979323846;
    Figure figure{"three_d", plot.title, "Re(" + plot.baseSymbol + ")",
                  "Im(" + plot.baseSymbol + ")", false};
    figure.zLabel = plot.logHeights ? "log10 " + plot.heightLabel : plot.heightLabel;
    // Azimuth -30, elevation 28: both base axes are clearly visible and the low elevation
    // keeps stem heights readable (at -90 deg azimuth an orthographic Im axis would
    // coincide with the height axis on screen; at 45 deg elevation heights shrink).
    figure.cameraAzimuthDegrees = -30.0;
    figure.cameraElevationDegrees = 28.0;
    figure.boxAxes3D = true;               // independent axis scaling, no IC marker
    double floor = 0.0;
    if (plot.logHeights) {
        double smallest = std::numeric_limits<double>::infinity();
        for (const double h : plot.heights) if (h > 0 && std::isfinite(h)) smallest = std::min(smallest, h);
        for (const auto* set : {&plot.overlays, &plot.curves})
            for (const auto& overlay : *set)
                for (const double h : overlay.heights) if (h > 0 && std::isfinite(h)) smallest = std::min(smallest, h);
        floor = std::isfinite(smallest) ? std::floor(std::log10(smallest)) : 0.0;
    }
    const double gapValue = std::numeric_limits<double>::quiet_NaN();
    const auto inside = [&](double re, double im) {
        if (!plot.referenceWindow) return true;
        const auto& w = *plot.referenceWindow;
        return re >= w[0] && re <= w[1] && im >= w[2] && im <= w[3];
    };
    if (plot.unitCircle) {
        Trace circle; circle.label = "unit circle"; circle.color = "#94a3b8";
        circle.lineStyle = "dash";
        // Fine sampling so a zoomed window still shows a smooth arc.
        const int points = plot.referenceWindow ? 20000 : 360;
        bool drawn = false;
        for (int point = 0; point <= points; ++point) {
            const double angle = 2.0 * pi * point / points;
            const double re = std::cos(angle), im = std::sin(angle);
            const bool in = inside(re, im);
            if (!in && (circle.x.empty() || std::isnan(circle.x.back()))) continue;  // one gap per break
            circle.x.push_back(in ? re : gapValue); circle.y.push_back(in ? im : gapValue);
            circle.z.push_back(floor); drawn = drawn || in;
        }
        if (drawn) figure.traces.push_back(std::move(circle));
    }
    if (plot.imaginaryAxis) {
        double low = 0.0, high = 0.0;
        const auto extend = [&](const std::vector<std::complex<double>>& bases) {
            for (const auto& z : bases) if (std::isfinite(z.imag())) { low = std::min(low, z.imag()); high = std::max(high, z.imag()); }
        };
        extend(plot.bases); for (const auto& overlay : plot.overlays) extend(overlay.bases);
        if (plot.referenceWindow) {
            low = std::max(low, (*plot.referenceWindow)[2]); high = std::min(high, (*plot.referenceWindow)[3]);
        }
        if (high > low && inside(0.0, 0.5 * (low + high))) {
            Trace axis; axis.label = "Re = 0"; axis.color = "#94a3b8"; axis.lineStyle = "dash";
            axis.x = {0.0, 0.0}; axis.y = {low, high}; axis.z = {floor, floor};
            figure.traces.push_back(std::move(axis));
        }
    }
    // One trace of base->tip pairs separated by NaN: the 3-D renderer skips any
    // segment with a non-finite end, so every stem is drawn on its own.
    const auto stemTrace = [&](const std::string& label, const std::string& color,
                               const std::vector<std::complex<double>>& bases,
                               const std::vector<double>& heights) {
        if (bases.size() != heights.size())
            throw std::invalid_argument("complex stem plot needs one height per base point");
        Trace stems; stems.label = label; stems.color = color;
        const double gap = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t index = 0; index < bases.size(); ++index) {
            const auto base = bases[index];
            const double raw = heights[index];
            if (!std::isfinite(base.real()) || !std::isfinite(base.imag()) || !std::isfinite(raw)) continue;
            if (plot.logHeights && !(raw > 0)) continue;
            const double tip = plot.logHeights ? std::log10(raw) : raw;
            stems.x.insert(stems.x.end(), {base.real(), base.real(), gap});
            stems.y.insert(stems.y.end(), {base.imag(), base.imag(), gap});
            stems.z.insert(stems.z.end(), {floor, tip, gap});
        }
        return stems;
    };
    figure.traces.push_back(stemTrace(plot.seriesLabel, plot.color, plot.bases, plot.heights));
    for (const auto& overlay : plot.overlays)
        figure.traces.push_back(stemTrace(overlay.label, overlay.color, overlay.bases, overlay.heights));
    for (const auto& curve : plot.curves) {
        if (curve.bases.size() != curve.heights.size())
            throw std::invalid_argument("complex stem plot curve needs one height per base point");
        Trace line; line.label = curve.label; line.color = curve.color;
        for (std::size_t index = 0; index < curve.bases.size(); ++index) {
            const double raw = curve.heights[index];
            const bool shown = std::isfinite(raw) && (!plot.logHeights || raw > 0);
            line.x.push_back(shown ? curve.bases[index].real() : gapValue);
            line.y.push_back(shown ? curve.bases[index].imag() : gapValue);
            line.z.push_back(shown ? (plot.logHeights ? std::log10(raw) : raw) : gapValue);
        }
        figure.traces.push_back(std::move(line));
    }
    return figure;
}

}  // namespace dsl

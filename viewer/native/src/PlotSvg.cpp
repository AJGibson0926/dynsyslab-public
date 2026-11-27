// Vector export of a Viewer payload page (PlotSvg.hpp). Draws the same traces, value scales and
// focus ranges as the canvas, as plain SVG 1.1 (W3C, Scalable Vector Graphics 1.1, 2nd ed., 2011)
// that browsers, Inkscape and LaTeX toolchains read. Heat maps use the viridis colour map (van der
// Walt & Smith, matplotlib, 2015); three-state views use an orthographic projection at the
// figure's camera.
#include "PlotSvg.hpp"

#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <sstream>

namespace dsl {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr const char* kPalette[] = {"#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
                                    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"};

std::string escaped(const std::string& text) {
    std::string out;
    for (const char c : text) {
        if (c == '&') out += "&amp;"; else if (c == '<') out += "&lt;"; else if (c == '>') out += "&gt;";
        else if (c == '"') out += "&quot;"; else if (c != '$') out += c;
    }
    return out;
}

// A plot label (the Viewer's TeX-like math text: _{sub}, ^{sup}, \frac{a}{b}, \psi, ...) as SVG
// text content with sub- and superscript tspans; fractions are written a/b.
std::string mathSpans(const std::vector<MathTextSpan>& spans) {
    std::string out;
    for (const auto& span : spans) {
        if (span.role == MathTextRole::fraction) {
            out += "(" + mathSpans(span.numerator) + ")/(" + mathSpans(span.denominator) + ")";
            continue;
        }
        auto text = escaped(narrow(span.text));
        if (!span.numerator.empty()) text += mathSpans(span.numerator);
        if (span.role == MathTextRole::subscript) out += "<tspan baseline-shift=\"sub\" font-size=\"75%\">" + text + "</tspan>";
        else if (span.role == MathTextRole::superscript) out += "<tspan baseline-shift=\"super\" font-size=\"75%\">" + text + "</tspan>";
        else out += text;
    }
    return out;
}

std::string label(const std::string& text) { return mathSpans(parseMathText(wide(text))); }

std::string plainSpans(const std::vector<MathTextSpan>& spans) {
    std::string out;
    for (const auto& span : spans) {
        if (span.role == MathTextRole::fraction) { out += "(" + plainSpans(span.numerator) + ")/(" + plainSpans(span.denominator) + ")"; continue; }
        out += (span.role == MathTextRole::subscript ? "_" : span.role == MathTextRole::superscript ? "^" : "") + narrow(span.text)
             + plainSpans(span.numerator);
    }
    return out;
}

std::string number(double value) {
    char text[32];
    std::snprintf(text, sizeof(text), "%.6g", value);
    return text;
}

std::string viridis(double t) {
    static constexpr std::array<std::array<double, 3>, 5> stops{{
        {68, 1, 84}, {59, 82, 139}, {33, 145, 140}, {94, 201, 98}, {253, 231, 37}}};
    t = std::clamp(std::isfinite(t) ? t : 0.0, 0.0, 1.0) * 4.0;
    const auto k = std::min<std::size_t>(3, static_cast<std::size_t>(t));
    const double f = t - static_cast<double>(k);
    char text[8];
    std::snprintf(text, sizeof(text), "#%02X%02X%02X",
                  static_cast<int>(stops[k][0] + f * (stops[k + 1][0] - stops[k][0])),
                  static_cast<int>(stops[k][1] + f * (stops[k + 1][1] - stops[k][1])),
                  static_cast<int>(stops[k][2] + f * (stops[k + 1][2] - stops[k][2])));
    return text;
}

std::vector<double> niceTicks(double low, double high) {
    if (!(high > low)) return {low};
    const double span = high - low, raw = span / 5.0;
    const double magnitude = std::pow(10.0, std::floor(std::log10(raw)));
    const double step = magnitude * (raw / magnitude < 1.5 ? 1 : raw / magnitude < 3.5 ? 2 : raw / magnitude < 7.5 ? 5 : 10);
    std::vector<double> ticks;
    for (double tick = std::ceil(low / step) * step; tick <= high + 1e-9 * span; tick += step)
        ticks.push_back(std::abs(tick) < 1e-12 * span ? 0.0 : tick);
    return ticks;
}

std::string dashes(const std::string& style) {
    if (style == "dash" || style == "time-dash") return " stroke-dasharray=\"7 4\"";
    if (style == "sparse-dash") return " stroke-dasharray=\"4 8\"";
    if (style == "dot") return " stroke-dasharray=\"1.5 3\"";
    if (style == "dash-dot") return " stroke-dasharray=\"7 3 1.5 3\"";
    return {};
}

}  // namespace

// The plane the figure is drawn in: (x, y), or (x, y, z) projected at the figure's camera.
Projected projectedTrace(const Figure& figure, const Trace& trace) {
    if (figure.kind != "three_d" || trace.z.size() != trace.x.size()) return {trace.x, trace.y};
    const double azimuth = figure.cameraAzimuthDegrees.value_or(42.0) * kPi / 180.0;
    const double elevation = figure.cameraElevationDegrees.value_or(24.0) * kPi / 180.0;
    Projected out;
    for (std::size_t i = 0; i < trace.x.size(); ++i) {
        const double u = trace.x[i] * std::cos(azimuth) - trace.y[i] * std::sin(azimuth);
        const double depth = trace.x[i] * std::sin(azimuth) + trace.y[i] * std::cos(azimuth);
        out.x.push_back(u);
        out.y.push_back(trace.z[i] * std::cos(elevation) - depth * std::sin(elevation));
    }
    return out;
}


std::string plainLabel(const std::string& text) { return plainSpans(parseMathText(wide(text))); }

std::vector<std::size_t> pageFigures(const Payload& payload, int page) {
    const bool paged = std::any_of(payload.figures.begin(), payload.figures.end(),
                                   [](const Figure& f) { return f.dashboardPage >= 0; });
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < payload.figures.size(); ++i)
        if (!paged || payload.figures[i].dashboardPage == page) indices.push_back(i);
    return indices;
}

int pageCount(const Payload& payload) {
    int pages = 1;
    for (const auto& figure : payload.figures) pages = std::max(pages, figure.dashboardPage + 1);
    return pages;
}

namespace {
// Legend entries: labelled line or marker traces, at most eight, in trace order.
struct SvgLegendEntry { std::string label, plain, color, style; };
std::vector<SvgLegendEntry> svgLegendEntries(const Figure& figure) {
    std::vector<SvgLegendEntry> entries;
    for (std::size_t t = 0; t < figure.traces.size() && entries.size() < 8; ++t) {
        const auto& trace = figure.traces[t];
        if (trace.label.empty() || trace.heatmap) continue;
        entries.push_back({trace.label, plainLabel(trace.label), trace.color.empty() ? kPalette[t % 10] : trace.color,
                           trace.marker ? "markers" : trace.lineStyle});
    }
    return entries;
}
double svgTextWidth(const std::string& plain, double size) {  // average advance of the sans-serif face
    std::size_t characters = 0;
    for (const unsigned char c : plain) characters += (c & 0xC0) != 0x80;  // UTF-8 code points
    return 0.56 * size * static_cast<double>(characters);
}
std::string svgLegendSample(const SvgLegendEntry& entry, double x, double y) {
    if (entry.style == "markers" || entry.style == "open-markers")
        return "<circle cx=\"" + number(x + 9) + "\" cy=\"" + number(y) + "\" r=\"3\" fill=\"" + entry.color + "\"/>";
    return "<line x1=\"" + number(x) + "\" x2=\"" + number(x + 18) + "\" y1=\"" + number(y) + "\" y2=\"" + number(y)
        + "\" stroke=\"" + entry.color + "\" stroke-width=\"2\"" + dashes(entry.style) + "/>";
}
}  // namespace

std::string figureSvgBody(const Figure& figure, double left, double top, double width, double height) {
    std::ostringstream out;
    const auto scale = ordinateScale(figure);
    out << "<g class=\"figure\"><rect x=\"" << number(left) << "\" y=\"" << number(top) << "\" width=\""
        << number(width) << "\" height=\"" << number(height) << "\" fill=\"#FFFFFF\"/>"
        << "<text x=\"" << number(left + width / 2) << "\" y=\"" << number(top + 20)
        << "\" text-anchor=\"middle\" font-size=\"14\" font-weight=\"600\">" << label(figure.title) << "</text>";
    // Data limits in plotted (scaled) coordinates.
    constexpr double inf = std::numeric_limits<double>::infinity();
    double x0 = inf, x1 = -inf, y0 = inf, y1 = -inf;
    const auto plotted = [&](double y) { return figure.kind == "three_d" ? y : scale.forward(y); };
    for (const auto& trace : figure.traces) {
        const auto p = projectedTrace(figure, trace);
        for (std::size_t i = 0; i < p.x.size() && i < p.y.size(); ++i) {
            const double y = plotted(p.y[i]);
            if (!std::isfinite(p.x[i]) || !std::isfinite(y)) continue;
            if (figure.focusX && (p.x[i] < figure.focusX->first || p.x[i] > figure.focusX->second)) continue;
            x0 = std::min(x0, p.x[i]); x1 = std::max(x1, p.x[i]); y0 = std::min(y0, y); y1 = std::max(y1, y);
        }
    }
    if (figure.focusX) { x0 = figure.focusX->first; x1 = figure.focusX->second; }
    if (!std::isfinite(x0)) { x0 = 0; x1 = 1; y0 = 0; y1 = 1; }
    if (x1 <= x0) { x0 -= 0.5; x1 += 0.5; }
    if (y1 <= y0) { y0 -= 0.5; y1 += 0.5; }
    const double padY = 0.05 * (y1 - y0); y0 -= padY; y1 += padY;
    // A right ordinate (fitted or given) takes a margin on the right; it applies to linear plots.
    const auto* right = figure.rightAxis && scale.kind == ValueScaleKind::linear && figure.kind != "three_d"
        ? &*figure.rightAxis : nullptr;
    const bool banded = right && right->leftBand && right->rightBand;
    const double plotLeft = left + 70, plotWidth = width - 90 - (right ? 60 : 0);
    double plotTop = top + 34, plotHeight = height - 82;
    // Legend: inside, in the first corner with no data under it; a time or spectrum plot with data
    // in every corner, or legendPlacement "top", gets one row under the title (as PlotLegend.cpp).
    const auto entries = svgLegendEntries(figure);
    double boxWidth = 0;
    for (const auto& entry : entries) boxWidth = std::max(boxWidth, 26 + svgTextWidth(entry.plain, 11));
    const double boxHeight = 16.0 * static_cast<double>(entries.size()) + 4;
    std::optional<std::pair<double, double>> corner;  // top-left of the legend box
    const auto sxAt = [&](double x, double l, double w) { return l + (x - x0) / (x1 - x0) * w; };
    const auto syAt = [&](double y, double t, double h) { return t + h - (y - y0) / (y1 - y0) * h; };
    bool row = !entries.empty() && figure.legendPlacement == "top";
    if (!entries.empty() && !row) {
        const std::array<std::pair<double, double>, 4> corners{{{plotLeft + plotWidth - boxWidth - 6, plotTop + 6},
            {plotLeft + 6, plotTop + 6}, {plotLeft + plotWidth - boxWidth - 6, plotTop + plotHeight - boxHeight - 6},
            {plotLeft + 6, plotTop + plotHeight - boxHeight - 6}}};
        const bool searched = figure.kind == "time" || figure.kind == "fourier";
        for (const auto& candidate : corners) {
            bool clear = true;
            for (const auto& trace : figure.traces) {
                if (!searched || trace.heatmap) continue;
                const auto p = projectedTrace(figure, trace);
                for (std::size_t i = 0; i < p.x.size() && i < p.y.size() && clear; ++i) {
                    const double px = sxAt(p.x[i], plotLeft, plotWidth), py = syAt(plotted(p.y[i]), plotTop, plotHeight);
                    clear = !(px >= candidate.first - 4 && px <= candidate.first + boxWidth + 4
                              && py >= candidate.second - 4 && py <= candidate.second + boxHeight + 4);
                }
                if (!clear) break;
            }
            if (clear || !searched) { corner = candidate; break; }
        }
        row = !corner && figure.legendPlacement.empty();
    }
    if (row) { plotTop += 20; plotHeight -= 20; }
    const auto sx = [&](double x) { return sxAt(x, plotLeft, plotWidth); };
    const auto sy = [&](double y) { return syAt(y, plotTop, plotHeight); };
    out << "<rect x=\"" << number(plotLeft) << "\" y=\"" << number(plotTop) << "\" width=\"" << number(plotWidth)
        << "\" height=\"" << number(plotHeight) << "\" fill=\"none\" stroke=\"#555555\"/>";
    for (const double tick : niceTicks(x0, x1))
        out << "<line x1=\"" << number(sx(tick)) << "\" x2=\"" << number(sx(tick)) << "\" y1=\"" << number(plotTop)
            << "\" y2=\"" << number(plotTop + plotHeight) << "\" stroke=\"#E5E5E5\"/><text x=\"" << number(sx(tick))
            << "\" y=\"" << number(plotTop + plotHeight + 16) << "\" text-anchor=\"middle\" font-size=\"11\">"
            << number(tick) << "</text>";
    const auto yTick = [&](double at, const std::string& text, bool onRight) {
        out << "<line x1=\"" << number(plotLeft) << "\" x2=\"" << number(plotLeft + plotWidth) << "\" y1=\""
            << number(sy(at)) << "\" y2=\"" << number(sy(at)) << "\" stroke=\"#E5E5E5\"/><text x=\""
            << number(onRight ? plotLeft + plotWidth + 6 : plotLeft - 6) << "\" y=\"" << number(sy(at) + 4)
            << "\" text-anchor=\"" << (onRight ? "start" : "end") << "\" font-size=\"11\">" << text << "</text>";
    };
    const auto yTitle = [&](const std::string& text, double low, double high, bool onRight) {
        if (text.empty()) return;
        out << "<text transform=\"translate(" << number(onRight ? plotLeft + plotWidth + 52 : left + 16) << ' '
            << number((sy(low) + sy(high)) / 2) << ") rotate(-90)\" text-anchor=\"middle\" font-size=\"12\">"
            << label(text) << "</text>";
    };
    if (banded) {  // stacked axes: each band its own ticks, in its own units
        for (const double tick : niceTicks(right->leftBand->first, right->leftBand->second)) yTick(tick, number(tick), false);
        for (const double value : niceTicks((right->rightBand->first - right->offset) / right->scale,
                                            (right->rightBand->second - right->offset) / right->scale))
            yTick(right->offset + right->scale * value, number(value), true);
        yTitle(figure.yLabel, right->leftBand->first, right->leftBand->second, false);
        yTitle(right->label, right->rightBand->first, right->rightBand->second, true);
    } else {
        for (const double tick : niceTicks(y0, y1)) {
            const double value = figure.kind == "three_d" ? tick : scale.inverse(tick);
            yTick(tick, number(value), false);
            if (right) out << "<text x=\"" << number(plotLeft + plotWidth + 6) << "\" y=\"" << number(sy(tick) + 4)
                           << "\" font-size=\"11\">" << number((value - right->offset) / right->scale) << "</text>";
        }
        yTitle(figure.yLabel, y0, y1, false);
        if (right) yTitle(right->label, y0, y1, true);
    }
    out << "<text x=\"" << number(plotLeft + plotWidth / 2) << "\" y=\"" << number(top + height - 10)
        << "\" text-anchor=\"middle\" font-size=\"12\">" << label(figure.xLabel) << "</text>";
    out << "<svg x=\"" << number(plotLeft) << "\" y=\"" << number(plotTop) << "\" width=\"" << number(plotWidth)
        << "\" height=\"" << number(plotHeight) << "\" viewBox=\"" << number(plotLeft) << ' ' << number(plotTop) << ' '
        << number(plotWidth) << ' ' << number(plotHeight) << "\" overflow=\"hidden\">";
    for (std::size_t t = 0; t < figure.traces.size(); ++t) {
        const auto& trace = figure.traces[t];
        const std::string color = trace.color.empty() ? kPalette[t % 10] : trace.color;
        if (trace.heatmap && trace.zShape.size() == 2 && !trace.z.empty()) {
            const auto rows = static_cast<std::size_t>(trace.zShape[0]), cols = static_cast<std::size_t>(trace.zShape[1]);
            const auto [lowIt, highIt] = std::minmax_element(trace.z.begin(), trace.z.end());
            const double low = trace.colorLimits ? trace.colorLimits->first : *lowIt;
            const double high = trace.colorLimits ? trace.colorLimits->second : *highIt;
            const double cw = plotWidth / static_cast<double>(cols), ch = plotHeight / static_cast<double>(rows);
            for (std::size_t r = 0; r < rows; ++r)
                for (std::size_t c = 0; c < cols && r * cols + c < trace.z.size(); ++c)
                    out << "<rect x=\"" << number(plotLeft + c * cw) << "\" y=\"" << number(plotTop + plotHeight - (r + 1) * ch)
                        << "\" width=\"" << number(cw + 0.3) << "\" height=\"" << number(ch + 0.3) << "\" fill=\""
                        << viridis(high > low ? (trace.z[r * cols + c] - low) / (high - low) : 0.5) << "\"/>";
            continue;
        }
        const auto p = projectedTrace(figure, trace);
        if (trace.marker || trace.lineStyle == "markers" || trace.lineStyle == "open-markers") {
            const bool open = trace.markerShape == "open-circle" || trace.lineStyle == "open-markers";
            for (std::size_t i = 0; i < p.x.size() && i < p.y.size(); ++i) {
                const double y = plotted(p.y[i]);
                if (!std::isfinite(p.x[i]) || !std::isfinite(y)) continue;
                out << "<circle cx=\"" << number(sx(p.x[i])) << "\" cy=\"" << number(sy(y)) << "\" r=\"2.6\" "
                    << (open ? "fill=\"none\" stroke=\"" + color + "\"" : "fill=\"" + color + "\"") << "/>";
            }
        } else {
            out << "<path fill=\"none\" stroke=\"" << color << "\" stroke-width=\""
                << number(trace.lineWidthPoints > 0 ? trace.lineWidthPoints : 1.4) << "\" stroke-opacity=\""
                << number(trace.opacity) << '"' << dashes(trace.lineStyle) << " d=\"";
            bool pen = false;
            for (std::size_t i = 0; i < p.x.size() && i < p.y.size(); ++i) {
                const double y = plotted(p.y[i]);
                if (!std::isfinite(p.x[i]) || !std::isfinite(y)) { pen = false; continue; }
                out << (pen ? 'L' : 'M') << number(sx(p.x[i])) << ',' << number(sy(y)) << ' ';
                pen = true;
            }
            out << "\"/>";
        }
    }
    out << "</svg>";
    if (row) {  // one centred row between the title and the plot
        double total = 0;
        for (const auto& entry : entries) total += 26 + svgTextWidth(entry.plain, 11) + 18;
        double x = left + std::max(8.0, (width - total + 18) / 2);
        for (const auto& entry : entries) {
            out << svgLegendSample(entry, x, top + 40) << "<text x=\"" << number(x + 24) << "\" y=\"" << number(top + 44)
                << "\" font-size=\"11\">" << label(entry.label) << "</text>";
            x += 26 + svgTextWidth(entry.plain, 11) + 18;
        }
    } else if (corner) {
        for (std::size_t k = 0; k < entries.size(); ++k) {
            const double y = corner->second + 10 + 16.0 * static_cast<double>(k);
            out << svgLegendSample(entries[k], corner->first + 2, y) << "<text x=\"" << number(corner->first + 26)
                << "\" y=\"" << number(y + 4) << "\" font-size=\"11\">" << label(entries[k].label) << "</text>";
        }
    }
    out << "</g>";
    return out.str();
}

std::string payloadSvg(const Payload& payload, int page, int width, int height) {
    const auto indices = pageFigures(payload, page);
    const int columns = std::max(1, std::min<int>(payload.columns, static_cast<int>(indices.size())));
    const int rows = std::max(1, static_cast<int>((indices.size() + columns - 1) / columns));
    const double cellWidth = static_cast<double>(width) / columns, cellHeight = static_cast<double>(height) / rows;
    std::ostringstream out;
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width << "\" height=\"" << height
        << "\" viewBox=\"0 0 " << width << ' ' << height << "\" font-family=\"Segoe UI, Helvetica, Arial, sans-serif\">"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#FFFFFF\"/>";
    for (std::size_t k = 0; k < indices.size(); ++k)
        out << figureSvgBody(payload.figures[indices[k]], (k % columns) * cellWidth, (k / columns) * cellHeight,
                             cellWidth, cellHeight);
    out << "</svg>";
    return out.str();
}

}  // namespace dsl

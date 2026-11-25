// Viewer payload for dynsyslab.grid-field/v1 archives (the domain_field job): a function drawn on
// its domain.
//   1-D grid   curves (Re and Im for complex values)
//   2-D grid   heat maps; complex scalars default to domain coloring with a magnitude panel and
//              switch in place to real/imaginary or magnitude/phase (ComplexHeatmap.hpp)
//   3-D grid   five staggered slices across the third axis as 2-D panels, plus a translucent
//              stacked-slice view in 3-D (kind "domain_field_3d")
// Points that are not valid (undefined, push-forward conflict, no preimage) are blank in the value
// panels and drawn by status in a separate panel; push-forwards also show the fibre spread.
// Several frames (e.g. the Koopman action g o Phi^{k dt}) become one dashboard page each.
#include "DomainFieldPayload.hpp"
#include "ComplexHeatmap.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include <dynsyslab/core/GridFieldArchive.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>

namespace dsl {
namespace c = dynsyslab::core;
namespace {

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

std::string number(double v) {
    char text[32];
    std::snprintf(text, sizeof text, "%.4g", v);
    return text;
}

// Identifiers such as far_field_acoustic_echo read as words; x_1 keeps its subscript.
std::string displayName(std::string name) {
    const auto underscore = name.find('_');
    if (underscore != std::string::npos && name.size() - underscore > 2)
        std::replace(name.begin(), name.end(), '_', ' ');
    return name;
}

std::string axisName(const c::DomainGrid& grid, std::size_t axis) {
    const auto coordinate = grid.axes.at(axis).coordinate;
    return coordinate < grid.domain.names.size() && !grid.domain.names[coordinate].empty()
        ? displayName(grid.domain.names[coordinate]) : "q_" + std::to_string(coordinate + 1);
}

Figure heatFigure(const std::string& title, const std::vector<double>& x, const std::vector<double>& y,
                  std::vector<double> z, const std::string& xLabel, const std::string& yLabel) {
    Figure f{"domain_field", title, xLabel, yLabel, false};
    Trace t; t.label = title; t.heatmap = true; t.x = x; t.y = y; t.z = std::move(z);
    t.zShape = {static_cast<int>(y.size()), static_cast<int>(x.size())};
    t.encoding = "scalar";
    f.colorLabel = title;
    f.traces.push_back(std::move(t));
    return f;
}

// One output on a 2-D cut (rows = second axis): value panels for a real or complex function.
void addValuePanels(std::vector<Figure>& out, const std::string& name, const std::vector<double>& x,
                    const std::vector<double>& y, const std::vector<std::complex<double>>& values, bool complex,
                    const std::string& xLabel, const std::string& yLabel) {
    if (!complex) {
        std::vector<double> z;
        for (const auto v : values) z.push_back(v.real());
        out.push_back(heatFigure(name, x, y, std::move(z), xLabel, yLabel));
        return;
    }
    for (unsigned part = 0; part < 2; ++part) {
        auto f = heatFigure(name, x, y, {}, xLabel, yLabel);
        auto& t = f.traces.front();
        t.complexValues = values; t.complexPart = static_cast<int>(part); t.complexName = name;
        encodeComplexHeatmap(f, kDefaultComplexDisplay);
        out.push_back(std::move(f));
    }
}

// Values, status and spread of one output on the cut of the grid where axis 2 (if any) is at `slice`.
struct Cut { std::vector<std::complex<double>> values; std::vector<double> status, spread; };
Cut cut(const c::GridField& field, std::size_t slice) {
    const auto& g = field.grid;
    const std::size_t nx = g.axes[0].points, ny = g.axes.size() > 1 ? g.axes[1].points : 1;
    Cut out;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const auto k = i + nx * (j + ny * slice);
            const bool valid = field.status[k] == c::GridStatus::valid;
            out.values.push_back(valid ? field.value(k) : std::complex<double>(nan, nan));
            out.status.push_back(double(static_cast<int>(field.status[k])));
            if (!field.spread.empty()) out.spread.push_back(field.spread[k]);
        }
    return out;
}

}  // namespace

Figure buildDomainField3D(const c::GridField& field, const std::vector<std::size_t>& slices, const std::string& titleSuffix) {
    const auto& g = field.grid;
    const std::size_t nx = g.axes[0].points, ny = g.axes[1].points;
    const auto name = !field.outputNames.empty() ? field.outputNames.front() : std::string("g");
    Figure f{"three_d", name + " — slices" + titleSuffix, axisName(g, 0), axisName(g, 1), false};
    f.zLabel = axisName(g, 2);
    f.boxAxes3D = true;
    f.colorLabel = field.complexValued ? "phase (rad); brightness = magnitude" : name;
    // One colour scale for every plane of the first output.
    double low = std::numeric_limits<double>::infinity(), high = -low, largest = 0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        if (field.status[k] != c::GridStatus::valid) continue;
        const auto v = field.value(k);
        low = std::min(low, v.real()); high = std::max(high, v.real()); largest = std::max(largest, std::abs(v));
    }
    if (!(high > low)) high = low + 1;
    if (!(largest > 0)) largest = 1;
    const double x0 = g.axes[0].lower, x1 = g.axes[0].upper, y0 = g.axes[1].lower, y1 = g.axes[1].upper;
    for (const auto s : slices) {
        const double z = g.axes[2].value(s);
        Trace t;
        t.label = axisName(g, 2) + " = " + number(z);
        t.heatmap = true;
        t.encoding = field.complexValued ? "slice-domain-color" : "slice";
        t.x = {x0, x1, x1, x0, x0}; t.y = {y0, y0, y1, y1, y0}; t.z = std::vector<double>(5, z);
        t.zShape = {static_cast<int>(ny), static_cast<int>(nx)};
        t.opacity = 0.45;
        t.color = "#94a3b8";
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i) {
                const auto k = i + nx * (j + ny * s);
                const bool valid = field.status[k] == c::GridStatus::valid;
                const auto v = field.value(k);
                if (field.complexValued) {
                    t.colorValues.push_back(valid ? std::arg(v) : nan);
                    t.magnitudes.push_back(valid ? std::abs(v) / largest : nan);
                } else {
                    t.colorValues.push_back(valid ? (v.real() - low) / (high - low) : nan);
                }
            }
        f.traces.push_back(std::move(t));
    }
    return f;
}

Payload buildDomainFieldPayload(const std::filesystem::path& source, const c::NpzRecord& record) {
    const auto series = c::readGridFieldArchive(record);
    const auto& meta = series.metadata;
    Payload payload{source.u8string(), "domain_field", meta.at("label").string(source.stem().u8string()), 2};
    const auto& grid = series.frames.front().grid;
    const auto dims = grid.axes.size();
    const auto operation = meta.at("domain_field_settings").at("operation").string("evaluate");
    const bool complex = series.frames.front().complexValued;
    if (complex && dims >= 2) payload.complexDisplay = kDefaultComplexDisplay;
    const auto x = grid.axisValues(0);
    const auto y = dims > 1 ? grid.axisValues(1) : std::vector<double>{};
    std::vector<std::size_t> slices{0};
    if (dims == 3) {
        const auto n = grid.axes[2].points;
        slices = {0, (n - 1) / 4, (n - 1) / 2, 3 * (n - 1) / 4, n - 1};
        slices.erase(std::unique(slices.begin(), slices.end()), slices.end());
    }
    if (dims > 3) payload.warnings.push_back("Grids of more than three axes are shown as their first three axes.");
    for (std::size_t frame = 0; frame < series.frames.size(); ++frame) {
        const auto& field = series.frames[frame];
        const auto when = series.frames.size() > 1 ? " — t = " + number(series.frameTime[frame]) : std::string();
        std::vector<Figure> figures;
        for (std::size_t r = 0; r < field.outputs; ++r) {
            const auto one = field.output(r);
            auto name = r < field.outputNames.size() && !field.outputNames[r].empty() ? displayName(field.outputNames[r]) : "g";
            if (operation == "koopman_action") name = "U^t " + name;
            if (operation == "push_forward") name = "h_* " + name;
            if (operation == "pull_back") name = "h^* " + name;
            if (operation == "round_trip") name = "|h^* h_* g - g|";
            if (operation == "fold_set") name = "volume factor of h";
            if (dims == 1) {
                const auto values = cut(one, 0).values;
                Figure f{"domain_field_line", name + when, axisName(grid, 0), name, false};
                for (unsigned part = 0; part < (complex ? 2u : 1u); ++part) {
                    Trace t; t.label = complex ? (part ? "Im" : "Re") : name; t.x = x;
                    for (const auto v : values) t.y.push_back(part ? v.imag() : v.real());
                    f.traces.push_back(std::move(t));
                }
                figures.push_back(std::move(f));
                continue;
            }
            for (const auto s : slices) {
                const auto where = dims == 3 ? " at " + axisName(grid, 2) + " = " + number(grid.axes[2].value(s)) : std::string();
                const auto values = cut(one, s);
                addValuePanels(figures, name + where + when, x, y, values.values, complex, axisName(grid, 0), axisName(grid, 1));
                if (!one.spread.empty()) {
                    auto f = heatFigure("fibre spread of " + name + where, x, y, values.spread, axisName(grid, 0), axisName(grid, 1));
                    f.colorLabel = "max |g - fibre mean|";
                    figures.push_back(std::move(f));
                }
            }
        }
        // Status: kept distinct from the values (blank there).
        if (std::any_of(field.status.begin(), field.status.end(), [](auto s) { return s != c::GridStatus::valid; }) && dims >= 2)
            for (const auto s : slices) {
                auto f = heatFigure("Status" + std::string(dims == 3 ? " at " + axisName(grid, 2) + " = " + number(grid.axes[2].value(s)) : "") + when,
                                    x, y, cut(field, s).status, axisName(grid, 0), axisName(grid, 1));
                f.colorLabel = "0 valid · 1 undefined · 2 conflict · 3 no preimage";
                figures.push_back(std::move(f));
            }
        if (dims == 3) figures.push_back(buildDomainField3D(field, slices, when));
        for (auto& f : figures) {
            f.dashboardPage = static_cast<int>(frame);
            payload.figures.push_back(std::move(f));
        }
        payload.pageLabels.push_back(series.frames.size() > 1 ? "t = " + number(series.frameTime[frame]) : "Domain");
    }
    // Several frames (e.g. the Koopman action) play as an animation, with one colour range per
    // real-valued panel across every frame so colours compare between frames.
    if (series.frames.size() > 1 && !complex) {
        std::map<std::size_t, std::pair<double, double>> range;  // panel position -> range
        std::map<int, std::size_t> position;
        std::vector<std::size_t> slot(payload.figures.size());
        for (std::size_t i = 0; i < payload.figures.size(); ++i) slot[i] = position[payload.figures[i].dashboardPage]++;
        for (std::size_t i = 0; i < payload.figures.size(); ++i)
            for (const auto& t : payload.figures[i].traces) {
                if (!t.heatmap || t.encoding != "scalar") continue;
                auto& r = range.try_emplace(slot[i], std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()).first->second;
                for (const auto z : t.z) if (std::isfinite(z)) { r.first = std::min(r.first, z); r.second = std::max(r.second, z); }
            }
        for (std::size_t i = 0; i < payload.figures.size(); ++i)
            for (auto& t : payload.figures[i].traces)
                if (t.heatmap && t.encoding == "scalar" && range.count(slot[i]) && range[slot[i]].second > range[slot[i]].first)
                    t.colorLimits = range[slot[i]];
    }
    if (series.frames.size() > 1) {
        Payload::Animation animation;
        const auto& settings = meta.at("animation");
        animation.frameSeconds = settings.at("frame_seconds").number(animation.frameSeconds);
        animation.loop = settings.at("loop").boolean(true);
        animation.autoplay = settings.at("autoplay").boolean(true);
        payload.animation = animation;
    }
    return payload;
}

}  // namespace dsl

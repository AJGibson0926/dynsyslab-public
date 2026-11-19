#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
namespace {
struct Range {
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    void include(double value) {
        if (std::isfinite(value)) { low = std::min(low, value); high = std::max(high, value); }
    }
    void normalize(bool padding = true) {
        if (!std::isfinite(low) || !std::isfinite(high)) { low = -1; high = 1; }
        if (std::abs(high - low) < 1e-14) { low -= 0.5; high += 0.5; }
        const double margin = padding ? (high - low) * 0.045 : 0;
        low -= margin; high += margin;
    }
};
}
bool plotReferenceVisible(const Figure& figure, const Figure::ReferencePoint& reference) {
    if (std::any_of(figure.traces.begin(),figure.traces.end(),
        [](const auto& trace) { return trace.heatmap; })) return false;
    return figure.kind == "phase" || (figure.kind == "three_d" && !reference.equilibrium);
}
AxisLimits automaticAxisLimits(const Figure& figure) {
    Range xr, yr;
    const auto scale = ordinateScale(figure);
    for (const Trace& trace : figure.traces) {
        if (trace.role.rfind("phase_",0)==0) continue;
        if (trace.heatmap) {
            for (double value : trace.x) xr.include(value);
            for (double value : trace.y) yr.include(value);
            continue;
        }
        const std::size_t count = std::min(trace.x.size(), trace.y.size());
        for (std::size_t index = 0; index < count; ++index) {
            double x = trace.x[index], y = trace.y[index];
            if (!trace.z.empty() && index < trace.z.size()) {
                x -= 0.32 * trace.z[index];
                y += 0.20 * trace.z[index];
            }
            xr.include(x); yr.include(scale.forward(y));
        }
    }
    for (const auto& point : figure.references) {
        if (!plotReferenceVisible(figure,point) || point.equilibrium) continue;
        if (figure.kind == "phase" || !point.equilibrium) xr.include(point.x);
        yr.include(scale.forward(point.y));
    }
    const bool padding = !std::any_of(figure.traces.begin(),figure.traces.end(),[](const auto& t) { return t.heatmap; });
    xr.normalize(padding); yr.normalize(padding);
    // A short off-scale burst moves to the automatic detail inset; the main panel then
    // keeps the meaningful ordinate range (see automaticDetailInset).
    if (!figure.inset && scale.kind == ValueScaleKind::linear)
        if (const auto detail = automaticDetailInset(figure); detail && detail->transientWindow)
            return {xr.low, xr.high, detail->mainYMin, detail->mainYMax};
    return {xr.low, xr.high, scale.inverse(yr.low), scale.inverse(yr.high)};
}


bool canDisplayCoordinates(const Figure& figure,const std::string& coordinates) {
    return figure.coordinates==coordinates ||
        (figure.alternateCoordinates&&figure.alternateCoordinates->name==coordinates);
}
void changeDisplayCoordinates(Figure& figure,const std::string& coordinates) {
    if(figure.coordinates==coordinates||!canDisplayCoordinates(figure,coordinates)) return;
    const auto target=*figure.alternateCoordinates;
    Figure::CoordinateAxes previous{figure.coordinates,figure.xLabel,figure.yLabel,
        1/target.xScale,1/target.yScale,1/target.zScale,1/target.timeScale};
    previous.colorLabel = figure.colorLabel;
    previous.zLabel = figure.zLabel;
    for(std::size_t i=0;i<figure.traces.size();++i) {
        auto& trace=figure.traces[i]; previous.traceLabels.push_back(trace.label);
        for(auto& x:trace.x) x*=target.xScale;
        for(auto& y:trace.y) y*=target.yScale;
        if(!trace.heatmap || figure.kind == "observable_heatmap") for(auto& z:trace.z) z*=target.zScale;
        if(trace.encoding=="time") for(auto& t:trace.colorValues) t*=target.timeScale;
        if(i<target.traceLabels.size()) trace.label=target.traceLabels[i];
    }
    for(auto& point:figure.references) { point.x*=target.xScale;point.y*=target.yScale; }
    if(figure.inset) {
        figure.inset->xMin*=target.xScale;figure.inset->xMax*=target.xScale;
        figure.inset->yMin*=target.yScale;figure.inset->yMax*=target.yScale;
    }
    if (figure.title == figure.yLabel + " versus time")
        figure.title = target.yLabel + " versus time";
    figure.xLabel=target.xLabel;figure.yLabel=target.yLabel;
    if (!target.zLabel.empty()) figure.zLabel=target.zLabel;
    figure.coordinates=coordinates;
    if (!target.colorLabel.empty()) figure.colorLabel = target.colorLabel;
    else if (figure.kind == "phase" && !figure.colorLabel.empty())
        figure.colorLabel = coordinates == "dimensionless" ? "\\tau" : "t\\; (s)";
    figure.alternateCoordinates=previous;
}
}

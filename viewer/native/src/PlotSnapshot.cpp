#include "dynsyslab_viewer/ViewerSession.hpp"
#include "dynsyslab_viewer/MathText.hpp"
#include "dynsyslab_viewer/PlotPresentation.hpp"
#include "dynsyslab_viewer/PlotScale.hpp"
#include "PlotSeries.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace dsl {
namespace {
using Value = json::Value;
Value strings(const std::vector<std::string>& values) {
    Value::Array result; for (const auto& value : values) result.emplace_back(value);
    return result;
}
Value numbers(const std::vector<double>& values) {
    Value::Array result; result.reserve(values.size());
    for (const auto value : values) result.push_back(std::isfinite(value) ? Value(value) : Value());
    return result;
}
Value range(const std::vector<double>& values) {
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (const auto value : values) if (std::isfinite(value)) {
        low = std::min(low, value); high = std::max(high, value);
    }
    return std::isfinite(low) ? numbers({low, high}) : Value();
}
Value figureValue(const Figure& figure, bool samples) {
    Value::Object result{{"kind",figure.kind},{"title",scientificMathLabel(figure.title)},
        {"x_label",scientificMathLabel(figure.xLabel)},{"y_label",scientificMathLabel(figure.yLabel)},
        {"coordinates",figure.coordinates},{"square",figure.square},{"color_label",figure.colorLabel},
        {"z_label",scientificMathLabel(figure.zLabel)},
        {"title_color",figure.titleColor},
        {"dashboard_page",static_cast<double>(figure.dashboardPage)},
        {"grid_span",static_cast<double>(figure.gridSpan)},{"row_span",static_cast<double>(figure.rowSpan)}};
    const auto scale = resolveValueScale(figure);
    result["value_scale"] = figure.valueScale;
    result["observable_identity"] = figure.observableIdentity;
    result["effective_value_scale"] = Value::Object{{"kind",scale.name()},
        {"target",supportsValueScale(figure) ? scalesColorValues(figure) ? "color" : "y" : "none"},
        {"linear_threshold",scale.linearThreshold},{"reason",scale.reason}};
    const auto limits = automaticAxisLimits(figure);
    result["axis_limits"] = Value::Object{{"x",numbers({limits.xMin,limits.xMax})},
        {"y",numbers({limits.yMin,limits.yMax})}};
    Value::Array traces;
    for (const auto& trace : figure.traces) {
        Value::Object entry{{"label",scientificMathLabel(trace.label)},{"role",trace.role},
            {"color",trace.color},{"line_style",trace.lineStyle},{"encoding",trace.encoding},
            {"opacity",trace.opacity},
            {"marker",trace.marker},{"heatmap",trace.heatmap},
            {"line_width_points",static_cast<double>(plotSeriesWidthPoints(trace))},
            {"selection_id",static_cast<double>(trace.selectionId)},{"selected",trace.selected},
            {"sample_count",static_cast<double>(std::min(trace.x.size(),trace.y.size()))},
            {"z_range",range(trace.z)},{"x_range",range(trace.x)},{"y_range",range(trace.y)},{"time_range",
                trace.encoding == "time" ? range(trace.colorValues) : Value()}};
        if (samples) {
            entry["x"] = numbers(trace.x); entry["y"] = numbers(trace.y); entry["z"] = numbers(trace.z);
            entry["color_values"] = numbers(trace.colorValues); entry["magnitudes"] = numbers(trace.magnitudes);
            Value::Array shape; for (int n : trace.zShape) shape.emplace_back(static_cast<double>(n));
            entry["z_shape"] = std::move(shape);
        }
        traces.emplace_back(std::move(entry));
    }
    result["traces"] = std::move(traces);
    Value::Array references;
    for (const auto& ref : figure.references) if (plotReferenceVisible(figure,ref))
        references.emplace_back(Value::Object{{"x",ref.x},{"y",ref.y},{"label",ref.label},
            {"equilibrium",ref.equilibrium},{"color","#ff0000"},
            {"marker",ref.equilibrium ? "filled-circle" : "open-circle"}});
    result["references"] = std::move(references);
    if (figure.alternateCoordinates) {
        const auto& alt = *figure.alternateCoordinates;
        result["alternate_coordinates"] = Value::Object{{"name",alt.name},{"x_label",alt.xLabel},
            {"y_label",alt.yLabel},{"x_scale",alt.xScale},{"y_scale",alt.yScale},
            {"color_label",alt.colorLabel},{"z_label",scientificMathLabel(alt.zLabel)},
            {"z_scale",alt.zScale},{"time_scale",alt.timeScale},{"trace_labels",strings(alt.traceLabels)}};
    }
    if (figure.inset) {
        const auto& inset = *figure.inset; Value::Array indices;
        for (const auto i : inset.traceIndices) indices.emplace_back(static_cast<double>(i));
        result["detail_inset"] = Value::Object{{"xlim",numbers({inset.xMin,inset.xMax})},
            {"ylim",numbers({inset.yMin,inset.yMax})},{"trace_indices",std::move(indices)},
            {"automatic",inset.automatic}};
    }
    return result;
}
}
json::Value inspectPlot(const Payload& payload, bool includeSamples) {
    Value::Array modes; for (const auto i : payload.selectedModes) modes.emplace_back(static_cast<double>(i));
    Value::Array figures;
    for (const auto& figure : payload.figures) figures.push_back(figureValue(figure,includeSamples));
    Value result = Value::Object{{"schema",includeSamples?"dynsyslab.viewer.payload/v1":"dynsyslab.viewer.plot-inspection/v1"},
        {"source",payload.source},{"view",payload.view},
        {"layout",payload.layout},{"complex_display",payload.complexDisplay},{"selected_modes",std::move(modes)},
        {"tile",Value::Object{{"name",payload.tileName},{"columns",static_cast<double>(payload.columns)}}},
        {"page_labels",strings(payload.pageLabels)},
        {"warnings",strings(payload.warnings)},{"figures",std::move(figures)},
        {"rendering",Value::Object{{"sampling","published archive samples; no spline interpolation"},
            {"plot_background","#ffffff"},{"antialias",true},{"supersampling",2.0},{"phase_colormap","jet"},{"scalar_colormap","jet"},
            {"phase_linewidth_points",static_cast<double>(phaseLineWidthPoints)},
            {"curve_linewidth_points",static_cast<double>(curveLineWidthPoints)},
            {"target_linewidth_points",static_cast<double>(targetLineWidthPoints)},
            {"target_dash_pattern",numbers({targetDashPattern[0],targetDashPattern[1]})},
            {"target_dash_units","stroke widths"},{"phase_time","start of each saved segment"}}}};
    result.object().at("rendering").object()["three_d"] = Value::Object{
        {"renderer","native camera"},{"default_projection","perspective"},
        {"depth_sorted",true},{"interactive",true}};
    return result;
}
}

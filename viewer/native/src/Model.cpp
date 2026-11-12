#include "dynsyslab_viewer/Model.hpp"
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/System.hpp>

#include <windows.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dsl {
namespace {

std::vector<std::string> strings(const json::Value& value) {
    std::vector<std::string> result;
    for (const auto& item : value.array()) result.push_back(item.string());
    return result;
}

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    result.reserve(value.array().size());
    for (const auto& item : value.array())
        result.push_back(item.isNull() ? std::numeric_limits<double>::quiet_NaN() : item.number());
    return result;
}

std::vector<int> integers(const json::Value& value) {
    std::vector<int> result;
    for (const auto& item : value.array()) result.push_back(static_cast<int>(item.number()));
    return result;
}

Trace traceFrom(const json::Value& value) {
    Trace result;
    result.label = value.at("label").string();
    result.role = value.at("role").string();
    result.color = value.at("color").string();
    result.lineStyle = value.at("line_style").string("solid");
    result.lineWidthPoints = value.at("line_width_points").number(0.);
    result.encoding = value.at("encoding").string();
    result.opacity = std::clamp(value.at("opacity").number(1.0), 0.0, 1.0);
    result.selectionId = static_cast<int>(value.at("selection_id").number(-1));
    result.selected = value.at("selected").boolean();
    result.marker = value.at("marker").boolean();
    result.heatmap = value.at("heatmap").boolean();
    result.x = numbers(value.at("x"));
    result.y = numbers(value.at("y"));
    result.z = numbers(value.at("z"));
    result.colorValues = numbers(value.at("color_values"));
    result.magnitudes = numbers(value.at("magnitudes"));
    result.zShape = integers(value.at("z_shape"));
    return result;
}

std::optional<DetailInset> insetFrom(const json::Value& value) {
    if (value.isNull()) return std::nullopt;
    const auto& x = value.at("xlim").array();
    const auto& y = value.at("ylim").array();
    if (x.size() != 2 || y.size() != 2) return std::nullopt;
    DetailInset result{x[0].number(), x[1].number(), y[0].number(), y[1].number()};
    for (int index : integers(value.at("trace_indices")))
        if (index >= 0) result.traceIndices.push_back(static_cast<std::size_t>(index));
    result.automatic = value.at("automatic").boolean(false);
    return result;
}

Figure figureFrom(const json::Value& value) {
    Figure result;
    result.kind = value.at("kind").string();
    result.valueScale = value.at("value_scale").string("auto");
    result.observableIdentity = value.at("observable_identity").string();
    result.title = value.at("title").string();
    result.xLabel = value.at("x_label").string();
    result.yLabel = value.at("y_label").string();
    result.zLabel = value.at("z_label").string("z");
    result.square = value.at("square").boolean();
    result.titleColor = value.at("title_color").string();
    result.dashboardPage = static_cast<int>(value.at("dashboard_page").number(-1.0));
    result.gridSpan = static_cast<int>(value.at("grid_span").number(1.0));
    result.rowSpan = static_cast<int>(value.at("row_span").number(1.0));
    result.inset = insetFrom(value.at("detail_inset"));
    result.coordinates = value.at("coordinates").string();
    result.colorLabel = value.at("color_label").string();
    const auto& alternate = value.at("alternate_coordinates");
    if (!alternate.isNull()) {
        Figure::CoordinateAxes axes;
        axes.name = alternate.at("name").string();
        axes.xLabel = alternate.at("x_label").string();
        axes.yLabel = alternate.at("y_label").string();
        axes.zLabel = alternate.at("z_label").string("z");
        axes.xScale = alternate.at("x_scale").number(1);
        axes.yScale = alternate.at("y_scale").number(1);
        axes.zScale = alternate.at("z_scale").number(1);
        axes.timeScale = alternate.at("time_scale").number(1);
        axes.traceLabels = strings(alternate.at("trace_labels"));
        axes.colorLabel = alternate.at("color_label").string();
        result.alternateCoordinates = std::move(axes);
    }
    for (const auto& ref : value.at("references").array())
        result.references.push_back({ref.at("x").number(),ref.at("y").number(),
            ref.at("label").string(),ref.at("equilibrium").boolean()});
    for (const auto& trace : value.at("traces").array()) result.traces.push_back(traceFrom(trace));
    return result;
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char item) {
        return static_cast<char>(std::tolower(item));
    });
    return value;
}

std::string systemToken(std::string value) {
    value = lower(std::move(value));
    std::string result;
    for (const unsigned char item : value) {
        if (std::isalnum(item)) result.push_back(static_cast<char>(item));
    }
    return result;
}

}  // namespace

std::optional<DetailInset> automaticDetailInset(
    const Figure& figure, double disparityThreshold) {
    const std::string kind = lower(figure.kind);
    if (kind != "time" && kind != "observable_time" && kind != "klqe")
        return std::nullopt;
    struct Record {
        std::size_t index = 0;
        std::vector<double> x;
        std::vector<double> y;
        double lower = 0.0;
        double upper = 0.0;
        double robustSpan = 0.0;
        double fullSpan = 0.0;
    };
    auto quantile = [](std::vector<double> values, double probability) {
        std::sort(values.begin(), values.end());
        const double position = probability * static_cast<double>(values.size() - 1);
        const auto below = static_cast<std::size_t>(std::floor(position));
        const auto above = static_cast<std::size_t>(std::ceil(position));
        const double fraction = position - static_cast<double>(below);
        return values[below] * (1.0 - fraction) + values[above] * fraction;
    };
    std::vector<Record> records;
    for (std::size_t traceIndex = 0; traceIndex < figure.traces.size(); ++traceIndex) {
        const Trace& trace = figure.traces[traceIndex];
        if (trace.marker || trace.heatmap) continue;
        Record record; record.index = traceIndex;
        const std::size_t count = std::min(trace.x.size(), trace.y.size());
        for (std::size_t index = 0; index < count; ++index) {
            if (!std::isfinite(trace.x[index]) || !std::isfinite(trace.y[index])) continue;
            record.x.push_back(trace.x[index]); record.y.push_back(trace.y[index]);
        }
        if (record.x.size() < 20) continue;
        record.lower = quantile(record.y, 0.01);
        record.upper = quantile(record.y, 0.99);
        record.robustSpan = record.upper - record.lower;
        const auto extremes = std::minmax_element(record.y.begin(), record.y.end());
        record.fullSpan = *extremes.second - *extremes.first;
        records.push_back(std::move(record));
    }
    if (records.empty()) return std::nullopt;
    const double epsilon = std::numeric_limits<double>::epsilon();
    auto focus = records.end();
    double bestRatio = disparityThreshold;
    for (auto item = records.begin(); item != records.end(); ++item) {
        if (item->robustSpan <= epsilon) continue;
        const double ratio = item->fullSpan / item->robustSpan;
        if (ratio >= bestRatio) { bestRatio = ratio; focus = item; }
    }
    DetailInset result; result.automatic = true;
    if (focus != records.end()) {
        result.yMin = records.front().lower; result.yMax = records.front().upper;
        for (const auto& record : records) {
            result.traceIndices.push_back(record.index);
            result.yMin = std::min(result.yMin, record.lower);
            result.yMax = std::max(result.yMax, record.upper);
        }
        bool found = false;
        for (std::size_t index = 0; index < focus->x.size(); ++index) {
            if (focus->y[index] < focus->lower || focus->y[index] > focus->upper) continue;
            if (!found) { result.xMin = result.xMax = focus->x[index]; found = true; }
            else { result.xMin = std::min(result.xMin, focus->x[index]);
                   result.xMax = std::max(result.xMax, focus->x[index]); }
        }
        if (!found) {
            const auto x = std::minmax_element(focus->x.begin(), focus->x.end());
            result.xMin = *x.first; result.xMax = *x.second;
        }
    } else {
        double largest = 0.0, smallest = std::numeric_limits<double>::infinity();
        for (const auto& record : records) if (record.robustSpan > epsilon) {
            largest = std::max(largest, record.robustSpan);
            smallest = std::min(smallest, record.robustSpan);
        }
        if (!std::isfinite(smallest) || largest / smallest < disparityThreshold)
            return std::nullopt;
        bool found = false;
        for (const auto& record : records) {
            if (record.robustSpan <= epsilon || record.robustSpan > largest / disparityThreshold) continue;
            result.traceIndices.push_back(record.index);
            const auto x = std::minmax_element(record.x.begin(), record.x.end());
            if (!found) {
                result.xMin = *x.first; result.xMax = *x.second;
                result.yMin = record.lower; result.yMax = record.upper; found = true;
            } else {
                result.xMin = std::min(result.xMin, *x.first);
                result.xMax = std::max(result.xMax, *x.second);
                result.yMin = std::min(result.yMin, record.lower);
                result.yMax = std::max(result.yMax, record.upper);
            }
        }
        if (!found) return std::nullopt;
    }
    if (!(result.xMax > result.xMin)) return std::nullopt;
    double span = result.yMax - result.yMin;
    if (span <= epsilon) span = std::max(std::abs(result.yMin), 1.0) * 0.1;
    result.yMin -= 0.08 * span; result.yMax += 0.08 * span;
    return result;
}

Catalog Catalog::load(const std::wstring& path) {
    return Catalog::fromJson(json::dump(json::load(path)));
}

int Catalog::indexOf(const std::string& id) const {
    const auto found = std::find_if(items.begin(), items.end(), [&](const CatalogItem& item) {
        return item.id == id;
    });
    return found == items.end() ? -1 : static_cast<int>(found - items.begin());
}

Payload Payload::load(const std::wstring& path) {
    const json::Value root = json::load(path);
    if (root.at("schema").string() != "dynsyslab.viewer.payload/v1")
        throw std::runtime_error("Unsupported Viewer payload schema");
    Payload result;
    result.source = root.at("source").string();
    result.view = root.at("view").string();
    result.tileName = root.at("tile").at("name").string();
    result.columns = static_cast<int>(root.at("tile").at("columns").number(1.0));
    result.pageLabels = strings(root.at("page_labels"));
    result.warnings = strings(root.at("warnings"));
    result.layout = root.at("layout").string();
    result.complexDisplay = root.at("complex_display").string();
    for (const auto v : integers(root.at("selected_modes"))) if (v >= 0) result.selectedModes.push_back(v);
    for (const auto& figure : root.at("figures").array()) result.figures.push_back(figureFrom(figure));
    return result;
}

std::wstring wide(const std::string& value) {
    return dynsyslab::core::wide(value);
}

std::string narrow(const std::wstring& value) {
    return dynsyslab::core::narrow(value);
}

std::string canonicalSystemLabel(const std::string& value) {
    if (value.empty()) return {};
    const auto token = systemToken(value);
    static const auto systems = dynsyslab::core::availableSystems();
    for (const auto& descriptor : systems) {
        if (token == systemToken(descriptor.identifier) ||
            token == systemToken(descriptor.label)) return descriptor.label;
    }
    if (token == "rpe" || token == "rayleighplesset") return "Rayleigh-Plesset";
    if (token == "reddyszeri") return "Reddy-Szeri";
    return value;
}

}  // namespace dsl

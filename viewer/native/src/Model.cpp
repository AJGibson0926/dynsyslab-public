#include "dynsyslab_viewer/Model.hpp"
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/System.hpp>

#include <windows.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

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
    result.boxAxes3D = value.at("box_axes").boolean(false);
    if (const auto& camera = value.at("camera"); !camera.isNull()) {
        result.cameraAzimuthDegrees = camera.at("azimuth_degrees").number(42.0);
        result.cameraElevationDegrees = camera.at("elevation_degrees").number(24.0);
    }
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

struct EpisodeRecord {
    std::size_t index = 0;
    std::vector<double> x;
    std::vector<double> y;
    double lower = 0.0;
    double upper = 0.0;
    double robustSpan = 0.0;
    double fullSpan = 0.0;
};

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

// The transient episode of a time figure as [start, end], if there is one. Off-scale
// samples (beyond the robust band by more than one band width) must all fall inside it, and
// it may cover at most 30% of the record.
std::optional<std::pair<double, double>> transientEpisode(
    const std::vector<EpisodeRecord>& records, double lower, double upper) {
    const double band = upper - lower;
    if (!(band > 0)) return std::nullopt;
    const double low = lower - band, high = upper + band;
    double allMin = std::numeric_limits<double>::infinity(), allMax = -allMin;
    double firstBurst = allMin, lastBurst = -allMin;
    for (const auto& record : records)
        for (std::size_t index = 0; index < record.x.size(); ++index) {
            allMin = std::min(allMin, record.x[index]); allMax = std::max(allMax, record.x[index]);
            if (record.y[index] < low || record.y[index] > high) {
                firstBurst = std::min(firstBurst, record.x[index]);
                lastBurst = std::max(lastBurst, record.x[index]);
            }
        }
    const double total = allMax - allMin;
    if (!(total > 0) || !std::isfinite(firstBurst)) return std::nullopt;
    // Disagreement between the traces, sampled on a common grid.
    constexpr std::size_t grid = 800;
    std::vector<double> spread(grid, 0.0);
    std::vector<bool> offscale(grid, false);
    for (std::size_t k = 0; k < grid; ++k) {
        const double x = allMin + total * static_cast<double>(k) / (grid - 1);
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        int seen = 0;
        for (const auto& record : records) {
            if (record.x.empty() || x < record.x.front() || x > record.x.back()) continue;
            const auto it = std::lower_bound(record.x.begin(), record.x.end(), x);
            if (it == record.x.end()) continue;
            const std::size_t j = static_cast<std::size_t>(it - record.x.begin());
            double y = record.y[j];
            if (j > 0 && record.x[j] > record.x[j - 1]) {
                const double w = (x - record.x[j - 1]) / (record.x[j] - record.x[j - 1]);
                y = record.y[j - 1] * (1 - w) + record.y[j] * w;
            }
            lo = std::min(lo, y); hi = std::max(hi, y); ++seen;
        }
        if (seen >= 2) spread[k] = hi - lo;
        const double cell = total / (grid - 1);
        offscale[k] = x >= firstBurst - cell && x <= lastBurst + cell;
    }
    std::vector<double> sorted = spread;
    std::nth_element(sorted.begin(), sorted.begin() + grid / 2, sorted.end());
    const double threshold = std::max(0.10 * band, 3.0 * sorted[grid / 2]);
    std::size_t start = grid, end = 0, lull = 0;
    const std::size_t tolerance = grid * 3 / 100;
    for (std::size_t k = 0; k < grid; ++k) {
        const bool active = offscale[k] || spread[k] > threshold;
        if (start == grid) { if (active) start = end = k; continue; }
        if (active) { end = k; lull = 0; }
        else if (++lull > tolerance) break;
    }
    if (start == grid) return std::nullopt;
    const double pad = 0.01 * total;
    const double x0 = std::max(allMin, allMin + total * start / (grid - 1) - pad);
    const double x1 = std::min(allMax, allMin + total * end / (grid - 1) + pad);
    if (lastBurst > x1 || firstBurst < x0 || x1 - x0 > 0.30 * total) return std::nullopt;
    return std::make_pair(x0, x1);
}

}  // namespace

std::optional<DetailInset> automaticDetailInset(
    const Figure& figure, double disparityThreshold) {
    const std::string kind = lower(figure.kind);
    if (kind != "time" && kind != "observable_time" && kind != "klqe")
        return std::nullopt;
    using Record = EpisodeRecord;
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
        // Transient episode: if the off-scale samples sit inside one short episode (a start-up
        // spike and the settling that follows it), the inset shows that whole episode and the
        // main panel keeps the robust range, so both panels carry visible information. The
        // episode runs from the first off-scale sample for as long as the traces disagree
        // with one another markedly (for example a controller ringing as it attaches to its
        // target), tolerating short lulls.
        if (const auto episode = transientEpisode(records, result.yMin, result.yMax)) {
            const double x0 = episode->first, x1 = episode->second;
            const double band = result.yMax - result.yMin;
            double y0 = std::numeric_limits<double>::infinity(), y1 = -y0;
            for (const auto& record : records)
                for (std::size_t index = 0; index < record.x.size(); ++index)
                    if (record.x[index] >= x0 && record.x[index] <= x1) {
                        y0 = std::min(y0, record.y[index]); y1 = std::max(y1, record.y[index]);
                    }
            const double settled = std::max(std::abs(result.yMin), std::abs(result.yMax));
            if (x1 > x0 && y1 > y0) {
                result.transientWindow = true;
                result.mainYMin = result.yMin - 0.08 * band;
                result.mainYMax = result.yMax + 0.08 * band;
                result.xMin = x0; result.xMax = x1;
                // Spikes spanning decades beyond the settled level read best on a symmetric
                // log ordinate whose linear core is the settled range.
                // The linear core holds the episode's on-scale behaviour (e.g. ringing); the
                // off-scale spike goes to compressed decades beyond it.
                double core = 0.0;
                for (const auto& record : records)
                    for (std::size_t index = 0; index < record.x.size(); ++index)
                        if (record.x[index] >= x0 && record.x[index] <= x1
                            && record.y[index] >= result.yMin - band && record.y[index] <= result.yMax + band)
                            core = std::max(core, std::abs(record.y[index]));
                if (!(core > epsilon)) core = settled;
                if (core > epsilon && std::max(std::abs(y0), std::abs(y1)) > 30.0 * core) {
                    constexpr double weight = 0.35;
                    result.symmetricLog = true; result.linearThreshold = core; result.decadeWeight = weight;
                    // Centre the core vertically: limits symmetric about zero in the
                    // transformed coordinate, with a little headroom.
                    auto forward = [core](double v) {
                        const double q = std::abs(v) / core;
                        return std::copysign(q <= 1 ? q : 1 + weight * std::log10(q), v);
                    };
                    auto inverse = [core](double u) {
                        const double q = std::abs(u);
                        return std::copysign(core * (q <= 1 ? q : std::pow(10., std::min((q - 1) / weight, 307.))), u);
                    };
                    const double reach = 1.06 * std::max(std::abs(forward(y0)), std::abs(forward(y1)));
                    result.yMin = inverse(-reach); result.yMax = inverse(reach);
                } else {
                    const double span = y1 - y0;
                    result.yMin = y0 - 0.06 * span; result.yMax = y1 + 0.06 * span;
                }
                return result;
            }
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

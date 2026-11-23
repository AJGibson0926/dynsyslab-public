// Post-processing method "control_score": scores controlled runs against a reference signal and
// optionally selects one by a declared rule.
//
// Settings: {"reference":PATH (a reference_signal archive),
//            "windows":{"startup":[0,10],"steady":[100,null]},
//            "selection":{"metric":"steady.radius_rmse","direction":"minimize",
//                         "tie_break":"source_order","require_complete":true}}
// Per window and run: tracking RMSE of state component 0, relative drive error, drive
// correlation, and peak drive. Arrays: one matrix per metric (runs x windows), run_end_time,
// run_completed, selected_index (-1 when nothing qualifies). metadata names the selected source.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;

struct Window { std::string name; double start, stop; };

double linear(const std::vector<double>& xp, const std::vector<double>& fp, double x) {
    if (x <= xp.front()) return fp.front();
    if (x >= xp.back()) return fp.back();
    const auto j = static_cast<std::size_t>(std::upper_bound(xp.begin(), xp.end(), x) - xp.begin()) - 1;
    return fp[j] + (fp[j + 1] - fp[j]) * (x - xp[j]) / (xp[j + 1] - xp[j]);
}

std::vector<double> firstColumn(const c::NpyArray& array) {
    const auto values = array.doubles();
    if (array.shape.size() < 2) return values;
    std::vector<double> out(array.shape[0]);
    for (std::size_t r = 0; r < out.size(); ++r) out[r] = values[r * array.shape[1]];
    return out;
}

}  // namespace

std::filesystem::path runControlScore(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                      const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("control_score needs controlled runs");
    const auto referencePath = settings.at("reference").string();
    if (referencePath.empty()) throw std::invalid_argument("control_score needs a reference_signal archive");
    const auto reference = c::readNpz(std::filesystem::u8path(referencePath));
    const auto referenceTime = reference.at("time").doubles();
    const auto referenceState = firstColumn(reference.at("target_values"));
    const auto referenceDrive = firstColumn(reference.at("inverse_drive"));
    std::vector<Window> windows;
    for (const auto& [name, range] : settings.at("windows").object()) {
        const auto& bounds = range.array();
        if (bounds.size() != 2) throw std::invalid_argument("control_score windows are [start, stop]");
        windows.push_back({name, bounds[0].number(-HUGE_VAL), bounds[1].isNull() ? HUGE_VAL : bounds[1].number()});
    }
    if (windows.empty()) windows.push_back({"all", -HUGE_VAL, HUGE_VAL});
    const std::vector<std::string> metrics{"radius_rmse", "drive_relative_error", "drive_correlation", "peak_drive"};
    const auto runs = sources.size(), columns = windows.size();
    std::vector<std::vector<double>> table(metrics.size(), std::vector<double>(runs * columns, std::numeric_limits<double>::quiet_NaN()));
    std::vector<double> endTime(runs);
    std::vector<std::int64_t> completed(runs);
    for (std::size_t run = 0; run < runs; ++run) {
        const auto record = c::readNpz(sources[run]);
        const auto time = record.at("time").doubles();
        const auto state = firstColumn(record.at("values"));
        const auto drive = firstColumn(record.at("control_values"));
        endTime[run] = time.empty() ? std::numeric_limits<double>::quiet_NaN() : time.back();
        const auto metadata = json::parse(record.metadataJson());
        completed[run] = metadata.at("integration_terminated").isNull() || !metadata.at("integration_terminated").boolean(false) ? 1 : 0;
        for (std::size_t w = 0; w < columns; ++w) {
            double squares = 0.0, driveError = 0.0, driveNorm = 0.0, peak = 0.0, su = 0.0, sr = 0.0, suu = 0.0, srr = 0.0, sur = 0.0;
            std::size_t count = 0;
            for (std::size_t k = 0; k < time.size(); ++k) {
                if (time[k] < windows[w].start || time[k] > windows[w].stop) continue;
                const double xr = linear(referenceTime, referenceState, time[k]), ur = linear(referenceTime, referenceDrive, time[k]);
                squares += (state[k] - xr) * (state[k] - xr);
                driveError += (drive[k] - ur) * (drive[k] - ur); driveNorm += ur * ur;
                peak = std::max(peak, std::abs(drive[k]));
                su += drive[k]; sr += ur; suu += drive[k] * drive[k]; srr += ur * ur; sur += drive[k] * ur; ++count;
            }
            if (!count) continue;
            const double n = static_cast<double>(count);
            table[0][run * columns + w] = std::sqrt(squares / n);
            table[1][run * columns + w] = driveNorm > 0.0 ? std::sqrt(driveError / driveNorm) : std::numeric_limits<double>::quiet_NaN();
            const double denominator = std::sqrt((n * suu - su * su) * (n * srr - sr * sr));
            table[2][run * columns + w] = denominator > 0.0 ? (n * sur - su * sr) / denominator : std::numeric_limits<double>::quiet_NaN();
            table[3][run * columns + w] = peak;
        }
    }
    std::int64_t selected = -1;
    json::Value::Object meta{{"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("control-score")},
                             {"label", json::Value(label)}, {"score_settings", settings}};
    const auto& selection = settings.at("selection");
    if (!selection.isNull()) {
        const auto rule = selection.at("metric").string();
        const auto dot = rule.find('.');
        if (dot == std::string::npos) throw std::invalid_argument("selection.metric is WINDOW.METRIC");
        const auto windowName = rule.substr(0, dot), metricName = rule.substr(dot + 1);
        const auto w = std::find_if(windows.begin(), windows.end(), [&](const Window& item) { return item.name == windowName; });
        const auto m = std::find(metrics.begin(), metrics.end(), metricName);
        if (w == windows.end() || m == metrics.end()) throw std::invalid_argument("selection.metric names an unknown window or metric");
        const auto column = static_cast<std::size_t>(w - windows.begin());
        const auto& values = table[static_cast<std::size_t>(m - metrics.begin())];
        const bool maximize = selection.at("direction").string("minimize") == "maximize";
        const bool requireComplete = selection.at("require_complete").boolean(true);
        if (selection.at("tie_break").string("source_order") != "source_order")
            throw std::invalid_argument("selection.tie_break supports source_order");
        for (std::size_t run = 0; run < runs; ++run) {
            const double value = values[run * columns + column];
            if (!std::isfinite(value) || (requireComplete && !completed[run])) continue;
            if (selected < 0) { selected = static_cast<std::int64_t>(run); continue; }
            const double best = values[static_cast<std::size_t>(selected) * columns + column];
            if (maximize ? value > best : value < best) selected = static_cast<std::int64_t>(run);
        }
        if (selected >= 0) meta["selected_source"] = json::Value(sources[static_cast<std::size_t>(selected)].generic_u8string());
    }
    c::NpzRecord record;
    for (std::size_t k = 0; k < metrics.size(); ++k) record.arrays[metrics[k]] = c::NpyArray::fromDoubles({runs, columns}, table[k]);
    record.arrays["run_end_time"] = c::NpyArray::fromDoubles({runs}, endTime);
    record.arrays["run_completed"] = c::NpyArray::fromInt64s({runs}, completed);
    record.arrays["selected_index"] = c::NpyArray::fromInt64s({1}, {selected});
    json::Value::Array windowNames;
    for (const auto& window : windows) windowNames.emplace_back(window.name);
    meta["windows"] = json::Value(std::move(windowNames));
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = root / L"derived" / L"control-score";
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    c::writeNpzAtomicNew(path, record);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail

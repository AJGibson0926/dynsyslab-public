// Post-processing method "control_score": scores controlled runs against a reference signal and
// optionally selects one by a declared rule.
//
// Settings: {"reference":PATH (a reference_signal archive; optional with "estimation"),
//            "windows":{"startup":[0,10],"steady":[100,null]},
//            "selection":{"metric":"steady.radius_rmse","direction":"minimize",
//                         "tie_break":"source_order","require_complete":true}}
// Per window and run: tracking RMSE of state component 0, relative drive error, drive
// correlation, and peak drive. Optional metrics (core/SignalMetrics.hpp):
//   "tone": {"frequency_hz": f, "harmonics": H}  tone_amplitude, tone_phase, tone_residual_rms,
//                                                harmonic_content of state component 0;
//   "settling_band": b                           settling_time of the tracking error;
//   "drive_lowpass_hz": fc                       lowpass_drive_relative_error (both drives filtered);
//   runs with a state estimator                  estimate_rmse (all state components) and
//                                                innovation_lag1_autocorrelation (whiteness);
//   "compare_to": run index                      run_difference_rmse (state component 0);
//   "estimation": {...}                          estimator and shadow-controller scores
//                                                (ControlScoreEstimation.hpp);
//   "run_values": {"NAME": [one per run]}         stored as array NAME (an abscissa for plots).
// Arrays: one matrix per metric (runs x windows), run_end_time, run_completed, selected_index
// (-1 when nothing qualifies). metadata names the selected source.
#include "PostprocessInternal.hpp"
#include "ControlScoreEstimation.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SignalMetrics.hpp>
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
    // The tracking metrics need a reference_signal archive; estimation scores alone do not.
    const auto referencePath = settings.at("reference").string();
    const bool hasReference = !referencePath.empty();
    if (!hasReference && settings.at("estimation").isNull())
        throw std::invalid_argument("control_score needs a reference_signal archive or estimation scores");
    std::vector<double> referenceTime, referenceState, referenceDrive;
    if (hasReference) {
        const auto reference = c::readNpz(std::filesystem::u8path(referencePath));
        referenceTime = reference.at("time").doubles();
        referenceState = firstColumn(reference.at("target_values"));
        referenceDrive = firstColumn(reference.at("inverse_drive"));
    }
    std::vector<Window> windows;
    for (const auto& [name, range] : settings.at("windows").object()) {
        const auto& bounds = range.array();
        if (bounds.size() != 2) throw std::invalid_argument("control_score windows are [start, stop]");
        windows.push_back({name, bounds[0].number(-HUGE_VAL), bounds[1].isNull() ? HUGE_VAL : bounds[1].number()});
    }
    if (windows.empty()) windows.push_back({"all", -HUGE_VAL, HUGE_VAL});
    std::vector<std::string> metrics{"radius_rmse", "drive_relative_error", "drive_correlation", "peak_drive"};
    const auto& tone = settings.at("tone");
    const double settlingBand = settings.at("settling_band").number(0.0);
    const double lowpassHz = settings.at("drive_lowpass_hz").number(0.0);
    const double compareTo = settings.at("compare_to").number(-1.0);
    if (compareTo >= double(sources.size())) throw std::invalid_argument("control_score compare_to is not a run index");
    if (!tone.isNull()) for (const char* name : {"tone_amplitude", "tone_phase", "tone_residual_rms", "harmonic_content"}) metrics.emplace_back(name);
    if (settlingBand > 0.0) metrics.emplace_back("settling_time");
    if (lowpassHz > 0.0) metrics.emplace_back("lowpass_drive_relative_error");
    metrics.emplace_back("estimate_rmse");
    metrics.emplace_back("innovation_lag1_autocorrelation");
    if (compareTo >= 0.0) metrics.emplace_back("run_difference_rmse");
    const auto& estimation = settings.at("estimation");  // ControlScoreEstimation.hpp
    for (const auto& name : estimationMetricNames(estimation)) metrics.push_back(name);
    const auto metricColumn = [&metrics](const char* name) {
        return static_cast<std::size_t>(std::find(metrics.begin(), metrics.end(), name) - metrics.begin()); };
    std::vector<double> compareTime, compareState;
    if (compareTo >= 0.0) {
        const auto other = c::readNpz(sources[static_cast<std::size_t>(compareTo)]);
        compareTime = other.at("time").doubles(); compareState = firstColumn(other.at("values"));
    }
    const auto runs = sources.size(), columns = windows.size();
    std::vector<std::vector<double>> table(metrics.size(), std::vector<double>(runs * columns, std::numeric_limits<double>::quiet_NaN()));
    std::vector<double> endTime(runs);
    std::vector<std::int64_t> completed(runs);
    for (std::size_t run = 0; run < runs; ++run) {
        const auto record = c::readNpz(sources[run]);
        const auto time = record.at("time").doubles();
        const auto state = firstColumn(record.at("values"));
        const auto drive = hasReference ? firstColumn(record.at("control_values")) : std::vector<double>{};
        endTime[run] = time.empty() ? std::numeric_limits<double>::quiet_NaN() : time.back();
        const auto metadata = json::parse(record.metadataJson());
        completed[run] = metadata.at("integration_terminated").isNull() || !metadata.at("integration_terminated").boolean(false) ? 1 : 0;
        for (std::size_t w = 0; w < columns; ++w) {
            for (const auto& [name, value] : scoreEstimation(record, windows[w].start, windows[w].stop, estimation))
                table[static_cast<std::size_t>(std::find(metrics.begin(), metrics.end(), name) - metrics.begin())][run * columns + w] = value;
            if (!hasReference) continue;
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
            // Window samples for the optional metrics.
            std::vector<double> t, x, e, u, ur;
            for (std::size_t k = 0; k < time.size(); ++k) {
                if (time[k] < windows[w].start || time[k] > windows[w].stop) continue;
                t.push_back(time[k]); x.push_back(state[k]);
                e.push_back(state[k] - linear(referenceTime, referenceState, time[k]));
                u.push_back(drive[k]); ur.push_back(linear(referenceTime, referenceDrive, time[k]));
            }
            const auto cell = [&](const char* name) -> double& { return table[metricColumn(name)][run * columns + w]; };
            if (!tone.isNull() && t.size() > 3) {
                const auto fit = c::fitTone(t, x, tone.at("frequency_hz").number(),
                    static_cast<std::size_t>(tone.at("harmonics").number(1)));
                cell("tone_amplitude") = fit.amplitudes.front(); cell("tone_phase") = fit.phases.front();
                cell("tone_residual_rms") = fit.residualRms; cell("harmonic_content") = fit.harmonicContent;
            }
            if (settlingBand > 0.0) cell("settling_time") = c::settlingTime(t, e, settlingBand);
            if (lowpassHz > 0.0 && !t.empty()) {
                const auto fu = c::lowPass(t, u, lowpassHz), fr = c::lowPass(t, ur, lowpassHz);
                double num = 0.0, den = 0.0;
                for (std::size_t k = 0; k < t.size(); ++k) { num += (fu[k] - fr[k]) * (fu[k] - fr[k]); den += fr[k] * fr[k]; }
                cell("lowpass_drive_relative_error") = den > 0.0 ? std::sqrt(num / den) : std::numeric_limits<double>::quiet_NaN();
            }
            if (record.contains("estimator_time") && record.contains("state_estimate")) {
                const auto et = record.at("estimator_time").doubles();
                const auto estimate = record.at("state_estimate").doubles();
                const auto values = record.at("values").doubles();
                const std::size_t dim = record.at("values").shape.size() > 1 ? record.at("values").shape[1] : 1;
                const std::size_t edim = et.empty() ? 0 : estimate.size() / et.size();
                double estimateSquares = 0.0; std::size_t estimateCount = 0;
                for (std::size_t k = 0; k < et.size() && edim == dim; ++k) {
                    if (et[k] < windows[w].start || et[k] > windows[w].stop) continue;
                    for (std::size_t j = 0; j < dim; ++j) {
                        std::vector<double> component(time.size());
                        for (std::size_t r = 0; r < time.size(); ++r) component[r] = values[r * dim + j];
                        const double truth = linear(time, component, et[k]);
                        estimateSquares += (estimate[k * dim + j] - truth) * (estimate[k * dim + j] - truth); ++estimateCount;
                    }
                }
                if (estimateCount) cell("estimate_rmse") = std::sqrt(estimateSquares / double(estimateCount));
                if (record.contains("innovation")) {
                    const auto innovation = firstColumn(record.at("innovation"));
                    std::vector<double> inWindow;
                    for (std::size_t k = 0; k < et.size() && k < innovation.size(); ++k)
                        if (et[k] >= windows[w].start && et[k] <= windows[w].stop) inWindow.push_back(innovation[k]);
                    cell("innovation_lag1_autocorrelation") = c::autocorrelation(inWindow, 1);
                }
            }
            if (compareTo >= 0.0 && !compareTime.empty() && !t.empty()) {
                double differenceSquares = 0.0;
                for (std::size_t k = 0; k < t.size(); ++k) {
                    const double d = x[k] - linear(compareTime, compareState, t[k]);
                    differenceSquares += d * d;
                }
                cell("run_difference_rmse") = std::sqrt(differenceSquares / double(t.size()));
            }
        }
    }
    std::int64_t selected = -1;
    json::Value::Object meta{{"schema", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("control-score")},
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
    for (const auto& [name, values] : settings.at("run_values").object()) {  // a per-run abscissa for plots
        if (values.array().size() != runs) throw std::invalid_argument("control_score run_values." + name + " needs one value per run");
        std::vector<double> column;
        for (const auto& value : values.array()) column.push_back(value.isNull() ? std::numeric_limits<double>::quiet_NaN() : value.number());
        record.arrays[name] = c::NpyArray::fromDoubles({runs}, column);
    }
    record.arrays["run_completed"] = c::NpyArray::fromInt64s({runs}, completed);
    record.arrays["selected_index"] = c::NpyArray::fromInt64s({1}, {selected});
    json::Value::Array windowNames;
    for (const auto& window : windows) windowNames.emplace_back(window.name);
    meta["windows"] = json::Value(std::move(windowNames));
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = kindFolder(root, dynsyslab::core::object_kind::report);
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(label + "__" + timestamp()) + ".npz");
    stampDataObject(record, dynsyslab::core::object_kind::report, "control-score", {});
    c::writeNpzAtomicNew(path, record);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail

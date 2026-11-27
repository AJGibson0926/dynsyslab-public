// koopman_validation extras (settings):
//   "windows": [[a, b], ...]   anchored relative L2 error of each eigenfunction over each window,
//                              re-anchored at the window's first sample:
//                              ||psi(x(t)) - exp(lambda (t - t_a)) psi(x(t_a))|| / ||psi(x(t))||;
//   "binned_defect": {"axis": i, "bins": k}  the one-step defect |psi(x_{k+1}) - mu psi(x_k)| /
//                              |psi(x_{k+1})| averaged in k equal bins of state component i;
//   "binned_defect": {"axes": [i, j], "bins": [kx, ky] | "x_edges"/"y_edges", "horizon": h}
//                              the 2-D finite-horizon defect over occupied cells: pointwise
//                              |psi(t+h) - exp(lambda h) psi(t)| / max(|psi(t+h)|, sqrt(eps) max|psi|),
//                              averaged per cell with weights |psi(t+h)|, stored as log10 (floored at
//                              1e-12; NaN in empty cells) with the occupied mask;
//   "segments": {"starts": [t0, t1, ...], "horizon": h}  the relative L2 error of
//                              psi(t+h) ~ exp(lambda h) psi(t) for t in each [t_i, t_i+1];
//   "history": true            psi_history (samples x modes) and psi_predicted, the anchored
//                              prediction exp(lambda t) psi_0, for display.
// The anchored test follows Koopman eigenfunction validation on trajectories (Mezic 2005,
// Nonlinear Dyn. 41:309-325).
#include "KromValidationWindows.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
using Complex = std::complex<double>;

// Whole samples in the horizon h (time units) on the uniform record clock.
std::size_t horizonShift(const json::Value& settings, const std::vector<double>& time, std::size_t samples) {
    if (samples < 2) throw std::invalid_argument("finite-horizon defects need at least two samples");
    const double dt = time[1] - time[0], h = settings.at("horizon").number(dt);
    const double shift = std::round(h / dt);
    if (!(dt > 0) || !(shift >= 1) || shift >= double(samples))
        throw std::invalid_argument("horizon must be a positive whole number of samples shorter than the record");
    return static_cast<std::size_t>(shift);
}

void appendSegmentErrors(dynsyslab::core::NpzRecord& archive, const json::Value& settings,
                         const std::vector<Complex>& rates, const std::vector<Complex>& psi,
                         const std::vector<double>& time, std::size_t samples, std::size_t modes) {
    const auto& starts = settings.at("starts").array();
    if (starts.size() < 2) throw std::invalid_argument("segments need at least two start times");
    const auto shift = horizonShift(settings, time, samples);
    const double h = time[shift] - time[0];
    const auto count = starts.size() - 1;
    std::vector<double> errors(modes * count, std::numeric_limits<double>::quiet_NaN()), bounds;
    for (const auto& t : starts) bounds.push_back(t.number());
    for (std::size_t j = 0; j < modes; ++j) {
        const auto growth = std::exp(rates[j] * h);
        for (std::size_t g = 0; g < count; ++g) {
            double squares = 0.0, norm = 0.0;
            for (std::size_t k = 0; k + shift < samples; ++k) {
                if (time[k] < bounds[g] || time[k] > bounds[g + 1]) continue;
                const auto actual = psi[(k + shift) * modes + j];
                squares += std::norm(actual - growth * psi[k * modes + j]);
                norm += std::norm(actual);
            }
            if (norm > 0.0) errors[j * count + g] = std::sqrt(squares / norm);
        }
    }
    archive.arrays["segment_horizon_relative_l2"] = dynsyslab::core::NpyArray::fromDoubles({modes, count}, errors);
    archive.arrays["segment_starts"] = dynsyslab::core::NpyArray::fromDoubles({bounds.size()}, bounds);
}

std::vector<double> cellEdges(const json::Value& settings, const char* key, std::size_t bins,
                              const std::vector<double>& values, std::size_t dimension, std::size_t axis) {
    std::vector<double> edges;
    for (const auto& e : settings.at(key).array()) edges.push_back(e.number());
    if (!edges.empty()) {
        if (edges.size() < 2 || !std::is_sorted(edges.begin(), edges.end()))
            throw std::invalid_argument("binned_defect edges must be increasing");
        return edges;
    }
    double low = HUGE_VAL, high = -HUGE_VAL;
    for (std::size_t k = 0; k * dimension < values.size(); ++k) {
        low = std::min(low, values[k * dimension + axis]);
        high = std::max(high, values[k * dimension + axis]);
    }
    for (std::size_t e = 0; e <= bins; ++e) edges.push_back(low + (high - low) * double(e) / double(bins));
    return edges;
}

std::size_t cellOf(double x, const std::vector<double>& edges) {
    const auto it = std::upper_bound(edges.begin(), edges.end(), x);
    const auto index = it == edges.begin() ? 0 : std::size_t(it - edges.begin()) - 1;
    return std::min(index, edges.size() - 2);
}

void appendCellDefect(dynsyslab::core::NpzRecord& archive, const json::Value& settings,
                      const std::vector<Complex>& rates, const std::vector<Complex>& psi,
                      const std::vector<double>& time, const std::vector<double>& values,
                      std::size_t samples, std::size_t dimension, std::size_t modes) {
    const auto& axes = settings.at("axes").array();
    const auto& bins = settings.at("bins");
    if (axes.size() != 2) throw std::invalid_argument("binned_defect axes are [i, j]");
    std::size_t axis[2], count[2];
    for (int a = 0; a < 2; ++a) {
        const double v = axes[a].number(-1), b = bins.isNull() ? 10.0 : bins.array().at(a).number(10);
        if (v < 0 || v >= double(dimension) || !(b >= 1))
            throw std::invalid_argument("binned_defect needs state axes and positive bin counts");
        axis[a] = std::size_t(v);
        count[a] = std::size_t(b);
    }
    const auto xEdges = cellEdges(settings, "x_edges", count[0], values, dimension, axis[0]);
    const auto yEdges = cellEdges(settings, "y_edges", count[1], values, dimension, axis[1]);
    const auto columns = xEdges.size() - 1, rows = yEdges.size() - 1;
    const auto shift = horizonShift(settings, time, samples);
    const double h = time[shift] - time[0], nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> field(modes * rows * columns, nan), occupied(rows * columns, 0.0);
    std::vector<std::size_t> cell(samples - shift);
    for (std::size_t k = 0; k + shift < samples; ++k)
        cell[k] = cellOf(values[k * dimension + axis[1]], yEdges) * columns + cellOf(values[k * dimension + axis[0]], xEdges);
    for (std::size_t j = 0; j < modes; ++j) {
        double scale = std::numeric_limits<double>::min();
        for (std::size_t k = shift; k < samples; ++k) scale = std::max(scale, std::abs(psi[k * modes + j]));
        const double floor = std::sqrt(std::numeric_limits<double>::epsilon()) * scale;
        const auto growth = std::exp(rates[j] * h);
        std::vector<double> numerator(rows * columns, 0.0), denominator(rows * columns, 0.0);
        for (std::size_t k = 0; k + shift < samples; ++k) {
            const auto actual = psi[(k + shift) * modes + j];
            const double weight = std::abs(actual);
            numerator[cell[k]] += weight * std::abs(actual - growth * psi[k * modes + j]) / std::max(weight, floor);
            denominator[cell[k]] += weight;
        }
        for (std::size_t c = 0; c < rows * columns; ++c)
            if (denominator[c] > 0.0) {
                field[j * rows * columns + c] = std::log10(std::max(numerator[c] / denominator[c], 1e-12));
                occupied[c] = 1.0;
            }
    }
    archive.arrays["cell_log10_defect"] = dynsyslab::core::NpyArray::fromDoubles({modes, rows, columns}, field);
    archive.arrays["cell_defect_occupied"] = dynsyslab::core::NpyArray::fromDoubles({rows, columns}, occupied);
    archive.arrays["cell_defect_x_edges"] = dynsyslab::core::NpyArray::fromDoubles({xEdges.size()}, xEdges);
    archive.arrays["cell_defect_y_edges"] = dynsyslab::core::NpyArray::fromDoubles({yEdges.size()}, yEdges);
}
}  // namespace

void appendValidationWindows(dynsyslab::core::NpzRecord& archive, json::Value& meta,
                             const dynsyslab::core::KoopmanModel& model, const std::vector<double>& time,
                             const std::vector<double>& values, std::size_t samples, std::size_t dimension,
                             const json::Value& settings) {
    const auto& windows = settings.at("windows").array();
    const auto& binned = settings.at("binned_defect");
    const auto& segments = settings.at("segments");
    const bool history = settings.at("history").boolean(false);
    if (windows.empty() && binned.isNull() && segments.isNull() && !history) return;
    if (!model.spectrum().hasContinuousRates())
        throw std::invalid_argument("validation windows need continuous eigenvalues");
    const auto rates = model.spectrum().continuousRates();
    const std::size_t modes = model.liftedDimension();
    std::vector<std::complex<double>> psi(samples * modes);
    for (std::size_t k = 0; k < samples; ++k) {
        const auto value = model.evaluate(values.data() + k * dimension);
        std::copy(value.begin(), value.end(), psi.begin() + k * modes);
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!windows.empty()) {
        std::vector<double> errors(modes * windows.size(), nan);
        json::Value::Array bounds;
        for (std::size_t w = 0; w < windows.size(); ++w) {
            const auto& range = windows[w].array();
            if (range.size() != 2) throw std::invalid_argument("validation windows are [start, stop]");
            const double a = range[0].number(), b = range[1].isNull() ? HUGE_VAL : range[1].number();
            bounds.push_back(windows[w]);
            std::size_t anchor = samples;
            for (std::size_t k = 0; k < samples; ++k) if (time[k] >= a) { anchor = k; break; }
            if (anchor == samples) continue;
            for (std::size_t j = 0; j < modes; ++j) {
                double squares = 0.0, norm = 0.0;
                for (std::size_t k = anchor; k < samples && time[k] <= b; ++k) {
                    const auto predicted = std::exp(rates[j] * (time[k] - time[anchor])) * psi[anchor * modes + j];
                    squares += std::norm(psi[k * modes + j] - predicted);
                    norm += std::norm(psi[k * modes + j]);
                }
                if (norm > 0.0) errors[j * windows.size() + w] = std::sqrt(squares / norm);
            }
        }
        archive.arrays["window_anchored_relative_l2"] =
            dynsyslab::core::NpyArray::fromDoubles({modes, windows.size()}, errors);
        meta.object()["validation_windows"] = json::Value(std::move(bounds));
    }
    if (!segments.isNull()) appendSegmentErrors(archive, segments, rates, psi, time, samples, modes);
    if (!binned.isNull() && !binned.at("axes").isNull())
        appendCellDefect(archive, binned, rates, psi, time, values, samples, dimension, modes);
    else if (!binned.isNull()) {
        const double axisValue = binned.at("axis").number(0), binsValue = binned.at("bins").number(10);
        if (axisValue < 0 || axisValue >= double(dimension) || !(binsValue >= 1))
            throw std::invalid_argument("binned_defect needs a state axis and a positive bin count");
        const auto axis = static_cast<std::size_t>(axisValue), bins = static_cast<std::size_t>(binsValue);
        double low = HUGE_VAL, high = -HUGE_VAL;
        for (std::size_t k = 0; k < samples; ++k) { low = std::min(low, values[k * dimension + axis]); high = std::max(high, values[k * dimension + axis]); }
        std::vector<double> sums(modes * bins, 0.0), counts(bins, 0.0), edges(bins + 1);
        for (std::size_t e = 0; e <= bins; ++e) edges[e] = low + (high - low) * double(e) / double(bins);
        for (std::size_t k = 0; k + 1 < samples; ++k) {
            const double x = values[k * dimension + axis];
            const auto bin = std::min(bins - 1, static_cast<std::size_t>(high > low ? (x - low) / (high - low) * double(bins) : 0));
            counts[bin] += 1.0;
            const double dt = time[k + 1] - time[k];
            for (std::size_t j = 0; j < modes; ++j) {
                const auto next = psi[(k + 1) * modes + j];
                const double scale = std::abs(next);
                if (scale > 0.0) sums[j * bins + bin] += std::abs(next - std::exp(rates[j] * dt) * psi[k * modes + j]) / scale;
            }
        }
        for (std::size_t j = 0; j < modes; ++j)
            for (std::size_t b = 0; b < bins; ++b) sums[j * bins + b] = counts[b] > 0 ? sums[j * bins + b] / counts[b] : nan;
        archive.arrays["binned_one_step_defect"] = dynsyslab::core::NpyArray::fromDoubles({modes, bins}, sums);
        archive.arrays["binned_defect_edges"] = dynsyslab::core::NpyArray::fromDoubles({bins + 1}, edges);
        archive.arrays["binned_defect_counts"] = dynsyslab::core::NpyArray::fromDoubles({bins}, counts);
    }
    if (history) {
        std::vector<std::complex<double>> predicted(samples * modes);
        for (std::size_t k = 0; k < samples; ++k)
            for (std::size_t j = 0; j < modes; ++j)
                predicted[k * modes + j] = std::exp(rates[j] * (time[k] - time[0])) * psi[j];
        archive.arrays["psi_history"] = dynsyslab::core::NpyArray::fromComplex({samples, modes}, psi);
        archive.arrays["psi_predicted"] = dynsyslab::core::NpyArray::fromComplex({samples, modes}, predicted);
        archive.arrays["psi_history_time"] = dynsyslab::core::NpyArray::fromDoubles({samples}, time);
    }
}

}  // namespace dynsyslab::workflow::postprocess_detail

// derivative_benchmark: how well each numerical differentiator recovers the time derivative of a
// measured channel under each kind and strength of measurement noise.
//
// settings:
//   "observable": "far_field_acoustic_echo",            the measured channel (catalog observable)
//   "reference_derivative": "far_field_acoustic_echo_dot",  its exact derivative (for scoring only)
//   "time_window": [t0, t1],                            optional sample window
//   "noises":  [{"label": "white 1%", "kind": "white", "level": 0.01}, ...]
//   "methods": [{"label": "SG", "method": "savitzky_golay", "order": 3, "cutoffs": [5, 10, 20]}, ...]
//   "seeds": 1,                                         noise realizations averaged per entry
//   "trace_window": [t0, t1]                            optional: store the best derivative traces there
// A method entry without "cutoffs" is run once with its own parameters. Scores are relative RMS
// errors over the interior (2% trimmed each end). Causal methods also report the lag that minimizes
// their error and the error after removing it (lag is the price of causality).
//
// Output: one archive with results_json (every run), best_json (best cutoff per noise x method),
// and the best derivative traces for plotting.
#include "PostprocessInternal.hpp"
#include "DerivativeSettings.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
namespace {
double relErr(const std::vector<double>& est, const std::vector<double>& ref, std::size_t skip, std::size_t lag = 0) {
    double e = 0, r = 0;
    for (std::size_t i = skip + lag; i + skip < ref.size(); ++i) { const double d = est[i] - ref[i - lag]; e += d * d; r += ref[i - lag] * ref[i - lag]; }
    return r > 0 ? std::sqrt(e / r) : 0.0;
}
}  // namespace

std::filesystem::path runDerivativeBenchmark(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                             const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("derivative_benchmark needs a source record");
    const auto source = loadSource(sources.front());
    const auto catalog = c::createSystem(systemConfiguration(source))->observableCatalog();
    const auto observable = settings.at("observable").string(), referenceName = settings.at("reference_derivative").string();
    if (observable.empty() || referenceName.empty()) throw std::invalid_argument("derivative_benchmark needs observable and reference_derivative");
    const auto& window = settings.at("time_window").array();
    std::size_t i0 = 0, i1 = source.samples;
    if (window.size() == 2) {
        while (i0 < source.samples && source.time[i0] < window[0].number()) ++i0;
        i1 = i0; while (i1 < source.samples && source.time[i1] <= window[1].number()) ++i1;
    }
    const std::size_t n = i1 - i0;
    if (n < 50) throw std::invalid_argument("derivative_benchmark window has too few samples");
    const double dt = source.time[i0 + 1] - source.time[i0];
    std::vector<double> y(n), ref(n), time(n);
    const auto& map = catalog.at(observable); const auto& refMap = catalog.at(referenceName);
    for (std::size_t i = 0; i < n; ++i) {
        const double* x = &source.values[(i0 + i) * source.dimension];
        y[i] = map.evaluateReal(x).front(); ref[i] = refMap.evaluateReal(x).front(); time[i] = source.time[i0 + i];
    }
    const std::size_t skip = n / 50;
    const auto seeds = static_cast<std::size_t>(std::max(1.0, settings.at("seeds").number(1)));
    const auto& noises = settings.at("noises").array(); const auto& methods = settings.at("methods").array();
    if (noises.empty() || methods.empty()) throw std::invalid_argument("derivative_benchmark needs noises and methods");
    // traces
    std::size_t t0 = 0, t1 = 0;
    if (const auto& tw = settings.at("trace_window").array(); tw.size() == 2) {
        while (t0 < n && time[t0] < tw[0].number()) ++t0;
        t1 = t0; while (t1 < n && time[t1] <= tw[1].number()) ++t1;
    }
    const std::size_t T = t1 - t0;
    std::vector<double> traceDerivative(noises.size() * methods.size() * T, 0.0), traceMeasured(noises.size() * T);
    json::Value::Array results, best;
    for (std::size_t a = 0; a < noises.size(); ++a) {
        const auto noiseSpec = parseNoiseSpec(noises[a]);
        const auto noiseLabel = noises[a].at("label").string(noiseSpec.kind);
        std::vector<std::vector<double>> noisy(seeds);
        for (std::size_t s = 0; s < seeds; ++s) { auto ns = noiseSpec; ns.seed += 7919u * s; noisy[s] = c::addMeasurementNoise(y, dt, ns); }
        for (std::size_t k = 0; k < T; ++k) traceMeasured[a * T + k] = noisy[0][t0 + k];
        for (std::size_t b = 0; b < methods.size(); ++b) {
            auto spec = parseDerivativeSpec(methods[b]);
            const auto methodLabel = methods[b].at("label").string(spec.method);
            const bool causal = c::isCausalDerivativeMethod(spec.method);
            std::vector<double> cutoffs;
            for (const auto& v : methods[b].at("cutoffs").array()) cutoffs.push_back(v.number());
            if (cutoffs.empty()) cutoffs.push_back(spec.cutoff);
            double bestErr = INFINITY, bestCut = 0, bestLagErr = 0, bestLag = 0; std::vector<double> bestTrace;
            for (const double cut : cutoffs) {
                spec.cutoff = cut;
                double err = 0, valErr = 0, lagErr = 0, lagTime = 0; std::vector<double> first; bool failed = false; std::string why;
                for (std::size_t s = 0; s < seeds; ++s) {
                    c::DerivativeResult out;
                    try { out = c::differentiate(noisy[s], dt, spec); }
                    catch (const std::exception& e) { failed = true; why = e.what(); break; }
                    err += relErr(out.derivative, ref, skip) / double(seeds);
                    valErr += relErr(out.values, y, skip) / double(seeds);
                    if (causal) {  // best lag up to 2% of the record
                        double bestL = INFINITY; std::size_t lagBest = 0;
                        for (std::size_t l = 0; l <= std::min<std::size_t>(skip, 400); ++l) { const double e = relErr(out.derivative, ref, skip, l); if (e < bestL) { bestL = e; lagBest = l; } }
                        lagErr += bestL / double(seeds); lagTime += double(lagBest) * dt / double(seeds);
                    }
                    if (s == 0) first = std::move(out.derivative);
                }
                json::Value::Object row{{"noise", json::Value(noiseLabel)}, {"method", json::Value(methodLabel)},
                                        {"cutoff", json::Value(cut)}, {"causal", json::Value(causal)}};
                if (failed) { row["error"] = json::Value(why); results.emplace_back(std::move(row)); continue; }
                row["derivative_error"] = json::Value(err); row["value_error"] = json::Value(valErr);
                if (causal) { row["lag_time"] = json::Value(lagTime); row["lag_compensated_error"] = json::Value(lagErr); }
                results.emplace_back(std::move(row));
                if (err < bestErr) { bestErr = err; bestCut = cut; bestLagErr = lagErr; bestLag = lagTime; bestTrace = std::move(first); }
            }
            json::Value::Object row{{"noise", json::Value(noiseLabel)}, {"method", json::Value(methodLabel)}, {"causal", json::Value(causal)},
                                    {"best_cutoff", json::Value(bestCut)}, {"derivative_error", json::Value(std::isfinite(bestErr) ? bestErr : -1.0)}};
            if (causal) { row["lag_time"] = json::Value(bestLag); row["lag_compensated_error"] = json::Value(bestLagErr); }
            best.emplace_back(std::move(row));
            if (!bestTrace.empty()) for (std::size_t k = 0; k < T; ++k) traceDerivative[(a * methods.size() + b) * T + k] = bestTrace[t0 + k];
        }
    }
    c::NpzRecord archive;
    archive.arrays["results_json"] = c::NpyArray::fromString(json::dump(json::Value(results)));
    archive.arrays["best_json"] = c::NpyArray::fromString(json::dump(json::Value(best)));
    if (T) {
        archive.arrays["trace_time"] = c::NpyArray::fromDoubles({T}, std::vector<double>(time.begin() + std::ptrdiff_t(t0), time.begin() + std::ptrdiff_t(t1)));
        archive.arrays["trace_reference"] = c::NpyArray::fromDoubles({T}, std::vector<double>(ref.begin() + std::ptrdiff_t(t0), ref.begin() + std::ptrdiff_t(t1)));
        archive.arrays["trace_clean"] = c::NpyArray::fromDoubles({T}, std::vector<double>(y.begin() + std::ptrdiff_t(t0), y.begin() + std::ptrdiff_t(t1)));
        archive.arrays["trace_measured"] = c::NpyArray::fromDoubles({noises.size(), T}, traceMeasured);
        archive.arrays["trace_derivative"] = c::NpyArray::fromDoubles({noises.size(), methods.size(), T}, traceDerivative);
    }
    json::Value::Object meta{{"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("derivative-benchmark")},
        {"label", json::Value(label)}, {"source_archive", json::Value(sources.front().generic_u8string())},
        {"observable", json::Value(observable)}, {"reference_derivative", json::Value(referenceName)},
        {"delta_t", json::Value(dt)}, {"samples", json::Value(double(n))}, {"settings", settings},
        {"independent_archive", json::Value(true)}};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = objectFolder(root, "derivative-benchmark");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(sources.front().stem().u8string() + "__derivative-benchmark__" + timestamp()) + ".npz");
    c::writeNpzAtomic(path, archive);
    return path;
}
}  // namespace dynsyslab::workflow::postprocess_detail

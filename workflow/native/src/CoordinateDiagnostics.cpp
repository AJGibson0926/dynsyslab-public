#include "CoordinateDiagnostics.hpp"
#include <dynsyslab/core/CoordinateClosure.hpp>
#include <dynsyslab/core/PathCrossings.hpp>
#include <dynsyslab/core/System.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;

// RMS distance of the rows of a (rows x width) matrix from their mean.
double rmsRadius(const std::vector<double>& v, std::size_t rows, std::size_t width) {
    std::vector<double> mean(width, 0.0);
    for (std::size_t k = 0; k < rows; ++k) for (std::size_t c = 0; c < width; ++c) mean[c] += v[k * width + c];
    for (auto& m : mean) m /= double(rows);
    double sum = 0.0;
    for (std::size_t k = 0; k < rows; ++k) for (std::size_t c = 0; c < width; ++c) {
        const double d = v[k * width + c] - mean[c]; sum += d * d; }
    return std::sqrt(sum / double(rows));
}

std::vector<double> at(const std::vector<double>& v, std::size_t width, std::size_t k, double f) {
    std::vector<double> out(width);
    for (std::size_t c = 0; c < width; ++c) out[c] = v[k * width + c] + f * (v[(k + 1) * width + c] - v[k * width + c]);
    return out;
}

double distance(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0; for (std::size_t c = 0; c < a.size(); ++c) s += (a[c] - b[c]) * (a[c] - b[c]); return std::sqrt(s);
}

// Coordinate velocities: the recorded derivatives, else centred differences.
std::vector<double> velocities(const c::CoordinateSeries& series) {
    const auto& z = series.coordinates;
    const std::size_t n = z.sampleCount(), d = z.componentCount();
    if (series.derivatives) return {series.derivatives->valuesData(), series.derivatives->valuesData() + n * d};
    std::vector<double> out(n * d, 0.0);
    const double* t = z.timeData(); const double* v = z.valuesData();
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t a = k ? k - 1 : 0, b = k + 1 < n ? k + 1 : k;
        if (b == a) continue;
        for (std::size_t c = 0; c < d; ++c) out[k * d + c] = (v[b * d + c] - v[a * d + c]) / (t[b] - t[a]);
    }
    return out;
}
}  // namespace

void appendCoordinateDiagnostics(const SourceData& source, const J& settings, const c::CoordinateSeries& series,
    const std::optional<std::vector<double>>& states, const c::Observable* analyticMap, J& meta, c::NpzRecord& archive) {
    const auto& requested = settings.at("diagnostics");
    if (requested.isNull() || (std::holds_alternative<bool>(requested.storage()) && !requested.boolean())) return;
    const auto& z = series.coordinates;
    const std::size_t n = z.sampleCount(), d = z.componentCount();
    if (n < 4 || d < 2) return;
    const std::vector<double> values(z.valuesData(), z.valuesData() + n * d), time(z.timeData(), z.timeData() + n);
    const std::size_t stateWidth = states ? states->size() / n : 0;
    J::Object report;

    // --- crossings in the plotted pair, lifted into all coordinates
    const auto& crossing = requested.at("crossings");
    const auto& axes = crossing.at("axes").array();
    const std::size_t ax = axes.size() == 2 ? std::size_t(axes[0].number()) : 0;
    const std::size_t ay = axes.size() == 2 ? std::size_t(axes[1].number()) : 1;
    if (ax >= d || ay >= d || ax == ay) throw std::invalid_argument("diagnostics crossings.axes must name two distinct coordinates");
    const double tolerance = crossing.at("tolerance").number(1e-3);
    const std::size_t keep = std::size_t(crossing.at("max_witnesses").number(20));
    std::vector<double> x(n), y(n);
    for (std::size_t k = 0; k < n; ++k) { x[k] = values[k * d + ax]; y[k] = values[k * d + ay]; }
    const auto found = c::findPathCrossings(x, y);
    const auto rates = velocities(series);
    const double zRadius = rmsRadius(values, n, d), xRadius = states ? rmsRadius(*states, n, stateWidth) : 0.0;
    double speed = 0.0; for (double r : rates) speed += r * r; speed = std::sqrt(speed / double(n));
    std::vector<double> segments, metrics;  // metrics: u, v, tau_a, tau_b, lifted, state, velocity
    struct Witness { double lifted, state, velocity, ta, tb; std::vector<double> xa, xb; };
    std::vector<Witness> witnesses;
    std::size_t genuine = 0, noninjective = 0, obstructions = 0;
    double minimumLifted = INFINITY;
    for (const auto& hit : found) {
        const auto za = at(values, d, hit.first, hit.u), zb = at(values, d, hit.second, hit.v);
        const auto va = at(rates, d, hit.first, hit.u), vb = at(rates, d, hit.second, hit.v);
        Witness w{distance(za, zb) / zRadius, NAN, speed > 0 ? distance(va, vb) / speed : NAN,
                  time[hit.first] + hit.u * (time[hit.first + 1] - time[hit.first]),
                  time[hit.second] + hit.v * (time[hit.second + 1] - time[hit.second]), {}, {}};
        if (states) {
            w.xa = at(*states, stateWidth, hit.first, hit.u); w.xb = at(*states, stateWidth, hit.second, hit.v);
            w.state = distance(w.xa, w.xb) / xRadius;
        }
        minimumLifted = std::min(minimumLifted, w.lifted);
        if (w.lifted <= tolerance) {
            ++genuine;
            if (states && w.state > tolerance) { ++noninjective; if (w.velocity > tolerance) ++obstructions; }
        }
        segments.push_back(double(hit.first)); segments.push_back(double(hit.second));
        for (double m : {hit.u, hit.v, w.ta, w.tb, w.lifted, w.state, w.velocity}) metrics.push_back(m);
        witnesses.push_back(std::move(w));
    }
    std::sort(witnesses.begin(), witnesses.end(), [](const Witness& a, const Witness& b) { return a.lifted < b.lifted; });
    J::Array listed;
    for (std::size_t k = 0; k < witnesses.size() && k < keep; ++k) {
        const auto& w = witnesses[k];
        J::Object item{{"tau", J::Array{w.ta, w.tb}}, {"lifted_separation", w.lifted}, {"velocity_mismatch", w.velocity}};
        if (states) { item["state_separation"] = w.state; item["states"] = J::Array{J::Array(w.xa.begin(), w.xa.end()), J::Array(w.xb.begin(), w.xb.end())}; }
        listed.emplace_back(std::move(item));
    }
    report["crossings"] = J::Object{{"axes", J::Array{double(ax), double(ay)}}, {"count", double(found.size())},
        {"tolerance", tolerance}, {"genuine", double(genuine)}, {"noninjective_witnesses", double(noninjective)},
        {"closure_obstructions", double(obstructions)},
        {"minimum_lifted_separation", found.empty() ? J{} : J(minimumLifted)},
        {"separations", "fractions of the RMS radius of the coordinates, of the physical state and of the coordinate speed"},
        {"velocity_source", series.derivatives ? "recorded derivatives" : "centred differences"},
        {"witnesses", std::move(listed)}};
    if (!found.empty()) {
        archive.arrays["crossing_segments"] = c::NpyArray::fromDoubles({found.size(), 2}, segments);
        archive.arrays["crossing_metrics"] = c::NpyArray::fromDoubles({found.size(), 7}, metrics);
    }

    // --- analytic closure search for a map of the physical state
    const auto& closure = requested.at("closure");
    const bool closureOff = std::holds_alternative<bool>(closure.storage()) && !closure.boolean();
    if (analyticMap && states && !closureOff) {
        const auto system = c::createSystem(systemConfiguration(source));
        c::ObservableClosureOptions options;
        options.highestDerivativeOrder = 0;
        options.sampleCount = std::size_t(closure.at("sample_count").number(512));
        options.stateMinimum.assign(stateWidth, INFINITY); options.stateMaximum.assign(stateWidth, -INFINITY);
        for (std::size_t k = 0; k < n; ++k) for (std::size_t c = 0; c < stateWidth; ++c) {
            options.stateMinimum[c] = std::min(options.stateMinimum[c], (*states)[k * stateWidth + c]);
            options.stateMaximum[c] = std::max(options.stateMaximum[c], (*states)[k * stateWidth + c]); }
        for (const auto& w : witnesses)  // seed with the trajectory's own coincidences
            if (w.lifted <= tolerance && !w.xa.empty() && options.pairSeeds.size() < 32 * stateWidth) {
                options.pairSeeds.insert(options.pairSeeds.end(), w.xa.begin(), w.xa.end());
                options.pairSeeds.insert(options.pairSeeds.end(), w.xb.begin(), w.xb.end()); }
        try {
            const auto result = c::analyzeObservableClosure(*system, *analyticMap, options);
            const auto& order = result.orders.front();
            auto witness = [](const c::ObservableCollisionWitness& w) {
                J::Object o{{"found", w.found}};
                if (w.found) { o["left_state"] = J::Array(w.leftState.begin(), w.leftState.end());
                    o["right_state"] = J::Array(w.rightState.begin(), w.rightState.end());
                    o["normalized_residual"] = w.normalizedResidual; o["normalized_state_separation"] = w.normalizedStateSeparation;
                    o["normalized_successor_difference"] = w.normalizedSuccessorDifference; }
                return J(std::move(o)); };
            report["closure"] = J::Object{{"method", "native observable closure search over the record's state box"},
                {"accepted_samples", double(result.acceptedSamples)}, {"minimum_sampled_rank", double(order.minimumSampledRank)},
                {"coordinate_collision", witness(order.coordinateCollision)},
                {"closure_obstruction", witness(order.closureObstruction)},
                {"scope", "finite search: a witness disproves injectivity/closure; none found is not a proof"}};
        } catch (const std::exception& error) {
            report["closure"] = J::Object{{"error", std::string(error.what())}};
        }
    }
    meta.object()["coordinate_diagnostics"] = J(std::move(report));
}
} // namespace dynsyslab::workflow::postprocess_detail

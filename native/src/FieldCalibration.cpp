// Fixed-support trajectory calibration (FieldModels.hpp). Ported from
// tools/experiments/miscellaneous_computes/calibration.py. The nonzero coefficients of a polynomial
// model are refitted against short free trajectories (multiple shooting starts; Bock and Plitt
// 1984, IFAC Proc. 17:1603-1608) by bounded damped Gauss-Newton (Nocedal and Wright 2006,
// chapter 10); an incomplete or worse candidate is rejected, keeping the uncalibrated model.
#include "FieldModelsInternal.hpp"

#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

struct Shot {
    std::vector<double> time, state, weight;  // weight: samples x n (waveform), empty for RMS scaling
};

}  // namespace

CalibrationResult calibrateField(const json::Value& definition, const std::vector<FieldRecord>& records, const CalibrationOptions& o) {
    if (definition.at("kind").string() != "polynomial_model") throw std::invalid_argument("calibration: a polynomial_model field");
    const auto model = SindyDynamicsModel::fromDefinition(definition.at("model"));
    const auto n = model.dimension();
    if (!model.isVectorField() || records.empty() || !std::isfinite(o.duration) || o.duration <= 0 || !o.startsPerRecord
        || !o.maximumEvaluations || o.coefficientPenalty < 0 || !(o.changeBound > 0))
        throw std::invalid_argument("calibration: invalid shooting settings");
    const auto& C = model.coefficients();
    std::vector<std::size_t> active;
    for (std::size_t i = 0; i < C.size(); ++i) if (C[i] != 0) active.push_back(i);
    if (active.empty()) throw std::invalid_argument("calibration: empty support");
    const bool waveform = !o.globalWeights.empty();
    std::vector<Shot> shots;
    for (const auto& r : records) {
        if (r.dimension != n || r.state.size() != r.samples() * n || r.samples() < 2) throw std::invalid_argument("calibration: record dimensions");
        const auto& t = r.time;
        std::size_t last = 0;
        bool any = false;
        for (std::size_t i = 0; i < t.size(); ++i) if (t[i] + o.duration <= t.back()) { last = i; any = true; }
        if (!any) throw std::invalid_argument("calibration: training record shorter than requested shot");
        std::vector<std::size_t> starts;
        for (std::size_t s = 0; s < o.startsPerRecord; ++s) {
            const auto i = o.startsPerRecord == 1 ? 0 : std::size_t(std::floor(double(last) * double(s) / double(o.startsPerRecord - 1)));
            if (starts.empty() || starts.back() != i) starts.push_back(i);
        }
        for (const auto i : starts) {
            const auto stop = std::size_t(std::upper_bound(t.begin(), t.end(), t[i] + o.duration) - t.begin());
            if (stop <= i + 1) throw std::invalid_argument("calibration: shot has fewer than two output samples");
            Shot shot{{t.begin() + std::ptrdiff_t(i), t.begin() + std::ptrdiff_t(stop)},
                      {r.state.begin() + std::ptrdiff_t(i * n), r.state.begin() + std::ptrdiff_t(stop * n)}, {}};
            if (waveform)
                shot.weight = waveformWeights(shot.time, shot.state, n, o.globalWeights, o.localWeights.empty() ? std::vector<double>{0} : o.localWeights,
                                              o.waveformDuration > 0 ? o.waveformDuration : o.duration, o.floors.empty() ? std::vector<double>{1e-12} : o.floors);
            shots.push_back(std::move(shot));
        }
    }
    std::vector<double> scale(n, 0);
    std::size_t count = 0;
    for (const auto& s : shots) { for (std::size_t k = 0; k < s.time.size(); ++k) for (std::size_t i = 0; i < n; ++i) scale[i] += s.state[k * n + i] * s.state[k * n + i]; count += s.time.size(); }
    for (auto& v : scale) { v = std::sqrt(v / double(count)); if (!(v > 0)) v = 1; }
    std::vector<double> c0, unit;
    for (const auto i : active) { c0.push_back(C[i]); unit.push_back(std::max(std::abs(C[i]), 1e-10)); }
    const auto candidate = [&](const std::vector<double>& change) {
        auto cc = C;
        for (std::size_t a = 0; a < active.size(); ++a) cc[active[a]] = c0[a] + unit[a] * change[a];
        return SindyDynamicsModel(model.library(), cc, n);
    };
    bool lastComplete = true;
    const auto residual = [&](const std::vector<double>& change) {
        const auto field = candidate(change).vectorField();
        std::vector<double> out;
        lastComplete = true;
        for (const auto& s : shots) {
            const auto size = s.state.size();
            const double norm = waveform ? 1.0 : 1.0 / std::sqrt(double(size));
            std::vector<double> part(size, 1000.0 * norm);
            try {
                const std::vector<double> x0(s.state.begin(), s.state.begin() + std::ptrdiff_t(n));
                const auto pred = rolloutField(field, x0, s.time, o.maximumStep);
                for (std::size_t e = 0; e < pred.state.size(); ++e) {
                    const double d = pred.state[e] - s.state[e];
                    part[e] = waveform ? d * s.weight[e] : d / scale[e % n] * norm;
                }
                lastComplete = lastComplete && pred.complete;
            } catch (const std::exception&) { lastComplete = false; }
            for (auto& v : part) if (!std::isfinite(v)) { v = 1000.0 * norm; lastComplete = false; }
            out.insert(out.end(), part.begin(), part.end());
        }
        for (const auto d : change) out.push_back(std::sqrt(o.coefficientPenalty) * d);
        return out;
    };
    const auto loss = [](const std::vector<double>& r) { double s = 0; for (const auto v : r) s += v * v; return s; };
    const std::vector<double> zero(active.size(), 0.0);
    CalibrationResult result;
    result.initialLoss = loss(residual(zero));
    numerics::NonlinearLeastSquaresOptions lo;
    lo.maximumIterations = o.maximumEvaluations;
    lo.finiteDifferenceStep = 1e-6;
    const auto solution = numerics::boundedNonlinearLeastSquares(residual, zero, std::vector<double>(active.size(), -o.changeBound),
                                                                 std::vector<double>(active.size(), o.changeBound), lo);
    result.candidateLoss = loss(residual(solution.parameters));
    result.accepted = lastComplete && result.candidateLoss <= result.initialLoss;
    result.shots = shots.size(); result.terms = active.size(); result.evaluations = solution.evaluations;
    result.definition = result.accepted
        ? json::Value(json::Value::Object{{"kind", "polynomial_model"}, {"model", candidate(solution.parameters).definition()},
                                          {"family", definition.at("family")}})
        : definition;
    return result;
}

}  // namespace dynsyslab::core

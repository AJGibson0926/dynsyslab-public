#include <dynsyslab/core/DelayRollout.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

DelayRollout rolloutDelayModel(const DelayRightHandSide& rhs, std::size_t dimension,
                               const std::vector<double>& history, const DelayRolloutOptions& options) {
    const double h = options.step;
    const auto& lags = options.lags;
    if (!rhs || !dimension || lags.empty() || !options.outputSamples || !std::isfinite(h) || h <= 0.0
        || lags.front() != 0.0 || history.size() % dimension)
        throw std::invalid_argument("delay rollout needs a model, a positive step, lags starting at 0 and output rows");
    const std::size_t historyRows = history.size() / dimension, taps = lags.size();
    if (historyRows < 2) throw std::invalid_argument("delay rollout needs at least two history rows");
    for (std::size_t j = 1; j < taps; ++j)
        if (!std::isfinite(lags[j]) || lags[j] < h || lags[j] <= lags[j - 1]
            || lags[j] / h > static_cast<double>(historyRows - 1))
            throw std::invalid_argument("delay rollout lags must increase, be at least one step and fit the history");
    if (!std::all_of(history.begin(), history.end(), [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument("delay rollout history must be finite");
    const auto& limits = options.magnitudeLimits;
    if (!limits.empty() && (limits.size() != dimension
        || std::any_of(limits.begin(), limits.end(), [](double value) { return !(value > 0.0); })))
        throw std::invalid_argument("delay rollout magnitude limits must be positive, one per component");

    // Lag j at RK phase p (0: t, 1: t + h/2, 2: t + h) lies `behind` whole samples plus `fraction`
    // before the newest stored sample.
    struct Offset { std::size_t behind; double fraction; };
    std::vector<Offset> offsets[3];
    for (int phase = 0; phase < 3; ++phase)
        for (std::size_t j = 1; j < taps; ++j) {
            const double distance = lags[j] / h - 0.5 * phase;
            const double whole = std::floor(distance);
            offsets[phase].push_back({static_cast<std::size_t>(whole), distance - whole});
        }
    DelayRollout result;
    result.values.reserve(options.outputSamples * dimension);
    const auto row = [&](std::size_t index) {
        return index < historyRows ? history.data() + index * dimension
                                   : result.values.data() + (index - historyRows + 1) * dimension;
    };
    const auto guard = [&](const double* x) {
        for (std::size_t v = 0; v < dimension; ++v) {
            if (!std::isfinite(x[v])) { result.status = DelayRolloutStatus::nonFinite; result.reason = "non-finite state"; return false; }
            if (!limits.empty() && std::abs(x[v]) >= limits[v]) {
                result.status = DelayRolloutStatus::magnitudeGuard;
                result.reason = "state component " + std::to_string(v) + " reached its magnitude limit";
                return false;
            }
        }
        return true;
    };
    const double* first = history.data() + (historyRows - 1) * dimension;
    if (!guard(first)) return result;
    result.values.insert(result.values.end(), first, first + dimension);
    result.samples = 1;
    std::vector<double> window(taps * dimension), stage(dimension), k1(dimension), k2(dimension),
        k3(dimension), k4(dimension), next(dimension);
    for (std::size_t i = 0; i + 1 < options.outputSamples; ++i) {
        const std::size_t last = historyRows - 1 + i;
        const auto evaluate = [&](const double* current, int phase, double* velocity) {
            if (!guard(current)) return false;
            std::copy_n(current, dimension, window.data());
            for (std::size_t j = 1; j < taps; ++j) {
                const auto offset = offsets[phase][j - 1];
                const std::size_t late = last - offset.behind;
                const double* b = row(late);
                const double* a = offset.fraction == 0.0 ? b : row(late - 1);
                for (std::size_t v = 0; v < dimension; ++v)
                    window[j * dimension + v] = b[v] + offset.fraction * (a[v] - b[v]);
            }
            rhs(window.data(), velocity);
            for (std::size_t v = 0; v < dimension; ++v)
                if (!std::isfinite(velocity[v])) {
                    result.status = DelayRolloutStatus::nonFinite; result.reason = "non-finite derivative";
                    return false;
                }
            return true;
        };
        const std::vector<double> x(result.values.end() - static_cast<std::ptrdiff_t>(dimension), result.values.end());
        if (!evaluate(x.data(), 0, k1.data())) return result;
        for (std::size_t v = 0; v < dimension; ++v) stage[v] = x[v] + 0.5 * h * k1[v];
        if (!evaluate(stage.data(), 1, k2.data())) return result;
        for (std::size_t v = 0; v < dimension; ++v) stage[v] = x[v] + 0.5 * h * k2[v];
        if (!evaluate(stage.data(), 1, k3.data())) return result;
        for (std::size_t v = 0; v < dimension; ++v) stage[v] = x[v] + h * k3[v];
        if (!evaluate(stage.data(), 2, k4.data())) return result;
        for (std::size_t v = 0; v < dimension; ++v)
            next[v] = x[v] + (h / 6.0) * (k1[v] + 2.0 * k2[v] + 2.0 * k3[v] + k4[v]);
        if (!guard(next.data())) return result;
        result.values.insert(result.values.end(), next.begin(), next.end());
        ++result.samples;
    }
    return result;
}

}  // namespace dynsyslab::core

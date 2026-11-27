// Delayed observations at physical-time lags (DelaySampling.hpp).
#include <dynsyslab/core/DelaySampling.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void validLags(const std::vector<double>& endpoints, const std::vector<double>& lags) {
    if (endpoints.empty() || lags.empty()) throw std::invalid_argument("delay sampling needs endpoints and lags");
    for (const double value : endpoints) if (!std::isfinite(value)) throw std::invalid_argument("delay endpoints must be finite");
    for (const double value : lags) if (!std::isfinite(value) || value < 0) throw std::invalid_argument("delay lags must be finite and nonnegative");
}

}  // namespace

std::vector<double> delaySamples(const TimeSeries& series, const std::vector<double>& endpoints, const std::vector<double>& lags) {
    validLags(endpoints, lags);
    const std::size_t n = series.sampleCount(), d = series.componentCount();
    if (n < 2) throw std::invalid_argument("delay sampling needs at least two samples");
    const double* t = series.timeData();
    double smallest = std::numeric_limits<double>::infinity();
    for (std::size_t k = 1; k < n; ++k) {
        if (!(t[k] > t[k-1])) throw std::invalid_argument("delay sampling needs an increasing clock");
        smallest = std::min(smallest, t[k] - t[k-1]);
    }
    std::vector<double> out;
    out.reserve(endpoints.size() * lags.size() * d);
    for (const double end : endpoints)
        for (const double lag : lags) {
            double q = end - lag;
            // Snap roundoff-sized differences to the nearest tick (including the record's start).
            const auto right = std::size_t(std::lower_bound(t, t + n, q) - t);
            const std::size_t near = right == 0 ? 0 : right == n ? n - 1
                : (std::abs(t[right] - q) < std::abs(t[right-1] - q) ? right : right - 1);
            const double eps = std::min(8 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(q)), 1e-6 * smallest);
            if (std::abs(t[near] - q) <= eps) q = t[near];
            if (q < t[0] || q > t[n-1]) throw std::invalid_argument("delay window reaches outside the record (incomplete causal history)");
            const auto hi = std::size_t(std::lower_bound(t, t + n, q) - t);  // first sample at or after q
            if (t[hi] == q) {  // on a tick: that sample exactly
                for (std::size_t c = 0; c < d; ++c) out.push_back(series.value(hi, c));
                continue;
            }
            if (t[hi] > end + eps) throw std::invalid_argument("a delay tap would read a sample after its endpoint");
            const std::size_t lo = hi - 1;
            const double alpha = (q - t[lo]) / (t[hi] - t[lo]);
            for (std::size_t c = 0; c < d; ++c)
                out.push_back(series.value(lo, c) + alpha * (series.value(hi, c) - series.value(lo, c)));
        }
    return out;
}

std::vector<double> sampleDelayJet(const TimeSeries& jet, const std::vector<double>& endpoints, const std::vector<double>& lags) {
    validLags(endpoints, lags);
    const std::size_t n = jet.sampleCount();
    if (jet.componentCount() != 3 || n < 2) throw std::invalid_argument("a delay jet holds value, rate and acceleration columns");
    const double t0 = jet.time(0), dt = (jet.time(n - 1) - t0) / double(n - 1);
    for (std::size_t k = 1; k < n; ++k)
        if (std::abs(jet.time(k) - (t0 + double(k) * dt)) > 1e-9 * dt * double(k + 1))
            throw std::invalid_argument("a delay jet needs a uniform clock");
    std::vector<double> out;
    out.reserve(endpoints.size() * lags.size() * 3);
    for (const double end : endpoints)
        for (const double lag : lags) {
            double ticks = (end - lag - t0) / dt;
            const double nearest = std::round(ticks);
            if (std::abs(ticks - nearest) <= 2e-8) ticks = nearest;
            if (ticks < 0 || ticks > double(n - 1)) throw std::invalid_argument("a delay jet tap needs source history it does not have");
            if (jet.time(std::size_t(std::ceil(ticks))) > end + dt * 2e-8)
                throw std::invalid_argument("fractional delay interpolation would use an unobserved future sample");
            const auto i = std::min(std::size_t(ticks), n - 2);
            const double u = ticks - double(i);
            const double p = jet.value(i, 0), v = jet.value(i, 1), a = jet.value(i, 2);
            const double q = jet.value(i + 1, 0), w = jet.value(i + 1, 1), b = jet.value(i + 1, 2);
            // y(u) = p + C1 u + C2 u^2 + c3 u^3 + c4 u^4 + c5 u^5 on u in [0, 1].
            const double C1 = dt * v, C2 = 0.5 * dt * dt * a;
            const double gap = (q - p) - C1 - C2, slope = dt * w - C1 - 2 * C2, curve = dt * dt * b - 2 * C2;
            const double c3 = 10 * gap - 4 * slope + 0.5 * curve, c4 = -15 * gap + 7 * slope - curve, c5 = 6 * gap - 3 * slope + 0.5 * curve;
            out.push_back(p + u * (C1 + u * (C2 + u * (c3 + u * (c4 + u * c5)))));
            out.push_back((C1 + u * (2 * C2 + u * (3 * c3 + u * (4 * c4 + 5 * u * c5)))) / dt);
            out.push_back((2 * C2 + u * (6 * c3 + u * (12 * c4 + 20 * u * c5))) / (dt * dt));
        }
    return out;
}

}  // namespace dynsyslab::core

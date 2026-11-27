// Derivative-consistent scalar readout (core/TangentReadout.hpp, which cites the method).
#include <dynsyslab/core/TangentReadout.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
std::vector<double> rmsScales(const std::vector<double>& b, std::size_t rows, std::vector<double> scales) {
    if (!rows || b.size() % rows) throw std::invalid_argument("tangent readout: features must align with labels");
    const auto d = b.size() / rows;
    if (!scales.empty()) return scales;
    scales.assign(d, 0.0);
    for (std::size_t i = 0; i < b.size(); ++i) scales[i % d] += b[i] * b[i];
    for (auto& s : scales) s = std::sqrt(s / double(rows));
    return scales;
}
}  // namespace

TangentReadout::TangentReadout(std::vector<double> features, const std::vector<double>& rates, std::vector<double> labels,
                               const std::vector<double>& labelRates, std::vector<double> scales, double floor)
    : labels_(std::move(labels)), scales_(rmsScales(features, labels_.size(), std::move(scales))),
      search_(features, std::vector<double>(labels_.size(), 0.0), scales_, 1) {
    const auto n = labels_.size(), d = scales_.size();
    if (rates.size() != features.size() || labelRates.size() != n || !(floor > 0))
        throw std::invalid_argument("tangent readout: expected one scalar output, aligned velocities and a positive floor");
    x_.resize(n * d); slopes_.resize(n * d); amplitude_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        double norm = 0;
        std::vector<double> dx(d);
        for (std::size_t j = 0; j < d; ++j) {
            const double b = features[i * d + j];
            x_[i * d + j] = std::asinh(b / scales_[j]);
            dx[j] = rates[i * d + j] / std::hypot(scales_[j], b);
            norm += dx[j] * dx[j];
        }
        for (std::size_t j = 0; j < d; ++j) slopes_[i * d + j] = labelRates[i] * dx[j] / std::max(norm, 1e-14);
        amplitude_[i] = std::max(floor, std::abs(labels_[i]));
    }
}

TangentReadout::Phase TangentReadout::phase(const double* b, const double* v, std::size_t k) const {
    const auto d = scales_.size();
    if (!b || !v || k < 1 || k + 1 > labels_.size()) throw std::invalid_argument("tangent readout: invalid query");
    std::vector<double> x(d), dx(d);
    for (std::size_t j = 0; j < d; ++j) {
        x[j] = std::asinh(b[j] / scales_[j]);
        dx[j] = v[j] / std::hypot(scales_[j], b[j]);
    }
    const auto hits = search_.predict(b, k + 1).rows;  // nearest first; the last sets the support radius
    std::vector<double> distance(k + 1), speed(k + 1);
    for (std::size_t m = 0; m <= k; ++m) {
        double squares = 0, along = 0;
        for (std::size_t j = 0; j < d; ++j) {
            const double off = x[j] - x_[hits[m] * d + j];
            squares += off * off; along += off * dx[j];
        }
        distance[m] = std::sqrt(squares);
        speed[m] = along / std::max(distance[m], 1e-14);  // d(distance)/dt along the query velocity
    }
    const double radius = std::max(distance[k], 1e-13);
    double total = 0, weighted = 0, dweighted = 0;
    std::vector<double> w(k), dw(k), g(k), dg(k);
    for (std::size_t m = 0; m < k; ++m) {
        const auto i = hits[m];
        const double r = distance[m] / radius, dr = (speed[m] - r * speed[k]) / radius;
        w[m] = std::pow(1 - r, 4) * (1 + 4 * r);
        dw[m] = -20 * r * std::pow(1 - r, 3) * dr;
        double phase = 0, slope = 0;
        for (std::size_t j = 0; j < d; ++j) {
            phase += slopes_[i * d + j] * (x[j] - x_[i * d + j]);
            slope += slopes_[i * d + j] * dx[j];
        }
        const double t = std::tanh(phase / amplitude_[i]);
        g[m] = labels_[i] + amplitude_[i] * t;
        dg[m] = (1 - t * t) * slope;
        total += w[m];
    }
    if (total < 1e-25) {  // every neighbour at the support radius: equal weights
        for (std::size_t m = 0; m < k; ++m) { w[m] = 1; dw[m] = 0; }
        total = double(k);
    }
    for (std::size_t m = 0; m < k; ++m) weighted += w[m] * g[m];
    const double value = weighted / total;
    for (std::size_t m = 0; m < k; ++m) dweighted += dw[m] * (g[m] - value) + w[m] * dg[m];
    return {value, dweighted / total};
}

}  // namespace dynsyslab::core

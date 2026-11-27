// Flow maps on grids, the Cauchy-Green tensor and the FTLE (LagrangianCoherence.hpp).
#include <dynsyslab/core/LagrangianCoherence.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// Classical RK4 of the (possibly time-dependent) field from t0 over T in equal steps <= h
// (Hairer, Norsett and Wanner, Solving ODEs I, 1993, section II.1).
bool flow(const ODESystem& system, double* x, double t0, double T, double maximumStep) {
    const std::size_t n = system.descriptor().stateDimension;
    const auto steps = std::max<std::size_t>(1, std::size_t(std::ceil(std::abs(T) / maximumStep)));
    const double h = T / double(steps);
    std::vector<double> k1(n), k2(n), k3(n), k4(n), y(n);
    // Unforced systems skip the driver's clock and input conversions on every evaluation.
    const bool forced = system.configuration().driver.amplitude != 0;
    const auto f = [&](double time, const double* state, double* out) {
        if (forced) system.rhs(time, state, out);
        else system.rhsWithInput(time, state, 0, out);
    };
    double t = t0;
    for (std::size_t s = 0; s < steps; ++s) {
        f(t, x, k1.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k1[i];
        f(t + h / 2, y.data(), k2.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k2[i];
        f(t + h / 2, y.data(), k3.data());
        for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h * k3[i];
        f(t + h, y.data(), k4.data());
        for (std::size_t i = 0; i < n; ++i) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        t = t0 + double(s + 1) * h;
        for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
        if (system.invalid(t, x)) return false;
    }
    return true;
}

void check(const ODESystem& system, const DomainGrid& grid, const FlowMapGridOptions& options) {
    grid.validate();
    if (grid.dimension() != system.descriptor().stateDimension)
        throw std::invalid_argument("FTLE: the grid must live in the system's state space");
    if (!std::isfinite(options.integrationTime) || options.integrationTime == 0 || !(options.maximumStep > 0)
        || !std::isfinite(options.initialTime) || !(options.auxiliarySpacing >= 0))
        throw std::invalid_argument("FTLE: needs a finite nonzero integration time, a positive step and spacing >= 0");
    for (const auto& axis : grid.axes)
        if (axis.points < 2) throw std::invalid_argument("FTLE: every gridded axis needs at least two nodes");
}

}  // namespace

std::vector<double> flowMapOnGrid(const ODESystem& system, const DomainGrid& grid, const FlowMapGridOptions& options) {
    check(system, grid, options);
    const auto n = grid.dimension();
    auto state = grid.points();
    for (std::size_t p = 0; p < grid.size(); ++p)
        if (!flow(system, state.data() + p * n, options.initialTime, options.integrationTime, options.maximumStep))
            std::fill(state.begin() + std::ptrdiff_t(p * n), state.begin() + std::ptrdiff_t((p + 1) * n), nan);
    return state;
}

void symmetricEigen(std::vector<double> a, std::size_t k, std::vector<double>& values, std::vector<double>& vectors) {
    if (a.size() != k * k) throw std::invalid_argument("symmetricEigen: matrix size");
    vectors.assign(k * k, 0);
    for (std::size_t i = 0; i < k; ++i) vectors[i * k + i] = 1;
    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = 0, scale = 0;
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t j = 0; j < k; ++j) (i == j ? scale : off) += a[i * k + j] * a[i * k + j];
        if (off <= 1e-30 * std::max(scale, 1e-300)) break;
        for (std::size_t p = 0; p + 1 < k; ++p)
            for (std::size_t q = p + 1; q < k; ++q) {
                const double apq = a[p * k + q];
                if (apq == 0) continue;
                const double theta = (a[q * k + q] - a[p * k + p]) / (2 * apq);
                const double t = (theta >= 0 ? 1 : -1) / (std::abs(theta) + std::sqrt(theta * theta + 1));
                const double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (std::size_t r = 0; r < k; ++r) {  // A <- A J (columns p, q)
                    const double arp = a[r * k + p], arq = a[r * k + q];
                    a[r * k + p] = c * arp - s * arq;
                    a[r * k + q] = s * arp + c * arq;
                }
                for (std::size_t r = 0; r < k; ++r) {  // A <- J^T A (rows p, q)
                    const double apr = a[p * k + r], aqr = a[q * k + r];
                    a[p * k + r] = c * apr - s * aqr;
                    a[q * k + r] = s * apr + c * aqr;
                }
                for (std::size_t r = 0; r < k; ++r) {
                    const double vrp = vectors[r * k + p], vrq = vectors[r * k + q];
                    vectors[r * k + p] = c * vrp - s * vrq;
                    vectors[r * k + q] = s * vrp + c * vrq;
                }
            }
    }
    std::vector<std::size_t> order(k);
    for (std::size_t i = 0; i < k; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](auto l, auto r) { return a[l * k + l] < a[r * k + r]; });
    values.resize(k);
    std::vector<double> sorted(k * k);
    for (std::size_t j = 0; j < k; ++j) {
        values[j] = a[order[j] * k + order[j]];
        for (std::size_t r = 0; r < k; ++r) sorted[r * k + j] = vectors[r * k + order[j]];
    }
    vectors = std::move(sorted);
}

CauchyGreenField cauchyGreenField(const ODESystem& system, const DomainGrid& grid, const FlowMapGridOptions& options) {
    check(system, grid, options);
    const auto n = grid.dimension(), k = grid.axes.size(), points = grid.size();
    CauchyGreenField out;
    out.grid = grid;
    out.integrationTime = options.integrationTime;
    out.axes = k;
    out.eigenvalues.assign(points * k, nan);
    out.eigenvectors.assign(points * k * k, nan);
    out.status.assign(points, GridStatus::valid);
    std::vector<std::size_t> stride(k, 1);
    for (std::size_t a = 1; a < k; ++a) stride[a] = stride[a - 1] * grid.axes[a - 1].points;

    // Each node's flow-map gradient columns dPhi/dx_a (n x k), by central differences.
    std::vector<double> gradient(points * n * k, nan);
    if (options.auxiliarySpacing > 0) {
        const double d = options.auxiliarySpacing;
        std::vector<double> x(n);
        for (std::size_t p = 0; p < points; ++p) {
            grid.point(p, x.data());
            for (std::size_t a = 0; a < k; ++a) {
                auto plus = x, minus = x;
                plus[grid.axes[a].coordinate] += d;
                minus[grid.axes[a].coordinate] -= d;
                const bool ok = flow(system, plus.data(), options.initialTime, options.integrationTime, options.maximumStep)
                             && flow(system, minus.data(), options.initialTime, options.integrationTime, options.maximumStep);
                for (std::size_t i = 0; i < n; ++i) gradient[(p * n + i) * k + a] = ok ? (plus[i] - minus[i]) / (2 * d) : nan;
            }
        }
    } else {
        const auto image = flowMapOnGrid(system, grid, options);
        for (std::size_t p = 0; p < points; ++p)
            for (std::size_t a = 0; a < k; ++a) {
                const auto m = grid.axes[a].points, i = (p / stride[a]) % m;
                const auto lo = i == 0 ? p : p - stride[a], hi = i + 1 == m ? p : p + stride[a];
                const double h = (grid.axes[a].upper - grid.axes[a].lower) / double(m - 1) * double((i == 0 || i + 1 == m) ? 1 : 2);
                for (std::size_t r = 0; r < n; ++r) gradient[(p * n + r) * k + a] = (image[hi * n + r] - image[lo * n + r]) / h;
            }
    }
    std::vector<double> c(k * k), values, vectors;
    for (std::size_t p = 0; p < points; ++p) {
        const double* g = gradient.data() + p * n * k;
        bool ok = true;
        for (std::size_t a = 0; a < k; ++a)
            for (std::size_t b = 0; b < k; ++b) {
                double s = 0;
                for (std::size_t r = 0; r < n; ++r) s += g[r * k + a] * g[r * k + b];
                c[a * k + b] = s;
                ok = ok && std::isfinite(s);
            }
        if (!ok) { out.status[p] = GridStatus::undefined; continue; }
        symmetricEigen(c, k, values, vectors);
        for (std::size_t j = 0; j < k; ++j) out.eigenvalues[p * k + j] = values[j];
        std::copy(vectors.begin(), vectors.end(), out.eigenvectors.begin() + std::ptrdiff_t(p * k * k));
    }
    return out;
}

std::vector<double> CauchyGreenField::ftle() const {
    std::vector<double> out(size(), nan);
    for (std::size_t p = 0; p < size(); ++p) {
        const double l = eigenvalues[p * axes + axes - 1];
        if (status[p] == GridStatus::valid && l > 0) out[p] = std::log(l) / (2 * std::abs(integrationTime));
    }
    return out;
}

}  // namespace dynsyslab::core

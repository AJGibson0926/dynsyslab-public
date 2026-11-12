#include "AnalyticQuadrature.hpp"
#include <algorithm>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::detail {
double analyticIntegral(const std::function<double(double)> &f, double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b))
        throw std::domain_error("Analytic quadrature requires finite endpoints");
    if (a == b)
        return 0;
    double error = 0, l1 = 0;
    // Normalize first: Boost 1.92 reports its embedded estimate on the internal
    // [-1,1] interval. Including the Jacobian here gives it physical integral units,
    // even for tiny clock-inversion brackets.
    const double half = (b - a) / 2, middle = a + half;
    if (!std::isfinite(half) || !std::isfinite(middle))
        throw std::domain_error("Analytic quadrature interval overflow");
    const double value = boost::math::quadrature::gauss_kronrod<double, 31>::integrate(
        [&](double z) {
            const double y = f(std::fma(half, z, middle)) * half;
            if (!std::isfinite(y))
                throw std::domain_error("Nonfinite analytic integrand");
            return y;
        },
        -1.0, 1.0, 14, 2e-12, &error, &l1);
    if (!std::isfinite(value) || error > 2e-10 * l1)
        throw std::runtime_error("Analytic quadrature did not converge");
    return value;
}
PeriodicPrimitive::PeriodicPrimitive(std::function<double(double)> rate, double period,
                                     bool positive)
    : rate_([rate, positive](double x) {
          const double y = rate(x);
          if (!std::isfinite(y) || (positive && !(y > 0)))
              throw std::domain_error("Analytic clock must have a finite positive rate");
          return y;
      }),
      period_(period), positive_(positive) {
    if (!(period > 0) || !std::isfinite(period))
        throw std::invalid_argument("Invalid analytic period");
    constexpr std::size_t cells = 512;
    table_.resize(cells + 1);
    for (std::size_t i = 1; i <= cells; ++i) {
        table_[i] =
            table_[i - 1] + analyticIntegral(rate_, period * (i - 1) / cells, period * i / cells);
        if (!std::isfinite(table_[i]) || (positive && !(table_[i] > table_[i - 1])))
            throw std::domain_error("Analytic clock is not numerically resolvable");
    }
}
double PeriodicPrimitive::value(double parameter) const {
    if (!std::isfinite(parameter) || std::abs(parameter / period_) > 1e12)
        throw std::domain_error("Analytic phase exceeds the winding resolution limit");
    const double cycles = std::floor(parameter / period_);
    const double phase = std::clamp(parameter - cycles * period_, 0.0, period_);
    const auto cells = table_.size() - 1;
    const auto cell = std::min(cells - 1, static_cast<std::size_t>(phase / period_ * cells));
    return cycles * table_.back() + table_[cell] +
           analyticIntegral(rate_, period_ * cell / cells, phase);
}
double PeriodicPrimitive::inverse(double value) const {
    if (!positive_)
        throw std::logic_error("Inversion requires a positive analytic rate");
    const double total = table_.back();
    if (!std::isfinite(value) || std::abs(value / total) > 1e12)
        throw std::domain_error("Analytic time exceeds the winding resolution limit");
    const double cycles = std::floor(value / total),
                 phase = std::clamp(value - cycles * total, 0.0, total);
    const auto cells = table_.size() - 1;
    const auto it = std::upper_bound(table_.begin(), table_.end(), phase);
    const auto cell = std::min(cells - 1, static_cast<std::size_t>(it - table_.begin() - 1));
    const double start = period_ * cell / cells;
    double lo = start, hi = period_ * (cell + 1) / cells;
    double s = lo + (hi - lo) * (phase - table_[cell]) / (table_[cell + 1] - table_[cell]);
    for (int k = 0; k < 40; ++k) {
        const double residual = table_[cell] + analyticIntegral(rate_, start, s) - phase;
        if (std::abs(residual) <= 2e-13 * total)
            return cycles * period_ + s;
        if (residual > 0)
            hi = s;
        else
            lo = s;
        const double next = s - residual / rate_(s);
        s = next > lo && next < hi ? next : (lo + hi) / 2;
    }
    throw std::runtime_error("Analytic clock inversion did not converge");
}
} // namespace dynsyslab::core::detail

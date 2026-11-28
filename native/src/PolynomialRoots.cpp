// Real polynomial roots by the Aberth-Ehrlich simultaneous iteration (Aberth 1973; Bini, Numer.
// Algorithms 13, 1996) with Newton polishing (PolynomialRoots.hpp).
#include <dynsyslab/core/PolynomialRoots.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace dynsyslab::core {
namespace {
using C = std::complex<long double>;

C evaluate(const std::vector<long double>& c, C x, C* derivative) {
    C p = c[0], d = 0;
    for (std::size_t i = 1; i < c.size(); ++i) { d = d * x + p; p = p * x + c[i]; }
    if (derivative) *derivative = d;
    return p;
}

}  // namespace

std::vector<double> realPolynomialRoots(const std::vector<double>& coefficients) {
    std::vector<long double> c(coefficients.begin(), coefficients.end());
    while (!c.empty() && c.front() == 0) c.erase(c.begin());
    for (const auto v : c)
        if (!std::isfinite(double(v))) return {};
    const std::size_t degree = c.empty() ? 0 : c.size() - 1;
    if (degree == 0) return {};
    const long double lead = c.front();
    for (auto& v : c) v /= lead;
    // Initial guesses on a circle of the Cauchy bound's radius, rotated off the real axis.
    long double radius = 0;
    for (std::size_t i = 1; i < c.size(); ++i) radius = std::max(radius, std::abs(c[i]));
    radius += 1;
    std::vector<C> z(degree);
    const long double pi = 3.14159265358979323846264338327950288L;
    for (std::size_t k = 0; k < degree; ++k) z[k] = std::polar(radius, 2 * pi * (long double(k) + 0.25L) / long double(degree));
    for (int iteration = 0; iteration < 500; ++iteration) {
        long double change = 0;
        for (std::size_t k = 0; k < degree; ++k) {
            C d;
            const C p = evaluate(c, z[k], &d);
            if (p == C(0)) continue;
            const C ratio = p / d;
            C sum = 0;
            for (std::size_t j = 0; j < degree; ++j)
                if (j != k) sum += C(1) / (z[k] - z[j]);
            const C step = ratio / (C(1) - ratio * sum);
            z[k] -= step;
            change = std::max(change, std::abs(step) / std::max(1.0L, std::abs(z[k])));
        }
        if (change < 1e-18L) break;
    }
    std::vector<double> real;
    for (const auto& root : z) {
        if (std::abs(root.imag()) > 1e-8L * std::max(1.0L, std::abs(root))) continue;
        long double x = root.real();
        for (int step = 0; step < 8; ++step) {  // Newton on the real polynomial
            C d;
            const C p = evaluate(c, C(x), &d);
            if (d.real() == 0) break;
            const long double next = x - p.real() / d.real();
            if (next == x) break;
            x = next;
        }
        real.push_back(double(x));
    }
    std::sort(real.begin(), real.end(), std::greater<double>());
    return real;
}

double realPolynomialRoot(std::size_t k, const std::vector<double>& coefficients) {
    const auto roots = realPolynomialRoots(coefficients);
    return k < roots.size() ? roots[k] : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace dynsyslab::core

// Finite-difference weights on arbitrary nodes and high-order central differences
// (NumericalDifferentiation.hpp). Fornberg's recursion: B. Fornberg, Generation of finite difference
// formulas on arbitrarily spaced grids, Math. Comp. 51:699-706, 1988.
#include <dynsyslab/core/NumericalDifferentiation.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<double> differenceWeights(const std::vector<double>& x, std::size_t m) {
    const auto n = x.size();
    if (m >= n) throw std::invalid_argument("difference weights need more nodes than the derivative order");
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (x[i] == x[j]) throw std::invalid_argument("difference weights need distinct nodes");
    // c[j][k]: weight of node j for the k-th derivative at 0 (Fornberg 1988, Table 1 algorithm).
    std::vector<std::vector<double>> c(n, std::vector<double>(m + 1, 0));
    double c1 = 1, c4 = x[0];
    c[0][0] = 1;
    for (std::size_t i = 1; i < n; ++i) {
        const auto mn = std::min(i, m);
        double c2 = 1;
        const double c5 = c4;
        c4 = x[i];
        for (std::size_t j = 0; j < i; ++j) {
            const double c3 = x[i] - x[j];
            c2 *= c3;
            if (j == i - 1)
                for (std::size_t k = mn; k >= 1; --k) c[i][k] = c1 * (double(k) * c[i - 1][k - 1] - c5 * c[i - 1][k]) / c2;
            if (j == i - 1) c[i][0] = -c1 * c5 * c[i - 1][0] / c2;
            for (std::size_t k = mn; k >= 1; --k) c[j][k] = (c4 * c[j][k] - double(k) * c[j][k - 1]) / c3;
            c[j][0] = c4 * c[j][0] / c3;
        }
        c1 = c2;
    }
    std::vector<double> out(n);
    for (std::size_t j = 0; j < n; ++j) out[j] = c[j][m];
    return out;
}

std::vector<double> highOrderCentralDifference(const std::vector<double>& y, double dt, std::size_t points) {
    if (points < 3 || points % 2 == 0) throw std::invalid_argument("central_high_order needs an odd stencil of at least 3 points");
    const auto n = y.size();
    if (n < points) throw std::invalid_argument("central_high_order needs at least as many samples as stencil points");
    const auto half = points / 2;
    std::vector<std::vector<double>> stencils(points);  // stencil s covers samples [s0, s0 + points)
    for (std::size_t shift = 0; shift < points; ++shift) {  // evaluation node index `shift` within the stencil
        std::vector<double> nodes(points);
        for (std::size_t j = 0; j < points; ++j) nodes[j] = double(j) - double(shift);
        stencils[shift] = differenceWeights(nodes, 1);
    }
    std::vector<double> out(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t start = k < half ? 0 : (k + half >= n ? n - points : k - half);
        const auto& w = stencils[k - start];
        double s = 0;
        for (std::size_t j = 0; j < points; ++j) s += w[j] * y[start + j];
        out[k] = s / dt;
    }
    return out;
}

}  // namespace dynsyslab::core

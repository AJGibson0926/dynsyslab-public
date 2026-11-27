// Affine equations and their exact propagation (AffineDynamics.hpp). Ported from the
// miscellaneous-computes side-car (integration.affine_equations, exact_affine).
#include <dynsyslab/core/AffineDynamics.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

AffineEquations affineEquations(const Observable& field, double scale, double tolerance) {
    const auto n = field.stateDimension();
    if (field.outputDimension() != n) throw std::invalid_argument("affine equations need a vector field (outputs = states)");
    AffineEquations out;
    out.dimension = n;
    const std::vector<double> zero(n, 0);
    out.c = field.evaluateReal(zero.data());
    out.A.assign(n * n, 0);
    for (std::size_t j = 0; j < n; ++j) {  // exact columns for an affine field: f(s e_j) - f(0) over s
        auto x = zero;
        x[j] = scale;
        const auto f = field.evaluateReal(x.data());
        for (std::size_t i = 0; i < n; ++i) out.A[i * n + j] = (f[i] - out.c[i]) / scale;
    }
    // Verify on probe points mixing every coordinate: any active nonlinear term shows here.
    double size = 0, worst = 0;
    for (int k = 1; k <= 8; ++k) {
        std::vector<double> x(n);
        for (std::size_t j = 0; j < n; ++j) x[j] = scale * std::sin(1.3 * double(k) + 0.7 * double(j) + 0.1 * double(k * j));
        const auto f = field.evaluateReal(x.data());
        for (std::size_t i = 0; i < n; ++i) {
            double affine = out.c[i];
            for (std::size_t j = 0; j < n; ++j) affine += out.A[i * n + j] * x[j];
            worst = std::max(worst, std::abs(f[i] - affine));
            size = std::max(size, std::abs(f[i]));
        }
    }
    if (worst > tolerance * std::max(1.0, size)) throw std::invalid_argument("the field has active nonlinear terms; it has no A, c form");
    return out;
}

std::vector<double> affinePropagation(const AffineEquations& e, const std::vector<double>& x0, const std::vector<double>& t) {
    const auto n = e.dimension, m = n + 1;
    if (x0.size() != n || t.empty()) throw std::invalid_argument("affine propagation: initial state and clock");
    for (std::size_t k = 1; k < t.size(); ++k) if (!(t[k] > t[k - 1])) throw std::invalid_argument("affine propagation: the clock must increase");
    std::vector<double> M(m * m, 0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) M[i * m + j] = e.A[i * n + j];
        M[i * m + n] = e.c[i];
    }
    std::vector<double> state(t.size() * m, 0);
    for (std::size_t i = 0; i < n; ++i) state[i] = x0[i];
    state[n] = 1;
    const auto apply = [&](const std::vector<double>& P, std::size_t from, std::size_t to) {
        for (std::size_t i = 0; i < m; ++i) {
            double s = 0;
            for (std::size_t j = 0; j < m; ++j) s += P[i * m + j] * state[from * m + j];
            state[to * m + i] = s;
        }
    };
    const double dt = t.size() > 1 ? (t.back() - t.front()) / double(t.size() - 1) : 0;
    bool uniform = true;
    for (std::size_t k = 0; k < t.size(); ++k) uniform = uniform && std::abs(t[k] - t[0] - double(k) * dt) <= 1e-9 * std::max(dt, 1e-300);
    if (uniform && t.size() > 1) {
        auto step = numerics::matrixExponential(M, m, dt);
        for (std::size_t k = 1; k < t.size(); k *= 2) {  // states k..2k-1 from 0..k-1 with exp(k dt M)
            for (std::size_t j = 0; j < k && k + j < t.size(); ++j) apply(step, j, k + j);
            if (2 * k < t.size()) {
                std::vector<double> square(m * m, 0);
                for (std::size_t a = 0; a < m; ++a)
                    for (std::size_t b = 0; b < m; ++b)
                        for (std::size_t c = 0; c < m; ++c) square[a * m + b] += step[a * m + c] * step[c * m + b];
                step = std::move(square);
            }
        }
    } else {
        for (std::size_t k = 1; k < t.size(); ++k) apply(numerics::matrixExponential(M, m, t[k] - t[0]), 0, k);
    }
    std::vector<double> out(t.size() * n);
    for (std::size_t k = 0; k < t.size(); ++k) std::copy(state.begin() + std::ptrdiff_t(k * m), state.begin() + std::ptrdiff_t(k * m + n), out.begin() + std::ptrdiff_t(k * n));
    return out;
}

}  // namespace dynsyslab::core

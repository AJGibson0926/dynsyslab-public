// Spectral submanifolds by the parameterization method (SpectralSubmanifold.hpp).
// Haller and Ponsioen 2016, Nonlinear Dyn. 86:1493-1534 (existence, spectral quotient);
// Cabre, Fontich and de la Llave 2003, Indiana Univ. Math. J. 52:283-328 (parameterization method);
// Jain and Haller 2022, Nonlinear Dyn. 107:1417-1450 (normal-form and graph styles, cohomological
// equations solved order by order).
#include <dynsyslab/core/SpectralSubmanifold.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>

#include "PowerSeries.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using detail::Complex;
using detail::MonomialTable;
using detail::Series;

std::vector<Complex> inverse(const std::vector<Complex>& matrix, std::size_t n) {
    std::vector<Complex> identity(n * n, 0);
    for (std::size_t i = 0; i < n; ++i) identity[i * n + i] = 1;
    return numerics::solveLeastSquares(matrix, n, n, identity, n);
}

// The values p^m of every monomial of the table at p.
std::vector<Complex> monomialValues(const MonomialTable& table, const Complex* p) {
    std::vector<Complex> out(table.size());
    for (std::size_t i = 0; i < table.size(); ++i) {
        Complex v = 1;
        for (std::size_t j = 0; j < table.variables(); ++j)
            for (std::size_t e = 0; e < table.exponents(i)[j]; ++e) v *= p[j];
        out[i] = v;
    }
    return out;
}

// d/dp_j of a series.
Series derivative(const Series& s, std::size_t j) {
    const auto& t = s.table();
    Series out(s.shared());
    std::vector<std::size_t> e(t.variables());
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!t.exponents(i)[j] || s.c[i] == Complex(0)) continue;
        std::copy(t.exponents(i), t.exponents(i) + t.variables(), e.begin());
        --e[j];
        out.c[t.index(e.data())] += double(t.exponents(i)[j]) * s.c[i];
    }
    return out;
}

std::shared_ptr<const MonomialTable> tableOf(const SpectralSubmanifold& s) {
    return std::make_shared<const MonomialTable>(s.dimension, s.order);
}

// Columns of `coefficients` (monomials x outputs, in p) as series in q with p = chart q.
std::vector<Series> inChart(const SpectralSubmanifold& s, const std::vector<Complex>& coefficients, std::size_t outputs) {
    const auto table = tableOf(s);
    const auto d = s.dimension;
    std::vector<Series> p;
    for (std::size_t j = 0; j < d; ++j) {
        Series pj(table);
        for (std::size_t l = 0; l < d; ++l) pj.c[1 + l] = s.chart[j * d + l];
        p.push_back(pj);
    }
    std::vector<Series> out(outputs, Series(table));
    for (std::size_t i = 0; i < table->size(); ++i) {
        Series term(table, 1);
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t e = 0; e < table->exponents(i)[j]; ++e) term = term * p[j];
        for (std::size_t r = 0; r < outputs; ++r)
            if (coefficients[i * outputs + r] != Complex(0)) out[r] = out[r] + term * coefficients[i * outputs + r];
    }
    return out;
}

// Real coefficients over polynomialExponentTable(d, order, true), features x outputs.
std::vector<double> canonical(const std::vector<Series>& series, std::size_t d, std::size_t order) {
    const auto exponents = polynomialExponentTable(d, order, true);
    const auto features = exponents.size() / d;
    std::vector<double> out(features * series.size());
    for (std::size_t f = 0; f < features; ++f) {
        const auto index = series.front().table().index(exponents.data() + f * d);
        for (std::size_t r = 0; r < series.size(); ++r) out[f * series.size() + r] = series[r].c[index].real();
    }
    return out;
}

}  // namespace

std::vector<Complex> SpectralSubmanifold::parameterization(const Complex* p) const {
    const auto table = tableOf(*this);
    const auto v = monomialValues(*table, p);
    std::vector<Complex> x(equilibrium.begin(), equilibrium.end());
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t r = 0; r < stateDimension; ++r) x[r] += W[i * stateDimension + r] * v[i];
    return x;
}

std::vector<Complex> SpectralSubmanifold::reducedDynamics(const Complex* p) const {
    const auto table = tableOf(*this);
    const auto v = monomialValues(*table, p);
    std::vector<Complex> out(dimension, 0);
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t r = 0; r < dimension; ++r) out[r] += R[i * dimension + r] * v[i];
    return out;
}

std::vector<double> SpectralSubmanifold::realParameterization() const {
    auto x = inChart(*this, W, stateDimension);
    for (std::size_t r = 0; r < stateDimension; ++r) x[r].c[0] += equilibrium[r];
    return canonical(x, dimension, order);
}

std::vector<double> SpectralSubmanifold::realReducedDynamics() const {
    const auto pdot = inChart(*this, R, dimension);  // p' as series in q
    const auto back = inverse(chart, dimension);       // q' = chart^{-1} p'
    std::vector<Series> qdot(dimension, Series(pdot.front().shared()));
    for (std::size_t a = 0; a < dimension; ++a)
        for (std::size_t j = 0; j < dimension; ++j) qdot[a] = qdot[a] + pdot[j] * back[a * dimension + j];
    return canonical(qdot, dimension, order);
}

double SpectralSubmanifold::invarianceError(const std::vector<AnalyticExpression>& field, double radius) const {
    const auto table = tableOf(*this);
    const auto n = stateDimension, d = dimension;
    double worst = 0, scale = 0;
    for (std::size_t k = 0; k < 64; ++k) {
        std::vector<double> q(d);  // a deterministic spread of directions on the sphere |q| = radius
        double norm = 0;
        for (std::size_t j = 0; j < d; ++j) { q[j] = std::cos(2.399963 * double(k * (j + 1)) + double(j)); norm += q[j] * q[j]; }
        std::vector<Complex> p(d, 0);
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t l = 0; l < d; ++l) p[j] += chart[j * d + l] * (radius * q[l] / std::sqrt(norm));
        const auto x = parameterization(p.data());
        std::vector<double> xr(n);
        for (std::size_t i = 0; i < n; ++i) xr[i] = x[i].real();
        const auto pdot = reducedDynamics(p.data());
        std::vector<Complex> tangent(n, 0);  // DW(p) R(p)
        for (std::size_t i = 0; i < table->size(); ++i)
            for (std::size_t j = 0; j < d; ++j) {
                const auto e = table->exponents(i)[j];
                if (!e) continue;
                Complex v = double(e) * pdot[j];
                for (std::size_t l = 0; l < d; ++l)
                    for (std::size_t m = 0; m < table->exponents(i)[l] - (l == j ? 1 : 0); ++m) v *= p[l];
                for (std::size_t r = 0; r < n; ++r) tangent[r] += W[i * n + r] * v;
            }
        for (std::size_t r = 0; r < n; ++r) {
            const double f = field[r].evaluate(xr.data());
            worst = std::max(worst, std::abs(f - tangent[r]));
            scale = std::max(scale, std::abs(f));
        }
    }
    return scale > 0 ? worst / scale : worst;
}

std::vector<double> equilibriumNear(const std::vector<AnalyticExpression>& field, std::vector<double> x) {
    const auto n = field.size();
    if (x.size() != n) throw std::invalid_argument("equilibrium: the guess must have one value per state");
    std::vector<std::vector<AnalyticExpression>> jacobian(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) jacobian[i].push_back(field[i].derivative(j));
    for (int iteration = 0; iteration < 60; ++iteration) {
        std::vector<double> f(n), J(n * n);
        double size = 0;
        for (std::size_t i = 0; i < n; ++i) {
            f[i] = field[i].evaluate(x.data());
            size = std::max(size, std::abs(f[i]));
            for (std::size_t j = 0; j < n; ++j) J[i * n + j] = jacobian[i][j].evaluate(x.data());
        }
        if (size < 1e-13) return x;
        const auto step = numerics::solveLeastSquares(J, n, n, f, 1);
        double change = 0;
        for (std::size_t i = 0; i < n; ++i) { x[i] -= step[i]; change = std::max(change, std::abs(step[i])); }
        if (!std::isfinite(change)) break;
        if (change < 1e-15 * (1 + *std::max_element(x.begin(), x.end(), [](double a, double b) { return std::abs(a) < std::abs(b); })))
            return x;
    }
    throw std::runtime_error("equilibrium: Newton's method did not converge");
}

std::vector<AnalyticExpression> polynomialField(std::size_t d, std::size_t degree, const std::vector<double>& c) {
    const auto exponents = polynomialExponentTable(d, degree, true);
    const auto features = exponents.size() / d;
    if (c.size() != features * d) throw std::invalid_argument("polynomialField: coefficients must be features x dimension");
    std::vector<AnalyticExpression> out(d, AnalyticExpression::constant(0, d));
    for (std::size_t f = 0; f < features; ++f) {
        auto monomial = AnalyticExpression::constant(1, d);
        for (std::size_t v = 0; v < d; ++v)
            if (exponents[f * d + v]) monomial = monomial * AnalyticExpression::state(d, v).power(double(exponents[f * d + v]));
        for (std::size_t r = 0; r < d; ++r)
            if (c[f * d + r] != 0) out[r] = out[r] + monomial * c[f * d + r];
    }
    return out;
}

SpectralSubmanifold spectralSubmanifold(const std::vector<AnalyticExpression>& field, const std::vector<double>& x0,
                                        const SpectralSubmanifoldOptions& options) {
    const auto n = field.size();
    if (!n || x0.size() != n) throw std::invalid_argument("SSM: one expression and one equilibrium value per state");
    for (const auto& f : field) if (f.stateDimension() != n) throw std::invalid_argument("SSM: expressions must be functions of the state");
    if (options.order < 1 || options.order > 40) throw std::invalid_argument("SSM: order must be between 1 and 40");
    SpectralSubmanifold out;
    out.stateDimension = n;
    out.order = options.order;
    out.equilibrium = x0;
    double residual = 0, magnitude = 0;
    std::vector<double> A(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        residual = std::max(residual, std::abs(field[i].evaluate(x0.data())));
        for (std::size_t j = 0; j < n; ++j) {
            A[i * n + j] = field[i].derivative(j).evaluate(x0.data());
            magnitude = std::max(magnitude, std::abs(A[i * n + j]));
        }
    }
    if (residual > 1e-8 * std::max(1.0, magnitude)) throw std::invalid_argument("SSM: the point is not an equilibrium");

    // Eigenvalues by decreasing real part; unit eigenvectors, conjugate pairs exactly conjugate.
    const auto eig = numerics::eigensystem(A, n);
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](auto a, auto b) {
        const auto la = eig.eigenvalues[a], lb = eig.eigenvalues[b];
        return std::abs(la.real() - lb.real()) > 1e-12 * (1 + magnitude) ? la.real() > lb.real() : la.imag() > lb.imag();
    });
    std::vector<Complex> lambda(n), V(n * n);
    for (std::size_t k = 0; k < n; ++k) {
        lambda[k] = eig.eigenvalues[order[k]];
        double norm = 0, top = -1;
        Complex phase = 1;
        for (std::size_t i = 0; i < n; ++i) {
            const auto v = eig.rightEigenvectors[i * n + order[k]];
            norm += std::norm(v);
            if (std::abs(v) > top + 1e-12) { top = std::abs(v); phase = std::conj(v) / std::abs(v); }
        }
        for (std::size_t i = 0; i < n; ++i) V[i * n + k] = eig.rightEigenvectors[i * n + order[k]] * phase / std::sqrt(norm);
    }
    const double tiny = 1e-9 * (1 + magnitude);
    std::vector<std::size_t> partner(n, n);
    for (std::size_t k = 0; k < n; ++k) {
        if (lambda[k].imag() <= tiny || partner[k] != n) continue;
        std::size_t best = n;
        for (std::size_t l = 0; l < n; ++l)
            if (l != k && partner[l] == n && lambda[l].imag() < -tiny
                && (best == n || std::abs(lambda[l] - std::conj(lambda[k])) < std::abs(lambda[best] - std::conj(lambda[k])))) best = l;
        if (best == n) continue;
        partner[k] = best; partner[best] = k;
        lambda[best] = std::conj(lambda[k]);
        for (std::size_t i = 0; i < n; ++i) V[i * n + best] = std::conj(V[i * n + k]);
    }
    out.eigenvalues = lambda;
    const auto U = inverse(V, n);  // rows: left eigenvectors with u_l^* v_k = delta

    // Master modes, each complex one followed by its conjugate.
    auto requested = options.masterModes;
    if (requested.empty()) requested = {0};
    for (const auto m : requested) {
        if (m >= n) throw std::invalid_argument("SSM: master mode index out of range");
        if (std::find(out.masterModes.begin(), out.masterModes.end(), m) != out.masterModes.end()) continue;
        if (std::abs(lambda[m].imag()) > tiny) {
            if (partner[m] == n) throw std::invalid_argument("SSM: a complex master mode without its conjugate");
            const auto first = lambda[m].imag() > 0 ? m : partner[m];
            out.masterModes.push_back(first);
            out.masterModes.push_back(partner[first]);
        } else {
            out.masterModes.push_back(m);
        }
    }
    const auto d = out.dimension = out.masterModes.size();
    out.chart.assign(d * d, 0);
    for (std::size_t j = 0; j < d; ++j) {
        const auto m = out.masterModes[j];
        if (std::abs(lambda[m].imag()) > tiny && lambda[m].imag() > 0) {
            out.chart[j * d + j] = 1; out.chart[j * d + j + 1] = Complex(0, 1);
            out.chart[(j + 1) * d + j] = 1; out.chart[(j + 1) * d + j + 1] = Complex(0, -1);
            ++j;
        } else {
            out.chart[j * d + j] = 1;
        }
    }
    double slowest = -std::numeric_limits<double>::infinity(), fastest = std::numeric_limits<double>::infinity(), scale = 0;
    for (std::size_t k = 0; k < n; ++k) {
        const bool master = std::find(out.masterModes.begin(), out.masterModes.end(), k) != out.masterModes.end();
        if (master) { slowest = std::max(slowest, lambda[k].real()); scale = std::max(scale, std::abs(lambda[k])); }
        else fastest = std::min(fastest, lambda[k].real());
    }
    if (!(slowest < 0)) {
        out.spectralQuotient = std::numeric_limits<double>::quiet_NaN();
        out.warnings.push_back("the master spectrum is not strictly stable: the manifold is a center-type (Lyapunov "
                               "subcenter) manifold, not the unique SSM of Haller and Ponsioen");
    } else {
        out.spectralQuotient = std::isfinite(fastest) ? std::floor(fastest / slowest) : 0;
        if (double(options.order) < out.spectralQuotient + 1)
            out.warnings.push_back("order below the spectral quotient + 1: the expansion is not of the unique smoothest SSM");
    }

    // The order-by-order solution.
    const auto table = std::make_shared<const MonomialTable>(d, options.order);
    std::vector<Series> Ws(n, Series(table)), Rs(d, Series(table));
    for (std::size_t j = 0; j < d; ++j) {
        for (std::size_t i = 0; i < n; ++i) Ws[i].c[1 + j] = V[i * n + out.masterModes[j]];
        Rs[j].c[1 + j] = lambda[out.masterModes[j]];
    }
    std::vector<json::Value> trees;
    for (const auto& f : field) trees.push_back(f.toJson());
    const double tolerance = options.resonanceTolerance * std::max(scale, 1e-300);
    for (std::size_t k = 2; k <= options.order; ++k) {
        std::vector<Series> x;
        for (std::size_t i = 0; i < n; ++i) x.push_back(Ws[i] + Series(table, x0[i]));
        std::vector<Series> F;
        for (const auto& tree : trees) F.push_back(detail::evaluate(tree, x));
        std::vector<Series> C(n, Series(table));  // [DW R] from the orders below k
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t i = 0; i < n; ++i) C[i] = C[i] + derivative(Ws[i], j) * Rs[j];
        for (std::size_t m = table->begin(k); m < table->end(k); ++m) {
            Complex mu = 0;
            for (std::size_t j = 0; j < d; ++j) mu += double(table->exponents(m)[j]) * lambda[out.masterModes[j]];
            std::vector<Complex> y(n, 0);  // U (C_m - N_m)
            for (std::size_t l = 0; l < n; ++l)
                for (std::size_t i = 0; i < n; ++i) y[l] += U[l * n + i] * (C[i].c[m] - F[i].c[m]);
            for (std::size_t l = 0; l < n; ++l) {
                const auto position = std::find(out.masterModes.begin(), out.masterModes.end(), l) - out.masterModes.begin();
                const bool master = std::size_t(position) < d;
                const auto gap = lambda[l] - mu;
                bool resonant = std::abs(gap) <= std::max(tolerance, 1e-13 * (1 + magnitude));
                if (master && options.style == SsmStyle::graph) resonant = true;
                if (master && options.imaginaryResonance && std::abs(lambda[l].imag()) > tiny
                    && std::abs(gap.imag()) <= tolerance) resonant = true;
                if (master && resonant) { Rs[std::size_t(position)].c[m] = -y[l]; y[l] = 0; continue; }
                if (!master && resonant) {
                    out.warnings.push_back("outer resonance at order " + std::to_string(k) + " with eigenvalue "
                                           + std::to_string(l) + ": the SSM does not exist at this order");
                    y[l] = 0;
                    continue;
                }
                y[l] /= gap;
            }
            for (std::size_t i = 0; i < n; ++i) {
                Complex w = 0;
                for (std::size_t l = 0; l < n; ++l) w += V[i * n + l] * y[l];
                Ws[i].c[m] = w;
            }
        }
    }
    out.exponents.reserve(table->size() * d);
    for (std::size_t i = 0; i < table->size(); ++i) out.exponents.insert(out.exponents.end(), table->exponents(i), table->exponents(i) + d);
    out.W.resize(table->size() * n);
    out.R.resize(table->size() * d);
    for (std::size_t i = 0; i < table->size(); ++i) {
        for (std::size_t r = 0; r < n; ++r) out.W[i * n + r] = Ws[r].c[i];
        for (std::size_t r = 0; r < d; ++r) out.R[i * d + r] = Rs[r].c[i];
    }
    return out;
}

}  // namespace dynsyslab::core

// Conserved quantities from run identity (RunInvariant.hpp).
#include <dynsyslab/core/RunInvariant.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::vector<double> features(const std::vector<std::size_t>& exponents, std::size_t d, const double* x) {
    std::vector<double> out(exponents.size() / d, 1);
    for (std::size_t f = 0; f < out.size(); ++f)
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t e = 0; e < exponents[f * d + j]; ++e) out[f] *= x[j];
    return out;
}

}  // namespace

double RunInvariant::evaluate(const double* x) const {
    const auto theta = features(polynomialExponentTable(dimension, degree, true), dimension, x);
    double s = 0;
    for (std::size_t f = 0; f < theta.size(); ++f) s += coefficients[f] * theta[f];
    return s;
}

RunInvariant fitRunInvariant(const std::vector<TimeSeries>& runs, std::size_t degree, std::size_t stride) {
    if (runs.size() < 2) throw std::invalid_argument("run invariant: at least two runs are needed");
    const auto d = runs.front().componentCount();
    if (!d || !degree || !stride) throw std::invalid_argument("run invariant: dimension, degree and stride must be positive");
    const auto exponents = polynomialExponentTable(d, degree, true);
    const auto F = exponents.size() / d, m = F - 1;  // feature 0 is the constant
    std::vector<std::vector<double>> means;
    std::vector<double> global(m, 0), within(m * m, 0), between(m * m, 0);
    std::size_t total = 0;
    for (const auto& run : runs) {
        if (run.componentCount() != d || run.sampleCount() < 2) throw std::invalid_argument("run invariant: runs differ in dimension or are too short");
        std::vector<std::vector<double>> theta;
        std::vector<double> mean(m, 0);
        for (std::size_t k = 0; k < run.sampleCount(); ++k) {
            theta.push_back(features(exponents, d, run.valuesData() + k * d));
            for (std::size_t a = 0; a < m; ++a) mean[a] += theta.back()[a + 1] / double(run.sampleCount());
        }
        for (const auto& t : theta)
            for (std::size_t a = 0; a < m; ++a)
                for (std::size_t b = 0; b < m; ++b) within[a * m + b] += (t[a + 1] - mean[a]) * (t[b + 1] - mean[b]);
        for (std::size_t a = 0; a < m; ++a) global[a] += double(run.sampleCount()) * mean[a];
        total += run.sampleCount();
        means.push_back(std::move(mean));
    }
    for (auto& g : global) g /= double(total);
    for (std::size_t r = 0; r < runs.size(); ++r)
        for (std::size_t a = 0; a < m; ++a)
            for (std::size_t b = 0; b < m; ++b)
                between[a * m + b] += double(runs[r].sampleCount()) * (means[r][a] - global[a]) * (means[r][b] - global[b]);
    for (auto& v : within) v /= double(total);
    for (auto& v : between) v /= double(total);
    double trace = 0;
    for (std::size_t a = 0; a < m; ++a) trace += within[a * m + a];
    for (std::size_t a = 0; a < m; ++a) within[a * m + a] += 1e-8 * trace / double(m);
    // Cholesky W = L L^T; the symmetric problem L^-1 B L^-T y = rho y with c = L^-T y.
    std::vector<double> L(m * m, 0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double s = within[i * m + j];
            for (std::size_t k = 0; k < j; ++k) s -= L[i * m + k] * L[j * m + k];
            if (i == j) {
                if (!(s > 0)) throw std::runtime_error("run invariant: the within-run scatter is not positive definite");
                L[i * m + i] = std::sqrt(s);
            } else {
                L[i * m + j] = s / L[j * m + j];
            }
        }
    const auto solveLower = [&](std::vector<double> v) {  // L^-1 v
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t k = 0; k < i; ++k) v[i] -= L[i * m + k] * v[k];
            v[i] /= L[i * m + i];
        }
        return v;
    };
    std::vector<double> half(m * m), C(m * m);  // half = L^-1 B, C = L^-1 half^T = L^-1 B L^-T
    for (std::size_t j = 0; j < m; ++j) {
        std::vector<double> column(m);
        for (std::size_t i = 0; i < m; ++i) column[i] = between[i * m + j];
        const auto y = solveLower(column);
        for (std::size_t i = 0; i < m; ++i) half[i * m + j] = y[i];
    }
    for (std::size_t j = 0; j < m; ++j) {
        std::vector<double> row(m);
        for (std::size_t i = 0; i < m; ++i) row[i] = half[j * m + i];
        const auto y = solveLower(row);
        for (std::size_t i = 0; i < m; ++i) C[i * m + j] = y[i];
    }
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < i; ++j) { const double s = (C[i * m + j] + C[j * m + i]) / 2; C[i * m + j] = C[j * m + i] = s; }
    const auto eig = numerics::symmetricEigensystem(C, m);
    const double rho = eig.eigenvalues.back();
    if (!(rho > 0)) throw std::runtime_error("run invariant: no feature separates the runs");
    std::vector<double> c(m, 0), y(m);
    for (std::size_t i = 0; i < m; ++i) y[i] = eig.eigenvectors[i * m + (m - 1)];
    for (std::size_t i = m; i-- > 0;) {  // L^T c = y
        double s = y[i];
        for (std::size_t k = i + 1; k < m; ++k) s -= L[k * m + i] * c[k];
        c[i] = s / L[i * m + i];
    }
    const auto pivot = std::size_t(std::max_element(c.begin(), c.end(), [](double a, double b) { return std::abs(a) < std::abs(b); }) - c.begin());
    if (c[pivot] < 0) for (auto& v : c) v = -v;
    RunInvariant out;
    out.dimension = d; out.degree = degree; out.separation = rho;
    out.coefficients.assign(F, 0);
    double offset = 0;
    for (std::size_t a = 0; a < m; ++a) { out.coefficients[a + 1] = c[a]; offset += global[a] * c[a]; }
    out.coefficients[0] = -offset;
    for (const auto& run : runs) {
        std::vector<double> v;
        for (std::size_t k = 0; k < run.sampleCount(); ++k) v.push_back(out.evaluate(run.valuesData() + k * d));
        RunInvariantRun r;
        double s = 0, ss = 0, num = 0, den = 0;
        for (const auto x : v) { s += x; ss += x * x; }
        r.mean = s / double(v.size());
        r.standardDeviation = std::sqrt(std::max(0.0, ss / double(v.size()) - r.mean * r.mean));
        const auto [lo, hi] = std::minmax_element(v.begin(), v.end());
        r.peakToPeak = *hi - *lo;
        for (std::size_t k = stride; k < v.size(); ++k) { num += std::pow(v[k] - v[k - stride], 2); den += v[k] * v[k]; }
        r.sampledInvariance = den > 0 ? std::sqrt(num / den) : 0;
        out.runs.push_back(r);
    }
    return out;
}

}  // namespace dynsyslab::core

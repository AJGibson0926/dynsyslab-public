#include <dynsyslab/core/EmpiricalObservability.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();

// Fixed-step RK4 with the input sampled at each stage; returns the output samples
// (samples x outputs) or an empty vector if the trajectory leaves the validity domain.
std::vector<double> outputPath(const ODESystem& system, const Observable& observable, std::vector<double> x,
                               double t0, const EmpiricalObservabilityOptions& o) {
    const std::size_t n = x.size(), p = observable.outputDimension();
    const double h = o.horizon / double((o.samples - 1) * o.substeps);
    std::vector<double> k1(n), k2(n), k3(n), k4(n), tmp(n), out;
    out.reserve(o.samples * p);
    const auto f = [&](double t, const double* s, double* d) {
        if (o.input) system.rhsWithInput(t, s, o.input(t), d); else system.rhs(t, s, d);
    };
    double t = t0;
    for (std::size_t sample = 0; sample < o.samples; ++sample) {
        if (sample) for (std::size_t step = 0; step < o.substeps; ++step) {
            f(t, x.data(), k1.data());
            for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + 0.5 * h * k1[i];
            f(t + 0.5 * h, tmp.data(), k2.data());
            for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + 0.5 * h * k2[i];
            f(t + 0.5 * h, tmp.data(), k3.data());
            for (std::size_t i = 0; i < n; ++i) tmp[i] = x[i] + h * k3[i];
            f(t + h, tmp.data(), k4.data());
            for (std::size_t i = 0; i < n; ++i) x[i] += h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            t += h;
        }
        for (double v : x) if (!std::isfinite(v)) return {};
        if (system.invalid(t, x.data())) return {};
        const auto y = observable.evaluateReal(x.data());
        for (double v : y) if (!std::isfinite(v)) return {};
        out.insert(out.end(), y.begin(), y.end());
    }
    return out;
}
// Eigenvalues descending, eigenvectors in matching columns, condition number.
void decompose(EmpiricalObservability& r) {
    const std::size_t n = r.states;
    const auto eig = numerics::symmetricEigensystem(r.gramian, n);  // ascending
    r.eigenvalues.resize(n); r.eigenvectors.assign(n * n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t src = n - 1 - k;
        r.eigenvalues[k] = eig.eigenvalues[src];
        for (std::size_t i = 0; i < n; ++i) r.eigenvectors[i * n + k] = eig.eigenvectors[i * n + src];
    }
    const double low = r.eigenvalues.back();
    r.conditionNumber = low > 0 ? r.eigenvalues.front() / low : kInf;
}
}  // namespace

EmpiricalObservability empiricalObservability(const ODESystem& system, const Observable& observable,
                                              const std::vector<double>& state, double time,
                                              const EmpiricalObservabilityOptions& o) {
    const std::size_t n = system.descriptor().stateDimension, p = observable.outputDimension();
    if (state.size() != n) throw std::invalid_argument("empirical observability: state has the wrong dimension");
    if (!(o.horizon > 0) || o.samples < 2 || o.substeps < 1)
        throw std::invalid_argument("empirical observability needs a positive horizon, >= 2 samples and >= 1 substep");
    std::vector<double> e = o.perturbation;
    if (e.empty()) for (double v : state) e.push_back(1e-6 * std::max(1.0, std::abs(v)));
    if (e.size() != n) throw std::invalid_argument("empirical observability: one perturbation per state");
    for (double v : e) if (!(v > 0) || !std::isfinite(v)) throw std::invalid_argument("perturbations must be positive");
    EmpiricalObservability r;
    r.states = n; r.outputs = p; r.valid = true;
    std::vector<std::vector<double>> sensitivity(n);
    r.outputJacobian.assign(p * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        auto plus = state, minus = state;
        plus[i] += e[i]; minus[i] -= e[i];
        const auto yp = observable.evaluateReal(plus.data()), ym = observable.evaluateReal(minus.data());
        for (std::size_t k = 0; k < p; ++k) r.outputJacobian[k * n + i] = (yp[k] - ym[k]) / (2 * e[i]);
        const auto a = outputPath(system, observable, plus, time, o), b = outputPath(system, observable, minus, time, o);
        if (a.empty() || b.empty()) { r.valid = false; sensitivity[i].assign(o.samples * p, 0.0); continue; }
        sensitivity[i].resize(a.size());
        for (std::size_t k = 0; k < a.size(); ++k) sensitivity[i][k] = (a[k] - b[k]) / (2 * e[i]);
    }
    const double dt = o.horizon / double(o.samples - 1);
    r.gramian.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i; j < n; ++j) {
            double sum = 0;
            for (std::size_t s = 0; s < o.samples; ++s) {
                double dot = 0;
                for (std::size_t k = 0; k < p; ++k) dot += sensitivity[i][s * p + k] * sensitivity[j][s * p + k];
                sum += (s == 0 || s + 1 == o.samples ? 0.5 : 1.0) * dot;
            }
            r.gramian[i * n + j] = r.gramian[j * n + i] = sum * dt;
        }
    decompose(r);
    return r;
}

EmpiricalObservability scaledObservability(const EmpiricalObservability& result, const std::vector<double>& scale) {
    const std::size_t n = result.states;
    if (scale.size() != n) throw std::invalid_argument("scaledObservability: one scale per state");
    auto r = result;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) r.gramian[i * n + j] *= scale[i] * scale[j];
    for (std::size_t k = 0; k < r.outputs; ++k)
        for (std::size_t i = 0; i < n; ++i) r.outputJacobian[k * n + i] *= scale[i];
    decompose(r);
    return r;
}

std::vector<double> estimationBoundStd(const EmpiricalObservability& r, double V) {
    const std::size_t n = r.states;
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double variance = 0;
        for (std::size_t k = 0; k < n; ++k) {
            const double v = r.eigenvectors[i * n + k];
            if (v == 0.0) continue;
            if (!(r.eigenvalues[k] > 0)) { variance = kInf; break; }
            variance += v * v / r.eigenvalues[k];
        }
        out[i] = std::sqrt(V * variance);
    }
    return out;
}

}  // namespace dynsyslab::core

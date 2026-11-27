// Backbone curves of two-dimensional SSMs, and data-driven SSMs (SpectralSubmanifold.hpp).
// Backbones: Szalai, Ehrhardt and Haller 2017, Proc. R. Soc. A 473:20160759; Breunung and Haller
// 2018, Proc. R. Soc. A 474:20180083. Data-driven SSMs: Cenedese, Axas, Bauerlein, Avila and Haller
// 2022, Nat. Commun. 13:872 (SSMLearn: SVD tangent space, polynomial graph, polynomial reduced
// dynamics, and its normal form).
#include <dynsyslab/core/SpectralSubmanifold.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Complex = std::complex<double>;

std::vector<double> features(const std::vector<std::size_t>& exponents, std::size_t d, const double* q) {
    const auto count = exponents.size() / d;
    std::vector<double> out(count, 1);
    for (std::size_t f = 0; f < count; ++f)
        for (std::size_t v = 0; v < d; ++v)
            for (std::size_t e = 0; e < exponents[f * d + v]; ++e) out[f] *= q[v];
    return out;
}

std::size_t degree(const std::vector<std::size_t>& exponents, std::size_t d, std::size_t f) {
    return std::accumulate(exponents.begin() + std::ptrdiff_t(f * d), exponents.begin() + std::ptrdiff_t((f + 1) * d), std::size_t(0));
}

std::vector<double> apply(const std::vector<std::size_t>& exponents, std::size_t d, const std::vector<double>& coefficients,
                          std::size_t outputs, const double* q) {
    const auto theta = features(exponents, d, q);
    std::vector<double> out(outputs, 0);
    for (std::size_t f = 0; f < theta.size(); ++f)
        for (std::size_t r = 0; r < outputs; ++r) out[r] += theta[f] * coefficients[f * outputs + r];
    return out;
}

// Least squares of Y (rows x outputs) on the columns of Theta (rows x columns) with a ridge.
std::vector<double> regress(const std::vector<double>& theta, std::size_t rows, std::size_t columns,
                            const std::vector<double>& y, std::size_t outputs, double ridge) {
    numerics::LeastSquaresOptions options;
    options.method = numerics::LeastSquaresMethod::divideAndConquerSvd;
    options.ridge = ridge;
    options.singularValueTolerance = 1e-12;
    return numerics::solveLeastSquares(theta, rows, columns, y, outputs, options);
}

double relativeRms(const std::vector<double>& error, const std::vector<double>& reference) {
    double e = 0, r = 0;
    for (std::size_t i = 0; i < error.size(); ++i) { e += error[i] * error[i]; r += reference[i] * reference[i]; }
    return r > 0 ? std::sqrt(e / r) : std::sqrt(e);
}

}  // namespace

BackboneCurve backboneCurve(const SpectralSubmanifold& ssm, const std::vector<double>& radii, const StateOutputs& outputs,
                            std::size_t outputCount, std::size_t angles) {
    if (ssm.dimension != 2 || std::abs(ssm.chart[1].imag()) != 1)
        throw std::invalid_argument("backbone curves need a two-dimensional SSM of a complex pair");
    const auto lambda = ssm.eigenvalues[ssm.masterModes[0]];
    // p' = lambda p + sum_j g_j p^{j+1} conj(p)^j: the other terms must vanish (normal-form style).
    std::vector<Complex> gamma;
    double other = 0, size = std::abs(lambda);
    for (std::size_t i = 0; i < ssm.monomials(); ++i) {
        const auto a = ssm.exponents[2 * i], b = ssm.exponents[2 * i + 1];
        const auto value = ssm.R[i * 2];
        if (a + b < 2) continue;
        if (a == b + 1) { gamma.resize(std::max(gamma.size(), b + 1)); gamma[b] = value; }
        else other = std::max(other, std::abs(value));
        size = std::max(size, std::abs(value));
    }
    if (other > 1e-8 * size) throw std::invalid_argument("backbone curves need the normal-form style reduced dynamics");
    BackboneCurve out;
    out.radius = radii;
    out.outputs = outputCount;
    const std::vector<double> x0(ssm.equilibrium);
    const auto base = outputs(x0.data());
    for (const double rho : radii) {
        double growth = lambda.real(), frequency = lambda.imag(), power = 1;
        for (std::size_t j = 1; j < gamma.size(); ++j) {
            power *= rho * rho;
            growth += gamma[j].real() * power;
            frequency += gamma[j].imag() * power;
        }
        out.frequency.push_back(frequency);
        out.damping.push_back(frequency != 0 ? -growth / std::abs(frequency) : std::numeric_limits<double>::quiet_NaN());
        std::vector<double> amplitude(outputCount, 0), x(ssm.stateDimension);
        for (std::size_t k = 0; k < angles; ++k) {
            const double theta = 2 * 3.141592653589793 * double(k) / double(angles);
            const Complex p[2] = {std::polar(rho, theta), std::polar(rho, -theta)};
            const auto w = ssm.parameterization(p);
            for (std::size_t i = 0; i < x.size(); ++i) x[i] = w[i].real();
            const auto value = outputs(x.data());
            for (std::size_t r = 0; r < outputCount; ++r) amplitude[r] = std::max(amplitude[r], std::abs(value[r] - base[r]));
        }
        out.amplitude.insert(out.amplitude.end(), amplitude.begin(), amplitude.end());
    }
    return out;
}

std::vector<double> SsmFit::reduce(const double* y) const {
    std::vector<double> eta(dimension, 0);
    for (std::size_t i = 0; i < observableDimension; ++i)
        for (std::size_t j = 0; j < dimension; ++j) eta[j] += tangent[i * dimension + j] * (y[i] - equilibrium[i]);
    return eta;
}

std::vector<double> SsmFit::lift(const double* eta) const {
    return apply(polynomialExponentTable(dimension, order, true), dimension, manifold, observableDimension, eta);
}

std::vector<double> SsmFit::rate(const double* eta) const {
    return apply(polynomialExponentTable(dimension, dynamicsOrder, true), dimension, dynamics, dimension, eta);
}

TimeSeries SsmFit::predict(const TimeSeries& trajectory, double maximumStep) const {
    const auto samples = trajectory.sampleCount();
    if (!samples || trajectory.componentCount() != observableDimension) throw std::invalid_argument("SSM prediction: trajectory shape");
    auto eta = reduce(trajectory.valuesData());
    std::vector<double> time(samples), values;
    for (std::size_t k = 0; k < samples; ++k) {
        time[k] = trajectory.time(k);
        if (k) {
            const double span = time[k] - time[k - 1];
            const auto steps = maximumStep > 0 ? std::max<std::size_t>(1, std::size_t(std::ceil(std::abs(span) / maximumStep))) : 4;
            const double h = span / double(steps);
            for (std::size_t s = 0; s < steps; ++s) {  // classical RK4
                const auto k1 = rate(eta.data());
                std::vector<double> y(dimension);
                for (std::size_t j = 0; j < dimension; ++j) y[j] = eta[j] + h / 2 * k1[j];
                const auto k2 = rate(y.data());
                for (std::size_t j = 0; j < dimension; ++j) y[j] = eta[j] + h / 2 * k2[j];
                const auto k3 = rate(y.data());
                for (std::size_t j = 0; j < dimension; ++j) y[j] = eta[j] + h * k3[j];
                const auto k4 = rate(y.data());
                for (std::size_t j = 0; j < dimension; ++j) eta[j] += h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            }
        }
        const auto y = lift(eta.data());
        values.insert(values.end(), y.begin(), y.end());
    }
    return TimeSeries(std::move(time), std::move(values), observableDimension);
}

SsmFit fitSpectralSubmanifold(const std::vector<TimeSeries>& trajectories, const SsmFitOptions& options) {
    if (trajectories.empty()) throw std::invalid_argument("SSM fit: no trajectories");
    const auto n = trajectories.front().componentCount(), d = options.dimension;
    if (!d || d > n || options.order < 1 || options.dynamicsOrder < 1)
        throw std::invalid_argument("SSM fit: need 1 <= dimension <= observables and positive orders");
    SsmFit fit;
    fit.observableDimension = n; fit.dimension = d; fit.order = options.order; fit.dynamicsOrder = options.dynamicsOrder;
    fit.equilibrium = options.equilibrium;
    if (fit.equilibrium.empty()) {
        fit.equilibrium.assign(n, 0);
        for (const auto& t : trajectories)
            for (std::size_t i = 0; i < n; ++i) fit.equilibrium[i] += t.value(t.sampleCount() - 1, i) / double(trajectories.size());
    }
    if (fit.equilibrium.size() != n) throw std::invalid_argument("SSM fit: equilibrium size");
    // Deviations Y and their derivatives, past each trajectory's transient.
    std::vector<double> Y, Ydot;
    for (const auto& t : trajectories) {
        if (t.componentCount() != n) throw std::invalid_argument("SSM fit: trajectories differ in dimension");
        if (t.sampleCount() < 5 || !t.uniformTimeSteps()) throw std::invalid_argument("SSM fit: each trajectory needs >= 5 uniform samples");
        std::size_t first = 0;
        while (first < t.sampleCount() && t.time(first) - t.time(0) < options.skipTime) ++first;
        if (t.sampleCount() - first < 5) continue;
        std::vector<std::vector<double>> rates(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<double> column(t.sampleCount());
            for (std::size_t k = 0; k < column.size(); ++k) column[k] = t.value(k, i);
            rates[i] = differentiate(column, t.deltaTime(), options.derivative).derivative;
        }
        for (std::size_t k = first; k < t.sampleCount(); ++k)
            for (std::size_t i = 0; i < n; ++i) { Y.push_back(t.value(k, i) - fit.equilibrium[i]); Ydot.push_back(rates[i][k]); }
    }
    const auto rows = Y.size() / n;
    // Tangent space: the leading right singular vectors of the data (rows x observables).
    const auto svd = numerics::singularValueDecomposition(Y, rows, n);
    if (svd.rank < d) throw std::invalid_argument("SSM fit: the data span fewer directions than the SSM dimension");
    fit.singularValues = svd.singularValues;
    fit.tangent.assign(n * d, 0);
    for (std::size_t j = 0; j < d; ++j)
        for (std::size_t i = 0; i < n; ++i) fit.tangent[i * d + j] = svd.rightTranspose[j * n + i];
    std::vector<double> eta(rows * d, 0), etaDot(rows * d, 0);
    for (std::size_t k = 0; k < rows; ++k)
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                eta[k * d + j] += fit.tangent[i * d + j] * Y[k * n + i];
                etaDot[k * d + j] += fit.tangent[i * d + j] * Ydot[k * n + i];
            }

    // The graph y - y0 = V eta + w(eta), w of degrees 2..order.
    const auto mExp = polynomialExponentTable(d, options.order, true);
    const auto mFeatures = mExp.size() / d;
    std::vector<std::size_t> nonlinear;
    for (std::size_t f = 0; f < mFeatures; ++f) if (degree(mExp, d, f) >= 2) nonlinear.push_back(f);
    fit.manifold.assign(mFeatures * n, 0);
    for (std::size_t f = 0; f < mFeatures; ++f) {
        const auto deg = degree(mExp, d, f);
        if (deg == 0) for (std::size_t i = 0; i < n; ++i) fit.manifold[f * n + i] = fit.equilibrium[i];
        if (deg == 1) {
            const auto j = std::size_t(std::find(mExp.begin() + std::ptrdiff_t(f * d), mExp.begin() + std::ptrdiff_t((f + 1) * d), 1)
                                       - (mExp.begin() + std::ptrdiff_t(f * d)));
            for (std::size_t i = 0; i < n; ++i) fit.manifold[f * n + i] = fit.tangent[i * d + j];
        }
    }
    std::vector<double> residual(rows * n);
    for (std::size_t k = 0; k < rows; ++k)
        for (std::size_t i = 0; i < n; ++i) {
            double linear = 0;
            for (std::size_t j = 0; j < d; ++j) linear += fit.tangent[i * d + j] * eta[k * d + j];
            residual[k * n + i] = Y[k * n + i] - linear;
        }
    if (!nonlinear.empty()) {
        std::vector<double> theta(rows * nonlinear.size());
        for (std::size_t k = 0; k < rows; ++k) {
            const auto all = features(mExp, d, eta.data() + k * d);
            for (std::size_t c = 0; c < nonlinear.size(); ++c) theta[k * nonlinear.size() + c] = all[nonlinear[c]];
        }
        const auto w = regress(theta, rows, nonlinear.size(), residual, n, options.ridge);
        for (std::size_t c = 0; c < nonlinear.size(); ++c)
            for (std::size_t i = 0; i < n; ++i) fit.manifold[nonlinear[c] * n + i] = w[c * n + i];
    }
    std::vector<double> error(rows * n);
    for (std::size_t k = 0; k < rows; ++k) {
        const auto y = fit.lift(eta.data() + k * d);
        for (std::size_t i = 0; i < n; ++i) error[k * n + i] = y[i] - fit.equilibrium[i] - Y[k * n + i];
    }
    fit.manifoldError = relativeRms(error, Y);

    // Reduced dynamics eta' = r(eta), degrees 1..dynamicsOrder (the equilibrium stays at eta = 0).
    const auto rExp = polynomialExponentTable(d, options.dynamicsOrder, true);
    const auto rFeatures = rExp.size() / d;
    std::vector<double> theta(rows * (rFeatures - 1));
    for (std::size_t k = 0; k < rows; ++k) {
        const auto all = features(rExp, d, eta.data() + k * d);
        std::copy(all.begin() + 1, all.end(), theta.begin() + std::ptrdiff_t(k * (rFeatures - 1)));
    }
    const auto r = regress(theta, rows, rFeatures - 1, etaDot, d, options.ridge);
    fit.dynamics.assign(rFeatures * d, 0);
    std::copy(r.begin(), r.end(), fit.dynamics.begin() + std::ptrdiff_t(d));
    for (std::size_t k = 0; k < rows; ++k) {
        const auto v = fit.rate(eta.data() + k * d);
        for (std::size_t j = 0; j < d; ++j) error[k * d + j] = v[j] - etaDot[k * d + j];
    }
    error.resize(rows * d);
    fit.dynamicsError = relativeRms(error, etaDot);

    // The normal form of the reduced dynamics: every mode is a master mode.
    SpectralSubmanifoldOptions normal;
    normal.order = options.normalFormOrder ? options.normalFormOrder : options.dynamicsOrder;
    for (std::size_t j = 0; j < d; ++j) normal.masterModes.push_back(j);
    fit.normalForm = spectralSubmanifold(polynomialField(d, options.dynamicsOrder, fit.dynamics), std::vector<double>(d, 0), normal);
    return fit;
}

}  // namespace dynsyslab::core

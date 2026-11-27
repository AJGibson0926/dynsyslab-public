// Spectral submanifolds (core/SpectralSubmanifold.hpp) on oscillators with known backbone curves:
// the Duffing oscillator x'' + 2 zeta x' + x + kappa x^3 = 0, whose conservative frequency at
// amplitude a is omega = pi sqrt(1 + kappa a^2) / (2 K(m)), m = kappa a^2 / (2 (1 + kappa a^2))
// (Kovacic and Brennan, The Duffing Equation, Wiley 2011, section 2.2), and a two-mode oscillator
// whose slow mode is that Duffing oscillator driven into a faster linear mode, so its slow SSM has
// the same backbone.
#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/SpectralSubmanifold.hpp>
#include <dynsyslab/core/SINDy.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
using E = c::AnalyticExpression;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

constexpr double kappa = 0.5;

double exactFrequency(double a) {
    const double s = kappa * a * a;
    return 3.141592653589793 * std::sqrt(1 + s) / (2 * c::ellipticK(s / (2 * (1 + s))));
}

std::vector<E> duffing(double zeta) {
    const auto x = E::state(2, 0), v = E::state(2, 1);
    return {v, -2 * zeta * v - x - kappa * x.power(3)};
}

// x1'' + 2 z1 x1' + x1 + kappa x1^3 = 0 and x2'' + 2 z2 w2 x2' + w2^2 x2 + g x1^2 = 0.
std::vector<E> twoMode(double z1, double z2, double w2, double g) {
    const auto x1 = E::state(4, 0), v1 = E::state(4, 1), x2 = E::state(4, 2), v2 = E::state(4, 3);
    return {v1, -2 * z1 * v1 - x1 - kappa * x1.power(3), v2, -2 * z2 * w2 * v2 - w2 * w2 * x2 - g * x1 * x1};
}

// omega at amplitude a (of state 0) by interpolation along the backbone.
double frequencyAt(const c::BackboneCurve& b, double a) {
    for (std::size_t k = 1; k < b.radius.size(); ++k) {
        const double a0 = b.amplitude[(k - 1) * b.outputs], a1 = b.amplitude[k * b.outputs];
        if (a0 <= a && a <= a1) return b.frequency[k - 1] + (b.frequency[k] - b.frequency[k - 1]) * (a - a0) / (a1 - a0);
    }
    throw std::runtime_error("amplitude outside the backbone");
}

std::vector<double> radii(double top, std::size_t n) {
    std::vector<double> out;
    for (std::size_t k = 0; k <= n; ++k) out.push_back(top * double(k) / double(n));
    return out;
}

const c::StateOutputs firstState = [](const double* x) { return std::vector<double>{x[0]}; };

void duffingTest() {
    c::SpectralSubmanifoldOptions options;
    options.order = 7;
    const auto field = duffing(1e-3);
    const auto ssm = c::spectralSubmanifold(field, {0, 0}, options);
    require(ssm.dimension == 2 && ssm.masterModes.size() == 2, "the oscillating pair spans the SSM");
    require(std::abs(ssm.eigenvalues[0] - std::complex<double>(-1e-3, std::sqrt(1 - 1e-6))) < 1e-12, "linear spectrum");
    const auto backbone = c::backboneCurve(ssm, radii(0.4, 200), firstState, 1);
    for (const double a : {0.1, 0.2, 0.3})
        require(std::abs(frequencyAt(backbone, a) - exactFrequency(a)) < 3e-6, "Duffing backbone matches the exact frequency");
    require(ssm.invarianceError(field, 0.2) < 1e-5, "the parameterization is invariant to its order");
    // The real chart: q-polynomials reproduce the complex parameterization.
    const auto Wq = ssm.realParameterization();
    const auto exponents = c::polynomialExponentTable(2, options.order, true);
    const double q[2] = {0.11, -0.07};
    const std::complex<double> p[2] = {{q[0], q[1]}, {q[0], -q[1]}};
    const auto x = ssm.parameterization(p);
    for (std::size_t r = 0; r < 2; ++r) {
        double value = 0;
        for (std::size_t f = 0; f < exponents.size() / 2; ++f)
            value += Wq[f * 2 + r] * std::pow(q[0], double(exponents[2 * f])) * std::pow(q[1], double(exponents[2 * f + 1]));
        require(std::abs(value - x[r].real()) < 1e-12 && std::abs(x[r].imag()) < 1e-12, "real chart of the SSM");
    }
    const auto Rq = ssm.realReducedDynamics();
    const auto pd = ssm.reducedDynamics(p);
    double qdot0 = 0;
    for (std::size_t f = 0; f < exponents.size() / 2; ++f)
        qdot0 += Rq[f * 2] * std::pow(q[0], double(exponents[2 * f])) * std::pow(q[1], double(exponents[2 * f + 1]));
    require(std::abs(qdot0 - pd[0].real()) < 1e-12, "real chart of the reduced dynamics");
}

void twoModeTest() {
    const auto field = twoMode(1e-3, 0.05, 2.7, 0.8);
    c::SpectralSubmanifoldOptions options;
    options.order = 7;
    const auto ssm = c::spectralSubmanifold(field, std::vector<double>(4, 0), options);
    require(ssm.dimension == 2 && std::abs(ssm.eigenvalues[ssm.masterModes[0]].imag() - 1) < 1e-5, "the slow mode is the master");
    require(ssm.spectralQuotient >= 100, "spectral quotient of the slow SSM");
    require(ssm.invarianceError(field, 0.2) < 1e-5, "two-mode SSM invariance");
    const auto backbone = c::backboneCurve(ssm, radii(0.4, 200), firstState, 1);
    for (const double a : {0.1, 0.2, 0.3})
        require(std::abs(frequencyAt(backbone, a) - exactFrequency(a)) < 3e-6, "two-mode slow SSM has the Duffing backbone");
    // Graph style describes the same manifold: its invariance error is equally small.
    options.style = c::SsmStyle::graph;
    const auto graph = c::spectralSubmanifold(field, std::vector<double>(4, 0), options);
    require(graph.invarianceError(field, 0.2) < 1e-5, "graph-style SSM invariance");
    bool threw = false;
    try { (void)c::backboneCurve(graph, {0.1}, firstState, 1); } catch (const std::invalid_argument&) { threw = true; }
    require(threw, "backbone curves need the normal form");
}

// RK4 trajectories of the two-mode oscillator.
c::TimeSeries simulate(const std::vector<E>& field, std::vector<double> x, double dt, std::size_t samples) {
    const auto n = x.size();
    std::vector<double> time, values;
    const auto f = [&](const std::vector<double>& y) {
        std::vector<double> out(n);
        for (std::size_t i = 0; i < n; ++i) out[i] = field[i].evaluate(y.data());
        return out;
    };
    for (std::size_t k = 0; k < samples; ++k) {
        time.push_back(double(k) * dt);
        values.insert(values.end(), x.begin(), x.end());
        for (int s = 0; s < 4; ++s) {
            const double h = dt / 4;
            const auto k1 = f(x);
            auto y = x;
            for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k1[i];
            const auto k2 = f(y);
            for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h / 2 * k2[i];
            const auto k3 = f(y);
            for (std::size_t i = 0; i < n; ++i) y[i] = x[i] + h * k3[i];
            const auto k4 = f(y);
            for (std::size_t i = 0; i < n; ++i) x[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
    }
    return c::TimeSeries(std::move(time), std::move(values), n);
}

void dataDrivenTest() {
    const auto field = twoMode(0.01, 0.05, 2.7, 0.8);
    std::vector<c::TimeSeries> data;
    for (int k = 0; k < 4; ++k) {
        const double phase = 1.3 * k;
        data.push_back(simulate(field, {0.45 * std::cos(phase), 0.45 * std::sin(phase), 0.05, 0}, 0.05, 5000));
    }
    c::SsmFitOptions options;
    options.order = 3;
    options.dynamicsOrder = 5;
    options.skipTime = 60;
    options.equilibrium = {0, 0, 0, 0};
    const auto fit = c::fitSpectralSubmanifold(data, options);
    require(fit.manifoldError < 1e-3, "data-driven SSM graph fits the data");
    require(fit.dynamicsError < 1e-3, "data-driven reduced dynamics fit the derivatives");
    const auto slow = fit.normalForm.eigenvalues[fit.normalForm.masterModes[0]];
    require(std::abs(slow.imag() - std::sqrt(1 - 1e-4)) < 2e-3 && std::abs(slow.real() + 0.01) < 2e-3, "identified slow spectrum");
    const auto lifted = [&fit](const double* eta) { return std::vector<double>{fit.lift(eta)[0]}; };
    const auto backbone = c::backboneCurve(fit.normalForm, radii(0.4, 200), lifted, 1);
    require(std::abs(frequencyAt(backbone, 0.2) - exactFrequency(0.2)) < 1e-3, "data-driven backbone matches the exact frequency");
    // Prediction of a held-out trajectory from its first post-transient sample.
    const auto test = simulate(field, {0.0, 0.4, -0.02, 0.01}, 0.05, 5000).cropTime(60, 200);
    const auto predicted = fit.predict(test);
    double error = 0, size = 0;
    for (std::size_t k = 0; k < test.sampleCount(); ++k) {
        error += std::pow(predicted.value(k, 0) - test.value(k, 0), 2);
        size += std::pow(test.value(k, 0), 2);
    }
    require(std::sqrt(error / size) < 0.05, "the data-driven SSM model predicts a held-out trajectory");
}

}  // namespace

int main() {
    try {
        duffingTest();
        twoModeTest();
        dataDrivenTest();
        std::cout << "spectral submanifold tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

// Phase D identification pieces in the core: delay taps at physical-time lags and the quintic delay
// jet, SINDy options (rescale, pooled trajectories, single-absolute sparsification, directional
// finite differences, a radial dictionary), vector-field SINDy ridge/weights/weak form, HAVOK
// derivative order, eigenvalue summaries, empirical spectral weight and the SVD driver; delay
// pairs, delay consistency, locked basis extension, the endpoint decoder and the tangent readout.
#include <dynsyslab/core/Decomposition.hpp>
#include <dynsyslab/core/DelayDiagnostics.hpp>
#include <dynsyslab/core/TangentReadout.hpp>
#include <dynsyslab/core/DelaySampling.hpp>
#include <dynsyslab/core/ForecastScore.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }

// A damped linear oscillator x' = A x, A = [[-0.1, 1], [-1, -0.1]]: lambda = -0.1 +- i.
std::vector<double> oscillator(double x0, double y0, std::size_t n, double dt, std::vector<double>* rates = nullptr) {
    std::vector<double> out;
    for (std::size_t k = 0; k < n; ++k) {
        const double t = dt * double(k), e = std::exp(-0.1 * t);
        const double x = e * (x0 * std::cos(t) + y0 * std::sin(t)), y = e * (y0 * std::cos(t) - x0 * std::sin(t));
        out.insert(out.end(), {x, y});
        if (rates) rates->insert(rates->end(), {-0.1 * x + y, -x - 0.1 * y});
    }
    return out;
}

void delayTaps() {
    std::vector<double> time, values;
    for (int k = 0; k < 11; ++k) { time.push_back(0.1 * k); values.push_back(2.0 * 0.1 * k + 1.0); }
    const c::TimeSeries line(time, values, 1);
    const auto taps = c::delaySamples(line, {0.9, 1.0}, {0.0, 0.25});
    require(std::abs(taps[0] - (2 * 0.9 + 1)) < 1e-12 && std::abs(taps[1] - (2 * 0.65 + 1)) < 1e-12
            && std::abs(taps[3] - (2 * 0.75 + 1)) < 1e-12, "linear taps reproduce a linear signal at fractional lags");
    bool refused = false;
    try { (void)c::delaySamples(line, {0.1}, {0.5}); } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "a window before the record is refused");
    // The quintic Hermite jet is exact for a quintic polynomial.
    std::vector<double> jet, clock;
    const auto p = [](double t) { return 1 + t - 2 * t * t + 0.5 * std::pow(t, 3) + 0.3 * std::pow(t, 4) - 0.1 * std::pow(t, 5); };
    const auto dp = [](double t) { return 1 - 4 * t + 1.5 * t * t + 1.2 * std::pow(t, 3) - 0.5 * std::pow(t, 4); };
    const auto ddp = [](double t) { return -4 + 3 * t + 3.6 * t * t - 2 * std::pow(t, 3); };
    for (int k = 0; k <= 10; ++k) { const double t = 0.2 * k; clock.push_back(t); jet.insert(jet.end(), {p(t), dp(t), ddp(t)}); }
    const auto out = c::sampleDelayJet(c::TimeSeries(clock, jet, 3), {1.8}, {0.33});
    const double q = 1.8 - 0.33;
    require(std::abs(out[0] - p(q)) < 1e-10 && std::abs(out[1] - dp(q)) < 1e-9 && std::abs(out[2] - ddp(q)) < 1e-8,
            "the C2 quintic delay jet is exact for a quintic");
}

void sindyOptions() {
    const std::size_t n = 400;
    const double dt = 0.02;
    std::vector<double> time;
    for (std::size_t k = 0; k < n; ++k) time.push_back(dt * double(k));
    std::vector<double> rates;
    const auto states = oscillator(1.0, 0.0, n, dt, &rates);
    const std::vector<std::complex<double>> lambda{{-0.1, 1.0}, {-0.1, -1.0}};
    const auto check = [&](const c::EigenfunctionFitResult& fit, const std::string& what, double tolerance) {
        require(fit.spectrum.residual(1).has_value(), what + ": two eigenfunctions");
        // psi(x(t)) e^{-lambda t} is constant along the trajectory.
        const auto first = fit.functions.evaluate(states.data());
        const auto later = fit.functions.evaluate(states.data() + 2 * 300);
        const auto predicted = first[0] * std::exp(lambda[0] * time[300]);
        require(std::abs(later[0] - predicted) < tolerance * std::abs(predicted), what + ": the eigenfunction advances by its eigenvalue");
    };
    c::KoopmanSINDyOptions base;
    base.polynomialDegree = 1;
    base.threshold = 0.0;
    base.includeConstant = false;
    base.normalizeStates = false;  // centering would make the linear features affine
    check(c::fitKoopmanEigenfunctions(states, n, 2, time, lambda, rates, base), "default", 1e-6);
    auto rescaled = base;
    rescaled.rescale = {2.0, 0.5};
    check(c::fitKoopmanEigenfunctions(states, n, 2, time, lambda, rates, rescaled), "rescale", 1e-6);
    auto directional = base;
    directional.directionalDerivative = "finite-difference";
    check(c::fitKoopmanEigenfunctions(states, n, 2, time, lambda, rates, directional), "directional finite difference", 1e-5);
    auto single = base;
    single.sparsification = "single-absolute";
    single.threshold = 1e-8;
    check(c::fitKoopmanEigenfunctions(states, n, 2, time, lambda, rates, single), "single-absolute", 1e-6);
    // Two trajectories pooled, finite differences within each (time restarts at the boundary).
    auto second = oscillator(0.0, 1.5, n, dt);
    auto pooledStates = states;
    pooledStates.insert(pooledStates.end(), second.begin(), second.end());
    auto pooledTime = time;
    pooledTime.insert(pooledTime.end(), time.begin(), time.end());
    auto pooled = base;
    pooled.segments = {n, n};
    check(c::fitKoopmanEigenfunctions(pooledStates, 2 * n, 2, pooledTime, lambda, {}, pooled), "pooled trajectories", 1e-2);
    bool refused = false;
    try { (void)c::fitKoopmanEigenfunctions(pooledStates, 2 * n, 2, pooledTime, lambda, {}, base); }
    catch (const std::invalid_argument&) { refused = true; }
    require(refused, "a restarting clock without segments is refused");
    // A dictionary of the linear coordinates written as general terms gives the same eigenfunctions.
    std::vector<c::DictionaryTerm> terms(2);
    terms[0].kind = "coordinate"; terms[0].index = 0;
    terms[1].kind = "coordinate"; terms[1].index = 1;
    auto dictionary = base;
    dictionary.dictionary = std::make_shared<const c::ExpansionDictionary>(2, terms);
    const auto general = c::fitKoopmanEigenfunctions(states, n, 2, time, lambda, rates, dictionary);
    require(!general.functions.monomial(), "a general dictionary is kept as given");
    check(general, "general dictionary", 1e-6);
}

void sindyDynamicsOptions() {
    const std::size_t n = 500;
    const double dt = 0.01;
    std::vector<double> time, rates;
    for (std::size_t k = 0; k < n; ++k) time.push_back(dt * double(k));
    const auto states = oscillator(1.0, 0.5, n, dt, &rates);
    const c::TimeSeries q(time, states, 2), dq(time, rates, 2);
    const auto library = c::polynomialDynamicsLibrary(2, 1);
    const auto coefficient = [](const c::SindyDynamicsModel& model, std::size_t feature, std::size_t output) {
        return model.coefficients()[feature * model.outputs() + output];
    };
    const auto exact = [&](const c::SindyDynamicsModel& model, double tolerance, const std::string& what) {
        // features: 1, x, y; xdot = -0.1 x + y, ydot = -x - 0.1 y
        require(std::abs(coefficient(model, 1, 0) + 0.1) < tolerance && std::abs(coefficient(model, 2, 0) - 1) < tolerance
                && std::abs(coefficient(model, 1, 1) + 1) < tolerance, what + " recovers the vector field");
    };
    c::SindyDynamicsOptions options;
    options.threshold = 1e-6;
    exact(c::fitSindyDynamics(q, dq, library, options), 1e-8, "plain least squares");
    auto weighted = options;
    weighted.sampleWeights = c::trapezoidWeights(q);
    exact(c::fitSindyDynamics(q, dq, library, weighted), 1e-8, "trapezoid-weighted");
    auto weak = options;
    weak.weakDurations = {0.2, 0.37};
    weak.weakStride = 5;
    exact(c::fitSindyDynamics(q, dq, library, weak), 2e-3, "weak (integral) form");
    auto ridge = options;
    ridge.ridge = 1e3;
    const auto shrunk = c::fitSindyDynamics(q, dq, library, ridge);
    require(std::abs(coefficient(shrunk, 2, 0)) < 0.9, "a large ridge shrinks the coefficients");
}

void spectraAndSvd() {
    const auto spectrum = c::Spectrum::fromContinuousRates({{-0.5, 2.0}, {-0.5, -2.0}, {-1.0, 0.0}});
    const auto value = spectrum.eigenvalue(0);
    require(value.growthRate() == -0.5 && value.angularFrequency() == 2.0 && std::abs(value.qualityFactor() - 2.0) < 1e-15
            && std::abs(value.cyclicFrequency() - 2.0 / (2 * 3.141592653589793)) < 1e-15
            && spectrum.eigenvalue(2).qualityFactor() == 0.0, "eigenvalue summaries: rate, omega, f, Q");
    c::EigenfunctionSamples samples{{0, 1, 2}, {0, 1, 2}, {{1, 0}, {0, 0}, {2, 0}, {0, 0}, {1, 0}, {0, 0}}, 2};
    const auto weight = c::empiricalSpectralWeight(samples);
    require(std::abs(weight[0] - 1.0 / 6) < 1e-15 && std::abs(weight[2] - 4.0 / 6) < 1e-15 && std::abs(weight[1] - 1.0 / 3) < 1e-15,
            "empirical spectral weight: |psi|^2 normalized per mode, uniform for a vanishing mode");
    const std::vector<double> matrix{3, 1, 0, 1, 2, 1, 0, 1, 4, 2, 0, 1};
    const auto dc = dynsyslab::numerics::singularValueDecomposition(matrix, 4, 3);
    std::vector<double> jacobi;
    {
        const dynsyslab::numerics::SvdDriverScope scope(dynsyslab::numerics::svdDriverFromName("gesvd"));
        require(dynsyslab::numerics::svdDriver() == dynsyslab::numerics::SvdDriver::jacobi, "the driver scope applies");
        jacobi = dynsyslab::numerics::singularValueDecomposition(matrix, 4, 3).singularValues;
    }
    require(dynsyslab::numerics::svdDriver() == dynsyslab::numerics::SvdDriver::divideAndConquer, "the driver scope restores");
    for (std::size_t k = 0; k < 3; ++k) require(std::abs(dc.singularValues[k] - jacobi[k]) < 1e-12, "both SVD drivers agree");
}

void forecastScoring() {
    std::vector<double> t, x, pt, p;
    for (int k = 0; k <= 1000; ++k) { t.push_back(0.01 * k); x.push_back(std::sin(2.0 * 0.01 * k)); }
    for (int k = 0; k <= 600; ++k) { pt.push_back(0.013 * k); p.push_back(std::sin(2.0 * (0.013 * k - 0.05))); }
    c::ForecastScoreOptions options;
    options.prominence = 0.5;
    const auto score = c::forecastScore(c::TimeSeries(t, x, 1), c::TimeSeries(pt, p, 1), options);
    require(std::abs(score.at("coverage").number() - 0.78) < 1e-9 && !score.at("complete").boolean(),
            "coverage of a shorter forecast on its own clock");
    const auto& matches = score.at("peaks").at("matches").array();
    require(score.at("peaks").at("reference_count").number() == 3 && std::abs(matches[0].at("time_error").number() - 0.05) < 0.01,
            "prominent peaks match with the forecast's delay");
    require(c::prominentPeaks({0, 2, 1, 1.5, 0, 3, 0}, 1.0) == std::vector<std::size_t>{1, 5}, "topographic prominence");
}

void havokOrders() {
    std::vector<double> time, x;
    for (int k = 0; k < 600; ++k) { time.push_back(0.01 * k); x.push_back(std::sin(1.3 * 0.01 * k) + 0.4 * std::sin(3.1 * 0.01 * k)); }
    c::HAVOKOptions options;
    options.delays = 40;
    options.rank = 5;
    for (const std::size_t order : {1u, 2u, 4u}) {
        options.derivativeOrder = order;
        const auto result = c::havokDecomposition(x, time, x.size(), 1, options);
        require(result.generator.size() == 16, "HAVOK fits with derivative order " + std::to_string(order));
    }
    options.derivativeOrder = 3;
    bool refused = false;
    try { (void)c::havokDecomposition(x, time, x.size(), 1, options); } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "an unsupported HAVOK derivative order is refused");
}

void delayDiagnostics() {
    // Delay pairs of p = t^2 at dt = 0.1: the secant is exact for a quadratic, curvature is 2.
    std::vector<double> p;
    for (int k = 0; k < 40; ++k) p.push_back(std::pow(0.1 * k, 2));
    const auto pairs = c::scalarDelayPairs(p, {20, 30}, {0, 5}, 0.1, 2);
    const double centre = 0.1 * (20 - 2);  // estimate time of end 20 at lag 0
    require(std::abs(pairs.state[1] - 2 * centre) < 1e-12 && std::abs(pairs.curvature[0] - 2) < 1e-9
            && std::abs(pairs.latency - 0.2) < 1e-15 && pairs.raw[1] == p[20] && pairs.raw[0] == p[16], "delay pair secants");
    bool refused = false;
    try { (void)c::scalarDelayPairs(p, {3}, {2}, 0.1, 1); } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "a pair outside its record is refused");

    // Delay consistency: taps of sin t at lags 0, 0.3 are exactly shifted copies.
    std::vector<double> time, taps;
    for (int k = 0; k < 400; ++k) {
        const double t = 0.01 * k;
        time.push_back(t);
        taps.insert(taps.end(), {std::sin(t), std::cos(t), std::sin(t - 0.3), std::cos(t - 0.3)});
    }
    const auto consistent = c::delayConsistency(time, taps, taps, {0.0, 0.3});
    require(consistent.value.size() == 1 && consistent.value[0] < 1e-7 && consistent.rate[0] < 1e-5, "consistent taps");
    auto shifted = taps;
    for (std::size_t k = 0; k < time.size(); ++k) shifted[k * 4 + 2] += 0.1;
    require(c::delayConsistency(time, shifted, taps, {0.0, 0.3}).value[0] > 0.05, "inconsistent taps are detected");

    // Locked basis: the leading eigenpair of a diagonal Gram stays, the next is appended.
    const std::vector<double> gram{9, 0, 0, 0, 4, 0, 0, 0, 1};
    const auto basis = c::extendLockedBasis(gram, 3, {1, 0, 0}, {3}, 2);
    require(basis.modes[0] == 1 && basis.modes[2] == 0 && std::abs(std::abs(basis.modes[3]) - 1) < 1e-12
            && basis.modes[3] > 0 && std::abs(basis.singularValues[1] - 2) < 1e-12 && basis.gramError < 1e-12, "locked prefix extension");
    refused = false;
    try { (void)c::extendLockedBasis(gram, 3, {0, 1, 0}, {3}, 2); } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "a prefix that is not an eigenpair is refused");

    // Endpoint decoder: y_end = d . xi on the native coordinates xi = U^T w / s.
    c::HankelLeftBasis hankel;
    hankel.rows = 3; hankel.rank = 2; hankel.leftSingularVectors = {0.6, 0.0, 0.0, 1.0, 0.8, 0.0}; hankel.singularValues = {5, 2};
    const auto d = c::historyEndpointDecoder(hankel, 2);
    const double w[] = {0.6 * 5, 2.0, 0.8 * 5};  // the window U diag(s) [1, 1]
    double xi[2] = {0, 0};
    for (int i = 0; i < 3; ++i) for (int m = 0; m < 2; ++m) xi[m] += hankel.leftSingularVectors[i * 2 + m] * w[i] / hankel.singularValues[m];
    require(std::abs(d[0] * xi[0] + d[1] * xi[1] - w[2]) < 1e-12, "endpoint decoder recovers the newest sample");
}

void tangentReadout() {
    // g(b) = b on a line b in [-2, 2] with unit velocity: value and rate are recovered between samples.
    std::vector<double> b, v, g, dg;
    for (int k = 0; k <= 80; ++k) { b.push_back(-2 + 0.05 * k); v.push_back(1); g.push_back(b.back()); dg.push_back(1); }
    const c::TangentReadout readout(b, v, g, dg, {1.0}, 1.0);
    const double query = 0.512, rate = 1.0;
    const auto phase = readout.phase(&query, &rate, 4);
    require(std::abs(phase.value - query) < 2e-3 && std::abs(phase.rate - 1) < 5e-2, "tangent readout value and rate");
}

}  // namespace

int main() {
    try {
        delayTaps();
        sindyOptions();
        sindyDynamicsOptions();
        spectraAndSvd();
        havokOrders();
        forecastScoring();
        delayDiagnostics();
        tangentReadout();
        std::cout << "identification core tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "identification core tests failed: " << error.what() << '\n';
        return 1;
    }
}

// Phase E core algorithms: signal metrics, the linear observer, reconstruction quality, path
// sampling, eigenfunction refit/gauge/prune, POD-Galerkin, harmonic envelopes, the symmetric-pulse
// input-map estimator and control-map path diagnostics.
#include <dynsyslab/core/ControlMapPath.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/EigenfunctionRefit.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/LinearObserver.hpp>
#include <dynsyslab/core/PathSampling.hpp>
#include <dynsyslab/core/PodGalerkin.hpp>
#include <dynsyslab/core/ReconstructionQuality.hpp>
#include <dynsyslab/core/SignalMetrics.hpp>
#include <dynsyslab/core/SymmetricPulseControlMap.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Targets.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace core = dynsyslab::core;
using Complex = std::complex<double>;
constexpr double pi = 3.14159265358979323846;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void signalMetrics() {
    std::vector<double> t, x, e;
    for (int k = 0; k < 2000; ++k) {
        const double time = 0.001 * k;
        t.push_back(time);
        x.push_back(0.5 + 2.0 * std::sin(2 * pi * 3.0 * time + 0.3) + 0.2 * std::sin(2 * pi * 6.0 * time));
        e.push_back(std::exp(-5.0 * time));
    }
    const auto fit = core::fitTone(t, x, 3.0, 2);
    require(std::abs(fit.offset - 0.5) < 1e-9 && std::abs(fit.amplitudes[0] - 2.0) < 1e-9
                && std::abs(fit.phases[0] - 0.3) < 1e-9 && std::abs(fit.harmonicContent - 0.1) < 1e-9,
            "tone fit did not recover offset, amplitude, phase and harmonic content");
    const double settle = core::settlingTime(t, e, 0.01);
    require(std::abs(settle - std::log(100.0) / 5.0) < 2e-3, "settling time is wrong");
    const auto filtered = core::lowPass(t, std::vector<double>(t.size(), 1.0), 10.0);
    require(std::abs(filtered.back() - 1.0) < 1e-12, "low-pass does not keep a constant");
    std::vector<double> alternating;
    for (int k = 0; k < 1000; ++k) alternating.push_back(k % 2 ? 1.0 : -1.0);
    require(core::autocorrelation(alternating, 1) < -0.99, "lag-1 autocorrelation of an alternating sequence");
}

void linearObserver() {
    // Harmonic oscillator x' = [0 1; -1 0] x, y = x0: the LQE estimate converges to the truth.
    const std::vector<double> A{0, 1, -1, 0}, C{1, 0};
    const auto model = core::designLinearQuadraticEstimator(A, {}, C, {1, 0, 0, 1}, {0.01}, 2, 0, 1);
    for (const auto& pole : model.observerPoles) require(pole.real() < 0, "LQE observer poles must be stable");
    core::LinearObserver observer(model, 0.01, {0.0, 0.0});
    std::vector<double> estimate;
    for (int k = 0; k <= 2000; ++k) {
        const double time = 0.01 * k;
        estimate = observer.observe(time, {std::cos(time)}, {});
    }
    const double time = 20.0;
    require(std::abs(estimate[0] - std::cos(time)) < 0.05 && std::abs(estimate[1] + std::sin(time)) < 0.05,
            "LQE estimate did not converge");
    core::StaticMapEstimator map{{2.0, 3.0}, {1.0, 0.0}, 2, 1};
    const auto x = map.estimate({2.0});
    require(x[0] == 5.0 && x[1] == 6.0, "static map estimate");
}

void reconstructionQuality() {
    std::vector<Complex> psi, y;
    for (int k = 0; k < 200; ++k) {
        const double a = std::sin(0.1 * k), b = std::cos(0.37 * k);
        psi.insert(psi.end(), {a, b});
        y.push_back(2 * a - b);
    }
    const auto good = core::assessReconstruction(psi, y, 200, 2, 1, 0.0, {});
    require(good.passed && good.r2[0] > 0.999999, "an exact readout must pass the held-out gate");
    core::ReconstructionQualityPolicy strict;
    strict.maximumConditionNumber = 1.0;
    const auto bad = core::assessReconstruction(psi, y, 200, 2, 1, 0.0, strict);
    require(!bad.passed && !bad.failures.empty(), "the condition-number threshold must fail");
}

void pathsAndRefit() {
    std::vector<double> path{0, 0, 1, 0, 1, 1};
    const auto resampled = core::resampleByArcLength(path, 3, 2, 5);
    require(resampled.size() == 10 && std::abs(resampled[4] - 1.0) < 1e-12 && std::abs(resampled[5]) < 1e-12
                && std::abs(resampled[9] - 1.0) < 1e-12, "arc-length resampling midpoint and end");
    const auto tube = core::tubeStates(path, 3, 2, 0.1, 4, 7);
    require(tube.size() == (3 + 12) * 2, "tube sample count");
    std::vector<Complex> values{1.0, 2.0, 3.0}, reference{Complex(0, 2), Complex(0, 4), Complex(0, 6)};
    const auto alpha = core::gaugeFactor(values, reference);
    require(std::abs(alpha - Complex(0, 2)) < 1e-12, "gauge factor");
    std::vector<Complex> coefficients{1.0, 1e-6, 0.5};
    require(core::pruneCoefficients(coefficients, 1e-3) == 1 && coefficients[1] == 0.0, "pruning");
    // theta = [x, x^2] on x_{k+1} = 0.5 x_k: psi = x is the multiplier-0.5 eigenfunction.
    std::vector<double> now, next;
    for (double x = 0.1; x < 1.0; x += 0.1) { now.insert(now.end(), {x, x * x}); next.insert(next.end(), {0.5 * x, 0.25 * x * x}); }
    double residual = 1;
    const auto c = core::refitAtMultiplier(now, next, now.size() / 2, 2, 0.5, &residual);
    require(std::abs(c[1]) < 1e-8 && std::abs(std::abs(c[0]) - 1.0) < 1e-8 && residual < 1e-8, "fixed-multiplier refit");
}

core::KoopmanModel duffingModel() {
    // psi = (x, v) on Duffing: the exact input map is grad psi . g = g.
    auto fit = core::polynomialEigenfunctions(2, {1, 0, 0, 1}, {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
        core::Spectrum::fromContinuousRates({{-1.0, 0.0}, {-2.0, 0.0}}), {0.0, 0.0}, {1.0, 1.0});
    core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "double_well_reference";
    return core::KoopmanModel(std::move(fit), std::move(system));
}

void podAndPulses() {
    core::SystemConfiguration configuration;
    configuration.identifier = "duffing"; configuration.preset = "double_well_reference";
    auto system = core::createSystem(configuration);
    std::vector<double> snapshots;
    for (int k = 0; k < 50; ++k) snapshots.insert(snapshots.end(), {std::cos(0.1 * k), std::sin(0.2 * k)});
    const auto pod = core::fitPodGalerkin(snapshots, 50, 2, 2, 0, "mean", *system);
    const auto x = std::vector<double>{0.3, -0.2};
    const auto back = pod.reconstruct(pod.project(x.data()).data());
    require(std::abs(back[0] - 0.3) < 1e-12 && std::abs(back[1] + 0.2) < 1e-12, "full-rank POD round trip");
    require(std::abs(pod.energyFraction - 1.0) < 1e-12 && pod.A.size() == 4, "POD-Galerkin operators");

    const auto model = duffingModel();
    std::vector<core::SymmetricPulseRun> runs;
    const double dt = 1e-4, u = 0.5;
    for (std::size_t pair = 0; pair < 3; ++pair) {
        const std::vector<double> anchor{0.2 * double(pair) - 0.1, 0.05 * double(pair)};
        for (int polarity : {-1, 1}) {
            std::vector<double> dx(2);
            system->rhsWithInput(0.0, anchor.data(), polarity * u, dx.data());
            runs.push_back({pair, polarity, anchor, {anchor[0] + dt * dx[0], anchor[1] + dt * dx[1]}, polarity * u, dt});
        }
    }
    const auto samples = core::estimateSymmetricPulseControlMap(model, runs);
    const auto exact = model.controlMap(runs.front().anchor.data(), "exact");
    for (std::size_t j = 0; j < 2; ++j)
        require(std::abs(samples.maps[j] - exact[j]) < 1e-9, "symmetric-pulse map equals the exact map to first order");
    const auto constant = core::fitConstantControlMapSamples(model, samples);
    require(constant.relativeResidual < 1e-9, "a state-independent map fits a constant exactly");

    core::KoopmanModel copy = model;
    copy.setControlMap("constant", constant.map);
    std::vector<double> states{0.1, 0.0, 0.2, 0.1, -0.1, 0.05};
    const auto path = core::analyzeControlMapPath(copy, states, 3, {});
    require(path.passed && path.relativeError[0] < 1e-9 && path.cosine[0] > 1 - 1e-12, "path diagnostics of an exact map");
}

void envelope() {
    const auto& target = core::rpeControlTarget("rpe_multisine_offset_4_mixed");
    const auto [low, high] = target.conservativeEnvelope();
    for (double time = 0; time < 20; time += 0.01) {
        const auto x = target.evaluate(time);
        require(x[0] >= low[0] && x[0] <= high[0] && x[1] >= low[1] && x[1] <= high[1], "envelope contains the target");
    }
}

}  // namespace

int main() {
    try {
        signalMetrics(); linearObserver(); reconstructionQuality(); pathsAndRefit(); podAndPulses(); envelope();
        std::cout << "dynsyslab_control_estimation_core_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_control_estimation_core_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}

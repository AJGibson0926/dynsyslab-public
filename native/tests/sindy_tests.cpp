#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <optional>
#include <stdexcept>

int main() {
    try {
        constexpr std::size_t samples = 101; const double dt = 0.01;
        std::vector<double> states, velocities, time;
        for (std::size_t index = 0; index < samples; ++index) {
            const double t = index * dt, x = std::exp(-t);
            time.push_back(t); states.push_back(x); velocities.push_back(-x);
        }
        dynsyslab::core::KoopmanSINDyOptions options;
        options.polynomialDegree = 2; options.threshold = 1e-8; options.maximumIterations = 2;
        options.normalizeStates = false;
        const auto result = dynsyslab::core::fitKoopmanEigenfunctions(
            states, samples, 1, time, {{-1.0, 0.0}}, velocities, options);
        if (result.functions.dictionary().size() != 3 || result.functions.coefficients().size() != 3)
            throw std::runtime_error("polynomial feature shape");
        const double point[]{0.5};
        const auto value = result.functions.evaluate(0, point);
        const auto gradient = result.functions.gradient(0, point);
        if (!std::isfinite(value.real()) || gradient.size() != 1 || result.spectrum.residual(0)->value > 1e-8)
            throw std::runtime_error("polynomial eigenfunction fit");
        std::vector<double> dynamicalStates, dynamicalVelocities;
        for (std::size_t index=0;index<samples;++index) {
            const double x=states[index];
            dynamicalStates.insert(dynamicalStates.end(),{x,0.25*x});
            dynamicalVelocities.insert(dynamicalVelocities.end(),{-x,-0.25*x});
        }
        dynsyslab::core::SystemConfiguration configuration;
        configuration.identifier="duffing"; configuration.preset="double_well_reference";
        const auto system=dynsyslab::core::createSystem(configuration);
        const auto coordinate=system->observableCatalog().at("position");
        const auto mapped=dynsyslab::core::fitKoopmanEigenfunctionsInCoordinates(
            dynamicalStates,samples,2,time,{{-1.0,0.0}},coordinate,
            dynamicalVelocities,options);
        if (mapped.functions.coordinates().dimension!=1
            || mapped.functions.coordinates().physicalStateDimension!=std::optional<std::size_t>(2)
            || mapped.functions.coordinates().names!=std::vector<std::string>{"position"}
            || mapped.spectrum.residual(0)->value>1e-8)
            throw std::runtime_error("observable-coordinate SINDy fit");

        // A damped linear oscillator xdot = A x, A = [[-s, w], [-w, -s]]: eigenvalues -s ± i w with
        // eigenfunctions x1 ∓ i x2 (phi = x1 - i x2 has lambda = -s + i w).
        const double s = 0.1, w = 2.0; constexpr std::size_t n = 1001;
        std::vector<double> oscillator, rates, clock;
        for (std::size_t k = 0; k < n; ++k) {
            const double t = 0.01 * k, decay = std::exp(-s * t);
            const double x1 = decay * std::cos(w * t), x2 = -decay * std::sin(w * t);
            clock.push_back(t); oscillator.insert(oscillator.end(), {x1, x2});
            rates.insert(rates.end(), {-s * x1 + w * x2, -w * x1 - s * x2});
        }
        const auto state = system->observableCatalog().at("state");
        dynsyslab::core::KoopmanSINDyOptions linear; linear.polynomialDegree = 2; linear.threshold = 1e-8;
        // 1. Eigenvalues from data alone (Kaiser, Kutz & Brunton 2021, Algorithm 1).
        const auto found = dynsyslab::core::estimateKoopmanEigenvaluesInCoordinates(
            oscillator, n, 2, clock, state, rates, linear);
        const auto near = [&](std::complex<double> target) {
            return std::any_of(found.eigenvalues.begin(), found.eigenvalues.end(),
                               [&](auto value) { return std::abs(value - target) < 1e-6; });
        };
        if (!near({-s, w}) || !near({-s, -w}) || near({0.0, 0.0}))
            throw std::runtime_error("eigenvalues from data: -s ± i w found, trivial zero dropped");
        // 2. An anchored fit reproduces reference values, including their scale and phase.
        const std::complex<double> factor = std::polar(3.0, 0.7);
        std::vector<double> anchorStates; std::vector<std::complex<double>> anchorValues;
        for (std::size_t k = 0; k < n; k += 10) {
            anchorStates.insert(anchorStates.end(), {oscillator[2 * k], oscillator[2 * k + 1]});
            anchorValues.push_back(factor * std::complex<double>(oscillator[2 * k], -oscillator[2 * k + 1]));
        }
        linear.anchorStates = &anchorStates; linear.anchorValues = &anchorValues;
        const auto anchored = dynsyslab::core::fitKoopmanEigenfunctionsInCoordinates(
            oscillator, n, 2, clock, {{-s, w}}, state, rates, linear);
        const double probe[]{0.3, -0.2};
        const auto expected = factor * std::complex<double>(0.3, 0.2);
        if (std::abs(anchored.functions.evaluate(0, probe) - expected) > 1e-6
            || anchored.spectrum.residual(0)->value > 1e-6)
            throw std::runtime_error("anchored SINDy fit does not reproduce the reference eigenfunction");
        std::cout << "dynsyslab_sindy_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_sindy_tests: FAIL: " << error.what() << '\n'; return 1;
    }
}

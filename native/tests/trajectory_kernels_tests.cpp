// Trajectory kernels with exact answers: Fornberg finite-difference weights, high-order central
// differences, affine equations and their exact propagation, backward prehistory of the harmonic
// oscillator, and clock comparison of a decimated record.
#include <dynsyslab/core/AffineDynamics.hpp>
#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/core/RunInvariant.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TrajectoryAudit.hpp>

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

}  // namespace

int main() {
    try {
        const auto w = c::differenceWeights({-1, 0, 1}, 1);
        require(std::abs(w[0] + 0.5) < 1e-15 && std::abs(w[1]) < 1e-15 && std::abs(w[2] - 0.5) < 1e-15, "three-point first derivative");
        const auto w2 = c::differenceWeights({-2, -1, 0, 1, 2}, 2);
        const double expected[5] = {-1.0 / 12, 4.0 / 3, -5.0 / 2, 4.0 / 3, -1.0 / 12};
        for (int k = 0; k < 5; ++k) require(std::abs(w2[k] - expected[k]) < 1e-13, "five-point second derivative");
        std::vector<double> y;
        const double dt = 0.01;
        for (int k = 0; k < 600; ++k) y.push_back(std::sin(double(k) * dt));
        const auto d = c::highOrderCentralDifference(y, dt, 9);
        double worst = 0;
        for (int k = 0; k < 600; ++k) worst = std::max(worst, std::abs(d[k] - std::cos(double(k) * dt)));
        require(worst < 1e-9, "eighth-order central difference, ends included");
        c::DerivativeSpec spec;
        spec.method = "central_high_order";
        require(std::abs(c::differentiate(y, dt, spec).derivative[300] - std::cos(3.0)) < 1e-12, "central_high_order method");

        // x' = -0.5 x + y + 1, y' = -x: affine, exact flow against a fine RK4 reference.
        const auto x = E::state(2, 0), v = E::state(2, 1);
        const auto field = c::Observable::fromExpressions({"f", "f", "", "native", {}, false}, {-0.5 * x + v + 1.0, -x});
        const auto eq = c::affineEquations(field);
        require(std::abs(eq.A[0] + 0.5) < 1e-12 && std::abs(eq.A[1] - 1) < 1e-12 && std::abs(eq.A[2] + 1) < 1e-12 && std::abs(eq.c[0] - 1) < 1e-12,
                "A and c of an affine field");
        std::vector<double> times;
        for (int k = 0; k <= 100; ++k) times.push_back(0.05 * k);
        const auto exact = c::affinePropagation(eq, {0.3, -0.2}, times);
        std::vector<double> state{0.3, -0.2};
        for (int k = 0; k < 5000; ++k) {  // RK4 to t = 5 with h = 1e-3
            const auto f = [&](const std::vector<double>& s) { return std::vector<double>{-0.5 * s[0] + s[1] + 1, -s[0]}; };
            const double h = 1e-3;
            const auto k1 = f(state);
            const auto k2 = f({state[0] + h / 2 * k1[0], state[1] + h / 2 * k1[1]});
            const auto k3 = f({state[0] + h / 2 * k2[0], state[1] + h / 2 * k2[1]});
            const auto k4 = f({state[0] + h * k3[0], state[1] + h * k3[1]});
            for (int i = 0; i < 2; ++i) state[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
        require(std::abs(exact[200] - state[0]) < 1e-10 && std::abs(exact[201] - state[1]) < 1e-10, "exact affine propagation by doubling");
        const auto uneven = c::affinePropagation(eq, {0.3, -0.2}, {0, 1.7, 5});
        require(std::abs(uneven[4] - state[0]) < 1e-10, "exact affine propagation on an uneven clock");
        bool refused = false;
        try { (void)c::affineEquations(c::Observable::fromExpressions({"g", "g", "", "native", {}, false}, {v, -x - x.power(3)})); }
        catch (const std::invalid_argument&) { refused = true; }
        require(refused, "a nonlinear field has no A, c form");

        // Prehistory of the undamped oscillator from (1, 0) at t0 = 0: x(t) = cos(t).
        c::SystemConfiguration config;
        config.identifier = "harmonic_oscillator";
        config.preset = "undamped";
        const auto system = c::createSystem(config);
        std::vector<double> before;
        for (int k = 20; k >= 1; --k) before.push_back(-0.1 * k);
        const auto history = c::backwardHistory(*system, {1, 0}, 0, before, 1e-3);
        for (std::size_t k = 0; k < before.size(); ++k)
            require(std::abs(history.states[2 * k] - std::cos(before[k])) < 1e-10, "backward prehistory of the oscillator");
        require(history.roundtripError < 1e-10, "the prehistory returns to the initial state");

        // A decimated record against its dense reference: linear interpolation error ~ h^2 / 8.
        std::vector<double> rt, ref, st, sampled;
        for (int k = 0; k <= 2000; ++k) { rt.push_back(0.001 * k); ref.push_back(std::sin(0.001 * k)); }
        for (int k = 0; k <= 40; ++k) { st.push_back(0.05 * k); sampled.push_back(std::sin(0.05 * k)); }
        const auto audit = c::compareClocks(rt, ref, st, sampled, 1);
        require(audit.maximumAbsolute[0] < 0.05 * 0.05 / 8 + 1e-9 && audit.maximumAbsolute[0] > 1e-5 && audit.referenceCount == 2001, "clock audit");
        // Runs of the unit oscillator at three energies: the learned invariant is x^2 + v^2 up to scale.
        std::vector<c::TimeSeries> runs;
        for (const double radius : {0.5, 1.0, 1.5}) {
            std::vector<double> rt, rv;
            for (int k = 0; k < 400; ++k) { const double t = 0.02 * k; rt.push_back(t); rv.push_back(radius * std::cos(t)); rv.push_back(-radius * std::sin(t)); }
            runs.emplace_back(rt, rv, 2);
        }
        const auto invariant = c::fitRunInvariant(runs, 2);
        // Table order: 1, x, v, x^2, x v, v^2 (graded); x^2 and v^2 carry equal weight, the rest none.
        const auto& iw = invariant.coefficients;
        require(std::abs(iw[3] - iw[5]) < 1e-6 * std::abs(iw[3]) && std::abs(iw[1]) + std::abs(iw[2]) + std::abs(iw[4]) < 1e-6 * std::abs(iw[3]),
                "the run invariant is the oscillator energy");
        for (const auto& r : invariant.runs) require(r.peakToPeak < 1e-6 * std::abs(iw[3]) && r.sampledInvariance < 1e-6, "constant along each run");
        std::cout << "trajectory kernels tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

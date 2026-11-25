// The catalog eigenfunction library (SystemEigenfunctionLibrary.cpp): every entry satisfies the
// generator identity grad(psi).f = lambda psi against the configured equations, the elliptic
// integrals match reference values, the exact expression gradients match finite differences, and a
// shell eigenfunction propagates as exp(i Omega(H) t) along an integrated orbit across its branch cut.
#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
using Complex = std::complex<double>;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::unique_ptr<c::ODESystem> make(const std::string& id, std::map<std::string, double> parameters) {
    c::SystemConfiguration config;
    config.identifier = id;
    config.parameters = std::move(parameters);
    return c::createSystem(config);
}

// Residual of each entry on the orbit through each state (shell eigenvalues taken on that orbit).
void check(const c::ODESystem& system, const std::vector<std::string>& ids, const std::vector<std::vector<double>>& states,
           double tolerance = 1e-8) {
    for (const auto& x : states) {
        const auto psi = c::libraryEigenfunctions(system, ids, x);
        const auto residual = c::eigenfunctionGeneratorResidual(system, psi, x);
        for (std::size_t j = 0; j < psi.size(); ++j)
            require(std::isfinite(residual[j]) && residual[j] < tolerance,
                    system.descriptor().identifier + "/" + ids[j] + " generator residual " + std::to_string(residual[j]));
    }
}

}  // namespace

int main() {
    try {
        // Reference values: K(1/2), E(1/2) (DLMF 19.2), F(pi/2|m) = K(m), continuation F(phi + pi) = F + 2K.
        require(std::abs(c::ellipticK(0.5) - 1.8540746773013719) < 1e-13, "K(1/2)");
        require(std::abs(c::ellipticE(0.5) - 1.3506438810476755) < 1e-13, "E(1/2)");
        require(std::abs(c::ellipticF(std::acos(-1.0) / 2, 0.7) - c::ellipticK(0.7)) < 1e-13, "F(pi/2|m) = K(m)");
        require(std::abs(c::ellipticF(2.0 + std::acos(-1.0), 0.3) - c::ellipticF(2.0, 0.3) - 2 * c::ellipticK(0.3)) < 1e-12,
                "F continuation");
        require(std::abs(c::ellipticEIncomplete(std::acos(-1.0) / 2, 0.4) - c::ellipticE(0.4)) < 1e-13, "E(pi/2|m) = E(m)");

        // Exact gradients of atan2 and the elliptic integrals against central differences.
        const auto expr = c::AnalyticExpression::parse("ellipf(atan2(y, x), 0.3*x^2) + ellipe(y, 0.2 + 0.1*y) + ellipk(0.4*x*y)",
                                                       {"x", "y"});
        const double p[]{0.7, 1.1};
        for (std::size_t axis = 0; axis < 2; ++axis) {
            double a[]{p[0], p[1]}, b[]{p[0], p[1]};
            a[axis] += 1e-6; b[axis] -= 1e-6;
            const double numeric = (expr.evaluate(a) - expr.evaluate(b)) / 2e-6;
            require(std::abs(expr.derivative(axis).evaluate(p) - numeric) < 1e-7, "elliptic expression gradient");
        }

        // Pendulum: libration and rotation phases (shell), energy, damped principal modes.
        const auto pendulum = make("nonlinear_pendulum", {{"damping_n_m_s", 0.0}});
        check(*pendulum, {"energy", "libration_phase"}, {{0.5, 0.3}, {-1.0, 1.2}, {2.0, -0.5}, {0.1, -3.0}});
        check(*pendulum, {"rotation_phase"}, {{0.3, 7.0}, {-2.0, -8.0}});
        const auto damped = make("nonlinear_pendulum", {{"damping_n_m_s", 0.4}});
        require(c::eigenfunctionGeneratorResidual(*damped, c::libraryEigenfunctions(*damped, {"principal", "principal_conjugate"}),
                                                  {1e-5, -2e-5})[0] < 1e-4, "pendulum principal mode near the origin");

        // Integrate a libration orbit across the atan2 branch cut: psi(x(t)) = exp(i Omega t) psi(x(0)).
        {
            const std::vector<double> x0{1.2, 0.0};
            const auto psi = c::libraryEigenfunctions(*pendulum, {"libration_phase"}, x0).front();
            const double dt = 1e-3;
            std::vector<double> x = x0, k1(2), k2(2), k3(2), k4(2), y(2);
            const Complex start = psi.evaluate(x0.data()).front();
            for (int step = 1; step <= 4000; ++step) {
                pendulum->rhs(0, x.data(), k1.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt / 2 * k1[i];
                pendulum->rhs(0, y.data(), k2.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt / 2 * k2[i];
                pendulum->rhs(0, y.data(), k3.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt * k3[i];
                pendulum->rhs(0, y.data(), k4.data());
                for (int i = 0; i < 2; ++i) x[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
                if (step % 250) continue;
                const Complex expected = psi.propagate(start, step * dt);
                require(std::abs(psi.evaluate(x.data()).front() - expected) < 1e-8, "libration phase propagates along the orbit");
            }
        }

        // Duffing double well: outer and single-well phases, energy.
        const auto duffing = make("duffing", {{"delta", 0.0}, {"alpha", 1.0}, {"beta", 1.0}});
        check(*duffing, {"energy", "outer_phase"}, {{0.0, 1.0}, {1.5, 0.3}, {-2.0, -0.5}});
        check(*duffing, {"well_phase"}, {{1.0, 0.1}, {-1.2, 0.2}, {0.6, -0.2}, {-0.7, -0.3}});

        // Linear oscillator modes are exact for any damping.
        const auto harmonic = make("harmonic_oscillator", {{"mass_kg", 1.0}, {"stiffness_n_m", 2.0}, {"damping_kg_s", 0.3}});
        check(*harmonic, {"mode", "mode_conjugate"}, {{0.4, -0.9}, {-2.0, 1.5}});

        check(*make("isotonic_oscillator", {}), {"energy", "isochronous_mode"}, {{1.4, 0.2}, {0.8, -0.5}});
        check(*make("morse_oscillator", {}), {"energy", "phase"}, {{0.3, 0.2}, {-0.2, 0.5}});
        check(*make("poschl_teller", {}), {"energy", "phase"}, {{0.3, 0.2}, {-0.5, -0.4}});
        check(*make("cosh_oscillator", {}), {"energy", "phase"}, {{0.5, 0.3}, {-1.0, 0.7}});
        check(*make("radial_twist", {}), {"action", "rotation"}, {{1.0, 0.5}, {-0.3, 0.8}});
        check(*make("kepler_planar", {}), {"energy", "angular_momentum", "mean_anomaly"},
              {{1.0, 0.0, 0.1, 0.8}, {0.5, 0.6, -0.7, 0.4}});
        check(*make("euler_top", {}), {"casimir", "energy"}, {{0.2, 0.3, 1.0}});
        check(*make("coupled_duffing_hamiltonian", {{"coupling", 0.0}}), {"energy", "energy_1", "energy_2"}, {{0.7, 0.4, 0.1, -0.2}});
        check(*make("triple_well_duffing", {}), {"energy"}, {{0.3, 0.1}, {1.1, -0.2}});
        check(*make("three_leaf_hamiltonian_clover", {}), {"energy"}, {{0.5, 0.1}});

        // Shell entries need a reference orbit; damping makes the conservative entries inapplicable.
        bool needsReference = false;
        try { (void)c::libraryEigenfunctions(*pendulum, {"libration_phase"}); } catch (const std::invalid_argument&) { needsReference = true; }
        require(needsReference, "shell entries require a reference state");
        bool rejected = false;
        try { (void)c::libraryEigenfunctions(*damped, {"energy"}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "damped pendulum has no energy invariant");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "system eigenfunction library tests passed\n";
    return 0;
}

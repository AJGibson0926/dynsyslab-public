// The system catalogue (native/catalog/systems/*.json): every entry satisfies the
// generator identity grad(psi).f = lambda psi against the configured equations, the elliptic
// integrals match reference values, the exact expression gradients match finite differences. Shell
// entries are spectral bands: their profile satisfies grad(phi).f = i level(x) phi, with level(x) =
// n Omega(H(x)), and propagates as exp(i level(x0) t) along an integrated orbit across its branch cut.
#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/PolynomialRoots.hpp>

#include <array>

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

// Each entry at each state: an eigenfunction's generator residual, or a band's level-set identity
// |grad(phi).f - i level(x) phi| / |phi| for a shell entry.
void check(const c::ODESystem& system, const std::vector<std::string>& ids, const std::vector<std::vector<double>>& states,
           double tolerance = 1e-8) {
    for (const auto& id : ids) {
        const auto bands = c::libraryBands(system, {id});
        for (const auto& x : states) {
            double residual = 0;
            if (bands.empty()) {
                residual = c::eigenfunctionGeneratorResidual(system, c::libraryEigenfunctions(system, {id}), x).front();
            } else {
                const auto& g = bands.front().eigenfunctions;
                require(g.levelSet() && bands.front().support == c::band_support::imaginaryAxis, id + " is a level-set band");
                std::vector<double> f(x.size());
                system.rhs(0.0, x.data(), f.data());
                const auto phi = g.functions->evaluate(x.data()).front();
                const auto gradient = g.functions->jacobian(x.data());
                Complex lie{};
                for (std::size_t i = 0; i < x.size(); ++i) lie += gradient[i] * f[i];
                const auto lambda = c::bandPoint(bands.front(), g.level->evaluate(x.data()).front().real());
                residual = std::abs(lie - lambda * phi) / std::max(std::abs(phi), 1e-300);
            }
            require(std::isfinite(residual) && residual < tolerance,
                    system.descriptor().identifier + "/" + id + " generator residual " + std::to_string(residual));
        }
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

        // The elliptic integral of the third kind against quadrature (finite-interval references, ~1e-15),
        // including the continuation past pi/2 and a full turn.
        for (const auto& [phi, n, m, value] : std::vector<std::array<double, 4>>{
                 {0.7, -0.4, 0.3, 0.6773346322830175}, {2.5, -1.7, 0.8, 1.936485818810892}, {-1.2, 0.5, 0.2, -1.579054926250551},
                 {3.0, -0.2, 0.5, 3.2170229595274775}, {std::acos(-1.0), -0.3, 0.6, 3.3695660523974014}})
            require(std::abs(c::ellipticPi(phi, n, m) - value) < 4e-15 * std::max(1.0, std::abs(value)), "Pi(n; phi|m) reference");
        require(std::abs(c::ellipticPi(0.9, 0, 0.4) - c::ellipticF(0.9, 0.4)) < 1e-15, "Pi(0; phi|m) = F(phi|m)");
        require(std::abs(c::carlsonRC(0.5, 1.5) - std::acos(std::sqrt(1.0 / 3))) < 1e-15, "R_C closed form");
        // Exact gradients of ellippi in phi, n and m, and of polyroot, against central differences;
        // named definitions share an intermediate.
        {
            const auto pi3 = c::AnalyticExpression::parse("ellippi(x + 0.3*y, -0.4 + 0.2*x*y, 0.3 + 0.1*y)", {"x", "y"});
            const auto root = c::AnalyticExpression::parse("c = 6 + y; polyroot(0, 1, 0, -7 - x, c)", {"x", "y"});
            const double q[]{0.8, 0.6};
            for (const auto* e : {&pi3, &root})
                for (std::size_t axis = 0; axis < 2; ++axis) {
                    double a[]{q[0], q[1]}, b[]{q[0], q[1]};
                    a[axis] += 1e-6; b[axis] -= 1e-6;
                    const double numeric = (e->evaluate(a) - e->evaluate(b)) / 2e-6;
                    require(std::abs(e->derivative(axis).evaluate(q) - numeric) < 1e-7, "exact gradient of ellippi / polyroot");
                }
            const double z[]{0.0, 0.0};  // x^3 - 7x + 6 = (x - 1)(x - 2)(x + 3)
            require(std::abs(root.evaluate(z) - 2) < 1e-14 && std::abs(c::realPolynomialRoot(2, {1, 0, -7, 6}) + 3) < 1e-14,
                    "polyroot: the largest and smallest real roots");
            require(std::isnan(c::realPolynomialRoot(0, {1, 0, 1})), "x^2 + 1 has no real root");
            const auto shared = c::AnalyticExpression::parse("s = x*x; t = s + 1; t*t - s", {"x"});
            const double x1[]{2.0};
            require(std::abs(shared.evaluate(x1) - 21) < 1e-14, "named definitions");
            require(std::abs(c::AnalyticExpression::fromJson(pi3.toJson()).evaluate(q) - pi3.evaluate(q)) < 1e-15
                        && std::abs(c::AnalyticExpression::fromJson(root.toJson()).evaluate(q) - root.evaluate(q)) < 1e-15,
                    "ellippi and polyroot round trip as JSON");
        }

        // Pendulum: libration and rotation phases (shell), energy, damped principal modes.
        const auto pendulum = make("nonlinear_pendulum", {{"damping_n_m_s", 0.0}});
        check(*pendulum, {"energy", "libration_phase"}, {{0.5, 0.3}, {-1.0, 1.2}, {2.0, -0.5}, {0.1, -3.0}});
        check(*pendulum, {"rotation_phase"}, {{0.3, 7.0}, {-2.0, -8.0}});
        const auto damped = make("nonlinear_pendulum", {{"damping_n_m_s", 0.4}});
        require(c::eigenfunctionGeneratorResidual(*damped, c::libraryEigenfunctions(*damped, {"principal", "principal_conjugate"}),
                                                  {1e-5, -2e-5})[0] < 1e-4, "pendulum principal mode near the origin");

        // Integrate a libration orbit across the atan2 branch cut: phi(x(t)) = exp(i level(x0) t) phi(x(0)).
        {
            const std::vector<double> x0{1.2, 0.0};
            const auto band = c::libraryBands(*pendulum, {"libration_phase"}).front();
            const auto& profile = *band.eigenfunctions.functions;
            const auto lambda = c::bandPoint(band, band.eigenfunctions.level->evaluate(x0.data()).front().real());
            const double dt = 1e-3;
            std::vector<double> x = x0, k1(2), k2(2), k3(2), k4(2), y(2);
            const Complex start = profile.evaluate(x0.data()).front();
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
                const Complex expected = std::exp(lambda * (step * dt)) * start;
                require(std::abs(profile.evaluate(x.data()).front() - expected) < 1e-8, "libration phase propagates along the orbit");
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

        // The inviscid RPE phase in elliptic integrals: the band's identity at states on small and large
        // orbits and on both branches, its frequency against the native analytic orbit's period, and the
        // phase carried along an integrated orbit through both turning radii.
        {
            const auto rpe = make("rpe", {{"mu_pa_s", 0.0}});
            check(*rpe, {"phase"}, {{1.08, 0.0}, {1.3, -0.4}, {0.95, 0.2}, {2.0, -1.1}, {0.6, 1.5}, {2.8, 0.0}}, 1e-9);
            const auto band = c::libraryBands(*rpe, {"phase"}).front();
            const double x0[]{1.08, 0.0}, x1[]{2.8, 0.0};
            const double pi = std::acos(-1.0);
            require(std::abs(band.eigenfunctions.level->evaluate(x0).front().real() - 2 * pi / 1.0021098121065386) < 1e-12
                        && std::abs(band.eigenfunctions.level->evaluate(x1).front().real() - 2 * pi / 1.7021651650443048) < 1e-12,
                    "the RPE frequency is the analytic orbit's");
            const auto& profile = *band.eigenfunctions.functions;
            std::vector<double> x{2.0, 0.0}, k1(2), k2(2), k3(2), k4(2), y(2);
            const Complex start = profile.evaluate(x.data()).front();
            const auto lambda = c::bandPoint(band, band.eigenfunctions.level->evaluate(x.data()).front().real());
            const double dt = 1e-4;
            for (int step = 1; step <= 20000; ++step) {
                rpe->rhs(0, x.data(), k1.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt / 2 * k1[i];
                rpe->rhs(0, y.data(), k2.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt / 2 * k2[i];
                rpe->rhs(0, y.data(), k3.data());
                for (int i = 0; i < 2; ++i) y[i] = x[i] + dt * k3[i];
                rpe->rhs(0, y.data(), k4.data());
                for (int i = 0; i < 2; ++i) x[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
                if (step % 1000) continue;
                require(std::abs(profile.evaluate(x.data()).front() - std::exp(lambda * (step * dt)) * start) < 1e-7,
                        "the RPE phase propagates along the orbit: t = " + std::to_string(step * dt));
            }
        }

        // Shell entries are bands, not eigenpairs; the library bank keeps them as its spectrum's bands.
        bool refused = false;
        try { (void)c::libraryEigenfunctions(*pendulum, {"libration_phase"}); } catch (const std::invalid_argument&) { refused = true; }
        require(refused, "a shell entry is not an eigenfunction");
        const auto bank = c::libraryEigenfunctionBank(*pendulum, {"energy", "libration_phase"});
        require(bank.spectrum.size() == 1 && bank.spectrum.bands().size() == 1 && bank.spectrum.band("libration_phase")
                    && bank.spectrum.band("libration_phase")->eigenfunctions.levelSet(),
                "the bank: the energy as an eigenpair, the libration phase as a band");
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

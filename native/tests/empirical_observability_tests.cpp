// Empirical (nonlinear-simulation) observability gramian: on a linear plant it must
// equal the finite-horizon observability gramian; scaling, Cramer-Rao bounds and the
// instantaneous output Jacobian are checked analytically.
#include <dynsyslab/core/EmpiricalObservability.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/EigenfunctionObservability.hpp>
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <iostream>
#include <string>

using namespace dynsyslab::core;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) { if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; } }
}  // namespace

int main() {
    try {
        // Linear Duffing: x'' = -delta x' + alpha x - beta x^3 with alpha = -1, beta = 0, delta = 0.5.
        SystemConfiguration config; config.identifier = "duffing";
        config.parameters = {{"alpha", -1.0}, {"beta", 0.0}, {"delta", 0.5}};
        const auto system = createSystem(config);
        const auto catalog = system->observableCatalog();
        const auto& position = catalog.at("position");
        EmpiricalObservabilityOptions o; o.horizon = 3.0; o.samples = 601; o.substeps = 10;
        const auto r = empiricalObservability(*system, position, {0.3, -0.2}, 0.0, o);
        const auto exact = finiteObservabilityGramian(makeStateSpace({0, 1, -1, -0.5}, {0, 1}, {1, 0}, {}, 2, 1, 1), 3.0);
        double err = 0, scale = 0;
        for (int k = 0; k < 4; ++k) { err = std::max(err, std::abs(r.gramian[k] - exact[k])); scale = std::max(scale, std::abs(exact[k])); }
        check(r.valid && err < 1e-5 * scale, "empirical gramian equals the finite-horizon gramian of a linear plant");
        check(r.eigenvalues[0] >= r.eigenvalues[1] && r.eigenvalues[1] > 0, "eigenvalues descending and positive");
        check(std::abs(r.outputJacobian[0] - 1.0) < 1e-9 && std::abs(r.outputJacobian[1]) < 1e-9, "dh/dx of position is (1, 0)");
        // Cramer-Rao bound equals sqrt(V (W^-1)_ii)
        const double det = r.gramian[0] * r.gramian[3] - r.gramian[1] * r.gramian[2];
        const auto b = estimationBoundStd(r, 0.04);
        check(std::abs(b[0] - std::sqrt(0.04 * r.gramian[3] / det)) < 1e-9 * b[0], "Cramer-Rao bound, first state");
        check(std::abs(b[1] - std::sqrt(0.04 * r.gramian[0] / det)) < 1e-9 * b[1], "Cramer-Rao bound, second state");
        // scaling x = S z multiplies W by S on both sides
        const auto sc = scaledObservability(r, {2.0, 0.5});
        check(std::abs(sc.gramian[1] - r.gramian[1]) < 1e-12 && std::abs(sc.gramian[0] - 4 * r.gramian[0]) < 1e-12,
              "scaled gramian is S W S");
        // nonlinear (beta = 1): still symmetric positive and finite
        config.parameters["beta"] = 1.0;
        const auto hard = createSystem(config);
        const auto rn = empiricalObservability(*hard, hard->observableCatalog().at("position"), {0.8, 0.0}, 0.0, o);
        check(rn.valid && rn.eigenvalues[1] > 0 && std::abs(rn.gramian[1] - rn.gramian[2]) < 1e-15, "nonlinear gramian symmetric positive");
        // Discrete modal obsv of a two-mode K-ROM: O = [c; c M], M = diag(exp(lambda dt)).
        {
            auto model = KoopmanModel::fromPolynomialCoefficients(config, 2, {1, 0, 0, 1}, {2., 1., -1., 1.}, {1., -2.});
            model.setReconstruction("sensor", {{1. / 3, -1. / 3}, 1, 2, 0, 0, 0});
            const double dt = 0.1;
            const auto m = modalObservabilityMatrix(model, "sensor", dt);
            const double a = 1. / 3, b = -1. / 3, m1 = std::exp(dt), m2 = std::exp(-2 * dt);
            // singular values of [[a, b], [a m1, b m2]] from its 2x2 Gram matrix
            const double g11 = a * a * (1 + m1 * m1), g22 = b * b * (1 + m2 * m2), g12 = a * b * (1 + m1 * m2);
            const double tr = g11 + g22, det = g11 * g22 - g12 * g12, disc = std::sqrt(tr * tr / 4 - det);
            check(m.modes == 2 && m.blocks == 2 && m.singularValues.size() == 2, "modal obsv shape (one coordinate per mode)");
            check(std::abs(m.singularValues[0] - std::sqrt(tr / 2 + disc)) < 1e-12 &&
                  std::abs(m.singularValues[1] - std::sqrt(tr / 2 - disc)) < 1e-12, "modal obsv singular values");
            check(std::abs(std::abs(m.matrix[2]) - a * m1) < 1e-12, "second block row is c M");
            const auto three = modalObservabilityMatrix(model, "sensor", dt, 3);
            check(three.blocks == 3 && three.matrix.size() == 6, "explicit block count");
            // Per-eigenfunction observability over T: G_ij = c_i c_j (exp((l_i + l_j) T) - 1) / (l_i + l_j).
            EigenfunctionObservabilityOptions eo; eo.horizon = 0.5; eo.noiseIntensity = 0.01; eo.processCovariance = 1.0;
            const auto e = eigenfunctionObservability(model, "sensor", eo);
            const double T = 0.5, l1 = 1.0, l2 = -2.0;
            const auto G = [&](double li, double lj, double ci, double cj) { return ci * cj * (std::exp((li + lj) * T) - 1) / (li + lj); };
            const double G11 = G(l1, l1, a, a), G22 = G(l2, l2, b, b), G12 = G(l1, l2, a, b), D = G11 * G22 - G12 * G12;
            check(std::abs(e.signal[0] - std::sqrt(G11)) < 1e-10 && std::abs(e.signal[1] - std::sqrt(G22)) < 1e-10, "signal = sqrt(G_nn)");
            check(std::abs(e.distinguishable[0] - std::sqrt(D / G22)) < 1e-8 * std::sqrt(D / G22) &&
                  std::abs(e.distinguishable[1] - std::sqrt(D / G11)) < 1e-8 * std::sqrt(D / G11), "distinguishable = 1/sqrt((G^-1)_nn)");
            check(std::abs(e.boundStd[0] - std::sqrt(0.01 * G22 / D)) < 1e-8 * e.boundStd[0], "Cramer-Rao bound");
            check(e.isolation[0] > 0 && e.isolation[0] < 1 && e.partner[0] == 0 && e.partner[1] == 1, "isolation in (0, 1), real modes self-paired");
            check(e.pairBoundStd[0] == e.boundStd[0] && e.kalmanStd.size() == 2 && e.kalmanStd[0] > 0, "pair bound and Kalman std");
        }
        // a run that leaves the validity domain is flagged rather than integrated through
        bool threw = false;
        try { (void)empiricalObservability(*system, position, {0.3}, 0.0, o); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "state dimension checked");
    } catch (const std::exception& error) {
        std::cerr << "empirical observability tests threw: " << error.what() << "\n";
        return 1;
    }
    if (failures) { std::cerr << failures << " empirical observability checks failed\n"; return 1; }
    std::cout << "empirical observability tests passed\n";
    return 0;
}

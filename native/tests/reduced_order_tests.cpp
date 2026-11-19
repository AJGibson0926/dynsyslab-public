// KROM -> LTI realization and the balanced-truncation estimation/control sidecars.
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/ReducedOrderControl.hpp>
#include <dynsyslab/core/ReducedOrderEstimation.hpp>
#include <dynsyslab/core/SystemNorms.hpp>

#include <cmath>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace dynsyslab::core;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) { if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; } }
}  // namespace

int main() {
    try {
        SystemConfiguration config; config.identifier = "lorenz";
        // Three linear eigenfunctions psi_i = x_i with eigenvalues -1, -8, -60.
        auto model = KoopmanModel::fromPolynomialCoefficients(config, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1},
            {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}},
            {{-1, 0}, {-8, 0}, {-60, 0}});
        model.setReconstruction("echo", {{1., 0.5, 0.05}, 1, 3});
        model.setControlMap("drive", ConstantControlMap{{1., 1., 0.2}, 3, 1, "test"});
        KoopmanStateSpaceOptions realization; realization.measurements = {"echo"}; realization.controlMap = "drive";
        const auto system = koopmanStateSpace(model, realization);
        check(system.states == 3 && system.inputs == 1 && system.outputs == 1, "KROM realization dimensions");
        check(std::abs(system.A[0] + 1) < 1e-14 && std::abs(system.A[4] + 8) < 1e-14 && std::abs(system.B[2] - 0.2) < 1e-14
              && std::abs(system.C[1] - 0.5) < 1e-14, "KROM realization entries");
        const double expectedDc = 1.0 / 1 + 0.5 / 8 + 0.05 * 0.2 / 60;
        check(std::abs(dcGain(system)[0] - expectedDc) < 1e-12, "KROM transfer function DC gain");
        auto discrete = realization; discrete.domain = TimeDomain::discrete; discrete.sampleTime = 0.01;
        const auto sampled = koopmanStateSpace(model, discrete);
        check(std::abs(sampled.A[0] - std::exp(-0.01)) < 1e-14, "sampled KROM transition");
        auto missing = realization; missing.controlMap = "exact";
        bool rejected = false;
        try { (void)koopmanStateSpace(model, missing); } catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "state-dependent map without an operating point is rejected");

        BalancedReductionOptions reduce; reduce.order = 1;
        const std::vector<double> Q{1, 0, 0, 0, 1, 0, 0, 0, 1};
        const auto lqr = designReducedKoopmanLQR(model, realization, reduce, Q, {0.1});
        check(lqr.reduction.order == 1 && lqr.liftedGain.size() == 3 && lqr.liftedClosedLoopStable, "reduced KLQR lifted gain stabilizes the full model");
        const auto lqe = designReducedKoopmanLQE(model, {{"echo"}, "", {}, {}, TimeDomain::continuous, 0.0}, reduce, Q, {0.01});
        check(lqe.reduction.order == 1 && lqe.liftedGain.size() == 3 && lqe.liftedObserverStable, "reduced KLQE lifted gain");
        check(lqe.reduction.errorBound >= lqe.reduction.lowerBound && lqe.reduction.lowerBound > 0, "reduction bounds reported");

        // KLQE stage option: the design's gain is replaced by the lifted reduced gain.
        const auto fullLqe = designKoopmanLQEFromMeasurements(model, "echo", Q, {0.01});
        auto reducedLqe = fullLqe;
        const auto applied = reduceKoopmanLQEDesign(reducedLqe, reduce, Q, {0.01});
        double gainGap = 0;
        for (std::size_t k = 0; k < 3; ++k) gainGap = std::max(gainGap, std::abs(reducedLqe.solution.gain[k] - lqe.liftedGain[k]));
        check(applied.reduction.order == 1 && gainGap < 1e-12, "reduceKoopmanLQEDesign writes the lifted gain");
        bool stableObserver = true;
        for (const auto& pole : reducedLqe.solution.observerPoles) stableObserver = stableObserver && pole.real() < 0;
        check(reducedLqe.solution.observerPoles.size() == 3 && stableObserver, "lifted KLQE observer poles are stable");
        BalancedReductionOptions keepAll; keepAll.order = 3;
        auto sameLqe = fullLqe;
        (void)reduceKoopmanLQEDesign(sameLqe, keepAll, Q, {0.01});
        gainGap = 0;
        for (std::size_t k = 0; k < 3; ++k) gainGap = std::max(gainGap, std::abs(sameLqe.solution.gain[k] - fullLqe.solution.gain[k]));
        check(gainGap < 1e-8 * (1 + std::abs(fullLqe.solution.gain[0])), "full-order balanced KLQE equals the KLQE gain");

        // KLQR option: the controller designs on (A, B, I) reduced and applies K_r P.
        const auto shared = std::make_shared<const KoopmanModel>(model);
        const double x[3] = {0.3, -0.2, 0.1};
        KoopmanLQRController full(shared, Q, {0.1}, {}, "drive");
        KoopmanLQRController balancedFull(shared, Q, {0.1}, {}, "drive");
        balancedFull.setBalancedReduction(keepAll);
        const double uFull = full.value(0.0, x).front(), uSame = balancedFull.value(0.0, x).front();
        check(std::abs(uFull - uSame) < 1e-8 * (1 + std::abs(uFull)), "full-order balanced KLQR equals KLQR");
        KoopmanLQRController reducedController(shared, Q, {0.1}, {}, "drive");
        reducedController.setBalancedReduction(reduce);
        (void)reducedController.value(0.0, x);
        KoopmanStateSpaceOptions performance; performance.controlMap = "drive";  // C_w = I
        const auto expected = designReducedKoopmanLQR(model, performance, reduce, Q, {0.1});
        gainGap = 0;
        for (std::size_t k = 0; k < 3; ++k)
            gainGap = std::max(gainGap, std::abs(reducedController.lastSolution().gain[k] - expected.liftedGain[k]));
        check(reducedController.lastReducedOrder() == 1 && gainGap < 1e-12, "KLQR controller uses the lifted reduced gain");
    } catch (const std::exception& error) {
        std::cerr << "reduced order tests threw: " << error.what() << "\n";
        return 1;
    }
    if (failures) { std::cerr << failures << " reduced order checks failed\n"; return 1; }
    std::cout << "reduced order tests passed\n";
    return 0;
}

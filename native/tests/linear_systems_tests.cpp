// Analytic checks for the classical linear-systems toolkit.
#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/core/SystemNorms.hpp>
#include <dynsyslab/core/SystemStructure.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace dynsyslab::core;
using C = std::complex<double>;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; }
}
bool near(double a, double b, double tol) { return std::abs(a - b) <= tol * std::max(1.0, std::abs(b)); }
bool hasPole(const std::vector<C>& poles, C p, double tol = 1e-8) {
    return std::any_of(poles.begin(), poles.end(), [&](C q) { return std::abs(q - p) <= tol * std::max(1.0, std::abs(p)); });
}
StateSpaceModel ss(std::vector<double> A, std::vector<double> B, std::vector<double> Cm,
                   std::size_t n, std::size_t m, std::size_t p, std::vector<double> D = {}) {
    return makeStateSpace(std::move(A), std::move(B), std::move(Cm), std::move(D), n, m, p);
}
// Three decoupled first-order lags with distinct time scales, rotated so that A is dense.
StateSpaceModel denseThreeMode() {
    const std::vector<double> T{1, 0.3, -0.2, 0.1, 1, 0.4, -0.3, 0.2, 1};
    const auto modal = ss({-1, 0, 0, 0, -5, 0, 0, 0, -30}, {1, 2, 3}, {1, 0.5, 0.2}, 3, 1, 1);
    // inverse of T by explicit solve through the toolkit's similarity check
    const double a = T[0], b = T[1], c = T[2], d = T[3], e = T[4], f = T[5], g = T[6], h = T[7], i = T[8];
    const double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    const std::vector<double> Ti{(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det,
                                 (f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det,
                                 (d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det};
    return similarityTransform(modal, Ti, T);
}
double maxResponseGap(const StateSpaceModel& a, const StateSpaceModel& b, const std::vector<double>& w) {
    const auto ra = frequencyResponse(a, w), rb = frequencyResponse(b, w);
    double gap = 0.0;
    for (std::size_t k = 0; k < ra.response.size(); ++k) gap = std::max(gap, std::abs(ra.response[k] - rb.response[k]));
    return gap;
}
}  // namespace

int main() {
    try {
        // --- Lyapunov and gramians ------------------------------------------------
        const auto lag = ss({-1, 0, 0, -2}, {1, 2}, {1, 2}, 2, 1, 1);
        const auto Wc = controllabilityGramian(lag);
        check(near(Wc[0], 0.5, 1e-12) && near(Wc[1], 2.0 / 3, 1e-12) && near(Wc[3], 1.0, 1e-12), "continuous Lyapunov closed form");
        auto dlag = makeStateSpace({0.5, 0, 0, 0.25}, {1, 2}, {1, 2}, {}, 2, 1, 1, TimeDomain::discrete, 0.1);
        const auto Wd = controllabilityGramian(dlag);
        check(near(Wd[0], 1 / 0.75, 1e-12) && near(Wd[1], 2 / 0.875, 1e-12) && near(Wd[3], 4 / 0.9375, 1e-12), "Stein closed form");
        const auto dense = denseThreeMode();
        const auto W0 = observabilityGramian(dense), Wfin = finiteObservabilityGramian(dense, 40.0);
        double gap = 0; for (std::size_t k = 0; k < 9; ++k) gap = std::max(gap, std::abs(W0[k] - Wfin[k]));
        check(gap < 1e-10, "finite Van Loan gramian converges to the Lyapunov gramian");
        const auto Wshort = finiteControllabilityGramian(lag, 0.7);
        check(near(Wshort[0], 0.5 * (1 - std::exp(-1.4)), 1e-10) && near(Wshort[1], 2.0 / 3 * (1 - std::exp(-2.1)), 1e-10), "finite gramian closed form");
        EmpiricalGramianOptions eo; eo.horizon = 0.7; eo.samples = 4001; eo.scales = {0.5, 2.0};
        const auto We = empiricalControllabilityGramian(lag, eo);
        check(near(We[0], Wshort[0], 1e-6) && near(We[1], Wshort[1], 1e-6), "empirical controllability gramian");
        const auto Woe = empiricalObservabilityGramian(lag, eo), Wof = finiteObservabilityGramian(lag, 0.7);
        check(near(Woe[3], Wof[3], 1e-6) && near(Woe[1], Wof[1], 1e-6), "empirical observability gramian");
        EmpiricalGramianOptions de; de.horizon = 200;
        const auto Wde = empiricalControllabilityGramian(dlag, de);
        check(near(Wde[1], Wd[1], 1e-12), "discrete empirical gramian");
        const auto sg = snapshotGramian({1, 2, 3, 4}, 2, 2, {1, 0.5});
        check(near(sg[0], 1 + 0.5 * 4, 1e-14) && near(sg[1], 3 + 0.5 * 8, 1e-14), "snapshot gramian weights");

        // --- Hankel singular values and balanced truncation ---------------------------
        const auto first = ss({-2}, {3}, {0.5}, 1, 1, 1);
        check(near(hankelSingularValues(first)[0], 1.5 * 0.5 / 2.0, 1e-12), "first-order HSV = |bc|/(2a)");
        const auto hsv = hankelSingularValues(dense);
        const auto bal = balancedRealization(dense);
        const auto bWc = controllabilityGramian(bal.balanced), bWo = observabilityGramian(bal.balanced);
        for (std::size_t i = 0; i < 3; ++i)
            check(near(bWc[i * 4], hsv[i], 1e-9) && near(bWo[i * 4], hsv[i], 1e-9), "balanced gramians equal diag(HSV)");
        check(std::abs(bWc[1]) < 1e-9 && std::abs(bWo[2]) < 1e-9, "balanced gramians are diagonal");
        for (auto method : {BalancedReductionMethod::truncation, BalancedReductionMethod::residualization}) {
            BalancedReductionOptions ro; ro.order = 1; ro.method = method;
            const auto red = balancedReduction(dense, ro);
            const auto error = parallelConnection(dense, red.reduced, -1.0);
            const double e = hInfinityNorm(error).value;
            check(e >= red.lowerBound * (1 - 1e-6) && e <= red.errorBound * (1 + 1e-6), std::string("balanced error bounds: ") + balancedReductionMethodName(method));
            if (method == BalancedReductionMethod::residualization)
                check(near(dcGain(red.reduced)[0], dcGain(dense)[0], 1e-10), "residualization preserves DC gain");
        }
        BalancedReductionOptions tol; tol.errorTolerance = 2.0 * (hsv[1] + hsv[2]) * 1.01;
        check(balancedReduction(dense, tol).order == 1, "order chosen from the error tolerance");
        auto unstable = ss({-1, 0, 0, 0, 0.5, 0, 0, 0, -20}, {1, 1, 1}, {1, 1, 0.01}, 3, 1, 1);
        BalancedReductionOptions keep; keep.order = 2;
        const auto ured = balancedReduction(unstable, keep);
        check(ured.unstableOrder == 1 && ured.stableOrder == 1 && hasPole(systemPoles(ured.reduced), {0.5, 0}), "unstable mode retained exactly");
        check(hasPole(systemPoles(ured.reduced), {-1, 0}, 1e-2), "dominant stable mode retained");
        const double uerr = hInfinityNorm(parallelConnection(ured.stablePart, ured.reducedStablePart, -1.0)).value;
        check(uerr >= ured.lowerBound * (1 - 1e-6) && uerr <= ured.errorBound * (1 + 1e-6), "unstable-plant reduction error bounds");
        const auto split = splitStableUnstable(unstable);
        check(split.stableStates == 2 && split.unstableStates == 1, "stable/unstable dimensions");

        // --- Data-driven realization --------------------------------------------------
        const auto plant = makeStateSpace({0.9, 0.2, 0, 0.5}, {1, 1}, {1, 0.3}, {0.1}, 2, 1, 1, TimeDomain::discrete, 0.01);
        const std::size_t K = 41;
        std::vector<double> markov(K);
        auto impulse = impulseResponse(plant, 0, 0, K);
        for (std::size_t k = 0; k < K; ++k) markov[k] = impulse.outputs[k];
        EraOptions era; era.sampleTime = 0.01;
        const auto realized = eigensystemRealization(markov, K, 1, 1, era);
        check(realized.order == 2 && hasPole(systemPoles(realized.model), {0.9, 0}, 1e-8) && hasPole(systemPoles(realized.model), {0.5, 0}, 1e-8), "ERA recovers poles and order");
        std::vector<double> time(400), u(400);
        for (std::size_t k = 0; k < 400; ++k) { time[k] = 0.01 * k; u[k] = std::sin(0.37 * k) + 0.5 * std::cos(1.91 * static_cast<double>((k * k) % 17)); }
        const auto record = simulateLinear(plant, time, u);
        const auto okid = observerKalmanMarkovParameters(u, record.outputs, 400, 1, 1, 6, K);
        double okidGap = 0; for (std::size_t k = 0; k < K; ++k) okidGap = std::max(okidGap, std::abs(okid[k] - markov[k]));
        check(okidGap < 1e-8, "OKID recovers Markov parameters from arbitrary input");
        {   // balanced POD with trapezoid-weighted impulse snapshots
            const std::size_t N = 6001; const double horizon = 12.0, dt = horizon / (N - 1);
            const auto primal = impulseResponse(dense, 0, horizon, N);
            auto adjointModel = makeStateSpace({dense.A[0], dense.A[3], dense.A[6], dense.A[1], dense.A[4], dense.A[7], dense.A[2], dense.A[5], dense.A[8]},
                                               dense.C, dense.B, {}, 3, 1, 1);
            const auto adjoint = impulseResponse(adjointModel, 0, horizon, N);
            std::vector<double> X(3 * N), Z(3 * N);
            for (std::size_t k = 0; k < N; ++k) {
                const double w = std::sqrt((k == 0 || k == N - 1) ? 0.5 * dt : dt);
                for (std::size_t i = 0; i < 3; ++i) { X[i * N + k] = w * primal.states[k * 3 + i]; Z[i * N + k] = w * adjoint.states[k * 3 + i]; }
            }
            const auto pod = balancedPod(X, Z, 3, N, N, 3);
            check(near(pod.hankelSingularValues[0], hsv[0], 1e-5) && near(pod.hankelSingularValues[2], hsv[2], 1e-4), "balanced POD approximates HSVs");
        }

        // --- Structure ----------------------------------------------------------------
        const auto mixed = ss({-1, 0, 0, 0, -2, 0, 0, 0, -3}, {1, 1, 0}, {1, 0, 1}, 3, 1, 1);
        const auto kd = kalmanDecomposition(mixed);
        check(kd.controllableObservable == 1 && kd.controllableUnobservable == 1 && kd.uncontrollableObservable == 1
              && kd.uncontrollableUnobservable == 0 && kd.structureResidual < 1e-10, "Kalman decomposition dimensions");
        const auto minimal = minimalRealization(mixed);
        check(minimal.minimal.states == 1 && hasPole(systemPoles(minimal.minimal), {-1, 0}), "minimal realization order");
        check(maxResponseGap(minimal.minimal, mixed, {0.0, 0.3, 2.0, 11.0}) < 1e-12, "minimal realization preserves the transfer function");
        const auto rep = analyzeStructure(mixed);
        check(!rep.controllable && !rep.observable && rep.stabilizable && rep.detectable && rep.pbhAgrees && rep.minimalOrder == 1, "structure report");
        const auto hidden = analyzeStructure(ss({-1, 0, 0, 2}, {1, 0}, {1, 1}, 2, 1, 1));
        check(!hidden.stabilizable && hidden.detectable && hidden.pbhAgrees && hasPole(hidden.uncontrollableModes, {2, 0}), "unstabilizable hidden mode");
        {   // A conjugate pair split twice (mode and partner as separate real 2x2 blocks): one
            // real output can observe at most two of the four coordinates.
            const double a = -0.06, b = 6.28;
            const auto doubled = ss({a, -b, 0, 0, b, a, 0, 0, 0, 0, a, b, 0, 0, -b, a}, {1, 0, 1, 0},
                                    {0.3, -0.2, 0.3, 0.2}, 4, 1, 1);
            const auto dr = analyzeStructure(doubled);
            check(dr.observableDimension == 2 && dr.controllableDimension == 2 && dr.minimalOrder <= 2 && dr.pbhAgrees,
                  "split-doubled conjugate pair is reduced to two states");
        }
        check(controllabilityStaircase(mixed).reachable == 2 && observabilityStaircase(mixed).reachable == 2, "staircase dimensions");

        // --- Design -------------------------------------------------------------------
        const auto integrator = ss({0, 1, 0, 0}, {0, 1}, {1, 0}, 2, 1, 1);
        const auto Kp = placePoles(integrator, {{-1, 1}, {-1, -1}});
        check(near(Kp[0], 2, 1e-10) && near(Kp[1], 2, 1e-10), "Ackermann double integrator");
        const auto Lp = placeObserverPoles(integrator, {-4, -5});
        check(near(Lp[0], 9, 1e-9) && near(Lp[1], 20, 1e-9), "observer placement");
        const auto mimo = ss({0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, -1, 2, -3, 0.5}, {0, 0, 1, 0, 0, 1, 1, 1}, {1, 0, 0, 0}, 4, 2, 1);
        (void)placePoles(mimo, {{-1, 0}, {-2, 0}, {-3, 1}, {-3, -1}});
        (void)placePoles(mimo, {{-1, 0}, {-1, 0}, {-2, 0}, {-2, 0}});  // repeated poles, multiplicity <= inputs
        const auto scalar = makeStateSpace({1}, {1}, {1}, {}, 1, 1, 1, TimeDomain::discrete, 1.0);
        const auto dlqr = designLQR(scalar, {1}, {1});
        const double phi = 0.5 * (1 + std::sqrt(5.0));
        check(near(dlqr.riccati[0], phi, 1e-12) && near(dlqr.gain[0], 1 / phi, 1e-12), "DARE golden ratio");
        const auto clqr = designLQR(ss({0}, {1}, {1}, 1, 1, 1), {1}, {1});
        check(near(clqr.riccati[0], 1, 1e-10) && near(clqr.gain[0], 1, 1e-10), "CARE scalar");
        const auto dlqe = designLQE(scalar, {1}, {1});
        check(near(dlqe.covariance[0], phi, 1e-12) && near(dlqe.gain[0], 1 / phi, 1e-12) && near(dlqe.filterGain[0], phi / (1 + phi), 1e-12), "discrete Kalman filter");
        const std::vector<double> I2{1, 0, 0, 1};
        const auto lqg = designLQG(integrator, I2, {1}, I2, {0.1});
        bool separated = lqg.closedLoopPoles.size() == 4;
        for (auto p : lqg.regulator.closedLoopPoles) separated &= hasPole(lqg.closedLoopPoles, p, 1e-7);
        for (auto p : lqg.estimator.observerPoles) separated &= hasPole(lqg.closedLoopPoles, p, 1e-7);
        check(separated, "LQG separation principle");
        check(isAsymptoticallyStable(stateFeedbackClosedLoop(integrator, Kp)), "state-feedback closed loop");

        // --- Norms, zeros, margins ---------------------------------------------------
        const auto lowpass = ss({-1}, {1}, {1}, 1, 1, 1);
        check(near(h2Norm(lowpass), std::sqrt(0.5), 1e-12), "H2 first order");
        const auto hinf1 = hInfinityNorm(lowpass);
        check(near(hinf1.value, 1.0, 1e-8) && hinf1.peakFrequency < 1e-6, "Hinf first order");
        const double zeta = 0.1;
        const auto resonant = ss({0, 1, -1, -2 * zeta}, {0, 1}, {1, 0}, 2, 1, 1);
        const auto hinf2 = hInfinityNorm(resonant);
        check(near(hinf2.value, 1 / (2 * zeta * std::sqrt(1 - zeta * zeta)), 1e-7) && near(hinf2.peakFrequency, std::sqrt(1 - 2 * zeta * zeta), 1e-5), "Hinf resonant peak");
        check(near(h2Norm(resonant), std::sqrt(1 / (4 * zeta)), 1e-10), "H2 resonant");
        const auto dres = discretize(resonant, 0.05);
        const auto hd = hInfinityNorm(dres);
        double gridPeak = 0; std::vector<double> w; for (int k = 0; k < 20001; ++k) w.push_back(0.9 + 0.2 * k / 20000.0);
        for (double s : frequencyResponse(dres, w).singularValues) gridPeak = std::max(gridPeak, s);
        check(near(hd.value, gridPeak, 1e-6), "discrete Hinf through the bilinear map");
        check(near(h2Norm(dres), std::sqrt([&] { double s = 0; auto r = impulseResponse(dres, 0, 0, 20000); for (double y : r.outputs) s += y * y; return s; }()), 1e-8), "discrete H2 equals impulse energy");
        const auto zeroed = ss({0, 1, -3, -4}, {0, 1}, {2, 1}, 2, 1, 1);
        const auto zeros = transmissionZeros(zeroed);
        check(zeros.size() == 1 && near(zeros[0].real(), -2, 1e-10), "transmission zero of (s+2)/((s+1)(s+3))");
        check(near(dcGain(zeroed)[0], 2.0 / 3.0, 1e-12), "DC gain");
        const auto damp = modalDamping(resonant);
        check(near(damp[0].dampingRatio, zeta, 1e-12) && near(damp[0].naturalFrequency, 1, 1e-12), "damping table");
        const auto loop1 = ss({0, 1, 0, -1}, {0, 1}, {1, 0}, 2, 1, 1);  // 1/(s(s+1))
        const auto m1 = stabilityMargins(loop1);
        check(near(m1.phaseMargin, 51.827, 1e-4) && std::isinf(m1.gainMargin), "phase margin of 1/(s(s+1))");
        const auto loop3 = ss({0, 1, 0, 0, 0, 1, -1, -3, -3}, {0, 0, 2}, {1, 0, 0}, 3, 1, 1);  // 2/(s+1)^3
        const auto m3 = stabilityMargins(loop3);
        check(near(m3.gainMargin, 4.0, 1e-8) && near(m3.phaseCrossover, std::sqrt(3.0), 1e-8), "gain margin of 2/(s+1)^3");

        // --- Discretization ------------------------------------------------------------
        const double h = 0.1;
        const auto tustin = discretize(resonant, h, DiscretizationMethod::tustin);
        check(maxResponseGap(continuousFromTustin(tustin), resonant, {0.1, 1.0, 3.0}) < 1e-12, "Tustin round trip");
        const double wd = 2.0;
        const auto gd = frequencyResponse(tustin, {wd}).response[0];
        const auto gc = frequencyResponse(resonant, {2 / h * std::tan(wd * h / 2)}).response[0];
        check(std::abs(gd - gc) < 1e-12, "Tustin frequency warping");
        const auto zoh = discretize(lowpass, h);
        check(near(zoh.A[0], std::exp(-h), 1e-14) && near(zoh.B[0], 1 - std::exp(-h), 1e-14), "zero-order hold");
        const auto step = stepResponse(lowpass, 0, 5.0, 51);
        check(near(step.outputs[50], 1 - std::exp(-5.0), 1e-12), "step response");
    } catch (const std::exception& error) {
        std::cerr << "linear systems tests threw: " << error.what() << "\n";
        return 1;
    }
    if (failures) { std::cerr << failures << " linear systems checks failed\n"; return 1; }
    std::cout << "linear systems tests passed\n";
    return 0;
}

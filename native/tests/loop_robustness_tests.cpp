// Analytic checks for classical loop robustness and transfer functions:
// LQR/Kalman guarantees, Doyle's LQG counterexample, loop transfer recovery,
// classical and disk margins, pole/zero/gain and Bode data.
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/LoopMargins.hpp>
#include <dynsyslab/core/StateSpace.hpp>
#include <dynsyslab/core/SystemNorms.hpp>
#include <dynsyslab/core/TransferFunction.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>

using namespace dynsyslab::core;
using Cx = std::complex<double>;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) { if (!ok) { ++failures; std::cerr << "FAIL: " << what << "\n"; } }
bool near(double a, double b, double tol) { return std::abs(a - b) <= tol * std::max(1.0, std::abs(b)); }
StateSpaceModel ss(std::vector<double> A, std::vector<double> B, std::vector<double> C, std::size_t n,
                   std::size_t m = 1, std::size_t p = 1) {
    return makeStateSpace(std::move(A), std::move(B), std::move(C), {}, n, m, p);
}
Cx at(const StateSpaceModel& m, double w) { return frequencyResponse(m, {w}).response.front(); }
std::vector<double> outer(double s, std::vector<double> v) {
    std::vector<double> o(v.size() * v.size());
    for (std::size_t i = 0; i < v.size(); ++i) for (std::size_t j = 0; j < v.size(); ++j) o[i * v.size() + j] = s * v[i] * v[j];
    return o;
}
}  // namespace

int main() {
    try {
        const std::vector<double> I2{1, 0, 0, 1};
        // --- LQR and Kalman loops: guaranteed |1 + L| >= 1 ------------------------------------
        const auto plant = ss({0, 1, 0, 0}, {0, 1}, {1, 0}, 2);  // double integrator, position sensor
        const auto lqr = designLQR(plant, I2, {1});
        const auto rl = loopRobustness(stateFeedbackInputLoop(plant, lqr.gain));
        check(rl.closedLoopStable && rl.minimumReturnDifference >= 1 - 1e-9 && rl.minimumReturnDifference <= 1 + 1e-9,
              "LQR return difference >= 1, reaching 1 as w -> inf");
        const auto hot = loopRobustness(stateFeedbackInputLoop(plant, designLQR(plant, {1e6, 0, 0, 1e6}, {1e-4}).gain));
        check(std::abs(hot.minimumReturnDifference - 1) <= 1e-6, "high-gain LQR loop: grid reaches roll-off");
        check(rl.diskGainLower <= 0.5 + 1e-9 && std::isinf(rl.diskGainUpper) && rl.diskPhaseMargin >= 60 - 1e-6,
              "LQR guaranteed margins [1/2, inf), 60 deg");
        check(rl.phaseMargin >= 60 - 1e-6 && std::isinf(rl.gainMarginUpper), "LQR classical margins");
        const auto lqe = designLQE(plant, I2, {1});
        const auto el = loopRobustness(estimatorOutputLoop(plant, lqe.gain));
        check(el.closedLoopStable && el.minimumReturnDifference >= 1 - 1e-9, "Kalman filter loop return difference >= 1");

        // --- Loop algebra: input and output loops equal -C_K G for SISO -------------------------
        const auto lqg = designLQG(plant, I2, {1}, I2, {1});
        const auto li = outputFeedbackInputLoop(plant, lqg.compensator);
        const auto lo = outputFeedbackOutputLoop(plant, lqg.compensator);
        for (double w : {0.3, 1.0, 4.0}) {
            const Cx direct = -at(lqg.compensator, w) * at(plant, w);
            check(std::abs(at(li, w) - direct) < 1e-10 * std::max(1.0, std::abs(direct)), "LQG input loop is -C_K G");
            check(std::abs(at(lo, w) - direct) < 1e-10 * std::max(1.0, std::abs(direct)), "SISO output loop equals input loop");
        }
        check(loopRobustness(li).closedLoopStable, "LQG loop nominally stable");

        // --- Doyle (1978): LQG has no guaranteed margin ----------------------------------------
        const auto doyle = ss({1, 1, 0, 1}, {0, 1}, {1, 0}, 2);
        const double q = 1000, sigma = 1000;
        const auto dl = designLQG(doyle, outer(q, {1, 1}), {1}, outer(sigma, {1, 1}), {1});
        const auto dr = loopRobustness(stateFeedbackInputLoop(doyle, dl.regulator.gain));
        const auto dg = loopRobustness(outputFeedbackInputLoop(doyle, dl.compensator));
        check(dr.minimumReturnDifference >= 1 - 1e-9, "Doyle example: LQR keeps its guarantee");
        check(dg.closedLoopStable && dg.minimumReturnDifference < 0.1 && dg.diskGainUpper < 1.11,
              "Doyle example: LQG input loop has almost no margin");

        // --- Loop transfer recovery: minimum-phase plant, LQG -> LQR at the input --------------
        const auto mp = ss({0, 1, -1, -1}, {0, 1}, {1, 1}, 2);  // zero at -1
        const auto mpLqr = designLQR(mp, I2, {1});
        const double target = loopRobustness(stateFeedbackInputLoop(mp, mpLqr.gain)).minimumReturnDifference;
        double previousGap = 1e300;
        for (double q2 : {0.0, 1e2, 1e4, 1e6}) {
            const auto design = designLQG(mp, I2, {1}, loopTransferRecoveryCovariance(mp, I2, q2), {1});
            const double alpha = loopRobustness(outputFeedbackInputLoop(mp, design.compensator)).minimumReturnDifference;
            const double gap = std::abs(alpha - target);
            check(gap <= previousGap + 1e-6, "LTR gap shrinks as q2 grows");
            previousGap = gap;
        }
        check(previousGap < 0.02, "LTR recovers the LQR return difference");

        // --- Classical margins both ways ---------------------------------------------------------
        const auto cubic = ss({0, 1, 0, 0, 0, 1, -1, -3, -3}, {0, 0, 2}, {1, 0, 0}, 3);  // 2/(s+1)^3
        const auto cm = loopRobustness(cubic);
        check(near(cm.gainMarginUpper, 4.0, 1e-8) && cm.gainMarginLower == 0.0 && near(cm.phaseCrossover, std::sqrt(3.0), 1e-8),
              "gain margin of 2/(s+1)^3");
        double alpha = 1e300, tpeak = 0;
        for (int k = 0; k <= 200000; ++k) {
            const double w = std::pow(10.0, -3 + 6.0 * k / 200000.0);
            const Cx l = at(cubic, w);
            alpha = std::min(alpha, std::abs(1.0 + l)); tpeak = std::max(tpeak, std::abs(l / (1.0 + l)));
        }
        check(near(cm.minimumReturnDifference, alpha, 1e-3) && near(cm.sensitivityPeak, 1 / alpha, 1e-3), "sensitivity peak");
        check(near(cm.complementarySensitivityPeak, tpeak, 1e-3), "complementary sensitivity peak");
        // a conditionally stable loop: unstable for both larger and smaller gain
        const auto cond = ss({0, 1, 0, 0, 0, 1, -0.5, -3, -3}, {0, 0, 1}, {3, 1.5, 0.2}, 3);
        const auto cc = loopRobustness(cond);
        check(!cc.closedLoopStable || cc.gainMarginLower >= 0.0, "lower gain margin is reported");

        // --- Transfer function: (s+2)/((s+1)(s+3)) ---------------------------------------------
        const auto zeroed = ss({0, 1, -3, -4}, {0, 1}, {2, 1}, 2);
        const auto tf = transferFunction(zeroed, 0, 0);
        check(tf.zeros.size() == 1 && near(tf.zeros[0].real(), -2, 1e-10) && near(tf.gain, 1.0, 1e-10), "zero and gain");
        check(tf.relativeDegree == 1 && tf.cancellations.empty(), "relative degree");
        check(tf.denominator.size() == 3 && near(tf.denominator[1], 4, 1e-10) && near(tf.denominator[2], 3, 1e-10), "denominator");
        check(tf.numerator.size() == 2 && near(tf.numerator[0], 1, 1e-10) && near(tf.numerator[1], 2, 1e-10), "numerator");
        const auto unobservable = ss({-1, 0, 0, -2}, {1, 1}, {1, 0}, 2);  // mode -2 unobservable
        const auto ut = transferFunction(unobservable, 0, 0);
        check(ut.cancellations.size() == 1, "pole/zero cancellation reported");
        const auto b = bode(cubic, 0, 0, {0.01, 1.0, 1000.0});
        check(near(b.magnitudeDb[0], 20 * std::log10(2.0), 1e-3) && near(b.phaseDeg[2], -270.0, 0.5), "Bode magnitude and unwrapped phase");
        const auto two = ss({-1, 0, 1, 0, -1, 0, 0, 0, -3}, {1, 0, 0, 0, 1, 0}, {1, 0, 0, 0, 0, 1}, 3, 2, 2);
        const auto ch = channel(two, 1, 1);
        check(ch.states == 3 && ch.inputs == 1 && ch.outputs == 1, "channel extraction");
    } catch (const std::exception& error) {
        std::cerr << "loop robustness tests threw: " << error.what() << "\n";
        return 1;
    }
    if (failures) { std::cerr << failures << " loop robustness checks failed\n"; return 1; }
    std::cout << "loop robustness tests passed\n";
    return 0;
}

#include <dynsyslab/core/LoopMargins.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::core {
namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kPi = 3.14159265358979323846;

StateSpaceModel finish(std::vector<double> A, std::vector<double> B, std::vector<double> C, double D,
                       std::size_t N, const StateSpaceModel& plant) {
    auto loop = makeStateSpace(std::move(A), std::move(B), std::move(C), {D}, N, 1, 1);
    loop.domain = plant.domain; loop.sampleTime = plant.sampleTime;
    return loop;
}
double entry(const std::vector<double>& v, std::size_t i) { return v.empty() ? 0.0 : v[i]; }
}  // namespace

StateSpaceModel stateFeedbackInputLoop(const StateSpaceModel& plant, const std::vector<double>& gain) {
    if (plant.inputs != 1) throw std::invalid_argument("state-feedback input loop needs a single-input plant");
    if (gain.size() != plant.states) throw std::invalid_argument("state-feedback gain must have one entry per state");
    return finish(plant.A, plant.B, gain, 0.0, plant.states, plant);
}

StateSpaceModel estimatorOutputLoop(const StateSpaceModel& plant, const std::vector<double>& gain) {
    if (plant.outputs != 1) throw std::invalid_argument("estimator output loop needs a single-output plant");
    if (gain.size() != plant.states) throw std::invalid_argument("estimator gain must have one entry per state");
    return finish(plant.A, gain, plant.C, 0.0, plant.states, plant);
}

StateSpaceModel outputFeedbackInputLoop(const StateSpaceModel& plant, const StateSpaceModel& compensator) {
    const std::size_t n = plant.states, k = compensator.states, p = plant.outputs, N = n + k;
    if (plant.inputs != 1) throw std::invalid_argument("input loop needs a single-input plant");
    if (compensator.inputs != p || compensator.outputs != 1)
        throw std::invalid_argument("compensator must map the plant outputs to the plant input");
    // x' = A x + B v, y = C x + D v; xc' = Ac xc + Bc y; L v = -(Cc xc + Dc y)
    std::vector<double> A(N * N, 0.0), B(N, 0.0), C(N, 0.0);
    double D = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) A[i * N + j] = plant.A[i * n + j];
        B[i] = plant.B[i];
    }
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) A[(n + i) * N + n + j] = compensator.A[i * k + j];
        for (std::size_t r = 0; r < p; ++r) {
            const double bc = compensator.B[i * p + r];
            for (std::size_t j = 0; j < n; ++j) A[(n + i) * N + j] += bc * plant.C[r * n + j];
            B[n + i] += bc * entry(plant.D, r);
        }
        C[n + i] = -compensator.C[i];
    }
    for (std::size_t r = 0; r < p; ++r) {
        const double dc = entry(compensator.D, r);
        for (std::size_t j = 0; j < n; ++j) C[j] -= dc * plant.C[r * n + j];
        D -= dc * entry(plant.D, r);
    }
    return finish(std::move(A), std::move(B), std::move(C), D, N, plant);
}

StateSpaceModel outputFeedbackOutputLoop(const StateSpaceModel& plant, const StateSpaceModel& compensator) {
    const std::size_t n = plant.states, k = compensator.states, m = plant.inputs, N = n + k;
    if (plant.outputs != 1) throw std::invalid_argument("output loop needs a single-output plant");
    if (compensator.inputs != 1 || compensator.outputs != m)
        throw std::invalid_argument("compensator must map the plant output to the plant inputs");
    // xc' = Ac xc + Bc w, u = Cc xc + Dc w; x' = A x + B u; L w = -(C x + D u)
    std::vector<double> A(N * N, 0.0), B(N, 0.0), C(N, 0.0);
    double D = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) A[i * N + j] = plant.A[i * n + j];
        for (std::size_t u = 0; u < m; ++u) {
            for (std::size_t j = 0; j < k; ++j) A[i * N + n + j] += plant.B[i * m + u] * compensator.C[u * k + j];
            B[i] += plant.B[i * m + u] * entry(compensator.D, u);
        }
        C[i] = -plant.C[i];
    }
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) A[(n + i) * N + n + j] = compensator.A[i * k + j];
        B[n + i] = compensator.B[i];
    }
    for (std::size_t u = 0; u < m; ++u) {
        for (std::size_t j = 0; j < k; ++j) C[n + j] -= entry(plant.D, u) * compensator.C[u * k + j];
        D -= entry(plant.D, u) * entry(compensator.D, u);
    }
    return finish(std::move(A), std::move(B), std::move(C), D, N, plant);
}

LoopRobustness loopRobustness(const StateSpaceModel& loop, std::size_t points) {
    if (loop.inputs != 1 || loop.outputs != 1) throw std::invalid_argument("loop robustness requires a SISO loop");
    LoopRobustness out;
    out.frequencies = logarithmicFrequencyGrid(loop, std::max<std::size_t>(points, 16));
    const auto value = [&](double w) { return frequencyResponse(loop, {w}).response.front(); };
    // High-gain loops are still large where the pole-based grid ends: extend it until
    // L has settled to its direct term (continuous) so crossings and |1 + L| minima are seen.
    if (loop.domain == TimeDomain::continuous) {
        const double settled = 1e-4 * std::max(1.0, std::abs(loop.D[0]));
        for (int decade = 0; decade < 12 && std::abs(value(out.frequencies.back()) - loop.D[0]) > settled; ++decade) {
            const double start = out.frequencies.back();
            for (int k = 1; k <= 50; ++k) out.frequencies.push_back(start * std::pow(10.0, k / 50.0));
        }
    }
    out.response = frequencyResponse(loop, out.frequencies).response;
    const auto refine = [&](double a, double b, auto&& f) {  // bisection in log frequency
        const bool fa = f(a) < 0.0;
        for (int i = 0; i < 60; ++i) { const double c = std::sqrt(a * b); if ((f(c) < 0.0) == fa) a = c; else b = c; }
        return std::sqrt(a * b);
    };
    const auto magnitude = [&](double w) { return std::abs(value(w)) - 1.0; };
    const auto imaginary = [&](double w) { return value(w).imag(); };
    out.gainMarginUpper = kInf; out.gainMarginLower = 0.0; out.phaseMargin = kInf;
    out.phaseCrossover = out.gainCrossover = kInf;
    double nearestUpper = kInf;
    for (std::size_t i = 0; i + 1 < out.frequencies.size(); ++i) {
        const double a = out.frequencies[i], b = out.frequencies[i + 1];
        const auto la = out.response[i], lb = out.response[i + 1];
        if ((std::abs(la) < 1.0) != (std::abs(lb) < 1.0)) {
            const double w = refine(a, b, magnitude);
            // angular distance of L(jw) from -1 on the unit circle, in [0, 180]
            const double phase = 180.0 - std::abs(std::arg(value(w))) * 180.0 / kPi;
            if (phase < out.phaseMargin) { out.phaseMargin = phase; out.gainCrossover = w; }
        }
        if ((la.imag() < 0.0) != (lb.imag() < 0.0)) {
            const double w = refine(a, b, imaginary);
            const auto l = value(w);
            if (l.real() >= 0.0 || std::abs(l) == 0.0) continue;
            const double factor = 1.0 / std::abs(l);
            if (factor > 1.0 && factor < nearestUpper) { nearestUpper = factor; out.gainMarginUpper = factor; out.phaseCrossover = w; }
            if (factor < 1.0 && factor > out.gainMarginLower) out.gainMarginLower = factor;
        }
    }
    out.minimumReturnDifference = kInf;
    out.complementarySensitivityPeak = 0.0;
    auto visit = [&](std::complex<double> l, double w) {
        const double d = std::abs(1.0 + l);
        if (d < out.minimumReturnDifference) { out.minimumReturnDifference = d; out.returnDifferenceFrequency = w; }
        if (d > 0.0) out.complementarySensitivityPeak = std::max(out.complementarySensitivityPeak, std::abs(l) / d);
    };
    for (std::size_t i = 0; i < out.response.size(); ++i) visit(out.response[i], out.frequencies[i]);
    if (loop.domain == TimeDomain::continuous) visit({loop.D[0], 0.0}, kInf);  // the curve ends at L(inf) = D
    if (loop.domain == TimeDomain::continuous) {  // the static gain closes the Nyquist curve
        try {  // integrating loops have no finite static gain (the curve starts at infinity)
            const auto dc = dcGain(loop);
            if (!dc.empty() && std::isfinite(dc[0])) visit({dc[0], 0.0}, 0.0);
        } catch (const std::invalid_argument&) {}
    }
    const double a = out.minimumReturnDifference;
    out.sensitivityPeak = a > 0.0 ? 1.0 / a : kInf;
    out.diskGainLower = 1.0 / (1.0 + a);
    out.diskGainUpper = a >= 1.0 ? kInf : 1.0 / (1.0 - a);
    out.diskPhaseMargin = 2.0 * std::asin(std::min(a, 2.0) / 2.0) * 180.0 / kPi;
    if (std::abs(1.0 + loop.D[0]) < 1e-14) throw std::invalid_argument("loop is ill-posed (1 + D = 0)");
    auto closed = loop;  // A - B C / (1 + D)
    const double scale = 1.0 / (1.0 + loop.D[0]);
    for (std::size_t i = 0; i < loop.states; ++i)
        for (std::size_t j = 0; j < loop.states; ++j) closed.A[i * loop.states + j] -= loop.B[i] * scale * loop.C[j];
    out.closedLoopStable = isAsymptoticallyStable(closed);
    return out;
}

std::vector<double> loopTransferRecoveryCovariance(const StateSpaceModel& plant,
                                                   const std::vector<double>& W, double q2) {
    const std::size_t n = plant.states, m = plant.inputs;
    if (W.size() != n * n) throw std::invalid_argument("process covariance must be states x states");
    if (!(q2 >= 0.0) || !std::isfinite(q2)) throw std::invalid_argument("recovery weight must be finite and nonnegative");
    auto out = W;
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j)
        for (std::size_t u = 0; u < m; ++u) out[i * n + j] += q2 * plant.B[i * m + u] * plant.B[j * m + u];
    return out;
}

}  // namespace dynsyslab::core

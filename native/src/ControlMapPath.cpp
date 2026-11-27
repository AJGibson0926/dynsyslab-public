#include <dynsyslab/core/ControlMapPath.hpp>

#include <dynsyslab/numerics/Control.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

ControlMapPathResult analyzeControlMapPath(const KoopmanModel& model, const std::vector<double>& states,
                                           std::size_t samples, const ControlMapPathOptions& options) {
    const std::size_t n = model.stateDimension(), split = model.splitDimension();
    if (states.size() != samples * n || !samples) throw std::invalid_argument("control-map path does not match the K-ROM state");
    if (model.inputDimension() != 1) throw std::invalid_argument("control-map path diagnostics support one input");
    const auto A = model.generatorSplit();
    const std::vector<double> Q = [&] { std::vector<double> q(split * split, 0.0);
        for (std::size_t i = 0; i < split; ++i) q[i * split + i] = options.statePenalty; return q; }();
    ControlMapPathResult result;
    std::size_t stable = 0;
    for (std::size_t k = 0; k < samples; k += std::max<std::size_t>(1, options.stride)) {
        const double* x = states.data() + k * n;
        const auto exact = model.controlMap(x, "exact");
        const auto fitted = model.controlMap(x, options.controlMap);
        double e = 0, a = 0, b = 0; std::complex<double> inner = 0;
        for (std::size_t j = 0; j < exact.size(); ++j) {
            e += std::norm(fitted[j] - exact[j]); a += std::norm(exact[j]); b += std::norm(fitted[j]);
            inner += std::conj(fitted[j]) * exact[j];
        }
        const double na = std::sqrt(a), nb = std::sqrt(b);
        result.relativeError.push_back(na > 0 ? std::sqrt(e) / na : HUGE_VAL);
        result.cosine.push_back(na > 0 && nb > 0 ? std::abs(inner) / (na * nb) : 0.0);
        result.normRatio.push_back(na > 0 ? nb / na : HUGE_VAL);
        const auto Bmap = model.controlMapSplit(x, options.controlMap);
        const auto Bexact = model.controlMapSplit(x, "exact");
        const auto lqr = numerics::solveContinuousLqr(A, Bmap, Q, {options.inputPenalty}, split, 1);
        auto closed = A;
        for (std::size_t i = 0; i < split; ++i)
            for (std::size_t j = 0; j < split; ++j) closed[i * split + j] -= Bexact[i] * lqr.gain[j];
        const auto poles = numerics::eigensystem(closed, split, false).eigenvalues;
        double worst = -HUGE_VAL;
        for (const auto& pole : poles) worst = std::max(worst, pole.real());
        result.maximumPoleReal.push_back(worst);
        stable += worst < 0.0;
        result.sampleIndices.push_back(k);
    }
    result.stableFraction = double(stable) / double(result.sampleIndices.size());
    const double worstError = *std::max_element(result.relativeError.begin(), result.relativeError.end());
    const double worstCosine = *std::min_element(result.cosine.begin(), result.cosine.end());
    if (worstError > options.maximumRelativeError) result.failures.push_back("relative error above maximum");
    if (worstCosine < options.minimumCosine) result.failures.push_back("direction cosine below minimum");
    if (result.stableFraction < options.minimumStableFraction) result.failures.push_back("truth-coupled poles unstable");
    result.passed = result.failures.empty();
    return result;
}

}  // namespace dynsyslab::core

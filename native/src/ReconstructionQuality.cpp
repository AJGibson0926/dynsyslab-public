#include <dynsyslab/core/ReconstructionQuality.hpp>

#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Complex = std::complex<double>;

std::vector<Complex> rows(const std::vector<Complex>& matrix, std::size_t width,
                          const std::vector<std::size_t>& selected) {
    std::vector<Complex> result;
    result.reserve(selected.size() * width);
    for (auto row : selected) result.insert(result.end(), matrix.begin() + row * width, matrix.begin() + (row + 1) * width);
    return result;
}

double condition(const std::vector<Complex>& psi, std::size_t samples, std::size_t modes) {
    std::vector<Complex> gram(modes * modes, 0.0);
    for (std::size_t s = 0; s < samples; ++s)
        for (std::size_t i = 0; i < modes; ++i)
            for (std::size_t j = 0; j < modes; ++j) gram[i * modes + j] += std::conj(psi[s * modes + i]) * psi[s * modes + j];
    const auto eigen = numerics::hermitianEigensystem(gram, modes);
    const double low = std::max(eigen.eigenvalues.front(), 0.0), high = eigen.eigenvalues.back();
    return low > 0.0 ? std::sqrt(high / low) : std::numeric_limits<double>::infinity();
}

}  // namespace

ReconstructionQuality assessReconstruction(const std::vector<Complex>& psi, const std::vector<Complex>& y,
                                           std::size_t samples, std::size_t modes, std::size_t outputs, double ridge,
                                           const ReconstructionQualityPolicy& policy,
                                           const std::vector<Complex>& validationPsi,
                                           const std::vector<Complex>& validationY, std::size_t validationSamples) {
    if (psi.size() != samples * modes || y.size() != samples * outputs)
        throw std::invalid_argument("reconstruction quality arrays are not sample aligned");
    std::vector<std::size_t> train, hold;
    if (validationSamples) {
        for (std::size_t s = 0; s < samples; ++s) train.push_back(s);
    } else {
        if (!(policy.validationFraction > 0.0 && policy.validationFraction < 1.0))
            throw std::invalid_argument("validation_fraction must lie in (0, 1)");
        const auto held = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(policy.validationFraction * samples)));
        if (held >= samples) throw std::invalid_argument("too few samples for a held-out split");
        if (policy.split == "tail") {
            for (std::size_t s = 0; s < samples; ++s) (s < samples - held ? train : hold).push_back(s);
        } else if (policy.split == "interleaved") {
            const double stride = double(samples) / double(held);
            std::vector<bool> isHeld(samples, false);
            for (std::size_t k = 0; k < held; ++k) isHeld[std::min(samples - 1, static_cast<std::size_t>((k + 0.5) * stride))] = true;
            for (std::size_t s = 0; s < samples; ++s) (isHeld[s] ? hold : train).push_back(s);
        } else throw std::invalid_argument("split must be interleaved or tail");
    }
    ReconstructionQuality result;
    const auto trainPsi = rows(psi, modes, train);
    result.fit = fitObservableReconstruction(trainPsi, train.size(), modes, rows(y, outputs, train), outputs, ridge);
    result.conditionNumber = condition(trainPsi, train.size(), modes);
    const auto scorePsi = validationSamples ? validationPsi : rows(psi, modes, hold);
    const auto scoreY = validationSamples ? validationY : rows(y, outputs, hold);
    const std::size_t scored = validationSamples ? validationSamples : hold.size();
    if (scorePsi.size() != scored * modes || scoreY.size() != scored * outputs)
        throw std::invalid_argument("validation reconstruction arrays are not sample aligned");
    const auto predicted = reconstructObservable(result.fit, scorePsi, scored);
    result.trainingSamples = train.size(); result.validationSamples = scored;
    result.r2.assign(outputs, 0.0); result.normalizedRmse.assign(outputs, 0.0);
    for (std::size_t o = 0; o < outputs; ++o) {
        double mean = 0.0;
        for (std::size_t s = 0; s < scored; ++s) mean += scoreY[s * outputs + o].real() / double(scored);
        double residual = 0.0, total = 0.0;
        for (std::size_t s = 0; s < scored; ++s) {
            const double truth = scoreY[s * outputs + o].real();
            const double fit = predicted[s * outputs + o].real();  // includes the fit's offset
            residual += (truth - fit) * (truth - fit); total += (truth - mean) * (truth - mean);
        }
        result.r2[o] = total > 0.0 ? 1.0 - residual / total : (residual == 0.0 ? 1.0 : -HUGE_VAL);
        result.normalizedRmse[o] = total > 0.0 ? std::sqrt(residual / total) : (residual == 0.0 ? 0.0 : HUGE_VAL);
        if (!(result.r2[o] >= policy.minimumR2)) result.failures.push_back("output " + std::to_string(o) + " R2 below minimum");
        if (!(result.normalizedRmse[o] <= policy.maximumNormalizedRmse))
            result.failures.push_back("output " + std::to_string(o) + " normalized RMSE above maximum");
    }
    if (!(result.conditionNumber <= policy.maximumConditionNumber)) result.failures.push_back("design condition number above maximum");
    result.passed = result.failures.empty();
    return result;
}

}  // namespace dynsyslab::core

#include <dynsyslab/core/KROM.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

KROMRollout rolloutKROM(const KROMRolloutRequest& request) {
    const std::size_t modeCount = request.eigenvalues.size();
    const std::size_t sampleCount = request.time.size();
    if (!modeCount || !sampleCount || !request.stateDimension ||
        request.initialCoordinates.size() != modeCount ||
        request.modes.size() != request.stateDimension * modeCount)
        throw std::invalid_argument("KROM array dimensions are inconsistent or empty");
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        if (!std::isfinite(request.time[sample]) ||
            (sample && request.time[sample] <= request.time[sample - 1]))
            throw std::invalid_argument("KROM time samples must be finite and strictly increasing");
    }
    KROMRollout result;
    result.samples = sampleCount;
    result.modes = modeCount;
    result.stateDimension = request.stateDimension;
    result.lifted.resize(sampleCount * modeCount);
    result.reconstructed.resize(sampleCount * request.stateDimension);
    const double initialTime = request.time.front();
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        for (std::size_t mode = 0; mode < modeCount; ++mode)
            result.lifted[sample * modeCount + mode] =
                std::exp(request.eigenvalues[mode] * (request.time[sample] - initialTime)) *
                request.initialCoordinates[mode];
        for (std::size_t row = 0; row < request.stateDimension; ++row) {
            std::complex<double> value{};
            for (std::size_t mode = 0; mode < modeCount; ++mode)
                value += request.modes[row * modeCount + mode] * result.lifted[sample * modeCount + mode];
            result.reconstructed[sample * request.stateDimension + row] = value;
        }
    }
    return result;
}

}  // namespace dynsyslab::core

#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double windowValue(FourierWindow window, std::size_t index, std::size_t count) {
    if (window == FourierWindow::rectangular || count == 1) return 1.0;
    const double angle = 2.0 * kPi * static_cast<double>(index)
        / static_cast<double>(count - 1);
    if (window == FourierWindow::hann) return 0.5 - 0.5 * std::cos(angle);
    if (window == FourierWindow::hamming) return 0.54 - 0.46 * std::cos(angle);
    return 0.42 - 0.5 * std::cos(angle) + 0.08 * std::cos(2.0 * angle);
}

std::shared_ptr<const FourierTransform> computeTransform(
    const TimeSeries& source, const FourierOptions& options) {
    if (source.sampleCount() < 2 || !source.uniformTimeSteps())
        throw std::invalid_argument("Fourier analysis requires at least two uniformly spaced samples");
    const std::size_t samples = source.sampleCount();
    const std::size_t components = source.componentCount();
    const std::size_t bins = options.oneSided ? samples / 2 + 1 : samples;
    auto result = std::make_shared<FourierTransform>();
    result->binCount = bins; result->componentCount = components;
    result->sourceSampleCount = samples; result->deltaTime = source.deltaTime();
    result->oneSided = options.oneSided; result->meanRemoved = options.removeMean;
    result->window = options.window; result->frequencies.resize(bins);
    result->coefficients.resize(bins * components);
    result->amplitudes.resize(bins * components);
    for (std::size_t bin = 0; bin < bins; ++bin) {
        const auto signedBin = !options.oneSided && bin > samples / 2
            ? static_cast<long long>(bin) - static_cast<long long>(samples)
            : static_cast<long long>(bin);
        result->frequencies[bin] = static_cast<double>(signedBin)
            / (static_cast<double>(samples) * result->deltaTime);
    }
    std::vector<double> weights(samples); double weightSum = 0.0;
    for (std::size_t sample = 0; sample < samples; ++sample) {
        weights[sample] = windowValue(options.window, sample, samples);
        weightSum += weights[sample];
    }
    const double coherentGain = weightSum / static_cast<double>(samples);
    if (!(coherentGain > std::numeric_limits<double>::epsilon()))
        throw std::invalid_argument("Fourier window has zero coherent gain");
    for (std::size_t component = 0; component < components; ++component) {
        double mean = 0.0;
        if (options.removeMean) {
            for (std::size_t sample = 0; sample < samples; ++sample) {
                const double value = source.value(sample, component);
                if (!std::isfinite(value))
                    throw std::invalid_argument("Fourier analysis requires finite samples");
                mean += value;
            }
            mean /= static_cast<double>(samples);
        }
        std::vector<std::complex<double>> input(samples);
        for (std::size_t sample = 0; sample < samples; ++sample) {
            const double value = source.value(sample, component);
            if (!std::isfinite(value))
                throw std::invalid_argument("Fourier analysis requires finite samples");
            input[sample] = (value - mean) * weights[sample];
        }
        const auto transformed = discreteFourierTransform(input);
        for (std::size_t bin = 0; bin < bins; ++bin) {
            const std::size_t target = bin * components + component;
            result->coefficients[target] = transformed[bin]
                / (static_cast<double>(samples) * coherentGain);
            double amplitude = std::abs(result->coefficients[target]);
            const bool nyquist = samples % 2 == 0 && bin == samples / 2;
            if (options.oneSided && bin != 0 && !nyquist) amplitude *= 2.0;
            result->amplitudes[target] = amplitude;
        }
    }
    return result;
}

}  // namespace

struct TimeSeries::FourierCache {
    FourierOptions options;
    std::shared_ptr<const FourierTransform> value;
};

std::complex<double> FourierTransform::coefficient(
    std::size_t bin, std::size_t component) const {
    if (bin >= binCount || component >= componentCount)
        throw std::out_of_range("Fourier coefficient index");
    return coefficients[bin * componentCount + component];
}

double FourierTransform::amplitude(std::size_t bin, std::size_t component) const {
    if (bin >= binCount || component >= componentCount)
        throw std::out_of_range("Fourier amplitude index");
    return amplitudes[bin * componentCount + component];
}

double FourierTransform::power(std::size_t bin, std::size_t component) const {
    const double value = amplitude(bin, component); return value * value;
}

double FourierTransform::phase(std::size_t bin, std::size_t component) const {
    return std::arg(coefficient(bin, component));
}

std::size_t FourierTransform::dominantBin(
    std::size_t component, bool includeZero) const {
    if (component >= componentCount || !binCount)
        throw std::out_of_range("Fourier component index");
    std::size_t best = includeZero ? 0 : std::min<std::size_t>(1, binCount - 1);
    for (std::size_t bin = best + 1; bin < binCount; ++bin)
        if (amplitude(bin, component) > amplitude(best, component)) best = bin;
    return best;
}

double FourierTransform::dominantFrequency(
    std::size_t component, bool includeZero) const {
    return frequencies[dominantBin(component, includeZero)];
}

std::shared_ptr<const FourierTransform> FourierSidecar::transform(
    FourierOptions options) const {
    if (!source_) throw std::logic_error("Fourier sidecar has no TimeSeries");
    if (source_->fourierCache_ && source_->fourierCache_->options == options)
        return source_->fourierCache_->value;
    auto value = computeTransform(*source_, options);
    auto cache = std::make_shared<TimeSeries::FourierCache>();
    cache->options = options; cache->value = value;
    source_->fourierCache_ = std::move(cache);
    return value;
}

bool FourierSidecar::cached(const FourierOptions& options) const noexcept {
    return source_ && source_->fourierCache_ && source_->fourierCache_->options == options;
}

void FourierSidecar::clear() const noexcept {
    if (source_) source_->fourierCache_.reset();
}

}  // namespace dynsyslab::core

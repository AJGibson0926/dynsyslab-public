#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

namespace dynsyslab::core {

enum class FourierWindow { rectangular, hann, hamming, blackman };

struct FourierOptions {
    bool oneSided{true};
    bool removeMean{false};
    FourierWindow window{FourierWindow::rectangular};

    [[nodiscard]] bool operator==(const FourierOptions& other) const noexcept;
    [[nodiscard]] bool operator!=(const FourierOptions& other) const noexcept {
        return !(*this == other);
    }
};

struct FourierTransform {
    std::vector<double> frequencies;
    std::vector<std::complex<double>> coefficients;
    std::vector<double> amplitudes;
    std::size_t binCount{};
    std::size_t componentCount{};
    std::size_t sourceSampleCount{};
    double deltaTime{};
    bool oneSided{};
    bool meanRemoved{};
    FourierWindow window{FourierWindow::rectangular};

    [[nodiscard]] std::complex<double> coefficient(
        std::size_t bin, std::size_t component) const;
    [[nodiscard]] double amplitude(std::size_t bin, std::size_t component) const;
    [[nodiscard]] double power(std::size_t bin, std::size_t component) const;
    [[nodiscard]] double phase(std::size_t bin, std::size_t component) const;
    [[nodiscard]] std::size_t dominantBin(
        std::size_t component, bool includeZero = false) const;
    [[nodiscard]] double dominantFrequency(
        std::size_t component, bool includeZero = false) const;
};

class FourierSidecar;

// Immutable numerical samples. Owning instances keep their buffers alive;
// view instances add no sample-sized allocation and require their source
// buffers to outlive the TimeSeries.
class TimeSeries {
public:
    TimeSeries(std::vector<double> time, std::vector<double> values,
               std::size_t componentCount);

    [[nodiscard]] static TimeSeries view(
        const std::vector<double>& time, const double* values,
        std::size_t sampleCount, std::size_t componentCount);

    [[nodiscard]] std::size_t sampleCount() const noexcept { return sampleCount_; }
    [[nodiscard]] std::size_t componentCount() const noexcept { return componentCount_; }
    [[nodiscard]] bool uniformTimeSteps() const noexcept { return uniformTimeSteps_; }
    [[nodiscard]] double deltaTime() const;
    [[nodiscard]] const double* timeData() const noexcept { return timeData_; }
    [[nodiscard]] const double* valuesData() const noexcept { return valuesData_; }
    [[nodiscard]] double time(std::size_t sample) const;
    [[nodiscard]] double value(std::size_t sample, std::size_t component) const;
    [[nodiscard]] FourierSidecar fourier() const;

private:
    struct FourierCache;
    TimeSeries() = default;
    void initialize(const double* time, const double* values,
                    std::size_t samples, std::size_t components);

    std::shared_ptr<const std::vector<double>> ownedTime_;
    std::shared_ptr<const std::vector<double>> ownedValues_;
    const double* timeData_{};
    const double* valuesData_{};
    std::size_t sampleCount_{};
    std::size_t componentCount_{};
    bool uniformTimeSteps_{};
    double deltaTime_{};
    mutable std::shared_ptr<const FourierCache> fourierCache_;

    friend class FourierSidecar;
};

// Lightweight facade over a TimeSeries. The latest requested transform is
// cached in process memory; requesting different options replaces that cache.
class FourierSidecar {
public:
    [[nodiscard]] std::shared_ptr<const FourierTransform> transform(
        FourierOptions options = {}) const;
    [[nodiscard]] bool cached(const FourierOptions& options = {}) const noexcept;
    void clear() const noexcept;

private:
    explicit FourierSidecar(const TimeSeries& source) : source_(&source) {}
    const TimeSeries* source_{};
    friend class TimeSeries;
};

}  // namespace dynsyslab::core

#pragma once
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <stdexcept>
#include <vector>
namespace dynsyslab::core::detail {
using Vector=Eigen::VectorXd;
inline std::size_t nextPowerOfTwo(std::size_t value) {
    std::size_t result = 1;
    while (result < value) result <<= 1U;
    return result;
}

struct CorrelationPlan {
    std::size_t fftSize{};
    std::size_t kernelSize{};
    std::size_t outputSize{};
    std::vector<std::complex<double>> signalSpectrum;

    CorrelationPlan(const std::vector<double>& signal, std::size_t kernel)
        : fftSize(nextPowerOfTwo(signal.size() + kernel - 1)),
          kernelSize(kernel), outputSize(signal.size() - kernel + 1) {
        std::vector<std::complex<double>> padded(fftSize);
        for (std::size_t index = 0; index < signal.size(); ++index)
            padded[index] = signal[index];
        Eigen::FFT<double> fft;
        fft.fwd(signalSpectrum, padded);
    }

    [[nodiscard]] Vector correlate(const Vector& kernel) const {
        if (static_cast<std::size_t>(kernel.size()) != kernelSize)
            throw std::invalid_argument("correlation kernel size changed");
        std::vector<std::complex<double>> padded(fftSize);
        for (std::size_t index = 0; index < kernelSize; ++index)
            padded[index] = kernel[static_cast<Eigen::Index>(kernelSize - 1 - index)];
        Eigen::FFT<double> fft;
        std::vector<std::complex<double>> spectrum, product(fftSize), convolution;
        fft.fwd(spectrum, padded);
        for (std::size_t index = 0; index < fftSize; ++index)
            product[index] = signalSpectrum[index] * spectrum[index];
        fft.inv(convolution, product);
        Vector result(static_cast<Eigen::Index>(outputSize));
        const std::size_t offset = kernelSize - 1;
        for (std::size_t index = 0; index < outputSize; ++index)
            result[static_cast<Eigen::Index>(index)] = convolution[offset + index].real();
        return result;
    }
};

} // namespace dynsyslab::core::detail

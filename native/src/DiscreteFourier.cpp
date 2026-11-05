#include <algorithm>
#include <cmath>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
bool powerOfTwo(std::size_t value) { return value && !(value & (value - 1)); }

void radixTwo(std::vector<std::complex<double>> &values, bool inverse) {
    const std::size_t count = values.size();
    for (std::size_t index = 1, reversed = 0; index < count; ++index) {
        std::size_t bit = count >> 1;
        for (; reversed & bit; bit >>= 1)
            reversed ^= bit;
        reversed ^= bit;
        if (index < reversed)
            std::swap(values[index], values[reversed]);
    }
    for (std::size_t length = 2; length <= count; length <<= 1) {
        const double angle = (inverse ? 2.0 : -2.0) * kPi / static_cast<double>(length);
        const std::complex<double> step = std::polar(1.0, angle);
        for (std::size_t offset = 0; offset < count; offset += length) {
            std::complex<double> weight{1.0, 0.0};
            for (std::size_t index = 0; index < length / 2; ++index) {
                const auto even = values[offset + index];
                const auto odd = values[offset + index + length / 2] * weight;
                values[offset + index] = even + odd;
                values[offset + index + length / 2] = even - odd;
                weight *= step;
            }
        }
    }
    if (inverse)
        for (auto &value : values)
            value /= static_cast<double>(count);
}

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &input) {
    if (powerOfTwo(input.size())) {
        auto result = input;
        radixTwo(result, false);
        return result;
    }
    if (input.size() > (std::numeric_limits<std::size_t>::max)() / 4)
        throw std::length_error("Fourier transform length is not addressable");
    std::size_t padded = 1;
    while (padded < input.size() * 2 - 1)
        padded <<= 1;
    std::vector<std::complex<double>> signal(padded), chirp(padded);
    for (std::size_t index = 0; index < input.size(); ++index) {
        const long double square = static_cast<long double>(index) * index;
        const double angle =
            kPi *
            static_cast<double>(std::fmod(square, static_cast<long double>(2 * input.size()))) /
            static_cast<double>(input.size());
        signal[index] = input[index] * std::polar(1.0, -angle);
        chirp[index] = std::polar(1.0, angle);
        if (index)
            chirp[padded - index] = chirp[index];
    }
    radixTwo(signal, false);
    radixTwo(chirp, false);
    for (std::size_t index = 0; index < padded; ++index)
        signal[index] *= chirp[index];
    radixTwo(signal, true);
    signal.resize(input.size());
    for (std::size_t index = 0; index < input.size(); ++index) {
        const long double square = static_cast<long double>(index) * index;
        const double angle =
            kPi *
            static_cast<double>(std::fmod(square, static_cast<long double>(2 * input.size()))) /
            static_cast<double>(input.size());
        signal[index] *= std::polar(1.0, -angle);
    }
    return signal;
}

} // namespace
std::vector<std::complex<double>>
discreteFourierTransform(const std::vector<std::complex<double>> &input, bool inverse) {
    if (input.empty())
        return {};
    for (const auto &x : input)
        if (!std::isfinite(x.real()) || !std::isfinite(x.imag()))
            throw std::invalid_argument("Fourier transform requires finite values");
    if (!inverse)
        return fft(input);
    auto result = input;
    if (powerOfTwo(result.size())) {
        radixTwo(result, true);
        return result;
    }
    for (auto &x : result)
        x = std::conj(x);
    result = fft(result);
    for (auto &x : result)
        x = std::conj(x) / static_cast<double>(result.size());
    return result;
}
} // namespace dynsyslab::core

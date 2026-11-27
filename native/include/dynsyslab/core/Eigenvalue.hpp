#pragma once
// One Koopman eigenvalue: along solutions its eigenfunction evolves as psi(x(t)) = exp(lambda t) psi(x(0)).
#include <complex>
#include <cstddef>
#include <limits>
#include <optional>

namespace dynsyslab::core {

struct Eigenvalue {
    std::size_t index{};
    std::complex<double> discrete;    // multiplier over one sample interval (exp(lambda dt))
    std::complex<double> continuous;  // the rate lambda (stored branch)
    double sampleInterval{};          // 0 for an eigenvalue given in continuous time
    std::optional<std::size_t> conjugate;  // index of the conjugate partner, if stored

    [[nodiscard]] static Eigenvalue fromRate(std::complex<double> rate, double sampleInterval = 0) {
        return {0, sampleInterval > 0 ? std::exp(rate * sampleInterval) : std::complex<double>(), rate, sampleInterval, {}};
    }
    // exp(lambda t): the factor an eigenfunction value gains over `elapsed` time.
    [[nodiscard]] std::complex<double> transition(double elapsed) const { return std::exp(continuous * elapsed); }
    // Modal summaries of lambda = sigma + i omega: the growth (positive) or decay (negative) rate
    // sigma, the signed angular frequency omega, the cyclic frequency omega / 2 pi, and the quality
    // factor Q = |omega| / (2 |sigma|) (0 for a real mode, infinite for an undamped oscillation).
    [[nodiscard]] double growthRate() const noexcept { return continuous.real(); }
    [[nodiscard]] double angularFrequency() const noexcept { return continuous.imag(); }
    [[nodiscard]] double cyclicFrequency() const noexcept { return continuous.imag() / (2.0 * 3.14159265358979323846); }
    [[nodiscard]] double qualityFactor() const noexcept {
        const double omega = std::abs(continuous.imag()), sigma = std::abs(continuous.real());
        if (omega == 0.0) return 0.0;
        return sigma == 0.0 ? std::numeric_limits<double>::infinity() : omega / (2.0 * sigma);
    }
};

}  // namespace dynsyslab::core

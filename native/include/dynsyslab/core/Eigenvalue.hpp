#pragma once
// One Koopman eigenvalue: along solutions its eigenfunction evolves as psi(x(t)) = exp(lambda t) psi(x(0)).
#include <complex>
#include <cstddef>
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
};

}  // namespace dynsyslab::core

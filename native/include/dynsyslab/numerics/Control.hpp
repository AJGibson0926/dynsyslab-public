#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

struct ContinuousLqrResult {
    std::vector<double> riccati; // state x state, row major
    std::vector<double> gain;    // input x state, row major
    std::vector<std::complex<double>> closedLoopPoles;
    double relativeResidual{};
};

struct FiniteHorizonLqrResult {
    std::vector<double> gain; // input x state, row major
    std::vector<double> cost; // state x state, row major
};

// Runtime-sized, row-major matrix values cross this boundary. The compiled
// backend owns the template-heavy factorization and eigensolver details.
[[nodiscard]] ContinuousLqrResult solveContinuousLqr(
    const std::vector<double>& generator, const std::vector<double>& input,
    const std::vector<double>& statePenalty,
    const std::vector<double>& inputPenalty,
    std::size_t states, std::size_t inputs);

[[nodiscard]] FiniteHorizonLqrResult solveFiniteHorizonLqr(
    const std::vector<double>& generator, const std::vector<double>& input,
    const std::vector<double>& statePenalty,
    const std::vector<double>& inputPenalty,
    std::size_t states, std::size_t inputs, std::size_t horizon,
    double sampleTime);

}  // namespace dynsyslab::numerics

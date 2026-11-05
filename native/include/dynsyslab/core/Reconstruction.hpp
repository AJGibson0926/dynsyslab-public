#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

class KoopmanModel;
class Observable;

struct ReconstructionResult {
    std::vector<std::complex<double>> matrix; // output x mode (C)
    std::size_t outputs{};
    std::size_t modes{};
    double residualRms{};
    double residualMaximum{};
    double relativeResidual{};
};

// Solve G ~= Psi C^T over samples. Both matrices are row major.
[[nodiscard]] ReconstructionResult fitObservableReconstruction(
    const std::vector<std::complex<double>>& eigenfunctionValues,
    std::size_t samples, std::size_t modes,
    const std::vector<std::complex<double>>& observableValues,
    std::size_t outputs, double ridge = 0.0);

[[nodiscard]] ReconstructionResult fitObservableReconstruction(
    const KoopmanModel& model, const std::vector<double>& dynamicalStates,
    std::size_t samples, const Observable& observable, double ridge = 0.0);

[[nodiscard]] std::vector<std::complex<double>> reconstructObservable(
    const ReconstructionResult& reconstruction,
    const std::vector<std::complex<double>>& eigenfunctionValues,
    std::size_t samples);

}  // namespace dynsyslab::core

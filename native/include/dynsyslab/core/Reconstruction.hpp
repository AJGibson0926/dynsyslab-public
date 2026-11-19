#pragma once

#include <complex>
#include <cstddef>
#include <vector>
#include <string>

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
    // Optional real intercept (one per output): h = offset + C Psi. Empty means zero.
    // An intercept lets a KROM drop a nearly-constant mode (e.g. the equilibrium
    // offset) that a measurement such as the echo cannot observe.
    std::vector<double> offset;
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

// Fit aligned measured outputs against supplied Koopman coordinates and attach
// the named map to this KROM. No equations of motion or state evaluator needed.
[[nodiscard]] ReconstructionResult fitAndSetObservableReconstruction(
    KoopmanModel& model, const std::string& name,
    const std::vector<std::complex<double>>& eigenfunctionValues,
    std::size_t samples, const std::vector<std::complex<double>>& observableValues,
    std::size_t outputs, double ridge = 0.0);

[[nodiscard]] ReconstructionResult fitAndSetObservableReconstruction(
    KoopmanModel& model, const std::string& name,
    const std::vector<double>& dynamicalStates, std::size_t samples,
    const Observable& observable, double ridge = 0.0);

[[nodiscard]] std::vector<std::complex<double>> reconstructObservable(
    const ReconstructionResult& reconstruction,
    const std::vector<std::complex<double>>& eigenfunctionValues,
    std::size_t samples);

}  // namespace dynsyslab::core

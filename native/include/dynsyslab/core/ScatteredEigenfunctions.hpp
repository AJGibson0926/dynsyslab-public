#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ScatteredEigenfunctionData {
    std::vector<std::complex<double>> continuousEigenvalues;
    std::vector<double> trainingStates;                 // sample x state
    std::vector<std::complex<double>> trainingValues;  // sample x mode
    std::size_t samples{};
    std::size_t stateDimension{};
    std::size_t neighbors{8};
    double finiteDifferenceRelativeStep{1.0e-6};
    std::string provenance;
};

// Immutable multi-output inverse-distance interpolant. One spatial search is
// shared by every eigenfunction, which keeps large H-DMD KROM evaluation
// substantially cheaper than one interpolator/tree per mode.
class ScatteredEigenfunctionSet {
public:
    explicit ScatteredEigenfunctionSet(ScatteredEigenfunctionData data);

    [[nodiscard]] const ScatteredEigenfunctionData& data() const noexcept;
    [[nodiscard]] std::size_t modes() const noexcept;
    [[nodiscard]] std::vector<std::complex<double>> evaluate(
        const double* state) const;
    // Row-major mode x state gradient.
    [[nodiscard]] std::vector<std::complex<double>> gradient(
        const double* state) const;
    [[nodiscard]] std::shared_ptr<const ScatteredEigenfunctionSet> extract(
        const std::vector<std::size_t>& modes) const;

private:
    struct Index;
    ScatteredEigenfunctionData data_;
    std::shared_ptr<const Index> index_;
};

}  // namespace dynsyslab::core

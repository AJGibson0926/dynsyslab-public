#pragma once

#include <cstddef>
#include <vector>

namespace dynsyslab::numerics {

struct SymmetricEigensystem {
    std::vector<double> eigenvalues;  // ascending
    std::vector<double> eigenvectors; // row-major, eigenvectors in columns
    std::size_t dimension{};
};

struct ZeroOrderHold {
    std::vector<double> transition; // states x states
    std::vector<double> response;   // states x inputs
    std::size_t states{};
    std::size_t inputs{};
};

[[nodiscard]] SymmetricEigensystem symmetricEigensystem(
    const std::vector<double>& matrix,std::size_t dimension);

[[nodiscard]] std::vector<double> matrixExponential(
    const std::vector<double>& matrix,std::size_t dimension,double scale=1.0);

[[nodiscard]] ZeroOrderHold zeroOrderHold(
    const std::vector<double>& generator,const std::vector<double>& input,
    std::size_t states,std::size_t inputs,double sampleTime);

[[nodiscard]] bool isDetectable(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs);

[[nodiscard]] std::vector<double> finiteObservabilityGramian(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs,double horizon);

[[nodiscard]] std::vector<double> infiniteObservabilityGramian(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs);

}  // namespace dynsyslab::numerics

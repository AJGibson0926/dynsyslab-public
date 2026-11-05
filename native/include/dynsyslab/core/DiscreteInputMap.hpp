#pragma once
#include <complex>
#include <cstddef>
#include <string>
#include <vector>
namespace dynsyslab::core {
// psi[k+1] = diag(mu) psi[k] + matrix u[k]. NOT a continuous B_psi.
struct DiscreteInputMap {
    std::vector<std::complex<double>> matrix; // eigenfunction x input
    std::size_t eigenfunctionCount{}, inputDimension{};
    double sampleInterval{};
    std::string provenance;
    std::string modelIdentity;
};
} // namespace dynsyslab::core

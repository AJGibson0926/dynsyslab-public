#pragma once
// How an Observable's (or an Eigenfunction's) values are represented on its domain q. Two kinds:
//   "expansion"      FunctionExpansion: g_r(q) = sum_k A_rk phi_k(q) in the universal dictionary
//                    (FunctionExpansion.hpp, foundation library). Formula observables are expansions
//                    with one expression term per output.
//   "interpolation"  FunctionInterpolation: interpolated from sampled values with a selectable
//                    interpolant (FunctionInterpolation.hpp, eigenfunctions library).
// A representation may hold several outputs that share one evaluation (one dictionary, one
// neighbour search); an Observable or Eigenfunction points at it and, if needed, at one output.
// Representations know nothing about eigenvalues.
#include <dynsyslab/core/FunctionDomain.hpp>

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionRepresentation {
public:
    virtual ~FunctionRepresentation() = default;
    [[nodiscard]] virtual std::string kind() const = 0;
    [[nodiscard]] virtual const FunctionDomain& coordinates() const noexcept = 0;
    [[nodiscard]] virtual std::size_t outputCount() const noexcept = 0;
    [[nodiscard]] virtual std::vector<std::complex<double>> evaluate(const double* q) const = 0;
    // outputs x dimension, row-major.
    [[nodiscard]] virtual std::vector<std::complex<double>> jacobian(const double* q) const = 0;
    // True when jacobian() is exact (closed form), false when it is a numerical derivative.
    [[nodiscard]] virtual bool exactJacobian() const noexcept = 0;
    [[nodiscard]] virtual std::shared_ptr<const FunctionRepresentation> selectOutputs(
        const std::vector<std::size_t>& outputs) const = 0;
    // The same function on another description of its domain (e.g. with its catalog map bound).
    [[nodiscard]] virtual std::shared_ptr<const FunctionRepresentation> onDomain(FunctionDomain coordinates) const = 0;
    [[nodiscard]] std::size_t coordinateDimension() const noexcept { return coordinates().dimension; }
};

}  // namespace dynsyslab::core

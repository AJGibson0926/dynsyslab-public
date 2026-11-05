#pragma once
#include <dynsyslab/core/EigenfunctionCoordinates.hpp>
#include <dynsyslab/core/EigenfunctionSpectrum.hpp>
#include <memory>

namespace dynsyslab::core {
struct PolynomialDictionary {
    std::size_t coordinateDimension{};
    std::vector<std::size_t> exponents; // feature x coordinate
    std::vector<double> center, scale;
    std::size_t featureCount() const noexcept;
    void validate() const;
};
class PolynomialEigenfunctionView;

// Immutable, algorithm-independent polynomial bank. No training counters,
// fit residuals, physical ODE, live history, or validation cache is required.
class PolynomialEigenfunctionSet {
  public:
    PolynomialEigenfunctionSet(PolynomialDictionary dictionary,
                               std::vector<std::complex<double>> coefficients,
                               EigenfunctionSpectrum spectrum,
                               EigenfunctionCoordinates coordinates);
    const PolynomialDictionary &dictionary() const noexcept;
    const std::vector<std::complex<double>> &coefficients() const noexcept;
    const EigenfunctionSpectrum &spectrum() const noexcept;
    const EigenfunctionCoordinates &coordinates() const noexcept;
    std::size_t coordinateDimension() const noexcept;
    std::size_t eigenfunctionCount() const noexcept;
    std::complex<double> evaluate(std::size_t mode, const double *coordinates) const;
    std::vector<std::complex<double>> evaluate(const double *coordinates) const;
    std::vector<std::complex<double>> gradient(std::size_t mode, const double *coordinates) const;
    PolynomialEigenfunctionSet extract(const std::vector<std::size_t> &modes) const;
    PolynomialEigenfunctionSet withCoordinates(EigenfunctionCoordinates coordinates) const;
    PolynomialEigenfunctionView eigenfunction(std::size_t mode) const;

  private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
// Shares the bank's immutable storage; remains valid after the original handle
// is destroyed. It does not duplicate the dictionary or coefficient matrix.
class PolynomialEigenfunctionView {
  public:
    std::size_t index() const noexcept { return index_; }
    std::complex<double> evaluate(const double *q) const { return bank_.evaluate(index_, q); }
    std::vector<std::complex<double>> gradient(const double *q) const {
        return bank_.gradient(index_, q);
    }
    std::optional<std::complex<double>> continuousRate() const {
        return bank_.spectrum().continuousRate(index_);
    }

  private:
    friend class PolynomialEigenfunctionSet;
    PolynomialEigenfunctionView(PolynomialEigenfunctionSet bank, std::size_t index)
        : bank_(std::move(bank)), index_(index) {}
    PolynomialEigenfunctionSet bank_;
    std::size_t index_;
};
} // namespace dynsyslab::core

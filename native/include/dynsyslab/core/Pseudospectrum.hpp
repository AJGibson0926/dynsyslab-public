#pragma once
// A sampled pseudospectrum: points z of the complex plane with the smallest residual
// min_g ||(K - z) g|| / ||g|| over the method's trial space at each. The epsilon-pseudospectrum
// sigma_eps = {z : residual(z) < eps} is a level set; sigma is contained in every sigma_eps and is
// their intersection over eps > 0 (Trefethen & Embree, Spectra and Pseudospectra, Princeton 2005).
// It is computed evidence about the spectrum, not part of it: ResDMD evaluates it from snapshot data
// as a finite quadrature, with the same definition and scope as its eigenpair residuals (Colbrook &
// Townsend, Comm. Pure Appl. Math. 77, 2024). For a sampled operator z is a multiplier (the plane of
// mu = exp(lambda dt)).
#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace dynsyslab::core {

class FunctionExpansion;

struct Pseudospectrum {
    std::string definition, scope;  // as EigenpairResidual (Eigenpair.hpp)
    std::vector<std::complex<double>> grid;               // the points z
    std::vector<double> minimum;                          // the smallest residual at each point
    // The minimizers g at each z (approximate eigenfunctions), one output per grid point, on the
    // operator's frozen dictionary; optional.
    std::shared_ptr<const FunctionExpansion> minimizers;

    [[nodiscard]] bool empty() const noexcept { return grid.empty(); }
    // The grid points in the sampled epsilon-pseudospectrum.
    [[nodiscard]] std::vector<std::size_t> levelSet(double epsilon) const {
        std::vector<std::size_t> inside;
        for (std::size_t k = 0; k < minimum.size(); ++k)
            if (minimum[k] < epsilon) inside.push_back(k);
        return inside;
    }
    void validate() const {
        if (minimum.size() != grid.size()) throw std::invalid_argument("pseudospectrum grid and minima do not match");
    }
};

}  // namespace dynsyslab::core

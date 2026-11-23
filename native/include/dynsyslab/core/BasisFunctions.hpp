#pragma once
#include <cstddef>
#include <string>
namespace dynsyslab::core::basis {
// Classical three-term recurrences: NIST DLMF chapters 18 and 10.
double chebyshev(std::size_t degree,double x);
double laguerre(std::size_t degree,double x);
double hermiteFunction(std::size_t degree,double x);
// Normalized Matern: 2^(1-nu)/Gamma(nu) z^nu K_nu(z), z=distance/bandwidth.
// With nu=2, amplitude=2, bandwidth=1/6 reproduces the SpecRKHS Duffing kernel.
double radialKernel(const std::string& kind,double distance,double bandwidth,
                    double smoothness=2,double amplitude=1);
bool isRadial(const std::string& kind);
bool isProduct(const std::string& kind);
} // namespace dynsyslab::core::basis

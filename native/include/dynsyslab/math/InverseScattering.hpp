#pragma once
#include <array>
#include <complex>
#include <vector>
namespace dynsyslab::math {
// Spatial operator L=-d_x^2+V. These are Lax/Schroedinger spectral data, not
// eigenfunctions of the Koopman operator on the space of potentials.
struct ScatteringCoefficient {
    std::complex<double> reflection,transmission;
    double fluxDefect{};
};
// Potential is constant on each cell, zero outside. RIGHT incidence:
// psi~exp(-ikx)+R exp(ikx) at +infinity; psi~T exp(-ikx) at -infinity.
// Exact for the supplied piecewise-constant potential; discretization/tails
// of an underlying smooth potential are the caller's responsibility.
ScatteringCoefficient rightScattering(const std::vector<double>&cellPotential,
    double left,double cellWidth,double wavenumber);
struct HillMonodromy {std::array<double,4> matrix;double discriminant{},determinantDefect{};};
// Periodic extension of the same cell potential; discriminant is FULL trace.
HillMonodromy hillMonodromy(const std::vector<double>&cellPotential,double cellWidth,double energy);
// Finite-box Dirichlet centered-difference eigenvalues below zero, using Sturm
// bisection. Potential values are at INTERNAL nodes, spacing length/(N+1).
// They approximate whole-line bound states only after box/grid convergence.
std::vector<double> dirichletNegativeSpectrum(const std::vector<double>&internalPotential,double length);
struct ReflectionlessData {
    std::vector<double> kappa,logNormingSquared;
};
// Convention: V_t-6 V V_x+V_xxx=0. Right norming coefficient c for an L2-unit
// bound state behaves psi~c exp(-kappa*x) as x->+infinity; stored log(c^2).
ReflectionlessData evolveKdV(const ReflectionlessData&,double elapsed);
std::complex<double> evolveKdVReflection(std::complex<double>,double wavenumber,double elapsed);
// Marchenko finite-rank determinant reconstruction, reflection identically zero.
double reflectionlessPotential(const ReflectionlessData&,double position);
}

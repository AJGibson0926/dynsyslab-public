#pragma once
// POD-Galerkin reduced-order model (Holmes, Lumley & Berkooz 1996, "Turbulence, Coherent
// Structures, Dynamical Systems and Symmetry", Cambridge, Ch. 3; Sirovich 1987, "Turbulence and
// the dynamics of coherent structures", Q. Appl. Math. 45:561-571). Snapshots x_k minus a
// reference x_r give orthonormal POD modes Phi (n x r) from the thin SVD; the state is
// x ~= x_r + Phi z. The Galerkin reduced dynamics are z' = Phi^T f(x_r + Phi z) + Phi^T B u, and
// the linear reduced operators at the reference A_r = Phi^T J(x_r) Phi, B_r = Phi^T B(x_r).
// Separate from the K-ROM: POD coordinates are not Koopman observables.
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct PodGalerkinModel {
    std::vector<double> modes;          // n x r, row-major, orthonormal columns
    std::vector<double> reference;      // n
    std::vector<double> singularValues; // all retained-candidate values, descending
    std::vector<double> A;              // r x r
    std::vector<double> B;              // r x m (empty without an input)
    std::size_t states{}, rank{}, inputs{};
    double energyFraction{};            // retained sum of s^2 / total

    [[nodiscard]] std::vector<double> project(const double* state) const;
    [[nodiscard]] std::vector<double> reconstruct(const double* coordinates) const;
    // Nonlinear Galerkin right-hand side Phi^T f(t, x_r + Phi z) (+ Phi^T B(x) u for one input).
    [[nodiscard]] std::vector<double> reducedRhs(const ODESystem& system, double time,
                                                 const double* coordinates, double input = 0.0) const;
};

// rank 0 keeps the smallest rank reaching `energy` of the snapshot energy. reference: "mean" or
// "zero". The linear operators come from the system at the reference.
[[nodiscard]] PodGalerkinModel fitPodGalerkin(const std::vector<double>& snapshots, std::size_t samples,
                                              std::size_t states, std::size_t rank, double energy,
                                              const std::string& reference, const ODESystem& system);

}  // namespace dynsyslab::core

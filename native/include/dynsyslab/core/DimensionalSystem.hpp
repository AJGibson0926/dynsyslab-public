#pragma once

// Dimensional integration of a nondimensionalised model.
//
// A model advanced in dimensionless coordinates (tau = t / T, x = X / S, u = U / P) declares
// its scales in CoordinateScaling (secondsPerModelTime = T, physicalStatePerModelState = S,
// physicalInputPerModelInput = P). The adapter integrates the same equations in physical
// coordinates by the exact change of variables
//     dX/dt = S f(t / T, X / S, U / P) / T,   dX/dt|_jac = diag(S) J diag(S)^-1 / T,
// so a solve can be archived in seconds and SI units without a second implementation of the
// physics (Buckingham, Phys. Rev. 4 (1914) 345 for the scaling; the bubble scales follow
// Gibson, Calvisi & Yee, SIAM J. Appl. Dyn. Syst. 21 (2022), doi:10.1137/21M1456078).

#include <dynsyslab/core/System.hpp>

#include <memory>
#include <string>

namespace dynsyslab::core {

// True when the system declares physical time, input and per-state scales.
[[nodiscard]] bool supportsDimensionalIntegration(const ODESystem& system);

// The model in physical coordinates (descriptor.preferredCoordinates = "dimensional").
[[nodiscard]] std::unique_ptr<ODESystem> dimensionalSystem(std::unique_ptr<ODESystem> model);
// The same, published under its own descriptor and identifier (a registered dimensional
// system such as rayleigh-plesset-dimensional-baseline); the labels are taken as given.
[[nodiscard]] std::unique_ptr<ODESystem> dimensionalSystem(std::unique_ptr<ODESystem> model,
                                                           const SystemDescriptor& descriptor);

// requested: "automatic" (or the system's own coordinates) returns the model unchanged;
// "dimensional" wraps a dimensionless model with declared scales; anything else throws.
[[nodiscard]] std::unique_ptr<ODESystem> systemInCoordinates(std::unique_ptr<ODESystem> model,
                                                             const std::string& requested);

}  // namespace dynsyslab::core

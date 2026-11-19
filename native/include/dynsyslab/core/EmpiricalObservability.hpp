#pragma once
// Local (empirical) observability of a nonlinear ODE system through a named
// observable, in the sense of Lall, Marsden and Glavaski: from a state x0 at
// time t0, perturb each state component by +/- e_i, simulate over a horizon with
// the same input u(t), and integrate the output sensitivities
//
//   W_ij = int_0^T  dy_i(t) . dy_j(t) dt,   dy_i = (y(x0 + e_i) - y(x0 - e_i)) / (2 e_i).
//
// For a smooth system and small e this is the finite-horizon observability
// gramian of the linearization along the trajectory through x0. Its eigenvalues
// measure how much output energy each state direction produces; with white
// measurement noise of intensity V, V W^-1 is the Cramer-Rao bound on the
// covariance of an unbiased estimate of x0 from the output over [t0, t0 + T].
#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace dynsyslab::core {

struct EmpiricalObservabilityOptions {
    double horizon{1.0};            // output window length (model time units)
    std::size_t samples{201};       // output samples on [t0, t0 + horizon] (trapezoid rule)
    std::size_t substeps{20};       // RK4 steps between samples
    std::vector<double> perturbation;  // e_i per state (absolute); empty = 1e-6 * max(1, |x0_i|)
    std::function<double(double)> input;  // u(t); empty = unforced
};

struct EmpiricalObservability {
    std::size_t states{}, outputs{};
    std::vector<double> gramian;        // states x states (row-major)
    std::vector<double> eigenvalues;    // descending
    std::vector<double> eigenvectors;   // states x states, eigenvector k in column k
    double conditionNumber{};           // lambda_max / lambda_min (infinite if lambda_min <= 0)
    std::vector<double> outputJacobian; // outputs x states: instantaneous dh/dx at x0 (central difference)
    bool valid{};                       // every perturbed trajectory stayed finite and inside the validity domain
};

[[nodiscard]] EmpiricalObservability empiricalObservability(const ODESystem& system, const Observable& observable,
    const std::vector<double>& state, double time, const EmpiricalObservabilityOptions& options);

// The same gramian in scaled coordinates x = S z (S = diag(scale)): S W S, re-decomposed,
// so eigenvalues of states with different units become comparable.
[[nodiscard]] EmpiricalObservability scaledObservability(const EmpiricalObservability& result,
                                                         const std::vector<double>& scale);

// Cramer-Rao standard deviations sqrt(V (W^-1)_ii) for measurement-noise intensity V
// (infinite for a direction with zero gramian eigenvalue).
[[nodiscard]] std::vector<double> estimationBoundStd(const EmpiricalObservability& result, double noiseIntensity);

}  // namespace dynsyslab::core

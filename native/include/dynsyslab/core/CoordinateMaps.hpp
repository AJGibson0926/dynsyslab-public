#pragma once
// Coordinate maps acting on functions. A coordinate map is an Observable h: X -> Q of a state space
// (an instantaneous coordinate transformation; history windows are measurement maps, not maps of
// the state). On functions it acts two ways:
//
//   pull-back     h^* G = G o h. Defined wherever h(x) lies in G's domain; this is how a function of
//                 q becomes a function of x (FunctionDomain's state map does exactly this).
//   push-forward  h_* g(q) = E[g | h = q], the fibre average over source samples. It is a function of
//                 q only where g is constant on each fibre h^{-1}(q); elsewhere the fibre spread
//                 measures how far it fails (status `conflict`), and points of Q that no source
//                 sample reaches have no preimage. (Conditional expectation; Lasota and Mackey,
//                 Chaos, Fractals, and Noise, 2nd ed. 1994, section 3.1.)
//
// The flow map Phi^t of a system is a coordinate map X -> X, and its pull-back is the Koopman
// operator, U^t g = g o Phi^t (Koopman 1931, PNAS 17:315-318); koopmanAction samples it on a grid.
#include <dynsyslab/core/DomainGrid.hpp>
#include <dynsyslab/core/System.hpp>

#include <vector>

namespace dynsyslab::core {

// G o h with the exact chain-rule Jacobian when both have one.
[[nodiscard]] Observable pullBack(const Observable& function, const Observable& map);

// Phi^t(x) of the unforced system by classical fourth-order Runge-Kutta with steps no longer than
// `maximumStep` (Hairer, Norsett and Wanner, Solving ODEs I, 1993, section II.1). Its Jacobian is
// a finite difference.
[[nodiscard]] Observable flowMap(std::shared_ptr<const ODESystem> system, double time, double maximumStep = 1e-2);

// g o Phi^{k dt} for k = 0..steps on the grid (the grid lives in the system's state space).
// Points whose orbit leaves the valid state set become undefined.
[[nodiscard]] std::vector<GridField> koopmanAction(const Observable& function, const ODESystem& system, const DomainGrid& grid,
                                                   double timeStep, std::size_t steps, double maximumStep = 1e-2);

struct PushForwardOptions {
    double tolerance{1e-3};  // relative to the largest |g| over the source samples
    // The source grid is refined by this factor on each axis before sampling, so that cells inside
    // the image of h are not missed where h stretches the source (a missed cell has no preimage).
    std::size_t oversampling{1};
};
// h_* g on `target` from the source samples `source`: per target cell the fibre mean (values), the
// fibre spread max |g - mean| (spread), the sample count (counts) and the status (valid, conflict,
// noPreimage).
[[nodiscard]] GridField pushForward(const Observable& function, const Observable& map, const DomainGrid& source,
                                    const DomainGrid& target, PushForwardOptions options = {});

// |h^* h_* g - g| on the source grid: where the push-forward loses g (status kept from the target
// cell: noPreimage, conflict). Zero exactly where g factors through h.
[[nodiscard]] GridField roundTripError(const Observable& function, const Observable& map, const DomainGrid& source,
                                       const DomainGrid& target, PushForwardOptions options = {});

// sqrt(det(Dh Dh^T)) on the grid: the local volume factor of h (area formula; Evans and Gariepy,
// Measure Theory and Fine Properties of Functions, 2015, section 3.3). Zero where Dh loses rank, so
// its zero set is the fold set where sheets of h meet.
[[nodiscard]] GridField mapVolumeFactor(const Observable& map, const DomainGrid& grid);

}  // namespace dynsyslab::core

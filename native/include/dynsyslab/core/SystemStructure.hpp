#pragma once
// Structural analysis: controllable/observable subspaces (staircase forms),
// the four-part Kalman decomposition, minimal realization and a combined
// controllability / observability / stabilizability / detectability report.
// Pbh.hpp supplies the mode-by-mode Hautus tests used as an independent check.
#include <dynsyslab/core/StateSpace.hpp>

#include <complex>
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct StaircaseForm {
    std::vector<double> transform; // n x n orthogonal; first `reachable` columns span the subspace
    std::size_t reachable{};       // controllable (or observable) dimension
    std::vector<std::size_t> blockSizes; // staircase block sizes
};
// Orthonormal Krylov staircase of range[B, AB, A^2B, ...].
[[nodiscard]] StaircaseForm controllabilityStaircase(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);
// Dual staircase of range[C^T, A^T C^T, ...]: the observable subspace.
[[nodiscard]] StaircaseForm observabilityStaircase(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);

// x = T z with z = [z_co; z_cno; z_nco; z_ncno] (c: controllable, o: observable).
struct KalmanDecomposition {
    StateSpaceModel transformed;
    std::vector<double> transform, inverse; // n x n
    std::size_t controllableObservable{}, controllableUnobservable{};
    std::size_t uncontrollableObservable{}, uncontrollableUnobservable{};
    double structureResidual{}; // relative norm of blocks that must vanish
};
[[nodiscard]] KalmanDecomposition kalmanDecomposition(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);

struct MinimalRealization {
    StateSpaceModel minimal;
    std::vector<double> projection; // r x n: z = P x
    std::vector<double> lift;       // n x r
    std::size_t removedStates{};
};
// The controllable-and-observable part; same transfer function, minimal order.
[[nodiscard]] MinimalRealization minimalRealization(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);

struct StructureReport {
    std::size_t controllableDimension{}, observableDimension{}, minimalOrder{};
    bool controllable{}, observable{}, stabilizable{}, detectable{};
    std::vector<std::complex<double>> uncontrollableModes, unobservableModes;
    bool pbhAgrees{}; // Hautus tests reach the same four verdicts
    // Hautus verdicts and smallest relative PBH margins sigma_min/sigma_max over
    // all eigenvalues. Near-repeated eigenvalues (e.g. a conjugate pair split
    // twice) give tiny margins, where the two numerical tests may disagree.
    bool pbhControllable{}, pbhObservable{}, pbhStabilizable{}, pbhDetectable{};
    double pbhControllabilityMargin{}, pbhObservabilityMargin{};
};
// Kalman-decomposition verdicts only (PBH fields left at their defaults).
[[nodiscard]] StructureReport kalmanStructureReport(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);
// Kalman verdicts cross-checked by the Hautus (PBH) tests.
[[nodiscard]] StructureReport analyzeStructure(const StateSpaceModel& model,
    double relativeTolerance = 1.0e-10);

}  // namespace dynsyslab::core

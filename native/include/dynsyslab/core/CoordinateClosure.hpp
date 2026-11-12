#pragma once

#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace dynsyslab::core {

// A scalar observable jet through order k is [h,L_f h,...,L_f^k h].
// Vector observables contribute every component at each derivative order.
[[nodiscard]] Observable observableJet(const Observable& observable,
    const std::vector<AnalyticExpression>& autonomousField,
    std::size_t highestDerivativeOrder);

struct ObservableClosureOptions {
    std::vector<double> stateMinimum;
    std::vector<double> stateMaximum;
    std::size_t highestDerivativeOrder{4};
    std::size_t sampleCount{512};
    std::size_t pairCandidateCount{64};
    std::size_t maximumRefinementIterations{80};
    std::uint64_t seed{1};
    double minimumNormalizedStateSeparation{0.05};
    double collisionTolerance{1.0e-9};
    double successorDifferenceTolerance{1.0e-8};
    double rankTolerance{1.0e-10};
    // Flattened [left state, right state] blocks, two states per seed.
    std::vector<double> pairSeeds;
};

struct ObservableCollisionWitness {
    bool found{};
    std::vector<double> leftState;
    std::vector<double> rightState;
    std::vector<double> leftCoordinates;
    std::vector<double> rightCoordinates;
    std::vector<double> successorDerivativeLeft;
    std::vector<double> successorDerivativeRight;
    double normalizedResidual{};
    double normalizedSuccessorDifference{};
    double normalizedStateSeparation{};
    std::size_t refinementIterations{};
};

struct ObservableClosureOrderReport {
    std::size_t highestDerivativeOrder{};
    std::size_t coordinateCount{};
    std::size_t minimumSampledRank{};
    double minimumSampledSingularValue{};
    std::vector<double> minimumRankState;
    // A coordinate collision disproves injectivity of the jet map. It only
    // obstructs closed jet dynamics when the successor derivative also differs.
    ObservableCollisionWitness coordinateCollision;
    ObservableCollisionWitness closureObstruction;
};

struct ObservableClosureReport {
    std::size_t stateDimension{};
    std::size_t observableDimension{};
    std::size_t requestedSamples{};
    std::size_t acceptedSamples{};
    std::vector<ObservableClosureOrderReport> orders;
};

// Deterministic numerical side-car. A closure obstruction (equal jets at
// distinct states with unequal successor derivatives) disproves closure. A
// coordinate collision merely disproves injectivity. Failure to find either is
// finite-search evidence, not a global proof.
[[nodiscard]] ObservableClosureReport analyzeObservableClosure(
    const ODESystem& system, const Observable& observable,
    const ObservableClosureOptions& options);

}  // namespace dynsyslab::core

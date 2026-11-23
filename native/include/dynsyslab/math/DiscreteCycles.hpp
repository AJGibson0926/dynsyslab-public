#pragma once
#include <dynsyslab/core/Json.hpp>
#include <cstddef>
#include <vector>
namespace dynsyslab::math {
struct CycleDefect {
    std::size_t period{}, comparisons{};
    double maximum{}, rms{};
};
struct SampledCycleResult {
    std::vector<CycleDefect> defects;
    std::size_t detectedPeriod{}; // zero means unresolved in the supplied candidate list
    std::vector<double> foldedValues; // phase-major; only for detected period
    double foldingRms{};
};
// Finite-sample recurrence diagnostic, not an existence/uniqueness certificate.
// Scales define the Euclidean norm; nonzero coordinatePeriods use circular distance.
SampledCycleResult diagnoseSampledCycles(const std::vector<double>& values,std::size_t dimension,
    const std::vector<std::size_t>& candidates,const std::vector<double>& scales,
    const std::vector<double>& coordinatePeriods,double tolerance,std::size_t minimumRepeats=4);
// Uniform measure on a supplied cycle; DFT coefficients and cyclic Koopman atoms.
// See Eisner et al. (2015), Operator Theoretic Aspects of Ergodic Theory, ch. 14/17.
core::json::Value discreteCycleSpectrum(const std::vector<double>& cycle,std::size_t dimension);
// Conditional expectation onto residues modulo quotientPeriod, which must
// divide the supplied cycle length. Standard finite-group Haar projection.
// Arithmetic coordinates only: unwrap angular coordinates before averaging.
core::json::Value projectDiscreteCycle(const std::vector<double>& cycle,
    std::size_t dimension,std::size_t quotientPeriod);
// Bounded exhaustive cyclic registration, preserving forward time order.
core::json::Value alignDiscreteCycles(const std::vector<double>& reference,
    const std::vector<double>& candidate,std::size_t dimension);
core::json::Value evaluateDiscreteCycleQuery(const core::json::Value&);
}

#pragma once
#include <dynsyslab/core/SindyDynamics.hpp>
namespace dynsyslab::core {
struct DynamicsLibraryOptions {
    std::size_t polynomialDegree{4};
    std::size_t inverseDegree{5};
    bool normalize{true};
    bool trigonometric{true},reciprocal{true},logarithmic{true},exponential{true},hyperbolic{true};
    bool mixedTerms{true};
    std::vector<double> frequencies{0.5,1,2};
    std::vector<double> exponentialRates{-1,-0.5,0.5,1};
};
struct DynamicsLibrary {
    Observable observable;
    std::vector<std::string> excluded;
};
// Build candidates from TRAINING coordinates only. All evaluation/derivative
// kernels belong to the shared Observable expression engine. Domain-restricted
// terms are excluded on a zero-crossing or nonpositive training interval, not
// clamped. The returned expressions are frozen for later data and rollouts.
DynamicsLibrary buildDynamicsLibrary(const TimeSeries& training,
    const DynamicsLibraryOptions& options={});
}

#pragma once
#include <complex>
#include <vector>
namespace dynsyslab::math {
struct InvariantReadoutFloor {
    std::complex<double> bestConstant;
    double meanSquaredFloor{}, totalWeight{};
};
// For observables with ZERO conditional mean on every invariant component,
// all linear combinations of their time translates have zero conditional mean.
// Orthogonal projection then bounds target reconstruction error from below by
// the squared component means (or their variance when a constant is permitted).
// Means and nonnegative component weights are supplied, not inferred. Weights
// are normalized. Sampled/approximate means produce a numerical estimate, not
// a certified continuous-measure bound. No nonlinear-readout restriction follows.
InvariantReadoutFloor zeroMeanReadoutFloor(
    const std::vector<std::complex<double>>& targetMeans,
    const std::vector<double>& componentWeights,bool allowConstant=true);
}

#pragma once
#include <dynsyslab/core/System.hpp>
namespace dynsyslab::core {
struct TorusSample {
    std::vector<std::size_t> shape;
    std::vector<double> angularFrequencies,states;
    std::vector<std::complex<double>> observables; // sample-major
    std::size_t stateDimension{},channels{};
    double maximumInvarianceResidual{},maximumWrapResidual{};
};
// Uniform product-Haar sampling of a supplied torus map; a known periodic
// trajectory can supply a one-angle map. Does not infer invariant tori from data.
// Invariance/wrap residuals are sampled numerical evidence, not a KAM certificate.
TorusSample sampleInvariantTorus(const ODESystem&,const AnalyticTrajectory&,
    const std::vector<Observable>&,const std::vector<std::size_t>&shape,double initialTime=0);
}

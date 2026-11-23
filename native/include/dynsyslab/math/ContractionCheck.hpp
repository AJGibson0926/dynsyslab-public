#pragma once
#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
#include <functional>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Logarithmic-norm test for a unique attracting steady response of x'=F(x,u(t)). If
// mu(Theta dF/dx Theta^{-1}) <= -c < 0 on a forward-invariant convex set, all solutions in it
// converge to each other: Lohmiller & Slotine, Automatica 34 (1998) 683-696; the Demidovich
// condition in Pavlov, van de Wouw & Nijmeijer, Uniform Output Regulation of Nonlinear Systems
// (2006), and their frequency-response use in IEEE Trans. Automat. Control 52 (2007) 1159-1165.
// Log norms: Dahlquist (1958); Soderlind, BIT 46 (2006) 631-652. Theta=diag(weights) here.
struct LogNorms {double two{},one{},infinity{};};
LogNorms logarithmicNorms(const ComplexMatrix& J,const std::vector<double>& weights={});
using JacobianFunction=std::function<ComplexMatrix(const std::vector<double>& point)>;
using BoxField=std::function<std::vector<double>(const std::vector<double>& point)>;
struct SampledContractionInput {
    std::size_t stateDimension{};
    std::vector<double> lower,upper;   // Box over (x..., declared inputs...).
    std::size_t samplesPerAxis{5};     // Tensor grid including the corners.
    std::vector<double> weights;       // Positive Theta diagonal; empty is the identity.
    JacobianFunction jacobian;         // dF/dx at a point of the box.
    BoxField field;                    // Optional F, to test forward invariance on the state faces.
};
struct SampledContractionReport {
    std::string label{"sampled"};      // A sampled result is evidence, never a proof.
    std::size_t samples{};
    LogNorms worst;
    std::vector<double> worstPoint;    // Of the best (smallest) of the three worst-case norms.
    std::string bestNorm;
    double rateEstimate{};             // Smallest worst-case norm over the three norms.
    bool contractiveAtSamples{false};
    std::optional<bool> invariantAtSamples;
    std::size_t invarianceViolations{};
    std::vector<std::string> warnings;
};
SampledContractionReport sampledContraction(const SampledContractionInput&);
// Entrywise enclosure lower<=J<=upper, supplied and trusted from the caller. Upper bounds:
// mu_inf, mu_1 by Gershgorin-type sums, and mu_2<=mu_2(mid)+||radius||_F.
struct IntervalContractionReport {
    std::string label{"caller_interval"};
    LogNorms upperBound;
    std::string bestNorm;
    double rateBound{};
    bool contractiveForInterval{false};
    std::string scope;
};
IntervalContractionReport intervalContraction(const std::vector<double>& lower,const std::vector<double>& upper,
    std::size_t dimension,const std::vector<double>& weights={});
}

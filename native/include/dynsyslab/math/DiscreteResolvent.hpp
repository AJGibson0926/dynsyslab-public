#pragma once
#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <optional>
#include <vector>
namespace dynsyslab::math {
// z-domain Koopman resolvent of the sampled flow U=exp(L h). The unilateral z-transform of
// y_k=(U^k f)(x0) is X(z)=[z (zI-U)^{-1} f](x0): Susuki, Mauroy & Drmac, Koopman Resolvents of
// Nonlinear Discrete-Time Systems, Proc. ECC (2024) 627-632. A pole of order q contributes
//   a_q z/(z-mu)^q  <->  a_q C(k,q-1) mu^{k-q+1}   (Jury, Theory and Application of the
// z-Transform Method, 1964). Converting Jordan coefficients uses Stirling numbers.
struct DiscretePole {
    std::complex<double> multiplier;
    ComplexVector laurent;   // laurent[q-1] multiplies z/(z-mu)^q. Order at most 8.
};
struct DiscretePoleModel {std::vector<DiscretePole> poles;double sampleInterval{1};};
void validateDiscreteModel(const DiscretePoleModel&);
std::complex<double> evaluateDiscreteResolvent(const DiscretePoleModel&,std::complex<double> z);
std::complex<double> evaluateDiscreteSequence(const DiscretePoleModel&,std::size_t k);
// Exact sampling y_k=y(k h) of a continuous pole model, with mu=exp(lambda h).
DiscretePoleModel discretizePoleModel(const PoleResidueModel&,double sampleInterval);
// Inverse on chosen branches: lambda=(log|mu|+i(arg mu+2 pi branch))/h. Empty means principal.
// Sampling cannot determine the branch; it must come from outside the samples.
PoleResidueModel continuousPoleModel(const DiscretePoleModel&,const std::vector<int>& branches={});
struct AliasEntry {
    std::size_t index{};
    std::complex<double> eigenvalue,multiplier,principal;
    int turns{};            // eigenvalue = principal + 2 pi i turns / h.
    bool aliased{false};    // |Im lambda| above the Nyquist rate pi/h.
    bool onNyquist{false};  // Im lambda = +-pi/h, where the principal branch is ambiguous.
};
struct AliasCollision {std::size_t first{},second{};double separation{};};
struct AliasReport {
    double nyquist{};
    std::vector<AliasEntry> entries;
    std::vector<AliasCollision> collisions;  // Distinct lambda, same mu: merged by sampling.
};
AliasReport aliasReport(const PoleResidueModel&,double sampleInterval,double tolerance=1e-10);
}

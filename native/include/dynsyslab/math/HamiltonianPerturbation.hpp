#pragma once
#include <dynsyslab/math/TorusFourier.hpp>

namespace dynsyslab::math {
struct HamiltonianFirstOrder {
    TorusFourierSeries generator, retained;
    std::vector<TorusFourierSeries> actionCorrection, angleCorrection;
    std::vector<TorusFourierSeries> retainedActionGradient;
    std::vector<double> meanFrequencyShift;
    double maximumInverseDivisor{};
};
// H(I,theta)=H0(I)+epsilon*H1(I,theta), at one specified action I0.
// Frequencies in H1 are grad H0(I0); hessian is row-major D^2 H0(I0).
// gradient[j] is partial H1/partial I_j, holding theta fixed.
// Old coordinates = (I,theta)+epsilon*(-partial_theta chi,partial_I chi)
// with omega.partial_theta chi = H1-retained. All results are coefficients
// of epsilon, NOT full nonlinear transformations. Small divisors are retained.
HamiltonianFirstOrder firstOrderHamiltonian(
    const TorusFourierSeries& perturbation,
    const std::vector<TorusFourierSeries>& actionGradient,
    const std::vector<double>& hessian, double divisorTolerance);

struct TwistDiagnostics {
    std::vector<double> kolmogorovSingularValues, isoenergeticSingularValues;
    double kolmogorovRelativeMinimum{}, isoenergeticRelativeMinimum{};
    bool kolmogorovNumericallyNondegenerate{}, isoenergeticNumericallyNondegenerate{};
};
bool twistDiagnosticsAvailable();
// Singular values of D^2 H0 and [D^2 H0,omega;omega^T,0]. Numerical diagnostics
// in the supplied action units, not interval certificates or a KAM theorem.
TwistDiagnostics hamiltonianTwist(const std::vector<double>& frequencies,
    const std::vector<double>& hessian, double relativeTolerance=1e-10);
}

#pragma once
#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
#include <functional>
#include <optional>
#include <string>
namespace dynsyslab::math {
// One-sided Laplace transform of sampled signals, Y(s)=int_0^inf e^{-st} y(t) dt, which on a
// Koopman orbit y(t)=f(S^t x0) equals [R(s;L)f](x0): Susuki, Mauroy & Mezic, SIAM J. Appl. Dyn.
// Syst. 20(4) (2021), eq. (4), arXiv:2009.11544.
//
// Oscillatory panels are integrated exactly against e^{-st} on a piecewise linear or cubic
// Hermite interpolant (the product-integration idea of Filon, Proc. R. Soc. Edinburgh 49
// (1930)), so accuracy does not degrade as |Im s| h grows.
enum class LaplaceInterpolant {Linear,CubicHermite};
struct SampledSignal {
    std::vector<double> times;  // Strictly increasing, nonnegative, need not be uniform.
    ComplexVector values;
    ComplexVector derivatives;  // Optional exact y'(t) for the cubic interpolant.
};
// Declared or fitted bound |y(t)| <= constant * exp(rate * t) on the tail.
struct GrowthBound {double rate{},constant{};bool declared{false};};
// Upper-envelope regression over blocks of samples. A finite window can never prove a growth
// bound: this is an ESTIMATE of the abscissa of convergence, not a certificate.
GrowthBound fitGrowthBound(const SampledSignal&,std::size_t blocks=8);
struct LaplaceTransformOptions {
    LaplaceInterpolant interpolant{LaplaceInterpolant::CubicHermite};
    std::optional<GrowthBound> growth;  // Empty: fitted from the samples.
};
struct LaplaceEstimate {
    std::complex<double> s,value;       // value is the integral over the sampled window.
    GrowthBound growth;
    bool insideRegion{false};           // Re s > growth.rate.
    std::optional<double> tailBound;    // |int_T^inf| <= C e^{(beta-sigma)T}/(sigma-beta).
};
LaplaceEstimate laplaceTransform(const SampledSignal&,std::complex<double> s,const LaplaceTransformOptions& ={});
// Window integral only, without any growth or tail estimate.
std::complex<double> laplaceWindowIntegral(const SampledSignal&,std::complex<double> s,LaplaceInterpolant);
// Node derivatives of the local 5-point interpolating polynomial on a nonuniform grid.
ComplexVector hermiteNodeDerivatives(const std::vector<double>& x,const ComplexVector& y);
// Unilateral z-transform X(z)=sum_k y_k z^{-k} for uniform samples.
std::complex<double> zTransform(const ComplexVector& values,std::complex<double> z);
// Exact map z=exp(s h), and its inverse on branch k: s=(log|z|+i(arg z+2 pi k))/h.
std::complex<double> zFromLaplace(std::complex<double> s,double sampleInterval);
std::complex<double> laplaceFromZ(std::complex<double> z,double sampleInterval,int branch=0);
using LaplaceFunction=std::function<std::complex<double>(std::complex<double>)>;
struct InverseLaplaceOptions {
    std::size_t nodes{32};                   // Even number of Talbot nodes.
    double shift{0};                         // Contour shift: all singularities need Re p < shift for Bromwich.
    ComplexVector singularities;             // Optional declared singularities, used to check scope.
    double dampingParameter{18.4};           // Bromwich/Euler A: discretization error about e^{-A}.
    std::size_t eulerTerms{11};              // Euler binomial averaging order.
    std::size_t seriesTerms{0};              // Bromwich base terms; 0 selects from the singularities.
};
struct InverseLaplaceResult {
    std::complex<double> value;
    bool scopeVerified{false};   // Every declared singularity lies inside the contour or left of the line.
    std::string method,scope;
};
// Modified Talbot contour with midpoint rule: Talbot, IMA J. Appl. Math. 23 (1979),
// parameters of Trefethen, Weideman & Schmelzer, BIT 46 (2006).
// For F analytic outside a sector around the negative axis.
// Poles off the negative axis converge more slowly than the 3.89^-N of the paper: use N>=48
// for oscillatory F. Not for poles with |Im p| t beyond about 0.83 N, nor singularities on the right.
InverseLaplaceResult talbotInverse(const LaplaceFunction&,double t,const InverseLaplaceOptions& ={});
// Bromwich trapezoid with Euler summation: Abate & Whitt, ORSA J. Comput. 7 (1995),
// Suited to poles near the imaginary axis at large |Im p| t. Needs f
// continuous at t and bounded growth below shift; discretization error is about e^{-A}.
InverseLaplaceResult bromwichInverse(const LaplaceFunction&,double t,const InverseLaplaceOptions& ={});
}

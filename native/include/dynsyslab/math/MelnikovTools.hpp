#pragma once
#include <dynsyslab/math/OrbitFamily.hpp>
#include <complex>
#include <functional>
#include <vector>
namespace dynsyslab::math {
// First-order Melnikov theory on one homoclinic or heteroclinic connection.
// Melnikov, Trans. Moscow Math. Soc. 12 (1963); Guckenheimer & Holmes, Nonlinear
// Oscillations, Dynamical Systems, and Bifurcations of Vector Fields (1983), section 4.5.
// The midpoint must lie on the connection; forward and backward flows are integrated to
// tailRadius from their saddles and the remainders are closed by exponential tails.
struct ConnectionProblem {
    PlanarField drift;PlanarPoint midpoint{},forwardSaddle{},backwardSaddle{};
    double forwardRate{},backwardRate{};  // positive hyperbolic rates at the two ends
    double tailRadius{1e-4},maximumTime{1e4};FlowOptions flow;
};
using ConnectionIntegrand=std::function<std::complex<double>(const PlanarPoint&,double)>;
struct ConnectionIntegral {
    std::complex<double> value,core,forwardTail,backwardTail;
    double forwardTime{},backwardTime{},forwardRadius{},backwardRadius{},forwardDecay{},backwardDecay{};
};
// int g(x(t),t) dt. oscillation is the angular frequency of any e^{i w t} factor in g, used
// by the tail g(T)/(r-i w). The decay r is measured from the integrand, bounded to
// [rate,4 rate]. The tail is a leading-order estimate, not a rigorous enclosure.
ConnectionIntegral connectionIntegral(const ConnectionProblem&,const ConnectionIntegrand&,double oscillation=0);
// int f0 ^ Y dt.
ConnectionIntegral melnikovWedge(const ConnectionProblem&,const PlanarField& perturbation);
// int Q dt for Hdot=-d Q. Equals -int f0^Y when f0 is Hamiltonian in canonical order.
ConnectionIntegral dampingIntegral(const ConnectionProblem&,const PlanarScalar& dissipation);
// A(w)=int f0^Y e^{i w t} dt for forcing F Y cos(w t+t0); its modulus is the amplitude.
ConnectionIntegral forcingIntegral(const ConnectionProblem&,const PlanarField& forcing,double w);
// Simple-zero threshold of M(t0)=-d D+F|A|cos(w t0+arg A): F/d=D/|A|.
double melnikovThreshold(double dampingIntegralValue,double forcingAmplitude);
struct LobeFlux {double ratio{},area{},fluxPerTime{};bool transverse{};};
// First-order lobe area int (M)+ dt0 over one forcing period and its mean rate.
// Rom-Kedar & Wiggins, Arch. Rational Mech. Anal. 109 (1990); MacKay, Meiss & Percival,
// Physica D 13 (1984). Zero when the separatrices do not intersect at first order.
LobeFlux lobeFlux(double d,double dampingIntegralValue,double F,double forcingAmplitude,double w);
// P_b=M_b/sum M for adiabatic separatrix crossing with nonnegative basin growth rates.
// Neishtadt, J. Appl. Math. Mech. 39 (1975); Henrard, Celestial Mechanics 27 (1982).
// Probabilities are for an ensemble with uniformly distributed phase at crossing.
std::vector<double> captureProbabilities(const std::vector<double>& lobeIntegrals);
struct ShellMelnikovRow {int harmonic{};double energy{},period{},coefficient{},product{},ratio{};};
// Compares T(E_n)|y_n(E_n)| on a shell ladder with the separatrix amplitude |A|.
std::vector<ShellMelnikovRow> shellMelnikovComparison(const std::vector<int>& harmonics,
    const std::vector<double>& energies,const std::vector<double>& periods,
    const std::vector<double>& coefficientMagnitudes,double amplitude);
}

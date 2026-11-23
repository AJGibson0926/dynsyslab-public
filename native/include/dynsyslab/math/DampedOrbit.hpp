#pragma once
#include <dynsyslab/math/OrbitFamily.hpp>
#include <complex>
#include <vector>
namespace dynsyslab::math {
// Weakly damped planar flow xdot=f0+d*Y with Hdot=-d*Q. First-order averaging:
// Arnold, Kozlov & Neishtadt, Mathematical Aspects of Classical and Celestial Mechanics,
// chapter 5. Valid away from the separatrix and for small d; no remainder bound is computed.
struct DampedAverage {double mean{};std::vector<double> chi;FourierSeries chiSeries;};
// Q at equal phases: D=<Q>, chi=(1/Omega) int (Q-D) dphi with zero mean.
DampedAverage dampedAverage(const std::vector<double>& qSamples,double frequency);
struct DampedTable {std::vector<double> energy,frequency,damping;};
struct DampedPhase {
    std::vector<double> energy,phaseFunction,phaseSlope,amplitude;
    double equilibriumFrequency{},realEigenvalue{},dampingScale{};
};
// S(E)=int_{Emin}^{E} (Omega-Omega0)/D dE, S'=(Omega-Omega0)/D, and a'=-Re(lambda)/(d D)
// with a=0 at the last table energy. The energy measure is what makes the WKB phase
// advance at Omega0 under Edot=-d Q. A first node with D=0 is treated as the equilibrium:
// S' is extrapolated there and a is -infinity for a stable eigenvalue.
DampedPhase dampedPhase(const DampedTable&,double equilibriumFrequency,double realEigenvalue,double dampingScale);
// exp(a(E)+i[phi+S(E)/d+S'(E) chi]); phi is measured from the section turning point.
std::complex<double> wkbEigenfunction(const DampedPhase&,double energy,double phase,double chi);
struct LinearMode {std::complex<double> eigenvalue;std::array<std::complex<double>,2> left;};
// Upper-half-plane (or slowest real) eigenvalue of [[a,b],[c,d]], left vector with w.v=1.
LinearMode linearMode(double a,double b,double c,double d);
struct LaplaceOptions {FlowOptions flow;double settleRadius{1e-6},maximumTime{1e7};};
struct LaplaceValue {std::complex<double> value;double time{},finalRadius{};};
// psi(x)=lim exp(-lambda t) w.(x(t)-x*). Mauroy & Mezic, IEEE Trans. Automat. Control
// 61 (2016), arXiv:1408.1379. Truncation error is O(settleRadius) relative.
LaplaceValue laplaceEigenfunction(const PlanarField&,const PlanarPoint& equilibrium,const LinearMode&,
    const PlanarPoint& x,const LaplaceOptions& ={});
inline double isochronPhase(std::complex<double> psi){return std::arg(psi);}
inline double isostableAmplitude(std::complex<double> psi){return std::abs(psi);}
// W=psi^n |psi|^s, s=(i w-n lambda)/Re(lambda): eigenvalue i w, and |W|=1 for real w.
// Defined on the punctured basin psi!=0.
std::complex<double> annularEigenfunction(std::complex<double> psi,std::complex<double> lambda,int n,double w);
// Stationary-phase weight sqrt(2 pi d D/|n Omega'|) at a shell n Omega(E)=w.
double shellWeight(double dampingScale,double damping,int n,double frequencySlope);
}

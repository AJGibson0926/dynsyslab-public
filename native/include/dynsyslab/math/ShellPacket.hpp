#pragma once
#include <dynsyslab/math/SpectralFilters.hpp>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Predicted spectral packets on a family of periodic orbits ("shells").
// Kernel: Colbrook, Drysdale & Horning, Rigged DMD, arXiv:2405.00782, reused
// through RationalCircleKernel, for which w(alpha;theta)=kernel(theta-alpha).
// Measure-preserving context: Colbrook, mpEDMD, SIAM J. Numer. Anal. 61 (2023).
struct OrbitShell {
    double energy{};
    FourierSeries series; // angularFrequency=Omega(E), coefficients g_n(E).
    double weight{1};     // Explicit nonnegative quadrature mass.
};
// Smallest delta>0 with kernel(delta)=kernel(0)/2, by scan and bisection.
// Returns pi when the kernel does not fall to half its peak on the circle.
double kernelHalfWidth(const RationalCircleKernel& kernel);
double shellPacketWeight(double alpha,double theta,const RationalCircleKernel& kernel);
// psi_n = w(n*Omega*dt;theta)*g_n*exp(i*n*Omega*targetDelay*dt).
FourierSeries predictedShellPacket(const FourierSeries& orbit,double sampleInterval,double theta,
                                   const RationalCircleKernel& kernel,double targetDelay=0);
struct ShellPacketPrediction {
    std::vector<FourierSeries> packets;
    std::vector<double> shellEnergy; // Mean square of each predicted packet.
    double density{},totalMass{},zeroHarmonicMass{},halfWidth{};
};
ShellPacketPrediction predictShellPackets(const std::vector<OrbitShell>& shells,double sampleInterval,
    double theta,const RationalCircleKernel& kernel,double targetDelay=0);
// Solutions of n*Omega(E)*dt=theta+2*pi*k on a piecewise-linear Omega(E) table.
struct ShellAlias {
    int harmonic{};
    long long turns{};
    double energy{},angularFrequency{},frequencySlope{};
    std::optional<double> ridgeWidth; // h/(|n|*dt*|Omega'|); empty on a flat segment.
};
std::vector<ShellAlias> shellAliasCatalogue(const std::vector<OrbitShell>& shells,double sampleInterval,
    double theta,int maximumHarmonic,double halfWidth);
struct ShellPacketScore {
    double relativeL2Error{};
    std::complex<double> complexCorrelation;
    double complexCorrelationMagnitude{};
    std::optional<double> magnitudeCorrelation,energyProfileCorrelation; // Empty if a variance vanishes.
};
// shellIndex groups samples for the energy profile; empty treats each sample as
// its own group. weights are nonnegative quadrature masses; empty means uniform.
ShellPacketScore scoreShellPacket(const std::vector<std::complex<double>>& computed,
    const std::vector<std::complex<double>>& predicted,
    const std::vector<std::size_t>& shellIndex={},const std::vector<double>& weights={});
struct DelayTargetDiagnostic {
    long long lagsBefore{},lagsAfter{},retainedSymmetricLag{};
    double asymmetry{},symmetricTailBound{};
    bool oneSided{false};
    std::string warning;
};
// Window holds lags firstLag..lastLag; the target sits at targetDelay. The
// symmetric kernel needs lags on both sides of the target.
DelayTargetDiagnostic diagnoseDelayTarget(long long firstLag,long long lastLag,long long targetDelay,
    const RationalCircleKernel& kernel,double asymmetryTolerance=.25);
struct ResolvableShellEstimate { double lastHarmonic{}; bool empirical{true}; std::string scope; };
// n_last=c*w/h. c is an EMPIRICAL caller constant and w a caller angular scale
// in radians per sample (for example Omega*dt). Not a theorem or a guarantee.
ResolvableShellEstimate resolvableShellEstimate(double empiricalConstant,double angularScale,double halfWidth);
}

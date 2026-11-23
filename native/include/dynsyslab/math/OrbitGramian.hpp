#pragma once
#include <array>
#include <complex>
#include <map>
#include <optional>
#include <utility>
#include <vector>
namespace dynsyslab::math {
using ComplexVector=std::vector<std::complex<double>>;
// (a) xdot=diag(lambda)x+Bu, y=Cx. Rows of inputs/outputs are channels over the modes.
// Gramians solve the Lyapunov equations entrywise (Cauchy structure). Moore, IEEE TAC 26
// (1981); Antoulas, Approximation of Large-Scale Dynamical Systems (2005), chapter 4.
struct ModalGramians {
    std::size_t dimension{};ComplexVector controllability,observability;  // row-major
    double h2Norm{};std::vector<double> hankelSingularValues;
};
// No horizon: infinite-time Gramians, all Re(lambda)<0. With tau: finite-horizon form with
// the exact limit tau on a vanishing denominator. Dense Jacobi; limited to 128 modes.
ModalGramians modalGramians(const ComplexVector& eigenvalues,const std::vector<ComplexVector>& inputs,
    const std::vector<ComplexVector>& outputs,std::optional<double> horizon={});
// (b) xdd+d xd+a x=u, y=x: sigma=(sqrt(4a+d^2) +/- d)/(4 a d), larger first.
std::array<double,2> oscillatorHankel(double stiffness,double damping);
// 2x2 real symmetric Gramians in (E,phi) or (I,theta) order.
struct OrbitGramianPair {std::array<double,3> controllability{},observability{};std::array<double,2> hankel{};};
// sqrt(eig(Wc Wo)) for symmetric positive semidefinite [ee,ep,pp] pairs.
std::array<double,2> hankelFromSymmetric(const std::array<double,3>& wc,const std::array<double,3>& wo);
struct OrbitGramianInput {
    double inputVariance{},outputVariance{},frequencySlope{},horizon{};
    double inputCross{},inputPhase{},outputCross{},outputEnergy{};  // optional lower-order terms
};
// (c) Undamped Jordan shear over horizon tau: W_EE=tau B, W_Ephi=Omega' B tau^2/2,
// W_phiphi=Omega'^2 B tau^3/3, plus lower-order cross and phase-channel terms.
OrbitGramianPair undampedOrbitGramians(const OrbitGramianInput&);
// Collocated law sigma=B |Omega_I| tau^2 (1/sqrt3 +/- 1/2), B=<W_theta^2>, larger first.
std::array<double,2> collocatedOrbitHankel(double B,double actionTwist,double horizon);
struct DampedGramianInput {
    std::vector<double> energy,inputVariance,outputVariance,damping,frequency;
    double dampingScale{},initialEnergy{},targetEnergy{},finalEnergy{};
};
struct DampedGramianResult {
    double cEE{},cEphi{},cPhiphi{},oPhiphi{},oEphi{},oEE{};
    double targetDamping{},targetFrequency{};std::array<double,2> hankel{};
};
// (d) Slow-passage quadratures between E0>Et>Ef on tabulated B,C,D,Omega. Panels are
// clipped at the limits; three-point Gauss on local cubic interpolants. D>0 is required.
DampedGramianResult dampedOrbitGramians(const DampedGramianInput&);
// (e) Collocated Hamiltonian forcing H+uW, W=sum W_n(I) e^{i n theta}.
struct JordanMode {int harmonic{};std::complex<double> amplitude,slope;};
// R(t)=sum_n [i n d_I|W_n|^2 - n^2 Omega_I t |W_n|^2] e^{i n Omega t}. Requires a real
// observable: W_{-n}=conj(W_n) for every supplied harmonic.
ComplexVector jordanImpulseResponse(const std::vector<JordanMode>&,double frequency,double actionTwist,
    const std::vector<double>& times);
// Frozen-action residues i n W_n dW_{-n}/dI, keyed by harmonic n.
std::map<int,std::complex<double>> frozenKromResidues(const std::vector<JordanMode>&);
// L_g e^{i n theta}=i n e^{i n theta} sum_m dW_m/dI e^{i m theta}; keyed by n+m.
std::map<int,std::complex<double>> exactInputLift(int harmonic,const std::vector<JordanMode>&);
}

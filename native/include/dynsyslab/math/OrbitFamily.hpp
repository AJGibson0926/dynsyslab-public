#pragma once
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/math/OrbitIntegrator.hpp>
#include <array>
#include <functional>
#include <vector>
namespace dynsyslab::math {
using PlanarPoint=std::array<double,2>;
using PlanarField=std::function<PlanarPoint(const PlanarPoint&)>;
using PlanarScalar=std::function<double(const PlanarPoint&)>;
// Period annulus around a centre. Orbits start on the ray centre+r*sectionDirection, so
// phase zero is that crossing; choose a reversibility axis to get turning-point phase.
// Action uses p dq for (coordinate,momentum), or, when density is set, the invariant area
// form rho dx^dy through the Poincare-lemma primitive around the centre.
// Arnold, Mathematical Methods of Classical Mechanics, section 50.
struct PlanarOrbitSystem {
    PlanarField drift;PlanarScalar hamiltonian;
    PlanarPoint centre{0,0},sectionDirection{1,0};
    double maximumRadius{};int coordinate{0},momentum{1};PlanarScalar density;
};
struct OrbitFamilyOptions {
    FlowOptions flow;std::size_t phaseSamples{256};int maximumHarmonic{16};
    double closureTolerance{1e-7},derivativeStep{1e-4};
};
struct PeriodicOrbit {
    double energy{},period{},frequency{},action{};PlanarPoint start{};
    int orientation{};double closureError{},energyDrift{};std::size_t steps{};
};
struct OrbitFourier {FourierSeries first,second;};
struct OrbitFamilyTable {
    std::vector<PeriodicOrbit> orbits;
    std::vector<double> energies,periods,frequencies,actions,frequencySlopes;
};
struct FrequencyBranch {std::size_t first{},last{};int direction{};double lower{},upper{};};
struct FrequencyBranchReport {std::vector<FrequencyBranch> branches;bool monotone{};std::size_t flatSegments{};};
struct FrequencyRoot {double energy{},frequency{},slope{};std::size_t branch{};bool refined{};};
struct ShellRoot {int harmonic{};FrequencyRoot root;};
struct SeparatrixFit {double lambda{},constant{},rmsResidual{};std::size_t samples{};};
using FrequencyEvaluator=std::function<double(double)>;

PlanarPoint orbitStart(const PlanarOrbitSystem&,double energy);
PeriodicOrbit periodicOrbit(const PlanarOrbitSystem&,double energy,const OrbitFamilyOptions& ={});
std::vector<PlanarPoint> orbitSamples(const PlanarOrbitSystem&,const PeriodicOrbit&,const OrbitFamilyOptions& ={});
OrbitFourier orbitFourier(const std::vector<PlanarPoint>& samples,double frequency,int maximumHarmonic);
// Phase in [0,2pi) of a point on its own energy level, measured from the section start.
double orbitPhase(const PlanarOrbitSystem&,const PlanarPoint&,const OrbitFamilyOptions& ={});
// Energies strictly increasing. dOmega/dE uses fourth-order central differences inside
// the supplied energy interval and second-order one-sided differences at its ends.
OrbitFamilyTable orbitFamily(const PlanarOrbitSystem&,const std::vector<double>& energies,const OrbitFamilyOptions& ={});
FrequencyBranchReport frequencyBranches(const std::vector<double>& energies,const std::vector<double>& frequencies);
// All solutions of Omega(E)=target on strictly monotone table branches. With an evaluator
// roots are refined by the Illinois method; otherwise table interpolation is linear.
std::vector<FrequencyRoot> invertFrequency(const std::vector<double>& energies,const std::vector<double>& frequencies,
    double target,const FrequencyEvaluator& evaluator={});
// Roots of n*Omega(E)=w for n=1..maximumHarmonic.
std::vector<ShellRoot> shellLadder(const std::vector<double>& energies,const std::vector<double>& frequencies,
    double w,int maximumHarmonic,const FrequencyEvaluator& evaluator={});
// Least-squares T=(1/lambda) log(C/delta), delta=|separatrixEnergy-E|. lambda is the fit
// coefficient: an orbit passing k saddles of exponent mu has lambda=mu/k.
SeparatrixFit fitSeparatrixPeriod(const std::vector<double>& energies,const std::vector<double>& periods,double separatrixEnergy);
}

#pragma once
#include <dynsyslab/math/LaplaceTransform.hpp>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Pole/residue models of the Koopman resolvent R(s)=(sI-L)^{-1} on a fixed observable and
// initial state: Susuki, Mauroy & Mezic, SIAM J. Appl. Dyn. Syst. 20(4) (2021), eqs. (9), (17),
// (22), (25). A pole of order q carries a Laurent principal part, which adds the Jordan
// (nilpotent) terms of Kato, Perturbation Theory for Linear Operators, sec. III.6.5:
//   R(s)f = sum_j sum_q a_{j,q}/(s-lambda_j)^q <-> y(t)=sum a_{j,q} t^{q-1}/(q-1)! e^{lambda_j t}.
struct ResolventPole {
    std::complex<double> eigenvalue;
    ComplexVector laurent; // laurent[q-1] multiplies (s-lambda)^{-q}.
};
struct PoleResidueModel {std::vector<ResolventPole> poles;};
void validatePoleModel(const PoleResidueModel&);
std::complex<double> evaluateResolvent(const PoleResidueModel&,std::complex<double> s);
std::complex<double> evaluateTimeResponse(const PoleResidueModel&,double t);
// Largest Re(lambda) over poles with a nonzero coefficient: the abscissa of convergence of
// the transform of this observable, which can lie left of the principal exponent when a
// residue vanishes (for example an orbit on a strong stable manifold).
double modelAbscissa(const PoleResidueModel&);
// Cauchy residue a_q=(1/2 pi i) oint F(s)(s-c)^{q-1} ds on a circle by the trapezoidal rule,
// which converges geometrically for F analytic on the circle: Trefethen & Weideman, SIAM
// Rev. 56 (2014). The projection reading is Kato, sec. III.6.4, and
// Remark 7 of Susuki, Mauroy & Mezic (2021).
struct CauchyOptions {
    std::size_t maximumOrder{2},initialNodes{32},maximumNodes{4096};
    double tolerance{1e-11};
    double orderTolerance{1e-9};  // Relative size below which a Laurent coefficient is zero.
};
struct LaurentEstimate {
    std::complex<double> center;
    double radius{};
    ComplexVector laurent;      // Indexed like ResolventPole::laurent, up to maximumOrder.
    std::size_t order{},nodes{};
    double relativeChange{};    // Between the last two node doublings.
    double higherOrderRatio{};  // |a_{maximumOrder+1}|/max|a_q|; large means the order is too low.
    bool converged{false};
};
LaurentEstimate cauchyLaurent(const LaplaceFunction&,std::complex<double> center,double radius,const CauchyOptions& ={});

// Pole lattices of the three expansion theorems. Multi-indices k are in N0.
struct LatticePole {
    std::complex<double> value;
    std::vector<int> k,m;   // Exponent orders, and harmonics (empty for an equilibrium).
};
struct LatticeCoincidence {std::size_t first{},second{};double separation{};};
struct PoleLattice {
    std::vector<LatticePole> poles;            // Sorted by descending Re, then |Im|.
    std::vector<LatticeCoincidence> coincidences; // Distinct indices at one value: resonance.
    std::optional<double> smallestHarmonicDivisor; // min nonzero |m.Omega| (torus small divisors).
    std::string scope;
};
// Theorem 8, eq. (17): k.lambda, Re(lambda_j)<0, distinct and non-resonant.
PoleLattice equilibriumPoleLattice(const ComplexVector& lambda,int maximumOrder,double coincidenceTolerance=1e-10);
// Theorem 10, eq. (22): k.nu + i m Omega. A complex exponent is defined only modulo i Omega;
// the lattice set does not depend on that choice.
PoleLattice limitCyclePoleLattice(const ComplexVector& nu,double Omega,int maximumOrder,int maximumHarmonic,
    double coincidenceTolerance=1e-10);
// Theorem 11, eq. (25): k.nu + i m.Omega. For two or more rationally independent Omega the
// imaginary lattice is dense, so no pole on Re s=0 is isolated; this returns a finite truncation.
PoleLattice torusPoleLattice(const ComplexVector& nu,const std::vector<double>& Omega,int maximumOrder,
    int maximumHarmonic,double coincidenceTolerance=1e-10);

// ||(sI-A)^{-1}||_2 = 1/sigma_min(sI-A), and the epsilon-pseudospectrum grid:
// Trefethen & Embree, Spectra and Pseudospectra (2005).
double resolventNorm(const ComplexMatrix& A,std::complex<double> s);
struct PseudospectrumGrid {
    std::vector<double> realAxis,imagAxis;
    std::vector<double> sigmaMin;  // Row-major: imagAxis rows by realAxis columns.
    ComplexVector eigenvalues;
};
PseudospectrumGrid pseudospectrumGrid(const ComplexMatrix& A,std::vector<double> realAxis,std::vector<double> imagAxis);
}

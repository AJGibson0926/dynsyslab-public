#pragma once
#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <dynsyslab/math/OrbitGramian.hpp>
#include <dynsyslab/math/ShellPacket.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>
#include <map>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Laplace-domain view of orbit tables on a period annulus. On one orbit the observable
// f=sum_n g_n e^{i n phi} has poles at i n Omega (Susuki, Mauroy & Mezic 2021, eq. (11)
// restricted to a periodic orbit). Over a family of shells, the autocorrelation pairing
//   Y(s)=sum_n int rho(E) |g_n(E)|^2 / (s - i n Omega(E)) dE
// has the boundary value (Plemelj-Sokhotski; Muskhelishvili, Singular Integral Equations, ch. 2)
//   Y(iw+0)=sum_n int rho|g_n|^2 [pi delta(w-n Omega) - i PV 1/(w-n Omega)] dE,
// so each simple root of n Omega(E)=w contributes pi rho|g_n|^2/|n Omega'(E)|.
struct ShellResolventTable {
    std::vector<double> energy;      // Strictly increasing nodes in one period annulus.
    std::vector<double> frequency;   // Omega(E) > 0.
    std::vector<double> density;     // rho(E) >= 0, per unit energy.
    std::map<int,ComplexVector> harmonics; // g_n(E) at the nodes.
};
void validateShellTable(const ShellResolventTable&);
// From ShellPacket shells: weights are quadrature masses, converted to densities by the
// trapezoid cell widths. Harmonics missing on a shell are zero.
ShellResolventTable shellTableFromOrbitShells(const std::vector<OrbitShell>&);
// One-orbit pole model: Y(s)=sum_n g_n e^{i n phase}/(s-i n Omega).
PoleResidueModel orbitHarmonicModel(const FourierSeries&,double phase=0);
struct ShellBoundaryOptions {
    double slopeTolerance{1e-10};    // |n Omega'| at or below this is a degenerate (non-simple) root.
    std::vector<double> damping;     // Optional D(E) at the nodes, for the stationary-phase weight.
    double dampingScale{0};
};
struct ShellBoundaryRoot {
    int harmonic{};
    double energy{},frequencySlope{},density{},weight{};  // weight=rho|g_n|^2/|n Omega'|.
    std::optional<double> stationaryPhaseWeight;          // DampedOrbit::shellWeight at the root.
};
struct ShellBoundaryValue {
    double frequency{};
    std::complex<double> value;  // absorptive + i dispersive.
    double absorptive{},dispersive{};
    std::vector<ShellBoundaryRoot> roots;
    double zeroHarmonicMass{};   // Mass of g_0: an atom at w=0, excluded from value when w=0.
    bool valid{true};            // False on a degenerate root or a root at a table edge.
    std::vector<std::string> warnings;
};
// Cubic Hermite interpolation of Omega and of rho|g_n|^2 on the table; the principal value uses
// singularity subtraction with Gauss-Legendre panels split at each root.
ShellBoundaryValue shellBoundaryValue(const ShellResolventTable&,double frequency,const ShellBoundaryOptions& ={});
// Laplace transform of R(t)=sum_n [i n d_I|W_n|^2 - n^2 Omega_I t |W_n|^2] e^{i n Omega t}
// (OrbitGramian::jordanImpulseResponse): second-order poles at i n Omega,
//   G(s)=sum_n [ i n d_I|W_n|^2/(s-i n Omega) - n^2 Omega_I |W_n|^2/(s-i n Omega)^2 ].
PoleResidueModel jordanOrbitModel(const std::vector<JordanMode>&,double frequency,double actionTwist);
// Damped-limit pairing Y_d(s)=sum_n int rho|g_n|^2/(s-i n Omega(E)+d gamma(E)) dE, Re s>=0, d>0.
// Each table panel is integrated in closed form for piecewise-linear Omega, rho|g|^2 and gamma,
// so the quadrature stays exact as d->0, where Y_d(iw) tends to shellBoundaryValue(w).
std::complex<double> dampedShellPairing(const ShellResolventTable&,std::complex<double> s,
    const std::vector<double>& dampingRate,double dampingScale);
// The Rigged DMD rational kernel as a resolvent combination of U=exp(L h): with
// z_j=e^{i theta - i eps p_j}, y_j=e^{i theta - i eps conj(p_j)} and residues a_j,
//   K = (1/4 pi) sum_j [ a_j (z_j+U)(z_j-U)^{-1} - conj(a_j) (y_j+U)(y_j-U)^{-1} ]
//     = constant I + sum_terms c_k z_k (z_k I-U)^{-1}.
// Colbrook, Drysdale & Horning, Rigged DMD, arXiv:2405.00782. For |mu|=1 this equals
// RationalCircleKernel(theta-arg mu).
struct RationalResolventTerm {std::complex<double> node,coefficient;};
struct RationalResolventCombination {std::complex<double> constant;std::vector<RationalResolventTerm> terms;};
RationalResolventCombination rationalKernelResolvent(const RationalCircleKernel&,double theta);
std::complex<double> applyRationalCombination(const RationalResolventCombination&,std::complex<double> mu);
// Applies the kernel to a pole model through g(U)=g(mu)I+g'(mu)(U-mu) on each Jordan block, for
// poles of order at most 2, with U=exp(L h) and U-mu=mu h N on a block with nilpotent part N.
PoleResidueModel rationalProjectedModel(const PoleResidueModel&,double sampleInterval,double theta,
    const RationalCircleKernel&);
}

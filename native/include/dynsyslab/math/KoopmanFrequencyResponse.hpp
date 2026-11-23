#pragma once
#include <dynsyslab/math/KoopmanResolvent.hpp>
#include <map>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Harmonic frequency response of a lifted bilinear model
//   z' = A z + (B u) + sum_j u_j N_j z,   y = C z,
// under u = d u0 e^{i w t}. Susuki, Katayama, Mauroy & Mezic, On Koopman Resolvents and Frequency
// Response of Nonlinear Systems, arXiv:2603.05771v2 (2026), Definition 4 and Theorem 1: H_n is the
// residue at i n w of the skew-product resolvent divided by u0^n. For this lift,
//   Z_1=(iwI-A)^{-1} B d,  Z_{n+1}=((n+1)iwI-A)^{-1} N_d Z_n,  H_n=C Z_n.
// The classical single-channel G(s) stays in core/TransferFunction.hpp and is not duplicated.
struct BilinearLiftedModel {
    ComplexMatrix A,B,C;
    std::vector<ComplexMatrix> N;  // One n-by-n matrix per input, or empty for a linear model.
};
void validateBilinearModel(const BilinearLiftedModel&);
struct HarmonicResponseOptions {
    int maximumOrder{4};
    ComplexVector direction;           // Input direction d; empty selects the first input.
    double singularityTolerance{1e-12}; // Relative sigma_min of (i n w I - A) that counts as resonant.
};
struct HarmonicResponse {
    double frequency{};
    std::vector<ComplexVector> outputs,states;  // H_n and Z_n; index n-1.
    std::vector<double> sigmaMin;               // sigma_min(i n w I - A) for each order.
    bool hurwitz{false};                        // A is Hurwitz, so a steady state attracts at u0->0.
    bool valid{false};                          // Hurwitz and no order is resonant.
    std::optional<double> seriesRadius;         // Ratio estimate of the u0 convergence radius.
    std::vector<std::string> warnings;
};
HarmonicResponse harmonicResponse(const BilinearLiftedModel&,double frequency,const HarmonicResponseOptions& ={});

// Symmetric Volterra transfer function of the bilinear model: Rugh, Nonlinear System Theory:
// The Volterra/Wiener Approach (1981), ch. 3,
//   H_n^sym(s_1..s_n) = (1/n!) sum_perm C (S_n I-A)^{-1} N (S_{n-1}I-A)^{-1} ... N (S_1 I-A)^{-1} B,
// with S_k the partial sums along the permutation. Evaluated by dynamic programming over
// subsets of {1..n}. A response to u=sum_k a_k e^{s_k t} has order-n part
// sum a_{k1}..a_{kn} H_n^sym(s_{k1},..,s_{kn}) e^{(s_{k1}+..+s_{kn})t}; so
// H_n^sym(iw,..,iw) equals the harmonic response above.
// channels[j] picks the input of argument j; empty means input 0 everywhere. n <= 8.
ComplexVector volterraKernel(const BilinearLiftedModel&,const ComplexVector& s,const std::vector<std::size_t>& channels={});

// Real forcing u=d a cos(wt) mixes harmonics: every Z_k couples to Z_{k-1} and Z_{k+1}. This solves
//   (i k w I-A) Z_k = (a/2)[B_d (delta_{k,1}+delta_{k,-1}) + N_d (Z_{k-1}+Z_{k+1})]
// on |k|<=K (harmonic balance; Hayashi, Nonlinear Oscillations in Physical Systems (1964)).
struct RealHarmonicBalance {
    int maximumHarmonic{};
    std::map<int,ComplexVector> states,outputs;
    double truncationResidual{}; // ||(a/2) N_d Z_{+-K}||, the term dropped at the truncation edge.
};
RealHarmonicBalance realHarmonicBalance(const BilinearLiftedModel&,double frequency,double amplitude,int maximumHarmonic,
    const ComplexVector& direction={});

// Floquet multipliers of z'=(A+N_d u(t))z over one forcing period (RK4). The periodic response
// exists and attracts for spectral radius < 1. Radius >= 1 marks a bifurcation or loss of the
// response at this amplitude, which invalidates H_n (arXiv:2603.05771v2, Remark 3).
struct BilinearFloquet {ComplexVector multipliers;double spectralRadius{};bool stable{false};};
BilinearFloquet bilinearFloquet(const BilinearLiftedModel&,double frequency,double amplitude,bool complexInput,
    const ComplexVector& direction={},std::size_t stepsPerPeriod=1024);

// Theorem 1 needs a pole of order one at i n w. Contour estimate of the Laurent part of F.
struct SimplePoleCheck {std::complex<double> residue,secondOrderCoefficient;bool simple{false},converged{false};};
SimplePoleCheck simplePoleCheck(const LaplaceFunction&,std::complex<double> pole,double radius,double tolerance=1e-8);
}

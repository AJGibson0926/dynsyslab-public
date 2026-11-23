#pragma once
#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
namespace dynsyslab::math {
// Frequency-domain Gramians of a stable complex realization (A,B,C):
//   W_c=(1/2pi) int R(iw) B B^H R(iw)^H dw,  W_o=(1/2pi) int R(iw)^H C^H C R(iw) dw,
// and ||G||_2^2=(1/2pi) int ||C R(iw) B||_F^2 dw=tr(C W_c C^H) (Parseval). Zhou, Doyle & Glover,
// Robust and Optimal Control (1996), ch. 4; Antoulas, Approximation of Large-Scale Dynamical
// Systems (2005), ch. 4. The map w=c tan(phi) with Gauss-Legendre panels graded toward each
// pole frequency handles lightly damped poles; panel doubling gives the error estimate.
struct FrequencyGramianOptions {
    double tolerance{1e-9};
    std::size_t maximumLevels{7};
};
struct FrequencyGramianResult {
    ComplexMatrix controllability,observability;
    double h2Norm{},relativeChange{};
    std::size_t evaluations{};
    bool converged{false};
};
FrequencyGramianResult frequencyGramians(const ComplexMatrix& A,const ComplexMatrix& B,const ComplexMatrix& C,
    const FrequencyGramianOptions& ={});
// sup_w sigma_max(C(iwI-A)^{-1}B) over all real w. Reuses core::hInfinityNorm
// (Boyd & Balakrishnan 1990; Bruinsma & Steinbuch 1990) on the real embedding [[Re,-Im],[Im,Re]],
// whose response is unitarily similar to diag(G(iw),conj(G(-iw))); the sign of the peak
// frequency is then resolved on the complex model. Needs the optional numerics library
// (DYNSYSLAB_MATH_NUMERICS); otherwise it throws std::runtime_error.
struct ResolventPeak {double value{},frequency{};std::size_t iterations{};};
ResolventPeak resolventHInfinity(const ComplexMatrix& A,const ComplexMatrix& B,const ComplexMatrix& C);
}

#pragma once
#include <dynsyslab/math/ResolventLinearAlgebra.hpp>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Pseudo-arclength continuation (Keller, in Applications of Bifurcation Theory, 1977) of periodic
// responses of x'=F(x,u,p), u=a cos(w t), by harmonic balance with alternating frequency-time
// evaluation (Nakhla & Vlach, IEEE Trans. Circuits Syst. 23 (1976) 85-91; Cameron & Griffin,
// J. Appl. Mech. 56 (1989) 149-154). Folds are sign changes of d(parameter)/d(arclength).
// Stability is by Floquet multipliers of the variational equation. A failed corrector, a
// non-finite field or a tangent jump ends the run with branchLost=true; it never continues.
using ForcedField=std::function<void(const double* x,double u,double p,double* dx)>;
enum class ContinuationVariable {Frequency,Amplitude,Parameter};
struct ContinuationProblem {
    std::size_t dimension{};
    ForcedField field;
    double frequency{1},amplitude{0},parameter{0};
    ContinuationVariable variable{ContinuationVariable::Frequency};
    double lower{},upper{};           // Continuation stops at these bounds.
    int harmonics{5};
    std::vector<double> initialCoefficients; // d(2K+1) as [a0,a1,b1,...] per state; empty is zero.
    double step{.01},minimumStep{1e-6},maximumStep{.05};
    std::size_t maximumSteps{1000};
    int direction{1};                 // Initial sign of the parameter change.
    std::size_t outputState{0};
    double tolerance{1e-9};
    bool computeStability{true};
};
struct ContinuationPoint {
    double value{},residual{};
    std::vector<double> coefficients;
    std::map<int,std::complex<double>> response; // c_n of e^{i n w t} for the output state, n>=0.
    std::optional<double> floquetRadius;
    std::optional<bool> stable;
};
struct ContinuationFold {
    double value{};                    // Quadratic estimate from the bracketing tangents.
    std::size_t afterPoint{};          // The fold lies between this point and the next.
    double responseMagnitude{};        // |c_1| interpolated at the fold.
};
struct ContinuationResult {
    std::vector<ContinuationPoint> points;
    std::vector<ContinuationFold> folds;
    std::string termination;           // parameter_bound, step_limit, corrector_failed, non_finite_field, tangent_jump, singular_tangent.
    bool branchLost{false};
    std::vector<std::string> warnings;
};
ContinuationResult continueHarmonicBalance(const ContinuationProblem&);
}

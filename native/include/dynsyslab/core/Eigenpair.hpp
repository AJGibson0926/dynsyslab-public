#pragma once
// One Koopman eigenpair: an eigenvalue, its eigenfunction, and optionally the residual of the pair.
//
// The residual is a property of the pair, not of the eigenvalue or the eigenfunction alone: it measures
// how well this function and this value satisfy the eigen-equation together. A method stores it under
// its own definition and scope; definitions are never mixed or converted.
//   "one-step"   ||psi(y) - mu psi(x)|| / ||psi(x)|| over snapshot pairs (x, y = F(x)), as a weighted
//                quadrature over the data: ResDMD (Colbrook & Townsend, Comm. Pure Appl. Math. 77,
//                2024), in the EDMD dictionary of Williams, Kevrekidis & Rowley (J. Nonlinear Sci. 25,
//                2015).
//   "generator"  the same idea in continuous time, RMS_k |grad psi(x_k) . f(x_k) - lambda psi(x_k)|
//                over the training rows (absolute, in the units of psi): the generator equation SINDy
//                solves for eigenfunctions (Kaiser, Kutz & Brunton, Mach. Learn.: Sci. Technol. 2,
//                035023, 2021).
//   "one-step-successor-mean"  the one-step residual normalized by ||psi(y)|| and averaged over
//                training runs: the pooled history EDMD's mode score (HistoryKoopman.hpp).
// Methods that compute none (H-DMD, DMD, EDMD, hand-given expansions, ...) leave it empty.
#include <dynsyslab/core/Eigenfunction.hpp>
#include <dynsyslab/core/Eigenvalue.hpp>

#include <optional>
#include <string>

namespace dynsyslab::core {

namespace residual_definition {
inline constexpr const char* oneStep = "one-step";
inline constexpr const char* generator = "generator";
inline constexpr const char* oneStepSuccessorMean = "one-step-successor-mean";
}  // namespace residual_definition

// Where a residual was evaluated: the training data of the fit, or data held out from it.
namespace residual_scope {
inline constexpr const char* empiricalTraining = "empirical-training";
inline constexpr const char* empiricalHeldOut = "empirical-held-out";
}  // namespace residual_scope

struct EigenpairResidual {
    double value{};
    std::string definition;  // residual_definition
    std::string scope;       // residual_scope
};

// An eigenvalue with its eigenfunction (when the holder has one: a Spectrum alone does not) and its
// residual (when the method computed one). KoopmanModel::eigenpair(j) assembles all three; the
// eigenfunction shares the model's representations rather than copying them.
struct Eigenpair {
    Eigenvalue eigenvalue;
    std::optional<Eigenfunction> eigenfunction;
    std::optional<EigenpairResidual> residual;
};

}  // namespace dynsyslab::core

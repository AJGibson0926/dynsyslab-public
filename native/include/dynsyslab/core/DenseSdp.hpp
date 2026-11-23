#pragma once
// Small dense conic solver for linear matrix inequalities, plus an independent
// a-posteriori verifier.
//
//   minimize c'x   subject to   F_b(x) = F_b0 + sum_i x_i F_bi >= 0  (PSD) for every block b
//
// Algorithm: operator splitting on the homogeneous self-dual embedding,
//   B. O'Donoghue, E. Chu, N. Parikh, S. Boyd, "Conic Optimization via Operator
//   Splitting and Homogeneous Self-Dual Embedding", J. Optim. Theory Appl. 169
//   (2016) 1042-1068 (the SCS method). This is a compact dense reimplementation,
//   not the SCS library. The embedding follows Ye, Todd & Mizuno, Math. Oper.
//   Res. 19 (1994) 53-67.
//
// Solver output is never a certificate. Use verifyLmiCertificate, which rebuilds
// every F_b(x) from the problem data and checks its minimum eigenvalue. That
// check is floating-point with a rounding allowance; it is NOT interval
// arithmetic and is not a machine-checked proof.
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

// All matrices are row-major, dimension x dimension, and symmetric.
struct LmiBlock {
    std::size_t dimension{};
    std::vector<double> constant;                  // F_0
    std::vector<std::vector<double>> coefficients; // F_i, one per variable; empty means zero
};

struct SdpProblem {
    std::size_t variables{};
    std::vector<double> objective;  // size variables; empty means pure feasibility
    std::vector<LmiBlock> blocks;
};

enum class SdpStatus { converged, primalInfeasible, dualInfeasible, iterationLimit };
[[nodiscard]] const char* sdpStatusName(SdpStatus status);

struct SdpOptions {
    std::size_t maxIterations{200000};
    std::size_t checkInterval{25};
    double tolerance{1e-7};              // relative primal, dual and gap residuals
    double infeasibilityTolerance{1e-7}; // normalized Farkas residual
    double relaxation{1.5};              // over-relaxation in (0, 2)
    bool equilibrate{true};
};

struct SdpResult {
    SdpStatus status{SdpStatus::iterationLimit};
    std::vector<double> x;                  // empty when no finite iterate exists
    std::vector<std::vector<double>> dual;  // one row-major multiplier per block
    double objective{};
    double primalResidual{};  // ||A x + s - b|| / (1 + ||b||)
    double dualResidual{};    // ||A'y + c|| / (1 + ||c||)
    double gap{};             // |c'x + b'y| / (1 + |c'x| + |b'y|)
    std::size_t iterations{};
};

[[nodiscard]] SdpResult solveDenseSdp(const SdpProblem& problem, const SdpOptions& options = {});

struct LmiCertificate {
    bool certified{false};
    double margin{};
    double worstMinEigenvalue{};              // smallest over blocks, before the allowance
    double roundingAllowance{};               // largest per-block allowance applied
    std::vector<double> blockMinEigenvalues;
    std::string reason;                       // empty when certified
    std::string scope;
};

// certified == true only if every reconstructed block is symmetric, finite, and
// has lambda_min(F_b(x)) - allowance >= margin. Independent of the solver path.
[[nodiscard]] LmiCertificate verifyLmiCertificate(const SdpProblem& problem,
    const std::vector<double>& x, double margin);

}  // namespace dynsyslab::core

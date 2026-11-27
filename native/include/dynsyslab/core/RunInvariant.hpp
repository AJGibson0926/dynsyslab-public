#pragma once
// A nonconstant conserved quantity (a Koopman eigenfunction at eigenvalue 0, multiplier 1; Mezic
// 2005, Nonlinear Dyn. 41:309-325) learned from run identity alone: the polynomial whose mean differs
// most between unforced runs while it varies least within each (the discriminant of Fisher 1936,
// Ann. Eugenics 7:179-188). With feature vectors theta (constant excluded), between-run scatter S_b
// and within-run scatter S_w, the coefficients maximize c^T S_b c / c^T S_w c: the leading
// generalized eigenvector of S_b c = rho (S_w + epsilon I) c, solved through the Cholesky factor of
// S_w + epsilon I (Golub and Van Loan, Matrix Computations, 4th ed. 2013, section 8.7). The constant
// term then centers the invariant on the pooled mean. Ported from tools/echo_sampled_koopman_diagnostics.py
// (fit_action_invariant, invariant_metrics); docs/architecture/TRAJECTORY_KERNELS.md.
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct RunInvariantRun {
    double mean{}, standardDeviation{}, peakToPeak{};
    double sampledInvariance{};  // l2 of I(x_{k+s}) - I(x_k) over l2 of I(x_{k+s}) along the run
};

struct RunInvariant {
    std::size_t dimension{}, degree{};
    std::vector<double> coefficients;  // over polynomialExponentTable(dimension, degree, true)
    double separation{};               // the generalized eigenvalue rho (between over within)
    std::vector<RunInvariantRun> runs;
    [[nodiscard]] double evaluate(const double* x) const;
};

[[nodiscard]] RunInvariant fitRunInvariant(const std::vector<TimeSeries>& runs, std::size_t degree = 3, std::size_t pairStride = 1);

}  // namespace dynsyslab::core

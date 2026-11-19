#pragma once
// KLQE process covariance W from an estimator spec.
//
//   "process_covariance": 1e4                      scalar * I (default 1)
//   "process_covariance": [w1, ..., wn] | n*n list diagonal or full
//   "process_covariance": {"from_records": [solution archives], "residual_horizon": h,
//                          "window": [t0, t1], "scale": 1, "floor": 1e-6, "add": 0}
//       W = scale * E[r r^T] / h + add * I, with the K-ROM's own prediction residual
//       r_k = psi(x(t_k + h)) - exp(A h) psi(x(t_k)) - Gamma(h) u(t_k)
//       on the records' states and inputs (zero-order hold); eigenvalues below
//       floor * max are raised to it so W stays positive definite. W then trusts the
//       K-ROM dynamics where they predict well and distrusts them where they do not.
//   "process_covariance_input_weight": q2          adds q2 * B B^T (constant map): model
//       error entering through the input channel (the estimation form of loop
//       transfer recovery).
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <vector>

namespace dynsyslab::workflow {

struct KlqeProcessCovariance {
    std::vector<double> matrix;  // split dimension squared, row-major
    core::json::Value evidence;  // how W was formed (null for a plain scalar/list)
};

[[nodiscard]] KlqeProcessCovariance klqeProcessCovariance(const core::json::Value& estimatorSpec,
                                                          const core::KoopmanModel& model);

}  // namespace dynsyslab::workflow

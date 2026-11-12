#pragma once
#include <cstddef>

// Autonomous RHS on one contiguous state. Return zero on success.
// Status: 0 complete, 1 state guard, 2 RHS/nonfinite failure, 3 substep budget,
// 4 invalid request, 5 allocation failure. Only the first *written rows are valid.
// A null limits vector disables the componentwise magnitude guard.
extern "C" {
using mc_rhs_v1 = int (*)(void*, const double*, double*);
int mc_rk4_v1(mc_rhs_v1 rhs, void* context, std::size_t dimension,
              const double* initial, const double* time, std::size_t count,
              double max_step, std::size_t max_steps, const double* limits,
              double* output, std::size_t* written, std::size_t* evaluations);
}

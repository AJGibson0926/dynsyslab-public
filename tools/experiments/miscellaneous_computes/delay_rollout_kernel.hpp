#pragma once
#include <cstddef>
// Evaluates the current state's derivative from [present, tau_1 delay, ..., tau_K delay] in row-major order.
using mc_delay_rhs_v1 = int (*)(void*, const double*, double*);
// History is uniformly sampled at the same step and ends at the initial state.
// Each positive lag must be at least one integration step. Output starts at the initial state.
// Status 0 complete, 1 magnitude guard, 2 RHS/nonfinite, 3 invalid input or allocation.
extern "C" int mc_delay_rk4_v1(mc_delay_rhs_v1 rhs, void* context, std::size_t dim,
    const double* lags, std::size_t taps, const double* history, std::size_t history_rows,
    double h, std::size_t output_rows, const double* limits, double* output, std::size_t* written);

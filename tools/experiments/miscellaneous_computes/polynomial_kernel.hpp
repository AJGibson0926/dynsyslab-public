#pragma once
#include <cstddef>
// Exponents are nonnegative, feature by input. Coefficients are output by
// feature. The factory owns a copy. eval is an autonomous mc_rhs_v1 for square
// fields and a single-row feature evaluator for rectangular readouts.
extern "C" {
void* mc_poly_new_v1(std::size_t n,std::size_t outputs,std::size_t features,
                     const int* powers,const double* coefficients);
int mc_poly_eval_v1(void* field,const double* state,double* result);
void mc_poly_drop_v1(void* field);
}

#ifndef DYNSYSLAB_LINALG_C_API_H
#define DYNSYSLAB_LINALG_C_API_H

#include <dynsyslab/c_api.h>

#if defined(_WIN32)
#  if defined(DYNSYSLAB_LINALG_BUILD)
#    define DSL_LINALG_API __declspec(dllexport)
#  else
#    define DSL_LINALG_API __declspec(dllimport)
#  endif
#else
#  define DSL_LINALG_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define DSL_LINALG_ABI_VERSION 1u

typedef struct dsl_complex64 {
    double real;
    double imag;
} dsl_complex64;

DSL_LINALG_API uint32_t dsl_linalg_abi_version(void);
DSL_LINALG_API dsl_status dsl_linalg_last_error(
    char* buffer, size_t capacity, size_t* required);

/* All matrices are contiguous row-major arrays. */
DSL_LINALG_API dsl_status dsl_linalg_svd_real(
    const double* matrix, size_t rows, size_t columns, size_t requested_rank,
    double* left, size_t left_capacity, double* singular_values,
    size_t singular_capacity, double* right_transpose,
    size_t right_capacity, size_t* output_rank);

DSL_LINALG_API dsl_status dsl_linalg_least_squares_real(
    const double* matrix, size_t rows, size_t columns,
    const double* right_hand_side, size_t outputs, double ridge,
    double* solution, size_t solution_capacity);

DSL_LINALG_API dsl_status dsl_linalg_least_squares_complex(
    const dsl_complex64* matrix, size_t rows, size_t columns,
    const dsl_complex64* right_hand_side, size_t outputs, double ridge,
    dsl_complex64* solution, size_t solution_capacity);

DSL_LINALG_API dsl_status dsl_linalg_eigen_real(
    const double* matrix, size_t dimension,
    dsl_complex64* eigenvalues, size_t eigenvalue_capacity,
    dsl_complex64* right_eigenvectors, size_t eigenvector_capacity);

#ifdef __cplusplus
}
#endif

#endif

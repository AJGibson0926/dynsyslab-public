#ifndef DYNSYSLAB_EXTENSION_API_H
#define DYNSYSLAB_EXTENSION_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DSL_EXTENSION_ABI_VERSION 1u
#define DSL_EXTENSION_ENTRY_NAME "dynsyslab_extension_v1"

typedef void (*dsl_extension_output_writer)(
    const char* data, size_t size, void* user_data);

typedef struct dsl_system_extension_v1 {
    uint32_t struct_size;
    const char* descriptor_json;
    void* (*create)(const char* configuration_json, char* error, size_t error_capacity);
    void (*destroy)(void* instance);
    const char* (*resolved_configuration_json)(const void* instance);
    void (*reset)(void* instance);
    int (*rhs)(const void* instance, double time, const double* state,
               size_t state_count, double input, double* derivative,
               size_t derivative_count, char* error, size_t error_capacity);
    int (*jacobian)(const void* instance, double time, const double* state,
                    size_t state_count, double* matrix, size_t matrix_count,
                    char* error, size_t error_capacity);
    int (*input_matrix)(const void* instance, const double* state,
                        size_t state_count, double* matrix, size_t matrix_count,
                        char* error, size_t error_capacity);
    int (*invalid)(const void* instance, double time, const double* state,
                   size_t state_count);
    size_t (*observable_count)(const void* instance);
    const char* (*observable_descriptor_json)(const void* instance, size_t index);
    int (*observable_evaluate)(const void* instance, size_t index,
                               const double* state, size_t state_count,
                               double* values, size_t value_count,
                               char* error, size_t error_capacity);
    int (*observable_jacobian)(const void* instance, size_t index,
                               const double* state, size_t state_count,
                               double* matrix, size_t matrix_count,
                               char* error, size_t error_capacity);
} dsl_system_extension_v1;

/* Optional non-ODE model families (docs/architecture/NATIVE_SOLVER.md). Each table carries an
 * opaque user_data pointer passed back to its callbacks. Callbacks return 0 on success.
 *
 * Delay model: x'(t) = F(x(t), x(t - tau_1), ..., x(t - tau_K)). descriptor_json:
 *   {"identifier","label","state_dimension","lags":[0,tau_1,...],"state_names":[...]}
 * rhs receives the (K + 1) x state_dimension row-major window, present state first. */
typedef struct dsl_delay_extension_v1 {
    uint32_t struct_size;
    const char* descriptor_json;
    void* user_data;
    int (*rhs)(void* user_data, const double* window, size_t window_count,
               double* derivative, size_t state_count, char* error, size_t error_capacity);
} dsl_delay_extension_v1;

/* Field (PDE) model, semi-discrete on a structured mesh (method of lines). descriptor_json:
 *   {"identifier","label","mesh":{"shape":[..],"spacing":[..],"origin":[..],
 *    "boundary":["periodic"|"edge"|"zero",...]},"field_shape":[..]} */
typedef struct dsl_field_extension_v1 {
    uint32_t struct_size;
    const char* descriptor_json;
    void* user_data;
    int (*initial_field)(void* user_data, const char* configuration_json, double* field,
                         size_t value_count, char* error, size_t error_capacity);
    int (*rhs)(void* user_data, double time, const double* field, size_t value_count,
               double* derivative, char* error, size_t error_capacity);
} dsl_field_extension_v1;

/* Lattice (Monte Carlo) model on integer site states. descriptor_json:
 *   {"identifier","label","sites":N,"observables":["magnetization",...]}
 * sweep advances one sweep in place and writes one value per observable. */
typedef struct dsl_lattice_extension_v1 {
    uint32_t struct_size;
    const char* descriptor_json;
    void* user_data;
    int (*initial)(void* user_data, uint64_t seed, int32_t* sites, size_t site_count,
                   char* error, size_t error_capacity);
    int (*sweep)(void* user_data, uint64_t seed, uint64_t sweep_index, int32_t* sites,
                 size_t site_count, double* observables, size_t observable_count,
                 char* error, size_t error_capacity);
} dsl_lattice_extension_v1;

typedef struct dsl_extension_v1 {
    uint32_t struct_size;
    uint32_t abi_version;
    const char* identifier;
    const char* version;
    const char* capabilities_json;
    size_t system_count;
    const dsl_system_extension_v1* systems;
    int (*execute_json)(const char* capability, const char* request_json,
                        dsl_extension_output_writer writer, void* writer_data,
                        char* error, size_t error_capacity);
    /* Optional families; read only when struct_size covers them (older v1 extensions end
     * at execute_json). */
    size_t delay_count;
    const dsl_delay_extension_v1* delays;
    size_t field_count;
    const dsl_field_extension_v1* fields;
    size_t lattice_count;
    const dsl_lattice_extension_v1* lattices;
} dsl_extension_v1;

#define DSL_EXTENSION_V1_CORE_SIZE (offsetof(dsl_extension_v1, delay_count))

typedef const dsl_extension_v1* (*dsl_extension_entry_v1)(void);

#if defined(_WIN32)
#  define DSL_EXTENSION_EXPORT __declspec(dllexport)
#else
#  define DSL_EXTENSION_EXPORT
#endif

#ifdef __cplusplus
}
#endif

#endif

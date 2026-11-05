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
} dsl_extension_v1;

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

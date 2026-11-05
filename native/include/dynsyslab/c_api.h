#ifndef DYNSYSLAB_C_API_H
#define DYNSYSLAB_C_API_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#  if defined(DYNSYSLAB_C_API_BUILD)
#    define DSL_API __declspec(dllexport)
#  else
#    define DSL_API __declspec(dllimport)
#  endif
#else
#  define DSL_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define DSL_ABI_VERSION 1u

typedef struct dsl_system_handle dsl_system_handle;

typedef enum dsl_status {
    DSL_STATUS_OK = 0,
    DSL_STATUS_INVALID_ARGUMENT = 1,
    DSL_STATUS_BUFFER_TOO_SMALL = 2,
    DSL_STATUS_NOT_FOUND = 3,
    DSL_STATUS_CANCELLED = 4,
    DSL_STATUS_RUNTIME_ERROR = 5,
    DSL_STATUS_RESOURCE_LIMIT = 6
} dsl_status;

DSL_API uint32_t dsl_abi_version(void);
DSL_API const char* dsl_library_version(void);

/* UTF-8 JSON/string outputs use caller-owned buffers. required includes the NUL. */
DSL_API dsl_status dsl_last_error(char* buffer, size_t capacity, size_t* required);
DSL_API dsl_status dsl_system_catalog_json(char* buffer, size_t capacity, size_t* required);

DSL_API dsl_status dsl_system_create_json(
    const char* configuration_json, dsl_system_handle** output);
DSL_API void dsl_system_destroy(dsl_system_handle* system);
DSL_API dsl_status dsl_system_dimension(
    const dsl_system_handle* system, size_t* dimension);
DSL_API dsl_status dsl_system_configuration_json(
    const dsl_system_handle* system, char* buffer, size_t capacity, size_t* required);
DSL_API dsl_status dsl_system_rhs(
    const dsl_system_handle* system, double time, const double* state,
    size_t state_count, double input, double* derivative, size_t derivative_count);

#ifdef __cplusplus
}
#endif

#endif

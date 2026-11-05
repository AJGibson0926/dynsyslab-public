#ifndef DYNSYSLAB_WORKFLOW_C_API_H
#define DYNSYSLAB_WORKFLOW_C_API_H

#include <dynsyslab/c_api.h>

#if defined(_WIN32)
#  if defined(DYNSYSLAB_WORKFLOW_C_API_BUILD)
#    define DSL_WORKFLOW_API __declspec(dllexport)
#  else
#    define DSL_WORKFLOW_API __declspec(dllimport)
#  endif
#else
#  define DSL_WORKFLOW_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define DSL_WORKFLOW_ABI_VERSION 4u

typedef struct dsl_workflow_result_handle dsl_workflow_result_handle;
typedef int (*dsl_cancellation_callback)(void* user_data);
typedef void (*dsl_progress_callback)(
    const char* stage, const char* message, double fraction,
    size_t completed, size_t planned, void* user_data);

typedef struct dsl_execution_callbacks {
    dsl_cancellation_callback cancellation;
    dsl_progress_callback progress;
    void* user_data;
    size_t maximum_threads;
    size_t maximum_memory_bytes;
} dsl_execution_callbacks;

DSL_WORKFLOW_API uint32_t dsl_workflow_abi_version(void);
/* Evaluate native observable values, Jacobians and autonomous flow derivatives. */
/* Explicit optional sidecar. Query mode lookup performs no writes or computation. */
DSL_WORKFLOW_API dsl_status dsl_workflow_phase_space_json(
    const char* query_json, const char* archive_root_utf8,
    const dsl_execution_callbacks* callbacks, dsl_workflow_result_handle** output);
/* List or add a system paper using dynsyslab.research-query/v1. */
DSL_WORKFLOW_API dsl_status dsl_workflow_research_json(
    const char* query_json, const char* archive_root_utf8, dsl_workflow_result_handle** output);
/* Read-only system-owned analytic and special-function queries (v1 schemas). */
DSL_WORKFLOW_API dsl_status dsl_workflow_analytic_json(
    const char* query_json, dsl_workflow_result_handle** output);
DSL_WORKFLOW_API dsl_status dsl_workflow_special_function_json(
    const char* query_json, dsl_workflow_result_handle** output);
DSL_WORKFLOW_API dsl_status dsl_workflow_observable_json(
    const char* query_json, dsl_workflow_result_handle** output);
DSL_WORKFLOW_API dsl_status dsl_workflow_plan_json(
    const char* job_json, dsl_workflow_result_handle** output);
DSL_WORKFLOW_API dsl_status dsl_workflow_generate_initial_conditions_json(
    const char* solve_job_json, dsl_workflow_result_handle** output);
DSL_WORKFLOW_API dsl_status dsl_workflow_run_json(
    const char* job_json, const char* archive_root_utf8,
    const char* progress_path_utf8, const dsl_execution_callbacks* callbacks,
    dsl_workflow_result_handle** output);
DSL_WORKFLOW_API void dsl_workflow_result_destroy(dsl_workflow_result_handle* result);
DSL_WORKFLOW_API dsl_status dsl_workflow_result_json(
    const dsl_workflow_result_handle* result, char* buffer,
    size_t capacity, size_t* required);
DSL_WORKFLOW_API dsl_status dsl_workflow_last_error(
    char* buffer, size_t capacity, size_t* required);

#ifdef __cplusplus
}
#endif

#endif

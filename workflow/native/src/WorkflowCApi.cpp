#include <dynsyslab/workflow_c_api.h>

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/workflow/KoopmanSpectralQuery.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/LieAlgebraQuery.hpp>
#include <dynsyslab/workflow/KromAuditQuery.hpp>
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow/Solve.hpp>

#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

struct dsl_workflow_result_handle { std::string json; };

namespace {
thread_local std::string lastError;

dsl_status copyText(const std::string& value, char* buffer, size_t capacity,
                    size_t* required) {
    const size_t needed = value.size() + 1;
    if (required) *required = needed;
    if (!buffer || capacity < needed) return DSL_STATUS_BUFFER_TOO_SMALL;
    std::memcpy(buffer, value.c_str(), needed);
    return DSL_STATUS_OK;
}

dynsyslab::core::json::Value paths(const std::vector<std::filesystem::path>& values) {
    dynsyslab::core::json::Value::Array result;
    for (const auto& value : values) result.emplace_back(value.generic_u8string());
    return dynsyslab::core::json::Value(std::move(result));
}

dynsyslab::core::ExecutionContext context(const dsl_execution_callbacks* callbacks) {
    dynsyslab::core::ExecutionContext result;
    if (!callbacks) return result;
    result.limits.maximumThreads = callbacks->maximum_threads;
    result.limits.maximumMemoryBytes = callbacks->maximum_memory_bytes;
    if (callbacks->cancellation) {
        result.cancellation = [callbacks] {
            return callbacks->cancellation(callbacks->user_data) != 0;
        };
    }
    if (callbacks->progress) {
        result.progress = [callbacks](const dynsyslab::core::ExecutionProgress& event) {
            callbacks->progress(event.stage.c_str(), event.message.c_str(), event.fraction,
                                event.completed, event.planned, callbacks->user_data);
        };
    }
    return result;
}

dsl_status failed(dsl_status status, const std::string& message) {
    lastError = message;
    return status;
}

dsl_workflow_result_handle* resultHandle(dynsyslab::core::json::Value value) {
    auto result = std::make_unique<dsl_workflow_result_handle>();
    result->json = dynsyslab::core::json::dump(value);
    return result.release();
}
}  // namespace

extern "C" {

uint32_t dsl_workflow_abi_version(void) { return DSL_WORKFLOW_ABI_VERSION; }

dsl_status dsl_workflow_koopman_spectral_json(const char* query_json,dsl_workflow_result_handle** output) {
    if(!query_json||!output)return failed(DSL_STATUS_INVALID_ARGUMENT,"Koopman query argument is null");
    *output=nullptr;
    try {
        *output=resultHandle(dynsyslab::workflow::evaluateKoopmanSpectralQuery(dynsyslab::core::json::parse(query_json)));
        lastError.clear();return DSL_STATUS_OK;
    } catch(const std::invalid_argument& e) {return failed(DSL_STATUS_INVALID_ARGUMENT,e.what());}
      catch(const std::exception& e) {return failed(DSL_STATUS_RUNTIME_ERROR,e.what());}
      catch(...) {return failed(DSL_STATUS_RUNTIME_ERROR,"Unknown Koopman spectral query error");}
}

dsl_status dsl_workflow_observable_json(const char* query_json,
                                      dsl_workflow_result_handle** output) {
    if (!query_json || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "observable query argument is null");
    *output = nullptr;
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateObservableQuery(
            dynsyslab::core::json::parse(query_json)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR, "unknown observable query exception");
    }
}

dsl_status dsl_workflow_krom_audit_json(const char* query_json,
                                      dsl_workflow_result_handle** output) {
    if (!output) return failed(DSL_STATUS_INVALID_ARGUMENT, "KROM audit output is null");
    *output = nullptr;
    if (!query_json) return failed(DSL_STATUS_INVALID_ARGUMENT, "KROM audit query is null");
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateKromAuditQuery(
            dynsyslab::core::json::parse(query_json)));
        lastError.clear(); return DSL_STATUS_OK;
    } catch (const std::invalid_argument& e) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, e.what());
    } catch (const std::exception& e) {
        return failed(DSL_STATUS_RUNTIME_ERROR, e.what());
    } catch (...) { return failed(DSL_STATUS_RUNTIME_ERROR, "Unknown KROM audit exception"); }
}

dsl_status dsl_workflow_lie_algebra_json(const char* query_json,
                                      dsl_workflow_result_handle** output) {
    if (!output) return failed(DSL_STATUS_INVALID_ARGUMENT, "Lie algebra output is null");
    *output = nullptr;
    if (!query_json) return failed(DSL_STATUS_INVALID_ARGUMENT, "Lie algebra query is null");
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateLieAlgebraQuery(
            dynsyslab::core::json::parse(query_json)));
        lastError.clear(); return DSL_STATUS_OK;
    } catch (const std::invalid_argument& e) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, e.what());
    } catch (const std::exception& e) {
        return failed(DSL_STATUS_RUNTIME_ERROR, e.what());
    } catch (...) { return failed(DSL_STATUS_RUNTIME_ERROR, "Unknown Lie algebra exception"); }
}

dsl_status dsl_workflow_analytic_json(const char* query_json,
                                      dsl_workflow_result_handle** output) {
    if (!query_json || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "analytic query argument is null");
    *output = nullptr;
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateAnalyticQuery(
            dynsyslab::core::json::parse(query_json)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR, "unknown analytic query exception");
    }
}

dsl_status dsl_workflow_special_function_json(const char* query_json,
                                      dsl_workflow_result_handle** output) {
    if (!query_json || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "special_function query argument is null");
    *output = nullptr;
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateSpecialFunctionQuery(
            dynsyslab::core::json::parse(query_json)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR, "unknown special_function query exception");
    }
}

dsl_status dsl_workflow_research_json(const char* query_json, const char* archive_root_utf8,
    dsl_workflow_result_handle** output) {
    if (!query_json || !archive_root_utf8 || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "research argument is null");
    *output = nullptr;
    try {
        *output = resultHandle(dynsyslab::workflow::evaluateResearchCatalogueQuery(
            dynsyslab::core::json::parse(query_json),
            std::filesystem::path(dynsyslab::core::wide(archive_root_utf8))));
        lastError.clear(); return DSL_STATUS_OK;
    } catch (const std::invalid_argument& e) {return failed(DSL_STATUS_INVALID_ARGUMENT, e.what());}
    catch (const std::exception& e) {return failed(DSL_STATUS_RUNTIME_ERROR, e.what());}
    catch (...) {return failed(DSL_STATUS_RUNTIME_ERROR, "unknown research exception");}
}

dsl_status dsl_workflow_phase_space_json(const char* query_json,const char* archive_root_utf8,
    const dsl_execution_callbacks* callbacks,dsl_workflow_result_handle** output) {
    if(!query_json||!archive_root_utf8||!output)return failed(DSL_STATUS_INVALID_ARGUMENT,"phase-space argument is null");
    *output=nullptr;
    try {
        auto execution=context(callbacks);
        *output=resultHandle(dynsyslab::workflow::evaluatePhaseSpaceQuery(
            dynsyslab::core::json::parse(query_json),std::filesystem::path(dynsyslab::core::wide(archive_root_utf8)),&execution));
        lastError.clear();return DSL_STATUS_OK;
    } catch(const std::invalid_argument& e){return failed(DSL_STATUS_INVALID_ARGUMENT,e.what());}
    catch(const std::exception& e){return failed(DSL_STATUS_RUNTIME_ERROR,e.what());}
    catch(...){return failed(DSL_STATUS_RUNTIME_ERROR,"unknown phase-space exception");}
}

dsl_status dsl_workflow_plan_json(const char* job_json,
                                  dsl_workflow_result_handle** output) {
    if (!job_json || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "workflow plan argument is null");
    *output = nullptr;
    try {
        *output = resultHandle(dynsyslab::workflow::planJob(
            dynsyslab::core::json::parse(job_json)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR, "unknown workflow planning exception");
    }
}

dsl_status dsl_workflow_generate_initial_conditions_json(
    const char* solve_job_json, dsl_workflow_result_handle** output) {
    if (!solve_job_json || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT,
            "initial-condition generation argument is null");
    *output = nullptr;
    try {
        const auto job = dynsyslab::core::json::parse(solve_job_json);
        if (job.at("kind").string() != "solve")
            throw std::invalid_argument(
                "initial-condition generation requires a solve job");
        *output = resultHandle(dynsyslab::workflow::odeInitialConditionSetJson(
            dynsyslab::workflow::resolveOdeInitialConditions(job)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR,
            "unknown initial-condition generation exception");
    }
}

dsl_status dsl_workflow_run_json(const char* job_json, const char* archive_root_utf8,
                                 const char* progress_path_utf8,
                                 const dsl_execution_callbacks* callbacks,
                                 dsl_workflow_result_handle** output) {
    if (!job_json || !archive_root_utf8 || !output)
        return failed(DSL_STATUS_INVALID_ARGUMENT, "workflow argument is null");
    *output = nullptr;
    try {
        const auto job = dynsyslab::core::json::parse(job_json);
        const auto archiveRoot = std::filesystem::path(dynsyslab::core::wide(archive_root_utf8));
        const auto progressPath = progress_path_utf8
            ? std::filesystem::path(dynsyslab::core::wide(progress_path_utf8))
            : std::filesystem::path{};
        auto execution = context(callbacks);
        const auto run = dynsyslab::workflow::executeJob(
            job, archiveRoot, progressPath, &execution);
        dynsyslab::core::json::Value::Object result;
        result["schema"] = "dynsyslab.workflow-result/v1";
        result["kind"] = run.kind;
        result["archives"] = paths(run.archives);
        if (!run.manifest.empty())
            result[run.kind == "campaign" ? "manifest" : "hierarchy_manifest"] =
                run.manifest.generic_u8string();
        *output = resultHandle(dynsyslab::core::json::Value(std::move(result)));
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const dynsyslab::core::ExecutionCancelled& error) {
        return failed(DSL_STATUS_CANCELLED, error.what());
    } catch (const dynsyslab::core::ExecutionLimitExceeded& error) {
        return failed(DSL_STATUS_RESOURCE_LIMIT, error.what());
    } catch (const std::invalid_argument& error) {
        return failed(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::exception& error) {
        return failed(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return failed(DSL_STATUS_RUNTIME_ERROR, "unknown workflow exception");
    }
}

void dsl_workflow_result_destroy(dsl_workflow_result_handle* result) { delete result; }

dsl_status dsl_workflow_result_json(const dsl_workflow_result_handle* result,
                                    char* buffer, size_t capacity, size_t* required) {
    if (!result) return failed(DSL_STATUS_INVALID_ARGUMENT, "result handle is null");
    return copyText(result->json, buffer, capacity, required);
}

dsl_status dsl_workflow_last_error(char* buffer, size_t capacity, size_t* required) {
    return copyText(lastError, buffer, capacity, required);
}

}  // extern "C"

#include <dynsyslab/workflow_c_api.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        require(dsl_workflow_abi_version() == DSL_WORKFLOW_ABI_VERSION,
                "workflow ABI version");
        dsl_workflow_result_handle* observable = nullptr;
        require(dsl_workflow_observable_json(nullptr, &observable) == DSL_STATUS_INVALID_ARGUMENT,
                "observable null argument");
        const char* query = R"({"schema":"dynsyslab.observable-query/v1","system":{"name":"rpe"},"state":[1.2,0.3],"observable":"far_field_acoustic_echo","flow_order":2})";
        require(dsl_workflow_observable_json(query, &observable) == DSL_STATUS_OK,
                "analytic observable through C ABI");
        size_t observableSize = 0;
        dsl_workflow_result_json(observable, nullptr, 0, &observableSize);
        std::vector<char> observableText(observableSize);
        require(dsl_workflow_result_json(observable, observableText.data(), observableSize, nullptr)
            == DSL_STATUS_OK, "observable query result");
        require(std::string(observableText.data()).find("analytic") != std::string::npos,
                "observable analytic result");
        dsl_workflow_result_destroy(observable);
        const auto root = std::filesystem::temp_directory_path() / "dynsyslab-workflow-c-api";
        std::error_code ignored; std::filesystem::remove_all(root, ignored);
        const std::string job = R"({
          "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"c-api-solve",
          "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
          "time":{"start":0,"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},
          "initial_conditions":{"mode":"default","count":1,"clusters":1,"seed":0,"values":[],"ranges":[]},
          "controller":{"kind":"none","settings":{},"target":{"kind":"constant","state":[0,0]}},
          "solver":{"method":"RK45","rtol":1e-8,"atol":1e-10},
          "archive":{"enabled":true}
        })";
        dsl_workflow_result_handle* result = nullptr;
        require(dsl_workflow_plan_json(job.c_str(), &result) == DSL_STATUS_OK,
                "workflow plan");
        size_t required = 0;
        require(dsl_workflow_result_json(result, nullptr, 0, &required)
            == DSL_STATUS_BUFFER_TOO_SMALL, "workflow plan size");
        std::vector<char> planned(required);
        require(dsl_workflow_result_json(result, planned.data(), planned.size(), &required)
            == DSL_STATUS_OK, "workflow plan output");
        require(std::string(planned.data()).find("workflow-plan") != std::string::npos,
                "workflow plan content");
        dsl_workflow_result_destroy(result);
        result = nullptr;
        require(dsl_workflow_generate_initial_conditions_json(job.c_str(), &result)
            == DSL_STATUS_OK, "workflow IC generation");
        required = 0;
        require(dsl_workflow_result_json(result, nullptr, 0, &required)
            == DSL_STATUS_BUFFER_TOO_SMALL, "workflow IC output size");
        std::vector<char> initialConditions(required);
        require(dsl_workflow_result_json(result, initialConditions.data(),
            initialConditions.size(), &required) == DSL_STATUS_OK,
            "workflow IC output");
        require(std::string(initialConditions.data()).find(
            "dynsyslab.initial-conditions/v1") != std::string::npos,
            "workflow IC response content");
        dsl_workflow_result_destroy(result);
        result = nullptr;
        const auto rootText = root.u8string();
        require(dsl_workflow_run_json(job.c_str(), rootText.c_str(), nullptr, nullptr, &result)
            == DSL_STATUS_OK, "workflow run");
        required = 0;
        require(dsl_workflow_result_json(result, nullptr, 0, &required)
            == DSL_STATUS_BUFFER_TOO_SMALL, "workflow result size");
        std::vector<char> text(required);
        require(dsl_workflow_result_json(result, text.data(), text.size(), &required)
            == DSL_STATUS_OK, "workflow result output");
        require(std::string(text.data()).find("hierarchy_manifest") != std::string::npos,
                "workflow result content");
        dsl_workflow_result_destroy(result);
        dsl_execution_callbacks limited{};
        limited.maximum_memory_bytes = 32;
        result = nullptr;
        require(dsl_workflow_run_json(job.c_str(), rootText.c_str(), nullptr,
            &limited, &result) == DSL_STATUS_RESOURCE_LIMIT,
            "workflow memory limit status");
        require(result == nullptr, "limited workflow unexpectedly returned a result");
        std::filesystem::remove_all(root, ignored);
        std::cout << "dynsyslab_workflow_c_api_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_workflow_c_api_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}

#include <dynsyslab/c_api.h>

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Lattice.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Version.hpp>

#include <algorithm>
#include <cstring>
#include <memory>
#include <new>
#include <string>

struct dsl_system_handle {
    std::unique_ptr<dynsyslab::core::ODESystem> value;
};

namespace {
thread_local std::string lastError;

dsl_status fail(dsl_status status, const std::string& message) {
    lastError = message;
    return status;
}

dsl_status copyText(const std::string& value, char* buffer, size_t capacity,
                    size_t* required) {
    const size_t needed = value.size() + 1;
    if (required) *required = needed;
    if (!buffer || capacity < needed) return DSL_STATUS_BUFFER_TOO_SMALL;
    std::memcpy(buffer, value.c_str(), needed);
    return DSL_STATUS_OK;
}

dynsyslab::core::json::Value stringMap(
    const std::map<std::string, std::string>& values) {
    dynsyslab::core::json::Value::Object result;
    for (const auto& item : values) result[item.first] = item.second;
    return dynsyslab::core::json::Value(std::move(result));
}

dynsyslab::core::json::Value numberMap(const std::map<std::string, double>& values) {
    dynsyslab::core::json::Value::Object result;
    for (const auto& item : values) result[item.first] = item.second;
    return dynsyslab::core::json::Value(std::move(result));
}

dynsyslab::core::json::Value configurationJson(
    const dynsyslab::core::ODESystem& system) {
    const auto& configuration = system.configuration();
    dynsyslab::core::json::Value::Object result;
    result["schema"] = "dynsyslab.resolved-system/v1";
    result["identifier"] = configuration.identifier;
    result["preset"] = configuration.preset;
    result["coordinates"] = system.descriptor().preferredCoordinates;
    result["equations_version"] = configuration.equationsVersion;
    result["equations_fingerprint"] = dynsyslab::core::equationsFingerprint(configuration);
    result["parameters"] = numberMap(configuration.parameters);
    result["parameter_units"] = stringMap(configuration.parameterUnits);
    result["dimensionless_groups"] = numberMap(configuration.dimensionlessGroups);
    result["dimensionless_group_units"] = stringMap(configuration.dimensionlessGroupUnits);
    dynsyslab::core::json::Value::Array initial;
    for (double value : system.initialState()) initial.emplace_back(value);
    result["initial_state"] = dynsyslab::core::json::Value(std::move(initial));
    return dynsyslab::core::json::Value(std::move(result));
}

dynsyslab::core::SystemConfiguration parseConfiguration(const char* text) {
    if (!text) throw std::invalid_argument("configuration_json is null");
    const auto root = dynsyslab::core::json::parse(text);
    dynsyslab::core::SystemConfiguration result;
    result.identifier = root.at("identifier").string(root.at("name").string());
    result.preset = root.at("preset").string();
    result.equationsVersion = root.at("equations_version").string();
    for (const auto& item : root.at("parameters").object())
        result.parameters[item.first] = item.second.number();
    for (const auto& item : root.at("initial_state").array())
        result.initialState.push_back(item.number());
    const auto& driver = root.at("driver");
    result.driver.amplitude = driver.at("amplitude").number();
    result.driver.frequencyHz = driver.at("frequency_hz").number();
    result.driver.phaseRadians = driver.at("phase_radians").number();
    if (result.identifier.empty()) throw std::invalid_argument("system identifier is empty");
    return result;
}

template<class Action>
dsl_status guarded(Action action) {
    try {
        action();
        lastError.clear();
        return DSL_STATUS_OK;
    } catch (const dynsyslab::core::ExecutionCancelled& error) {
        return fail(DSL_STATUS_CANCELLED, error.what());
    } catch (const std::invalid_argument& error) {
        return fail(DSL_STATUS_INVALID_ARGUMENT, error.what());
    } catch (const std::out_of_range& error) {
        return fail(DSL_STATUS_NOT_FOUND, error.what());
    } catch (const std::exception& error) {
        return fail(DSL_STATUS_RUNTIME_ERROR, error.what());
    } catch (...) {
        return fail(DSL_STATUS_RUNTIME_ERROR, "unknown native exception");
    }
}
}  // namespace

extern "C" {

uint32_t dsl_abi_version(void) { return DSL_ABI_VERSION; }

const char* dsl_library_version(void) {
    return dynsyslab::core::version.data();
}

dsl_status dsl_last_error(char* buffer, size_t capacity, size_t* required) {
    return copyText(lastError, buffer, capacity, required);
}

dsl_status dsl_system_catalog_json(char* buffer, size_t capacity, size_t* required) {
    std::string output;
    const auto status = guarded([&] {
        static_cast<void>(dynsyslab::core::mackeyGlassDescriptor());
        static_cast<void>(dynsyslab::core::ising2DDescriptor());
        static_cast<void>(dynsyslab::core::heatEquationDescriptor());
        dynsyslab::core::json::Value::Object root;
        root["schema"] = "dynsyslab.system-catalog/v1";
        dynsyslab::core::json::Value::Array systems;
        for (const auto& descriptor : dynsyslab::core::availableSystems()) {
            dynsyslab::core::json::Value::Object item;
            item["identifier"] = descriptor.identifier;
            item["label"] = descriptor.label;
            item["family"] = descriptor.family;
            item["state_dimension"] = static_cast<double>(descriptor.stateDimension);
            item["preferred_coordinates"] = descriptor.preferredCoordinates;
            systems.emplace_back(std::move(item));
        }
        root["systems"] = dynsyslab::core::json::Value(std::move(systems));
        output = dynsyslab::core::json::dump(dynsyslab::core::json::Value(std::move(root)));
    });
    return status == DSL_STATUS_OK ? copyText(output, buffer, capacity, required) : status;
}

dsl_status dsl_system_create_json(const char* configuration_json,
                                  dsl_system_handle** output) {
    if (!output) return fail(DSL_STATUS_INVALID_ARGUMENT, "output handle is null");
    *output = nullptr;
    return guarded([&] {
        auto handle = std::make_unique<dsl_system_handle>();
        handle->value = dynsyslab::core::createSystem(parseConfiguration(configuration_json));
        *output = handle.release();
    });
}

void dsl_system_destroy(dsl_system_handle* system) { delete system; }

dsl_status dsl_system_dimension(const dsl_system_handle* system, size_t* dimension) {
    if (!system || !system->value || !dimension)
        return fail(DSL_STATUS_INVALID_ARGUMENT, "invalid system or dimension pointer");
    *dimension = system->value->descriptor().stateDimension;
    return DSL_STATUS_OK;
}

dsl_status dsl_system_configuration_json(const dsl_system_handle* system,
                                         char* buffer, size_t capacity,
                                         size_t* required) {
    if (!system || !system->value)
        return fail(DSL_STATUS_INVALID_ARGUMENT, "invalid system handle");
    std::string output;
    const auto status = guarded([&] {
        output = dynsyslab::core::json::dump(configurationJson(*system->value));
    });
    return status == DSL_STATUS_OK ? copyText(output, buffer, capacity, required) : status;
}

dsl_status dsl_system_rhs(const dsl_system_handle* system, double time,
                          const double* state, size_t state_count, double input,
                          double* derivative, size_t derivative_count) {
    if (!system || !system->value || !state || !derivative)
        return fail(DSL_STATUS_INVALID_ARGUMENT, "invalid RHS argument");
    const size_t dimension = system->value->descriptor().stateDimension;
    if (state_count != dimension || derivative_count != dimension)
        return fail(DSL_STATUS_INVALID_ARGUMENT, "RHS dimension mismatch");
    return guarded([&] {
        system->value->rhsWithInput(time, state, input, derivative);
    });
}

}  // extern "C"

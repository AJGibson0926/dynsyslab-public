#include <dynsyslab/core/Extension.hpp>

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/extension_api.h>
#include "ExtensionModels.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include <utility>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace dynsyslab::core {
namespace {

struct Module {
    std::filesystem::path path;
    const dsl_extension_v1* extension{};
#ifdef _WIN32
    HMODULE handle{};
#else
    void* handle{};
#endif
};

void closeModule(
#ifdef _WIN32
    HMODULE handle
#else
    void* handle
#endif
) {
    if (!handle) return;
#ifdef _WIN32
    FreeLibrary(handle);
#else
    dlclose(handle);
#endif
}

std::string configurationJson(const SystemConfiguration& configuration) {
    json::Value::Object root;
    root["identifier"] = configuration.identifier;
    root["preset"] = configuration.preset;
    root["equations_version"] = configuration.equationsVersion;
    json::Value::Object parameters;
    for (const auto& item : configuration.parameters) parameters[item.first] = item.second;
    root["parameters"] = json::Value(std::move(parameters));
    json::Value::Array initial;
    for (double value : configuration.initialState) initial.emplace_back(value);
    root["initial_state"] = json::Value(std::move(initial));
    return json::dump(json::Value(std::move(root)));
}

void requireCallback(int status, const std::array<char, 1024>& error,
                     const char* operation) {
    if (status == 0) return;
    throw std::runtime_error(std::string("extension ") + operation + " failed: "
                             + (error[0] ? error.data() : "unspecified error"));
}

SystemDescriptor parseDescriptor(const dsl_system_extension_v1& api) {
    if (!api.descriptor_json) throw std::runtime_error("extension system descriptor is null");
    const auto value = json::parse(api.descriptor_json);
    SystemDescriptor result;
    result.identifier = value.at("identifier").string();
    result.label = value.at("label").string(result.identifier);
    result.stateDimension = static_cast<std::size_t>(value.at("state_dimension").number());
    result.preferredCoordinates = value.at("preferred_coordinates").string("dimensional");
    result.family = value.at("family").string("ode");
    for (const auto& item : value.at("state_names").array()) result.stateNames.push_back(item.string());
    for (const auto& item : value.at("state_labels").array()) result.stateLabels.push_back(item.string());
    for (const auto& item : value.at("presets").array())
        result.presets.push_back({item.at("key").string(), item.at("description").string()});
    if (result.identifier.empty() || result.stateDimension == 0)
        throw std::runtime_error("extension system descriptor is incomplete");
    return result;
}

class ExtensionSystem final : public ODESystem {
public:
    ExtensionSystem(const dsl_system_extension_v1& api,
                    SystemConfiguration configuration)
        : api_(api), descriptor_(parseDescriptor(api)), configuration_(std::move(configuration)) {
        if (configuration_.equationsVersion.empty())
            configuration_.equationsVersion = "extension." + descriptor_.identifier + "/1";
        std::array<char, 1024> error{};
        const std::string request = configurationJson(configuration_);
        instance_ = api_.create ? api_.create(request.c_str(), error.data(), error.size()) : nullptr;
        if (!instance_)
            throw std::runtime_error(std::string("extension system creation failed: ")
                + (error[0] ? error.data() : "null instance"));
        if (api_.resolved_configuration_json) {
            const char* text = api_.resolved_configuration_json(instance_);
            if (text) applyResolved(json::parse(text));
        }
    }

    ~ExtensionSystem() override { if (instance_ && api_.destroy) api_.destroy(instance_); }
    void resetRuntimeState() override { if (api_.reset) api_.reset(instance_); }
    const SystemDescriptor& descriptor() const override { return descriptor_; }
    const SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override {
        if (!configuration_.initialState.empty()) return configuration_.initialState;
        return std::vector<double>(descriptor_.stateDimension, 0.0);
    }
    void rhs(double time, const double* state, double* derivative) const override {
        rhsWithInput(time, state, configuration_.driver.value(time), derivative);
    }
    void rhsWithInput(double time, const double* state, double input,
                      double* derivative) const override {
        std::array<char, 1024> error{};
        requireCallback(api_.rhs(instance_, time, state, descriptor_.stateDimension,
            input, derivative, descriptor_.stateDimension, error.data(), error.size()),
            error, "RHS");
    }
    std::vector<double> jacobian(double time, const double* state) const override {
        std::vector<double> result(descriptor_.stateDimension * descriptor_.stateDimension);
        if (!api_.jacobian) return finiteDifference(time, state);
        std::array<char, 1024> error{};
        requireCallback(api_.jacobian(instance_, time, state, descriptor_.stateDimension,
            result.data(), result.size(), error.data(), error.size()), error, "Jacobian");
        return result;
    }
    std::vector<double> inputMatrix(const double* state) const override {
        std::vector<double> result(descriptor_.stateDimension);
        if (!api_.input_matrix) return result;
        std::array<char, 1024> error{};
        requireCallback(api_.input_matrix(instance_, state, descriptor_.stateDimension,
            result.data(), result.size(), error.data(), error.size()), error, "input matrix");
        return result;
    }
    bool invalid(double time, const double* state) const override {
        return api_.invalid && api_.invalid(instance_, time, state, descriptor_.stateDimension) != 0;
    }
    std::vector<Observable> modelObservables() const override {
        std::vector<Observable> result;
        if (!api_.observable_count || !api_.observable_descriptor_json
            || !api_.observable_evaluate) return result;
        const std::size_t count = api_.observable_count(instance_);
        result.reserve(count);
        for (std::size_t index=0; index<count; ++index) {
            const char* text=api_.observable_descriptor_json(instance_,index);
            if (!text) throw std::runtime_error("extension observable descriptor is null");
            const auto source=json::parse(text);
            const std::size_t output=static_cast<std::size_t>(
                source.at("output_dimension").number(1.0));
            ObservableDescriptor description{
                source.at("identifier").string(), source.at("label").string(),
                source.at("units").string("dimensionless"),
                source.at("coordinate_system").string(descriptor_.preferredCoordinates),
                {}, false};
            Observable::Evaluator evaluator=[this,index,output](const double* state) {
                std::vector<double> values(output); std::array<char,1024> error{};
                requireCallback(api_.observable_evaluate(instance_,index,state,
                    descriptor_.stateDimension,values.data(),values.size(),error.data(),
                    error.size()),error,"observable evaluation");
                Observable::Values result; result.reserve(output);
                for (double value:values) result.emplace_back(value,0.0);
                return result;
            };
            Observable::Jacobian jacobian;
            if (api_.observable_jacobian) {
                jacobian=[this,index,output](const double* state) {
                    std::vector<double> values(output*descriptor_.stateDimension);
                    std::array<char,1024> error{};
                    requireCallback(api_.observable_jacobian(instance_,index,state,
                        descriptor_.stateDimension,values.data(),values.size(),error.data(),
                        error.size()),error,"observable Jacobian");
                    Observable::Values result; result.reserve(values.size());
                    for (double value:values) result.emplace_back(value,0.0);
                    return result;
                };
            }
            result.emplace_back(std::move(description),descriptor_.stateDimension,output,
                                std::move(evaluator),std::move(jacobian));
        }
        return result;
    }

private:
    void applyResolved(const json::Value& value) {
        configuration_.preset = value.at("preset").string(configuration_.preset);
        configuration_.equationsVersion = value.at("equations_version").string(
            configuration_.equationsVersion);
        for (const auto& item : value.at("parameters").object())
            configuration_.parameters[item.first] = item.second.number();
        for (const auto& item : value.at("parameter_units").object())
            configuration_.parameterUnits[item.first] = item.second.string();
        for (const auto& item : value.at("dimensionless_groups").object())
            configuration_.dimensionlessGroups[item.first] = item.second.number();
        for (const auto& item : value.at("dimensionless_group_units").object())
            configuration_.dimensionlessGroupUnits[item.first] = item.second.string();
        if (!value.at("initial_state").array().empty()) {
            configuration_.initialState.clear();
            for (const auto& item : value.at("initial_state").array())
                configuration_.initialState.push_back(item.number());
        }
    }
    std::vector<double> finiteDifference(double time, const double* state) const {
        const std::size_t n = descriptor_.stateDimension;
        std::vector<double> result(n*n), plus(state,state+n), minus(plus), fp(n), fm(n);
        for (std::size_t column=0; column<n; ++column) {
            const double step=1.0e-6*(std::max)(1.0,std::abs(state[column]));
            plus[column]+=step; minus[column]-=step;
            rhs(time,plus.data(),fp.data()); rhs(time,minus.data(),fm.data());
            for (std::size_t row=0;row<n;++row)
                result[row*n+column]=(fp[row]-fm[row])/(2.0*step);
            plus[column]=state[column]; minus[column]=state[column];
        }
        return result;
    }

    const dsl_system_extension_v1& api_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
    void* instance_{};
};

}  // namespace

struct ExtensionRegistry::Impl {
    mutable std::mutex mutex;
    std::vector<Module> modules;
    ~Impl() {
        for (auto& module : modules) closeModule(module.handle);
    }
};

ExtensionRegistry::ExtensionRegistry() : impl_(std::make_unique<Impl>()) {}
ExtensionRegistry::~ExtensionRegistry() = default;
ExtensionRegistry& ExtensionRegistry::instance() {
    static ExtensionRegistry registry;
    return registry;
}

std::string ExtensionRegistry::load(const std::filesystem::path& path) {
#ifdef _WIN32
    HMODULE handle = LoadLibraryW(path.c_str());
    if (!handle) throw std::runtime_error("unable to load DynSysLab extension");
    auto entry = reinterpret_cast<dsl_extension_entry_v1>(
        GetProcAddress(handle, DSL_EXTENSION_ENTRY_NAME));
#else
    void* handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!handle) throw std::runtime_error(dlerror());
    auto entry = reinterpret_cast<dsl_extension_entry_v1>(dlsym(handle, DSL_EXTENSION_ENTRY_NAME));
#endif
    if (!entry) {
        closeModule(handle);
        throw std::runtime_error("DynSysLab extension entry point is missing");
    }
    const dsl_extension_v1* extension = entry();
    if (!extension || extension->abi_version != DSL_EXTENSION_ABI_VERSION
        || extension->struct_size < DSL_EXTENSION_V1_CORE_SIZE
        || !extension->identifier || !extension->identifier[0])
    {
        closeModule(handle);
        throw std::runtime_error("DynSysLab extension ABI descriptor is invalid");
    }
    std::lock_guard<std::mutex> lock(impl_->mutex);
    for (const auto& module : impl_->modules) {
        if (module.extension->identifier == std::string(extension->identifier)) {
            std::error_code error;
            const bool same = std::filesystem::equivalent(module.path, path, error);
            closeModule(handle);
            if (same && !error) return extension->identifier;
            throw std::runtime_error("DynSysLab extension identifier is already loaded");
        }
    }
    impl_->modules.push_back({path, extension, handle});
    return extension->identifier;
}

std::vector<ExtensionDescriptor> ExtensionRegistry::extensions() const {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    std::vector<ExtensionDescriptor> result;
    for (const auto& module : impl_->modules)
        result.push_back({module.extension->identifier,
            module.extension->version ? module.extension->version : "",
            module.extension->capabilities_json ? module.extension->capabilities_json : "{}",
            module.path});
    return result;
}

std::vector<SystemDescriptor> ExtensionRegistry::systems() const {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    std::vector<SystemDescriptor> result;
    for (const auto& module : impl_->modules)
        for (std::size_t index = 0; index < module.extension->system_count; ++index)
            result.push_back(parseDescriptor(module.extension->systems[index]));
    return result;
}

std::unique_ptr<ODESystem> ExtensionRegistry::createSystem(
    const SystemConfiguration& configuration) const {
    const dsl_system_extension_v1* selected = nullptr;
    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        for (const auto& module : impl_->modules)
            for (std::size_t index=0; index<module.extension->system_count; ++index) {
                const auto& api=module.extension->systems[index];
                if (parseDescriptor(api).identifier==configuration.identifier)
                    selected = &api;
            }
    }
    if (selected) return std::make_unique<ExtensionSystem>(*selected,configuration);
    return {};
}

std::vector<ExtensionModel> ExtensionRegistry::models() const {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    std::vector<ExtensionModel> result;
    for (const auto& module : impl_->modules) {
        auto models = detail::extensionModels(*module.extension);
        result.insert(result.end(), models.begin(), models.end());
    }
    return result;
}

std::optional<ExtensionModel> ExtensionRegistry::model(const std::string& identifier) const {
    for (auto& candidate : models())
        if (candidate.identifier == identifier) return candidate;
    return std::nullopt;
}

std::string ExtensionRegistry::execute(const std::string& extension,
                                       const std::string& capability,
                                       const std::string& requestJson) const {
    const dsl_extension_v1* selected = nullptr;
    {
        std::lock_guard<std::mutex> lock(impl_->mutex);
        for (const auto& module : impl_->modules)
            if (module.extension->identifier == extension) selected = module.extension;
    }
    if (!selected)
        throw std::invalid_argument("DynSysLab extension is not loaded: " + extension);
    if (!selected->execute_json)
        throw std::runtime_error("extension has no generic operation entry point");
    std::string output; std::array<char,1024> error{};
    const int status=selected->execute_json(capability.c_str(),requestJson.c_str(),
        [](const char* data,size_t size,void* raw) {
            static_cast<std::string*>(raw)->append(data,size);
        },&output,error.data(),error.size());
    requireCallback(status,error,"operation");
    return output;
}

}  // namespace dynsyslab::core

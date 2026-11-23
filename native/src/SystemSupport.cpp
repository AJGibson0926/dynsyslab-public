#include "SystemSupport.hpp"
#include <dynsyslab/core/SystemLabels.hpp>
#include "AnalyticSystems.hpp"
#include "DynamicalKnowledgeInternal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core::detail {

FunctionalSystem::FunctionalSystem(SystemDescriptor descriptor, SystemConfiguration configuration,
                                   std::vector<double> initialState, RHSFunction rhs,
                                   JacobianFunction jacobian, InputFunction input,
                                   InvalidFunction invalid, std::vector<double> controlledInitialState,
                                   std::vector<Observable> observables,
                                   std::map<std::string, std::string> observableAliases,
                                   std::map<std::string, std::vector<std::string>> coordinatePresets,
                                   std::vector<AnalyticExpression> drift,
                                   std::vector<AnalyticExpression> inputFields)
    : descriptor_(std::move(descriptor)), configuration_(std::move(configuration)),
      initialState_(std::move(initialState)), controlledInitialState_(std::move(controlledInitialState)),
      drift_(std::move(drift)), inputFields_(std::move(inputFields)), rhs_(std::move(rhs)),
      jacobian_(std::move(jacobian)), input_(std::move(input)), invalid_(std::move(invalid)),
      observables_(std::move(observables)), observableAliases_(std::move(observableAliases)),
      coordinatePresets_(std::move(coordinatePresets)) {
    if (initialState_.size() != descriptor_.stateDimension)
        throw std::invalid_argument("native system initial-state dimension mismatch");
    if (!controlledInitialState_.empty() && controlledInitialState_.size() != descriptor_.stateDimension)
        throw std::invalid_argument("native system controlled initial-state dimension mismatch");
    dynamicalKnowledge_ = builtInDynamicalKnowledge(descriptor_, configuration_);
    analyticSolutions_ = oscillatorAnalyticSolutions(configuration_);
    auto bubble = rpeAnalyticSolutions(configuration_);
    auto other = otherAnalyticSolutions(descriptor_, configuration_);
    auto integrable = integrableAnalyticSolutions(configuration_);
    analyticSolutions_.insert(analyticSolutions_.end(), bubble.begin(), bubble.end());
    analyticSolutions_.insert(analyticSolutions_.end(), other.begin(), other.end());
    analyticSolutions_.insert(analyticSolutions_.end(), integrable.begin(), integrable.end());
    for (auto& family : analyticSolutions_) {
        auto instantiate = family.instantiate;
        auto dimension = descriptor_.stateDimension;
        family.instantiate = [instantiate, dimension](const AnalyticState& initial, double t0) {
            if (initial.size() != dimension || !std::isfinite(t0)
                || !std::all_of(initial.begin(), initial.end(), [](double x) { return std::isfinite(x); }))
                return AnalyticTrajectory{{}, "Initial state/time must be finite and match the system dimension"};
            return instantiate(initial, t0);
        };
    }
}

void FunctionalSystem::rhs(double time, const double* state, double* derivative) const {
    const auto scaling = coordinates();
    const double clock = scaling.physicalTimeAvailable
        ? scaling.physicalTime(time) : time;
    const double physical = configuration_.driver.value(clock);
    const double input = scaling.physicalInputAvailable
        ? scaling.modelInput(physical) : physical;
    rhs_(time, state, input, derivative);
}
void FunctionalSystem::rhsWithInput(double time, const double* state, double input,
                                    double* derivative) const { rhs_(time, state, input, derivative); }
std::vector<double> FunctionalSystem::jacobian(double time, const double* state) const {
    return jacobian_(time, state);
}
std::vector<double> FunctionalSystem::inputMatrix(const double* state) const { return input_(state); }
bool FunctionalSystem::invalid(double time, const double* state) const { return invalid_(time, state); }

double parameter(const SystemConfiguration& configuration, const char* name, double fallback) {
    const auto found = configuration.parameters.find(name);
    return found == configuration.parameters.end() ? fallback : found->second;
}

namespace {
std::string inferredUnit(const std::string& name) {
    const auto endsWith = [&](const char* suffix) {
        const std::string value(suffix);
        return name.size() >= value.size()
            && name.compare(name.size() - value.size(), value.size(), value) == 0;
    };
    if (endsWith("_kg_m3")) return "kg/m^3";
    if (endsWith("_pa_s")) return "Pa s";
    if (endsWith("_n_m")) return "N/m";
    if (endsWith("_kg_s")) return "kg/s";
    if (endsWith("_m_s")) return "m/s";
    if (endsWith("_m_s2")) return "m/s^2";
    if (endsWith("_pa")) return "Pa";
    if (endsWith("_hz")) return "Hz";
    if (endsWith("_kg")) return "kg";
    if (endsWith("_m")) return "m";
    if (name.find("current") != std::string::npos) return "model current";
    return "dimensionless";
}
}  // namespace

double resolvedParameter(SystemConfiguration& configuration, const char* name,
                         double fallback, const char* unit) {
    const double value = parameter(configuration, name, fallback);
    configuration.parameters[name] = value;
    configuration.parameterUnits[name] = unit ? unit : inferredUnit(name);
    return value;
}

std::size_t integerParameter(const SystemConfiguration& configuration,
                             const char* name, std::size_t fallback) {
    const double value = parameter(configuration, name, static_cast<double>(fallback));
    if (!std::isfinite(value) || value < 1.0 || std::floor(value) != value)
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(value);
}

std::size_t resolvedIntegerParameter(SystemConfiguration& configuration,
                                     const char* name, std::size_t fallback) {
    const auto value = integerParameter(configuration, name, fallback);
    configuration.parameters[name] = static_cast<double>(value);
    configuration.parameterUnits[name] = "count";
    return value;
}

std::vector<double> selectedInitial(const SystemConfiguration& configuration,
                                    std::vector<double> fallback) {
    return configuration.initialState.empty() ? fallback : configuration.initialState;
}

std::vector<double> finiteDifferenceJacobian(const ODESystem& system, double time,
                                             const double* state) {
    const std::size_t dimension = system.descriptor().stateDimension;
    std::vector<double> matrix(dimension * dimension);
    std::vector<double> plus(state, state + dimension), minus(plus), fPlus(dimension), fMinus(dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
        const double step = 1.0e-6 * std::max(1.0, std::abs(state[column]));
        plus[column] += step;
        minus[column] -= step;
        system.rhs(time, plus.data(), fPlus.data());
        system.rhs(time, minus.data(), fMinus.data());
        for (std::size_t row = 0; row < dimension; ++row)
            matrix[row * dimension + column] = (fPlus[row] - fMinus[row]) / (2.0 * step);
        plus[column] = state[column];
        minus[column] = state[column];
    }
    return matrix;
}

bool genericInvalid(const double* state, std::size_t dimension, double bound) {
    double norm2 = 0.0;
    for (std::size_t index = 0; index < dimension; ++index) {
        if (!std::isfinite(state[index])) return true;
        norm2 += state[index] * state[index];
    }
    return std::sqrt(norm2) > bound;
}

SystemDescriptor descriptorFor(const std::string& identifier, std::size_t dimensionOverride) {
    const auto& catalog = systemCatalog();
    const auto found = std::find_if(catalog.begin(), catalog.end(), [&](const auto& item) {
        return item.identifier == identifier;
    });
    if (found == catalog.end()) throw std::invalid_argument("unknown native system descriptor");
    auto result = *found;
    if (dimensionOverride) {
        result.stateDimension = dimensionOverride;
        if (identifier == "kuramoto") result.stateLabels = systemStateLabels(result, dimensionOverride);
    }
    return result;
}

}  // namespace dynsyslab::core::detail

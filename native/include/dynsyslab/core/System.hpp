#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <dynsyslab/core/Observable.hpp>
#include <dynsyslab/core/AnalyticSolution.hpp>
#include <dynsyslab/core/DynamicalKnowledge.hpp>

namespace dynsyslab::core {

struct SinusoidalDriver {
    // Physical amplitude/frequency when the system declares physical
    // coordinate scaling; otherwise native model input/time values.
    double amplitude{};
    double frequencyHz{};
    double phaseRadians{};

    [[nodiscard]] double value(double time) const;
};

struct PresetDescriptor {
    std::string key;
    std::string description;
    // Recommended external driver (amplitude 0: free run), the source the preset reproduces,
    // and per-method processing settings and provenance as a JSON object (PresetMetadata.cpp).
    // A solve applies the driver only when the job asks for it (system.driver = "preset").
    SinusoidalDriver driver{};
    std::string reference;
    std::string processingJson{"{}"};
};

struct SystemConfiguration;

struct SystemDescriptor {
    std::string identifier;
    std::string label;
    std::size_t stateDimension{};
    std::string preferredCoordinates{"dimensional"};
    std::vector<std::string> stateNames;
    std::vector<std::string> stateLabels;
    std::vector<PresetDescriptor> presets;
    std::string family{"ode"};
    // System-level research is independent of individual trajectories and parameter runs.
    DynamicalKnowledge research;
    // Family modules bind their closed-form evaluators to resolved parameters.
    std::function<std::vector<AnalyticSolutionFamily>(const SystemConfiguration&)> analyticSolutions;
    ObservableDerivativeLinks observableTimeDerivatives;
};

// Conversion between the coordinates advanced by the native model and a
// physical clock/input when the model has a physical nondimensionalization.
// Abstract dimensionless models leave the availability flags false and use
// identity scales; their prescribed inputs therefore remain model-native.
struct CoordinateScaling {
    double secondsPerModelTime{1.0};
    double physicalInputPerModelInput{1.0};
    bool physicalTimeAvailable{false};
    bool physicalInputAvailable{false};
    std::string modelTimeUnits{"system time"};
    std::string modelInputUnits{"system input"};
    std::string physicalInputUnits{"system input"};
    // Physical value of one model unit of each state (for example R0 and R0/T for a bubble's
    // radius and wall velocity) and its unit. Empty: the state has no declared physical
    // scale and the system cannot be integrated in dimensional coordinates.
    std::vector<double> physicalStatePerModelState;
    std::vector<std::string> physicalStateUnits;

    [[nodiscard]] double physicalTime(double modelTime) const;
    [[nodiscard]] double modelTime(double physicalTime) const;
    [[nodiscard]] double physicalInput(double modelInput) const;
    [[nodiscard]] double modelInput(double physicalInput) const;
};

struct SystemConfiguration {
    std::string identifier;
    std::string preset;
    std::map<std::string, double> parameters;
    std::map<std::string, std::string> parameterUnits;
    std::map<std::string, double> dimensionlessGroups;
    std::map<std::string, std::string> dimensionlessGroupUnits;
    std::string equationsVersion;
    std::vector<double> initialState;
    SinusoidalDriver driver;
    CoordinateScaling coordinateScaling;
};

class ODESystem {
public:
    virtual ~ODESystem() = default;
    virtual void resetRuntimeState() {}
    [[nodiscard]] virtual const SystemDescriptor& descriptor() const = 0;
    [[nodiscard]] virtual const SystemConfiguration& configuration() const = 0;
    [[nodiscard]] CoordinateScaling coordinates() const;
    [[nodiscard]] virtual std::vector<double> initialState() const = 0;
    // Default used only when a controlled solve does not explicitly provide an IC.
    // Systems with a scientifically distinguished equilibrium override this value.
    [[nodiscard]] virtual std::vector<double> controlledInitialState() const {
        return initialState();
    }
    virtual void rhs(double time, const double* state, double* derivative) const = 0;
    virtual void rhsWithInput(double time, const double* state, double input,
                              double* derivative) const = 0;
    [[nodiscard]] virtual std::vector<double> jacobian(double time,
                                                       const double* state) const = 0;
    [[nodiscard]] virtual std::vector<double> inputMatrix(const double* state) const = 0;
    [[nodiscard]] virtual bool invalid(double time, const double* state) const = 0;
    // Autonomous, zero-input drift in native model coordinates; empty when unavailable.
    [[nodiscard]] virtual std::vector<AnalyticExpression> driftExpressions() const { return {}; }
    [[nodiscard]] ObservableCatalog observableCatalog() const;
    [[nodiscard]] virtual std::vector<Observable> modelObservables() const { return {}; }
    [[nodiscard]] virtual std::map<std::string, std::string> observableAliases() const {
        return {};
    }
    [[nodiscard]] virtual std::map<std::string, std::vector<std::string>>
        observableCoordinatePresets() const { return {}; }
    // Row-major state x input analytic fields, matching inputMatrix; empty when unavailable.
    [[nodiscard]] virtual std::vector<AnalyticExpression> inputExpressions() const { return {}; }
    // First-class system properties. Extensions may supply their own researched families/guides.
    [[nodiscard]] virtual std::vector<AnalyticSolutionFamily> analyticSolutions() const {
        return descriptor().analyticSolutions ? descriptor().analyticSolutions(configuration())
                                              : std::vector<AnalyticSolutionFamily>{};
    }
    [[nodiscard]] virtual DynamicalKnowledge dynamicalKnowledge() const { return descriptor().research; }
};

[[nodiscard]] const std::vector<SystemDescriptor>& systemCatalog();
using SystemCatalogProvider = std::vector<SystemDescriptor> (*)();
void registerSystemCatalogProvider(SystemCatalogProvider provider);
[[nodiscard]] std::vector<SystemDescriptor> availableSystems();
[[nodiscard]] std::unique_ptr<ODESystem> createSystem(const SystemConfiguration& configuration);
[[nodiscard]] std::string equationsFingerprint(const SystemConfiguration& configuration);

}  // namespace dynsyslab::core

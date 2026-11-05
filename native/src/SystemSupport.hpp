#pragma once

#include <dynsyslab/core/System.hpp>

#include <functional>
#include <map>

namespace dynsyslab::core::detail {

using RHSFunction = std::function<void(double, const double*, double, double*)>;
using JacobianFunction = std::function<std::vector<double>(double, const double*)>;
using InputFunction = std::function<std::vector<double>(const double*)>;
using InvalidFunction = std::function<bool(double, const double*)>;

class FunctionalSystem final : public ODESystem {
public:
    FunctionalSystem(SystemDescriptor descriptor, SystemConfiguration configuration,
                     std::vector<double> initialState, RHSFunction rhs,
                     JacobianFunction jacobian, InputFunction input,
                     InvalidFunction invalid, std::vector<double> controlledInitialState = {},
                     std::vector<Observable> observables = {},
                     std::map<std::string, std::string> observableAliases = {},
                     std::map<std::string, std::vector<std::string>> coordinatePresets = {},
                     std::vector<AnalyticExpression> drift = {},
                     std::vector<AnalyticExpression> inputFields = {});
    const SystemDescriptor& descriptor() const override { return descriptor_; }
    const SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return initialState_; }
    std::vector<double> controlledInitialState() const override {
        return controlledInitialState_.empty() ? initialState_ : controlledInitialState_;
    }
    void rhs(double time, const double* state, double* derivative) const override;
    void rhsWithInput(double time, const double* state, double input,
                      double* derivative) const override;
    std::vector<double> jacobian(double time, const double* state) const override;
    std::vector<double> inputMatrix(const double* state) const override;
    bool invalid(double time, const double* state) const override;
    DynamicalKnowledge dynamicalKnowledge() const override { return dynamicalKnowledge_; }
    std::vector<AnalyticSolutionFamily> analyticSolutions() const override { return analyticSolutions_; }
    std::vector<AnalyticExpression> driftExpressions() const override { return drift_; }
    std::vector<AnalyticExpression> inputExpressions() const override { return inputFields_; }
    std::vector<Observable> modelObservables() const override { return observables_; }
    std::map<std::string, std::string> observableAliases() const override {
        return observableAliases_;
    }
    std::map<std::string, std::vector<std::string>>
        observableCoordinatePresets() const override { return coordinatePresets_; }

private:
    DynamicalKnowledge dynamicalKnowledge_;
    std::vector<AnalyticSolutionFamily> analyticSolutions_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
    std::vector<double> initialState_;
    std::vector<double> controlledInitialState_;
    std::vector<AnalyticExpression> drift_, inputFields_;
    RHSFunction rhs_;
    JacobianFunction jacobian_;
    InputFunction input_;
    InvalidFunction invalid_;
    std::vector<Observable> observables_;
    std::map<std::string, std::string> observableAliases_;
    std::map<std::string, std::vector<std::string>> coordinatePresets_;
};

double parameter(const SystemConfiguration& configuration, const char* name,
                 double fallback);
double resolvedParameter(SystemConfiguration& configuration, const char* name,
                         double fallback, const char* unit = nullptr);
std::size_t integerParameter(const SystemConfiguration& configuration,
                             const char* name, std::size_t fallback);
std::size_t resolvedIntegerParameter(SystemConfiguration& configuration,
                                     const char* name, std::size_t fallback);
std::vector<double> selectedInitial(const SystemConfiguration& configuration,
                                    std::vector<double> fallback);
std::vector<double> finiteDifferenceJacobian(const ODESystem& system, double time,
                                             const double* state);
bool genericInvalid(const double* state, std::size_t dimension, double bound = 1.0e8);
SystemDescriptor descriptorFor(const std::string& identifier,
                               std::size_t dimensionOverride = 0);

std::unique_ptr<ODESystem> makeHarmonicSystem(SystemConfiguration configuration);
std::unique_ptr<ODESystem> makeBasicSystem(const SystemConfiguration& configuration);
std::unique_ptr<ODESystem> makeMechanicalSystem(const SystemConfiguration& configuration);
std::unique_ptr<ODESystem> makeBiologicalSystem(const SystemConfiguration& configuration);
std::unique_ptr<ODESystem> makeBubbleSystem(const SystemConfiguration& configuration);

}  // namespace dynsyslab::core::detail

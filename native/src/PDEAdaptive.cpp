#include <dynsyslab/core/PDE.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

// Adaptive method of lines: the spatially discretised field u' = F(t, u) is an ODE system, so
// it is integrated by the native ODE integrator with error control (Schiesser 1991; method
// references in Integrator.hpp).

namespace dynsyslab::core {
namespace {

class FieldSystem final : public ODESystem {
public:
    FieldSystem(const StructuredMesh& mesh, std::size_t values, const PDERightHandSide& rhs)
        : mesh_(mesh), values_(values), rhs_(rhs) {
        descriptor_.identifier = "pde_field";
        descriptor_.label = "Method-of-lines field";
        descriptor_.stateDimension = values;
        descriptor_.family = "pde";
    }
    const SystemDescriptor& descriptor() const override { return descriptor_; }
    const SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return std::vector<double>(values_, 0.0); }
    void rhs(double time, const double* state, double* derivative) const override {
        rhs_(time, mesh_, state, values_, derivative);
    }
    void rhsWithInput(double time, const double* state, double, double* derivative) const override {
        rhs(time, state, derivative);
    }
    // Forward-difference Jacobian (used only by the LSODA stiffness test).
    std::vector<double> jacobian(double time, const double* state) const override {
        std::vector<double> matrix(values_ * values_), base(values_), shifted(state, state + values_), f(values_);
        rhs(time, state, base.data());
        for (std::size_t c = 0; c < values_; ++c) {
            const double delta = 1.0e-7 * std::max(1.0, std::abs(state[c]));
            shifted[c] = state[c] + delta;
            rhs(time, shifted.data(), f.data());
            for (std::size_t r = 0; r < values_; ++r) matrix[r * values_ + c] = (f[r] - base[r]) / delta;
            shifted[c] = state[c];
        }
        return matrix;
    }
    std::vector<double> inputMatrix(const double*) const override { return {}; }
    bool invalid(double, const double* state) const override {
        return !std::all_of(state, state + values_, [](double value) { return std::isfinite(value); });
    }

private:
    const StructuredMesh& mesh_;
    std::size_t values_;
    const PDERightHandSide& rhs_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
};

}  // namespace

PDEIntegrationResult integratePDEAdaptive(
    const StructuredMesh& mesh, const std::vector<std::size_t>& fieldShape,
    const std::vector<double>& initialField, const PDEIntegrationOptions& options,
    const PDERightHandSide& rhs, const ExecutionContext* execution) {
    std::size_t values = 1;
    for (const auto extent : fieldShape) values *= extent;
    if (!values || values != initialField.size())
        throw std::invalid_argument("PDE initial field does not match its field shape");
    FieldSystem system(mesh, values, rhs);
    IntegrationOptions ode;
    ode.initialTime = options.initialTime;
    ode.finalTime = options.finalTime;
    ode.outputStep = options.outputStep;
    ode.method = options.adaptiveMethod;
    ode.relativeTolerance = options.relativeTolerance;
    ode.absoluteTolerance = options.absoluteTolerance;
    ode.maximumStep = options.maximumStep;
    // The output step caps RK45 steps unless a maximum is given; MOL fields are smooth in time.
    if (ode.maximumStep <= 0.0 && ode.method == IntegrationMethod::dormandPrince45)
        ode.maximumStep = options.outputStep;
    const auto solution = integrate(system, initialField, 1, ode, {}, {}, execution);
    PDEIntegrationResult result;
    result.time = solution.time;
    result.values = solution.values;
    result.fieldShape = fieldShape;
    result.sampleCount = solution.sampleCount;
    result.fieldValueCount = values;
    result.acceptedSteps = solution.statistics.acceptedSteps;
    return result;
}

}  // namespace dynsyslab::core

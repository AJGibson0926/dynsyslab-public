#include <dynsyslab/core/DimensionalSystem.hpp>

#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
namespace {

class DimensionalSystem final : public ODESystem {
public:
    explicit DimensionalSystem(std::unique_ptr<ODESystem> model, const SystemDescriptor* published = nullptr)
        : model_(std::move(model)), scale_(model_->coordinates()),
          descriptor_(published ? *published : model_->descriptor()), configuration_(model_->configuration()) {
        const auto& states = scale_.physicalStatePerModelState;
        descriptor_.preferredCoordinates = "dimensional";
        if (published) configuration_.identifier = published->identifier;
        else
            for (std::size_t index = 0; index < descriptor_.stateLabels.size() && index < states.size(); ++index)
                descriptor_.stateLabels[index] += " [" + scale_.physicalStateUnits[index] + "]";
        configuration_.coordinateScaling = {1.0, 1.0, true, true, "s", scale_.physicalInputUnits,
                                            scale_.physicalInputUnits, {}, {}};
        configuration_.coordinateScaling.physicalStatePerModelState.assign(states.size(), 1.0);
        configuration_.coordinateScaling.physicalStateUnits = scale_.physicalStateUnits;
    }

    void resetRuntimeState() override { model_->resetRuntimeState(); }
    const SystemDescriptor& descriptor() const override { return descriptor_; }
    const SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return physical(model_->initialState()); }
    std::vector<double> controlledInitialState() const override {
        return physical(model_->controlledInitialState());
    }
    void rhs(double time, const double* state, double* derivative) const override {
        const auto x = model(state);
        model_->rhs(time / scale_.secondsPerModelTime, x.data(), derivative);
        rate(derivative);
    }
    void rhsWithInput(double time, const double* state, double input, double* derivative) const override {
        const auto x = model(state);
        model_->rhsWithInput(time / scale_.secondsPerModelTime, x.data(),
                             input / scale_.physicalInputPerModelInput, derivative);
        rate(derivative);
    }
    std::vector<double> jacobian(double time, const double* state) const override {
        const auto x = model(state);
        auto matrix = model_->jacobian(time / scale_.secondsPerModelTime, x.data());
        const auto& s = scale_.physicalStatePerModelState;
        const std::size_t n = s.size();
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = 0; c < n; ++c)
                matrix[r * n + c] *= s[r] / (s[c] * scale_.secondsPerModelTime);
        return matrix;
    }
    std::vector<double> inputMatrix(const double* state) const override {
        const auto x = model(state);
        auto matrix = model_->inputMatrix(x.data());
        const auto& s = scale_.physicalStatePerModelState;
        const std::size_t n = s.size(), m = n ? matrix.size() / n : 0;
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = 0; c < m; ++c)
                matrix[r * m + c] *= s[r] / (scale_.secondsPerModelTime * scale_.physicalInputPerModelInput);
        return matrix;
    }
    bool invalid(double time, const double* state) const override {
        const auto x = model(state);
        return model_->invalid(time / scale_.secondsPerModelTime, x.data());
    }
    DynamicalKnowledge dynamicalKnowledge() const override { return model_->dynamicalKnowledge(); }
    std::vector<AnalyticSolutionFamily> analyticSolutions() const override { return {}; }

private:
    std::vector<double> model(const double* state) const {
        const auto& s = scale_.physicalStatePerModelState;
        std::vector<double> result(s.size());
        for (std::size_t index = 0; index < s.size(); ++index) result[index] = state[index] / s[index];
        return result;
    }
    std::vector<double> physical(std::vector<double> state) const {
        const auto& s = scale_.physicalStatePerModelState;
        for (std::size_t index = 0; index < s.size() && index < state.size(); ++index) state[index] *= s[index];
        return state;
    }
    void rate(double* derivative) const {
        const auto& s = scale_.physicalStatePerModelState;
        for (std::size_t index = 0; index < s.size(); ++index)
            derivative[index] *= s[index] / scale_.secondsPerModelTime;
    }

    std::unique_ptr<ODESystem> model_;
    CoordinateScaling scale_;
    SystemDescriptor descriptor_;
    SystemConfiguration configuration_;
};

}  // namespace

bool supportsDimensionalIntegration(const ODESystem& system) {
    const auto scale = system.coordinates();
    return scale.physicalTimeAvailable && scale.physicalInputAvailable
        && scale.physicalStatePerModelState.size() == system.descriptor().stateDimension
        && scale.physicalStateUnits.size() == scale.physicalStatePerModelState.size();
}

std::unique_ptr<ODESystem> dimensionalSystem(std::unique_ptr<ODESystem> model) {
    if (!model || !supportsDimensionalIntegration(*model))
        throw std::invalid_argument("this system declares no physical scales for dimensional integration");
    return std::make_unique<DimensionalSystem>(std::move(model));
}

std::unique_ptr<ODESystem> dimensionalSystem(std::unique_ptr<ODESystem> model,
                                             const SystemDescriptor& descriptor) {
    if (!model || !supportsDimensionalIntegration(*model))
        throw std::invalid_argument("this system declares no physical scales for dimensional integration");
    if (descriptor.stateDimension != model->descriptor().stateDimension)
        throw std::invalid_argument("dimensional descriptor does not match the model's state dimension");
    return std::make_unique<DimensionalSystem>(std::move(model), &descriptor);
}

std::unique_ptr<ODESystem> systemInCoordinates(std::unique_ptr<ODESystem> model, const std::string& requested) {
    const auto& own = model->descriptor().preferredCoordinates;
    if (requested.empty() || requested == "automatic" || requested == own) return model;
    if (requested == "dimensional" && own == "dimensionless") return dimensionalSystem(std::move(model));
    throw std::invalid_argument("native coordinate conversion from " + own + " to " + requested
                                + " is not available for this system");
}

}  // namespace dynsyslab::core

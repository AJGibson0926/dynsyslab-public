#pragma once

#include <dynsyslab/core/EstimationRuntime.hpp>
#include <dynsyslab/core/LQR.hpp>

namespace dynsyslab::workflow {

struct KoopmanFeedbackSample {
    core::KoopmanEstimate estimate;
    std::vector<double> measured;
    std::vector<double> previousAppliedInput;
    std::vector<double> commandedInput;
};

// Composition only: Estimator owns prediction/correction; Controller owns
// feedback. This object never evaluates an eigenfunction at a plant state.
// A sensor/actuator adapter supplies y and the previous actually applied u.
class KoopmanFeedback {
public:
    KoopmanFeedback(const core::KoopmanModel& model, core::KoopmanLQEDesign design,
        std::string controlRepresentation, std::vector<double> initialLifted,
        std::vector<double> Q, std::vector<double> R,
        core::StateTargetFunction target = {},
        core::KoopmanTargetSpace targetSpace = core::KoopmanTargetSpace::realLiftedValues,
        std::vector<double> measurementOffset = {}, std::vector<double> outputOffset = {});
    KoopmanFeedbackSample sample(double time, const std::vector<double>& measurement,
        const std::vector<double>& previousAppliedInput);
    const core::KoopmanModel& model() const noexcept { return *model_; }
    void setMeasurementFeedthrough(std::vector<double> feedthrough, int scaleOutput = -1) {
        observer_.setMeasurementFeedthrough(std::move(feedthrough), scaleOutput);
    }
    // Reduced-order KLQR gain (controller.settings.balanced_truncation).
    void setBalancedReduction(core::BalancedReductionOptions options) {
        controller_.setBalancedReduction(options);
    }
    const core::KoopmanLQRController& controller() const noexcept { return controller_; }
    core::KoopmanLQRController& controller() noexcept { return controller_; }

private:
    std::shared_ptr<const core::KoopmanModel> model_;
    core::KoopmanLQEObserver observer_;
    core::KoopmanLQRController controller_;
};
} // namespace dynsyslab::workflow

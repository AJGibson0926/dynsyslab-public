#include <dynsyslab/workflow/KoopmanFeedback.hpp>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow {
KoopmanFeedback::KoopmanFeedback(const core::KoopmanModel& model,
    core::KoopmanLQEDesign design,std::string representation,
    std::vector<double> initial,std::vector<double> Q,std::vector<double> R,
    core::StateTargetFunction target,core::KoopmanTargetSpace targetSpace,
    std::vector<double> measurementOffset,std::vector<double> outputOffset)
    :model_(std::make_shared<const core::KoopmanModel>(model)),
     observer_(*model_,std::move(design),representation,std::move(initial),
        std::move(measurementOffset),std::move(outputOffset)),
     controller_(model_,std::move(Q),std::move(R),std::move(target),
        std::move(representation),targetSpace) {}

KoopmanFeedbackSample KoopmanFeedback::sample(double time,
    const std::vector<double>& measurement,const std::vector<double>& previousApplied) {
    KoopmanFeedbackSample result;
    result.estimate=observer_.observe(time,measurement,previousApplied);
    result.measured=measurement;
    result.previousAppliedInput=previousApplied;
    result.commandedInput=controller_.valueFromLifted(
        time,model_->valuesFromSplit(result.estimate.lifted));
    for(double u:result.commandedInput)if(!std::isfinite(u))
        throw std::runtime_error("KLQR produced a nonfinite input");
    return result;
}
} // namespace dynsyslab::workflow

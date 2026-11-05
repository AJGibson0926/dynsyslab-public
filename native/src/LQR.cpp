#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <dynsyslab/numerics/Control.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
std::vector<double> referenceValues(const KoopmanModel& model,const StateTargetFunction& target,
    KoopmanTargetSpace space,double time) {
    if(!std::isfinite(time))throw std::invalid_argument("Controller time must be finite");
    if(!target)return std::vector<double>(model.splitDimension());
    const auto values=target(time);
    const auto dimension=space==KoopmanTargetSpace::realLiftedValues?model.splitDimension():
        space==KoopmanTargetSpace::observableCoordinates?model.coordinateDimension():model.stateDimension();
    if(values.size()!=dimension)throw std::invalid_argument("Koopman target dimension differs from declared target space");
    for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Koopman target must be finite");
    if(space==KoopmanTargetSpace::realLiftedValues)return values;
    return space==KoopmanTargetSpace::observableCoordinates?model.evaluateCoordinatesSplit(values.data()):
        model.evaluateSplit(values.data());
}
}

ContinuousLQRSolution solveContinuousLQR(
    const std::vector<double>& generator,const std::vector<double>& input,
    const std::vector<double>& statePenalty,
    const std::vector<double>& inputPenalty,
    std::size_t states,std::size_t inputs) {
    auto solved=numerics::solveContinuousLqr(
        generator,input,statePenalty,inputPenalty,states,inputs);
    ContinuousLQRSolution result;
    result.riccati=std::move(solved.riccati);
    result.gain=std::move(solved.gain);
    result.closedLoopPoles=std::move(solved.closedLoopPoles);
    result.relativeResidual=solved.relativeResidual;
    result.stateDimension=states;
    result.inputDimension=inputs;
    return result;
}

KoopmanLQRController::KoopmanLQRController(
    std::shared_ptr<const KoopmanModel> model,std::vector<double> statePenalty,
    std::vector<double> inputPenalty,StateTargetFunction target,
    std::string controlRepresentation,KoopmanTargetSpace targetSpace)
    : model_(std::move(model)),statePenalty_(std::move(statePenalty)),
      inputPenalty_(std::move(inputPenalty)),target_(std::move(target)),
      targetSpace_(targetSpace),representation_(std::move(controlRepresentation)) {
    if(!model_)throw std::invalid_argument("KLQR requires a KROM");
    requireContinuousControlModel(*model_,representation_);
    if(statePenalty_.size()!=model_->splitDimension()*model_->splitDimension()
        ||inputPenalty_.size()!=model_->inputDimension()*model_->inputDimension())
        throw std::invalid_argument("KLQR Q/R dimensions do not match KROM");
}

std::vector<double> KoopmanLQRController::value(
    double time,const double* modelState) const {
    const auto lifted=model_->evaluateSplit(modelState);
    return feedback(time,lifted,model_->controlMapSplit(modelState,representation_));
}

std::vector<double> KoopmanLQRController::valueFromCoordinates(double time,const double* q) const {
    return valueFromLifted(time,model_->evaluateCoordinates(q));
}
std::vector<double> KoopmanLQRController::valueFromLifted(double time,
    const std::vector<std::complex<double>>& values) const {
    const auto lifted=model_->splitValues(values);
    return feedback(time,lifted,model_->realControlMap(model_->controlMapFromLifted(values,representation_)));
}
std::vector<double> KoopmanLQRController::feedback(double time,const std::vector<double>& lifted,
    const std::vector<double>& inputMap) const {
    const auto reference=referenceValues(*model_,target_,targetSpace_,time);
    const bool constant=representation_!="exact"&&
        std::holds_alternative<ConstantControlMap>(model_->storedControlMap(representation_));
    // Compare B as well as its representation: a mutable alias to the model may
    // replace a constant map after construction. Never reuse a stale gain.
    if(!constant||last_.gain.empty()||cachedInputMap_!=inputMap) {
        last_=solveContinuousLQR(model_->generatorSplit(),inputMap,
            statePenalty_,inputPenalty_,model_->splitDimension(),model_->inputDimension());
        cachedInputMap_=constant?inputMap:std::vector<double>{};
        ++gainSolveCount_;
    }
    std::vector<double> result(model_->inputDimension());
    for(std::size_t input=0;input<model_->inputDimension();++input)
        for(std::size_t state=0;state<lifted.size();++state)
            result[input]-=last_.gain[input*lifted.size()+state]
                *(lifted[state]-reference[state]);
    return result;
}

const ContinuousLQRSolution& KoopmanLQRController::lastSolution() const noexcept {
    return last_;
}

const std::string& KoopmanLQRController::controlRepresentation() const noexcept {
    return representation_;
}

KoopmanMPCController::KoopmanMPCController(
    std::shared_ptr<const KoopmanModel> model,std::vector<double> statePenalty,
    std::vector<double> inputPenalty,StateTargetFunction target,std::size_t horizon,
    double sampleTime,std::string representation,std::vector<double> limits,KoopmanTargetSpace targetSpace)
    : model_(std::move(model)),statePenalty_(std::move(statePenalty)),
      inputPenalty_(std::move(inputPenalty)),target_(std::move(target)),
      targetSpace_(targetSpace),horizon_(horizon),sampleTime_(sampleTime),
      representation_(std::move(representation)),controlLimits_(std::move(limits)) {
    if(!model_||!model_->hasControlMap(representation_)||!horizon_
        ||!(sampleTime_>0.0)||!std::isfinite(sampleTime_))
        throw std::invalid_argument("invalid Koopman MPC configuration");
    requireContinuousControlModel(*model_,representation_);
    const std::size_t states=model_->splitDimension();
    const std::size_t inputs=model_->inputDimension();
    if(statePenalty_.size()!=states*states||inputPenalty_.size()!=inputs*inputs
        ||(!controlLimits_.empty()&&controlLimits_.size()!=1
            &&controlLimits_.size()!=inputs))
        throw std::invalid_argument("Koopman MPC matrix dimensions are inconsistent");
    for(const double limit:controlLimits_)
        if(!(limit>0.0)||!std::isfinite(limit))
            throw std::invalid_argument("Koopman MPC limits must be positive");
}

std::vector<double> KoopmanMPCController::value(
    double time,const double* state) const {
    return feedback(time,model_->evaluateSplit(state),model_->controlMapSplit(state,representation_));
}
std::vector<double> KoopmanMPCController::valueFromCoordinates(double time,const double* q) const {
    return valueFromLifted(time,model_->evaluateCoordinates(q));
}
std::vector<double> KoopmanMPCController::valueFromLifted(double time,
    const std::vector<std::complex<double>>& values) const {
    const auto lifted=model_->splitValues(values);
    return feedback(time,lifted,model_->realControlMap(model_->controlMapFromLifted(values,representation_)));
}
std::vector<double> KoopmanMPCController::feedback(double time,const std::vector<double>& lifted,
    const std::vector<double>& inputMap) const {
    const std::size_t states=model_->splitDimension();
    const std::size_t inputs=model_->inputDimension();
    auto solved=numerics::solveFiniteHorizonLqr(
        model_->generatorSplit(),inputMap,
        statePenalty_,inputPenalty_,states,inputs,horizon_,sampleTime_);
    lastGain_=std::move(solved.gain);
    const auto reference=referenceValues(*model_,target_,targetSpace_,time);
    std::vector<double> result(inputs);
    for(std::size_t input=0;input<inputs;++input) {
        for(std::size_t index=0;index<states;++index)
            result[input]-=lastGain_[input*states+index]
                *(lifted[index]-reference[index]);
        if(!controlLimits_.empty()) {
            const double limit=controlLimits_.size()==1
                ?controlLimits_.front():controlLimits_[input];
            result[input]=std::clamp(result[input],-limit,limit);
        }
    }
    return result;
}

const std::vector<double>& KoopmanMPCController::lastGain() const noexcept {
    return lastGain_;
}

}  // namespace dynsyslab::core

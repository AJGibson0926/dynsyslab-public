#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/ReducedOrderControl.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <dynsyslab/numerics/Control.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

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
    // One gain path for every LQR in DynSysLab: ControlDesign (CARE/DARE).
    auto solved=designLQR(makeStateSpace(generator,input,{},{},states,inputs,0),
        statePenalty,inputPenalty);
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
    if(representation_=="external") {
        if(!model_->spectrum().hasContinuousRates())
            throw std::invalid_argument("Continuous controller requires finite continuous rates");
    } else requireContinuousControlModel(*model_,representation_);
    if(statePenalty_.size()!=model_->splitDimension()*model_->splitDimension()
        ||inputPenalty_.size()!=model_->inputDimension()*model_->inputDimension())
        throw std::invalid_argument("KLQR Q/R dimensions do not match KROM");
}

std::vector<double> KoopmanLQRController::valueFromLiftedWithMap(double time,
    const std::vector<std::complex<double>>& values,
    const std::vector<std::complex<double>>& inputMap) const {
    if(inputMap.size()!=model_->liftedDimension()*model_->inputDimension())
        throw std::invalid_argument("External KLQR input map has the wrong dimension");
    return feedback(time,model_->splitValues(values),model_->realControlMap(inputMap));
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
ContinuousLQRSolution solveContinuousLQI(
    const std::vector<double>& A,const std::vector<double>& B,const std::vector<double>& C,
    const std::vector<double>& Q,const std::vector<double>& Qi,const std::vector<double>& R,
    std::size_t n,std::size_t m,std::size_t p) {
    if(A.size()!=n*n||B.size()!=n*m||C.size()!=p*n||Q.size()!=n*n||Qi.size()!=p*p||R.size()!=m*m||!p)
        throw std::invalid_argument("LQI dimensions are inconsistent");
    const std::size_t N=n+p;
    std::vector<double> Aa(N*N,0.0),Ba(N*m,0.0),Qa(N*N,0.0);
    for(std::size_t i=0;i<n;++i){
        for(std::size_t j=0;j<n;++j){Aa[i*N+j]=A[i*n+j];Qa[i*N+j]=Q[i*n+j];}
        for(std::size_t a=0;a<m;++a)Ba[i*m+a]=B[i*m+a];
    }
    for(std::size_t r=0;r<p;++r){
        for(std::size_t j=0;j<n;++j)Aa[(n+r)*N+j]=C[r*n+j];  // q' = C x (+ reference terms)
        for(std::size_t s=0;s<p;++s)Qa[(n+r)*N+n+s]=Qi[r*p+s];
    }
    return solveContinuousLQR(Aa,Ba,Qa,R,N,m);
}

namespace {
void checkChannel(const LiftedOutputChannel& c,std::size_t n,std::size_t targetDimension,const char* name) {
    const auto p=c.targetComponents.size();
    if(!p||c.outputMatrix.size()!=p*n||c.outputOffset.size()!=p)
        throw std::invalid_argument(std::string(name)+": output channel dimensions are inconsistent");
    for(auto k:c.targetComponents)if(k>=targetDimension)
        throw std::invalid_argument(std::string(name)+": target component is out of range");
}
std::vector<double> channelError(const LiftedOutputChannel& c,const std::vector<double>& lifted,
    const std::vector<double>& target) {  // y(z) - y_target
    const auto p=c.targetComponents.size(),n=lifted.size();
    std::vector<double> e(p);
    for(std::size_t r=0;r<p;++r){
        double y=c.outputOffset[r];
        for(std::size_t j=0;j<n;++j)y+=c.outputMatrix[r*n+j]*lifted[j];
        e[r]=y-target[c.targetComponents[r]];
    }
    return e;
}
}

void KoopmanLQRController::setIntegralAction(KoopmanIntegralAction action) {
    if(certified_)throw std::invalid_argument("KLQR integral action cannot be added to a certified gain");
    if(targetSpace_!=KoopmanTargetSpace::dynamicalState||!target_)
        throw std::invalid_argument("KLQR integral action needs a dynamical-state target");
    if(reduction_)throw std::invalid_argument("KLQR integral action and balanced_truncation cannot be combined");
    const auto n=model_->splitDimension(),p=action.channel.targetComponents.size();
    checkChannel(action.channel,n,model_->stateDimension(),"KLQR integral action");
    if(action.penalty.size()!=p*p)throw std::invalid_argument("KLQR integral action penalty must be outputs x outputs");
    if(!action.integratorLimit.empty()&&action.integratorLimit.size()!=p)
        throw std::invalid_argument("KLQR integral action limit needs one value per output");
    integral_=std::move(action);integralState_.assign(p,0.0);integralStarted_=false;
    last_={};cachedInputMap_.clear();
}

void KoopmanLQRController::setPIDCorrection(KoopmanPIDCorrection c) {
    if(certified_)throw std::invalid_argument("KLQR PID correction cannot be added to a certified gain");
    if(targetSpace_!=KoopmanTargetSpace::dynamicalState||!target_)
        throw std::invalid_argument("KLQR PID correction needs a dynamical-state target");
    const auto n=model_->splitDimension(),p=c.channel.targetComponents.size(),m=model_->inputDimension();
    checkChannel(c.channel,n,model_->stateDimension(),"KLQR PID correction");
    PIDOptions options;
    options.derivativeFilterTime=c.derivativeFilterTime;options.integratorLimit=c.integratorLimit;
    options.outputMinimum=c.outputMinimum;options.outputMaximum=c.outputMaximum;
    options.conditionalIntegration=c.conditionalIntegration;
    pid_.emplace(c.proportional,c.integral,c.derivative,p,m,options);
    pidChannel_=std::move(c.channel);
}

std::vector<double> KoopmanLQRController::feedback(double time,const std::vector<double>& lifted,
    const std::vector<double>& inputMap) const {
    const auto reference=referenceValues(*model_,target_,targetSpace_,time);
    if(certified_) {
        // u = u* + K (z - z*). A target that moves away from z* leaves the certificate.
        const auto& c=*certified_;
        const auto n=lifted.size(),m=model_->inputDimension();
        if(target_)for(std::size_t i=0;i<n;++i)
            if(std::abs(reference[i]-c.targetLifted[i])>1e-9*std::max(1.0,std::abs(c.targetLifted[i])))
                throw std::invalid_argument("KLQR target differs from the certified bilinear target");
        auto result=c.targetInput;
        for(std::size_t input=0;input<m;++input)for(std::size_t i=0;i<n;++i)
            result[input]+=c.gain[input*n+i]*(lifted[i]-c.targetLifted[i]);
        return result;
    }
    const bool constant=representation_!="exact"&&representation_!="external"&&
        std::holds_alternative<ConstantControlMap>(model_->storedControlMap(representation_));
    // Compare B as well as its representation: a mutable alias to the model may
    // replace a constant map after construction. Never reuse a stale gain.
    if(!constant||last_.gain.empty()||cachedInputMap_!=inputMap) {
        solveGain(inputMap);
        cachedInputMap_=constant?inputMap:std::vector<double>{};
        ++gainSolveCount_;
    }
    std::vector<double> result(model_->inputDimension());
    const std::size_t columns=lifted.size()+(integral_?integralState_.size():0);
    for(std::size_t input=0;input<model_->inputDimension();++input)
        for(std::size_t state=0;state<lifted.size();++state)
            result[input]-=last_.gain[input*columns+state]
                *(lifted[state]-reference[state]);
    if(integral_||pid_) {
        const auto target=target_(time);
        if(integral_) {
            // q' = y(z) - y_target, rectangular rule on the sampled clock.
            const auto e=channelError(integral_->channel,lifted,target);
            if(integralStarted_) {
                if(time<integralTime_)throw std::invalid_argument("KLQR integral action needs monotonic sample times");
                const double dt=time-integralTime_;
                for(std::size_t r=0;r<e.size();++r) {
                    integralState_[r]+=dt*e[r];
                    if(!integral_->integratorLimit.empty())integralState_[r]=std::clamp(integralState_[r],
                        -integral_->integratorLimit[r],integral_->integratorLimit[r]);
                }
            }
            integralTime_=time;integralStarted_=true;
            for(std::size_t input=0;input<model_->inputDimension();++input)
                for(std::size_t r=0;r<integralState_.size();++r)
                    result[input]-=last_.gain[input*columns+lifted.size()+r]*integralState_[r];
        }
        if(pid_) {
            auto e=channelError(*pidChannel_,lifted,target);
            for(double& v:e)v=-v;  // PID acts on y_target - y
            const auto correction=pid_->update(time,e);
            for(std::size_t input=0;input<result.size();++input)result[input]+=correction[input];
        }
    }
    return result;
}

void KoopmanLQRController::setCertifiedFeedback(CertifiedKoopmanFeedback feedback) {
    const auto n=model_->splitDimension(),m=model_->inputDimension();
    if(integral_||pid_||reduction_)
        throw std::invalid_argument("A certified gain cannot be combined with integral action, PID or balanced_truncation");
    if(feedback.gain.size()!=m*n||feedback.targetLifted.size()!=n||feedback.targetInput.size()!=m
        ||!(feedback.sampleTime>0.0)||!std::isfinite(feedback.sampleTime))
        throw std::invalid_argument("Certified KLQR feedback has inconsistent dimensions or sample time");
    const auto finite=[](const std::vector<double>& v){
        return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);});};
    if(!finite(feedback.gain)||!finite(feedback.targetLifted)||!finite(feedback.targetInput))
        throw std::invalid_argument("Certified KLQR feedback must be finite");
    certified_=std::move(feedback);
}

void KoopmanLQRController::setBalancedReduction(BalancedReductionOptions options) {
    if(certified_)throw std::invalid_argument("balanced_truncation cannot replace a certified gain");
    if(integral_)throw std::invalid_argument("KLQR integral action and balanced_truncation cannot be combined");
    reduction_=options;
    last_={};cachedInputMap_.clear();
}

void KoopmanLQRController::solveGain(const std::vector<double>& inputMap) const {
    const auto n=model_->splitDimension(),m=model_->inputDimension();
    if(integral_) {
        const auto p=integral_->channel.targetComponents.size();
        last_=solveContinuousLQI(model_->generatorSplit(),inputMap,integral_->channel.outputMatrix,
            statePenalty_,integral_->penalty,inputPenalty_,n,m,p);
        return;
    }
    if(!reduction_) {
        last_=solveContinuousLQR(model_->generatorSplit(),inputMap,statePenalty_,inputPenalty_,n,m);
        return;
    }
    std::vector<double> identity(n*n,0.0);
    for(std::size_t i=0;i<n;++i)identity[i*n+i]=1.0;
    const auto plant=makeStateSpace(model_->generatorSplit(),inputMap,identity,{},n,m,n);
    auto reduced=designReducedLQR(plant,*reduction_,statePenalty_,inputPenalty_);
    if(!reduced.liftedClosedLoopStable)
        throw std::runtime_error("balanced_truncation: the order-"+std::to_string(reduced.reduction.order)
            +" KLQR gain does not stabilize the full KROM; raise balanced_truncation.order");
    ContinuousLQRSolution result;
    result.gain=std::move(reduced.liftedGain);
    result.closedLoopPoles=std::move(reduced.liftedClosedLoopPoles);
    result.relativeResidual=reduced.regulator.relativeResidual;
    result.stateDimension=n;result.inputDimension=m;
    last_=std::move(result);
    lastReducedOrder_=reduced.reduction.order;
    lastReductionErrorBound_=reduced.reduction.errorBound;
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

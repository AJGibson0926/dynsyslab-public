#include <dynsyslab/core/OnlineShooting.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
OnlineShootingObserver::OnlineShootingObserver(ODESystem& model,std::vector<double> variances,
    std::vector<ShootingStart> starts,InitialStateMap initialize,InputAwareReadout read,
    OnlineShootingOptions options,IntegrationOptions solver,numerics::NonlinearLeastSquaresOptions optimizer,
    ShootingStartProvider provider)
    :model_(model),variances_(std::move(variances)),starts_(std::move(starts)),original_(starts_),
     initialize_(std::move(initialize)),read_(std::move(read)),options_(options),solver_(solver),
     optimizer_(optimizer),provider_(std::move(provider)){
    if(!initialize_||!read_||(!provider_&&starts_.empty())||variances_.empty()||
       !std::isfinite(options_.sampleTime)||options_.sampleTime<=0||options_.firstFitIndex<1||
       options_.fitEvery<1||options_.maximumSamples<2||options_.maximumSamples>1000000||
       options_.firstFitIndex>=options_.maximumSamples||variances_.size()>10000000/options_.maximumSamples)
        throw std::invalid_argument("Invalid online shooting configuration");
    for(double v:variances_)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Positive online measurement variances required");
    solver_.truncateOnFailure=false;
}
OnlineShootingSample OnlineShootingObserver::observe(double time,const std::vector<double>& y,double previous){
    const auto k=time_.size();
    if(!std::isfinite(time)||!std::isfinite(previous)||y.size()!=variances_.size()||k>=options_.maximumSamples)
        throw std::invalid_argument("Invalid online sensor sample or capacity exceeded");
    for(double v:y)if(!std::isfinite(v))throw std::invalid_argument("Finite online measurements required");
    if(k){
        const double expected=time_.front()+k*options_.sampleTime;
        if(time<=time_.back()||std::abs(time-expected)>1e-10*std::max(1.,std::abs(expected)))
            throw std::invalid_argument("Online observer requires one sample per uniform clock tick");
        if(!state_.empty()){
            auto options=solver_;options.initialTime=time_.back();options.finalTime=time;options.outputStep=time-time_.back();
            const auto prediction=integrateKnownInput(model_,state_,options,[previous](double){return previous;});
            const auto x=prediction.state(0,prediction.sampleCount-1);state_.assign(x,x+prediction.stateDimension);
        }
        commands_.back()=previous; // command after the preceding sensor sample
    }else beforeFirst_=previous;
    time_.push_back(time);observed_.insert(observed_.end(),y.begin(),y.end());commands_.push_back(previous);
    OnlineShootingSample sample;sample.index=k;sample.time=time;
    if(k>=options_.firstFitIndex&&(k-options_.firstFitIndex)%options_.fitEvery==0){
        if(!informationTime_&&provider_){starts_=provider_(time_,observed_);original_=starts_;}
        const auto input=KnownScalarInput::recorded(time_,commands_,beforeFirst_);
        auto fit=identifyByShooting(model_,time_,observed_,variances_,starts_,initialize_,read_,input,solver_,optimizer_);
        auto warm=starts_.at(fit.bestStart);warm.initial=fit.best.parameters;
        auto options=solver_;options.initialTime=time_.front();options.finalTime=time;options.outputStep=time-time_.front();
        const auto fitted=integrateKnownInput(model_,initialize_(fit.best.parameters),options,input);
        const auto x=fitted.state(0,fitted.sampleCount-1);state_.assign(x,x+fitted.stateDimension);
        starts_=options_.retainOriginalStarts?original_:std::vector<ShootingStart>{};starts_.insert(starts_.begin(),std::move(warm));
        if(starts_.size()>256)starts_.resize(256);
        informationTime_=time;sample.fit=std::move(fit);
    }
    sample.ready=informationTime_.has_value();sample.informationTime=informationTime_;sample.state=state_;return sample;
}
NativeSolution OnlineShootingObserver::forecast(double finalTime,double input,std::size_t intervals) const {
    if(!informationTime_||time_.empty()||!std::isfinite(finalTime)||finalTime<=time_.back()||
       !std::isfinite(input)||intervals<1||intervals>1000000)
        throw std::invalid_argument("Forecast needs an acquired estimate, finite held input and later endpoint");
    auto options=solver_;options.initialTime=time_.back();options.finalTime=finalTime;
    options.outputStep=std::nextafter((finalTime-time_.back())/double(intervals),std::numeric_limits<double>::infinity());
    return integrateKnownInput(model_,state_,options,[input](double){return input;});
}
}

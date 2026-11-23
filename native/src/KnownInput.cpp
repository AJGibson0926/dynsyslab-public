#include <dynsyslab/core/KnownInput.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
namespace dynsyslab::core {
KnownScalarInput KnownScalarInput::recorded(std::vector<double> time,std::vector<double> values,
    double beforeFirst,bool beforeCommand){
    if(time.size()<2||time.size()>1000000||time.size()!=values.size()||!std::isfinite(beforeFirst))
        throw std::invalid_argument("Recorded input requires matching bounded time/value arrays and finite initial input");
    for(std::size_t j=0;j<time.size();++j)
        if(!std::isfinite(time[j])||!std::isfinite(values[j])||(j&&time[j]<=time[j-1]))
            throw std::invalid_argument("Recorded input time must increase and all values must be finite");
    KnownScalarInput result;result.breaks_.push_back(time.front());
    for(std::size_t k=1;k+1<time.size();++k)if(values[k]!=values[k-1])result.breaks_.push_back(time[k]);
    result.breaks_.push_back(time.back());
    const auto t=std::make_shared<const std::vector<double>>(std::move(time));
    const auto u=std::make_shared<const std::vector<double>>(std::move(values));
    const auto index=[t](double x,bool left){
        if(!std::isfinite(x)||x<t->front()||x>t->back())throw std::out_of_range("Input replay outside recorded interval");
        auto it=left?std::lower_bound(t->begin(),t->end(),x):std::upper_bound(t->begin(),t->end(),x);
        return std::size_t(it-t->begin());
    };
    result.applied_=[u,index](double x){return u->at(index(x,false)-1);};
    result.measured_=[u,index,beforeFirst,beforeCommand](double x){const auto i=index(x,beforeCommand);return i?u->at(i-1):beforeFirst;};
    return result;
}
double KnownScalarInput::operator()(double time) const {
    const double v=applied_(time);if(!std::isfinite(v))throw std::domain_error("Nonfinite applied input");return v;
}
double KnownScalarInput::measurement(double time) const {
    const double v=measured_(time);if(!std::isfinite(v))throw std::domain_error("Nonfinite sensor input");return v;
}
NativeSolution integrateKnownInput(ODESystem& system,const std::vector<double>& initial,
    const IntegrationOptions& options,const KnownScalarInput& input){
    if(!input)throw std::invalid_argument("Known input required");
    if(input.breakpoints().empty())return integrate(system,initial,1,options,
        [&](std::size_t,double t,const double*){return input(t);});
    if(options.truncateOnFailure)throw std::invalid_argument("Recorded-input replay requires complete integration");
    if(!std::isfinite(options.initialTime)||!std::isfinite(options.finalTime)||options.finalTime<=options.initialTime||
       !std::isfinite(options.outputStep)||options.outputStep<=0||initial.size()!=system.descriptor().stateDimension)
        throw std::invalid_argument("Invalid replay clock or initial state");
    (void)input(options.initialTime);(void)input(options.finalTime);
    const long double intervals=(static_cast<long double>(options.finalTime)-options.initialTime)/options.outputStep;
    const long double nearest=std::round(intervals);
    const long double adjusted=std::abs(intervals-nearest)<=16*std::numeric_limits<double>::epsilon()*std::max(1.L,intervals)?nearest:intervals;
    const long double count=std::ceil(adjusted)+1;
    const auto dimension=initial.size();
    if(count>10000000||dimension==0||count*dimension>10000000)throw std::invalid_argument("Replay output resource bound exceeded");
    NativeSolution result;result.trajectoryCount=1;result.stateDimension=dimension;result.sampleCount=std::size_t(count);
    for(std::size_t k=0;k<result.sampleCount;++k)result.time.push_back(std::min(options.initialTime+k*options.outputStep,options.finalTime));
    result.time.back()=options.finalTime; // decimal clock roundoff must not duplicate the endpoint
    result.values=initial;result.controlValues.push_back(input(options.initialTime));
    std::vector<double> ends;for(double t:input.breakpoints())if(t>options.initialTime&&t<options.finalTime)ends.push_back(t);
    ends.push_back(options.finalTime);
    auto state=initial;double current=options.initialTime;std::size_t next=1;
    auto advance=[&](double stop,double step,double held){
        auto local=options;local.initialTime=current;local.finalTime=stop;local.outputStep=step;
        const auto piece=integrate(system,state,1,local,[held](std::size_t,double,const double*){return held;});
        result.coordinateSystem=piece.coordinateSystem;
        result.statistics.acceptedSteps+=piece.statistics.acceptedSteps;result.statistics.rejectedSteps+=piece.statistics.rejectedSteps;
        result.statistics.rightHandSideCalls+=piece.statistics.rightHandSideCalls;result.statistics.controllerCalls+=piece.statistics.controllerCalls;
        for(std::size_t j=1;j<piece.sampleCount&&next<result.time.size();++j){
            const auto t=piece.time[j],target=result.time[next];
            const double tolerance=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(t),std::abs(target)});
            if(std::abs(t-target)<=tolerance){const auto x=piece.state(0,j);result.values.insert(result.values.end(),x,x+dimension);
                result.controlValues.push_back(input(target));++next;}
        }
        const auto last=piece.state(0,piece.sampleCount-1);state.assign(last,last+dimension);current=stop;
    };
    for(double end:ends){
        const double held=input(current);
        // Align the local output grid after an off-grid switch; then emit the
        // remaining global samples in one solver call for this held interval.
        if(next<result.time.size()&&result.time[next]<end){
            const auto first=result.time[next];advance(first,first-current,held);
        }
        if(current<end)advance(end,options.outputStep,held);
    }
    if(next!=result.sampleCount)throw std::runtime_error("Recorded-input replay output clock mismatch");
    result.completedSamples={result.sampleCount};result.terminationReasons={""};return result;
}
}

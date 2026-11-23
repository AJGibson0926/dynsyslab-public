#include <dynsyslab/core/SequentialShooting.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
SequentialShootingResult sequentialShooting(ODESystem& system,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<double>& variances,const std::vector<ShootingStart>& original,
    const InitialStateMap& initialize,const InputAwareReadout& read,const KnownScalarInput& input,
    const SequentialShootingOptions& o,IntegrationOptions solver,
    const numerics::NonlinearLeastSquaresOptions& optimizer){
    const auto nm=variances.size(),dimension=system.descriptor().stateDimension;
    if(times.size()<2||times.size()>1000000||nm==0||nm>10000000/times.size()||
       measurements.size()!=times.size()*nm||o.updateIndices.empty()||o.updateIndices.size()>10000||
       o.predictionSubsteps==0||o.predictionSubsteps>10000||!dimension||
       times.size()>10000000/o.predictionSubsteps/dimension)
        throw std::invalid_argument("Invalid sequential shooting dimensions/resource bound");
    const double dt=(times.back()-times.front())/double(times.size()-1);
    if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("Invalid sequential observation step");
    for(std::size_t k=0;k<times.size();++k){
        const double expected=times.front()+k*dt;
        if(!std::isfinite(times[k])||(k&&times[k]<=times[k-1])||std::abs(times[k]-expected)>1e-10*std::max(1.,std::abs(expected)))
            throw std::invalid_argument("Sequential shooting needs a uniform increasing clock");
    }
    for(double y:measurements)if(!std::isfinite(y))throw std::invalid_argument("Finite sequential observations required");
    for(std::size_t j=0;j<o.updateIndices.size();++j){const auto k=o.updateIndices[j];
        if(k==0||k>=times.size()||(j&&k<=o.updateIndices[j-1]))throw std::invalid_argument("Update indices must increase and reference at least two past samples");
    }
    solver.truncateOnFailure=false;SequentialShootingResult result;result.stateDimension=dimension;
    auto starts=original;
    for(std::size_t j=0;j<o.updateIndices.size();++j){
        const auto k=o.updateIndices[j],end=j+1<o.updateIndices.size()?o.updateIndices[j+1]:times.size()-1;
        std::vector<double> past(times.begin(),times.begin()+k+1);
        std::vector<double> observed(measurements.begin(),measurements.begin()+(k+1)*nm);
        auto fit=identifyByShooting(system,past,observed,variances,starts,initialize,read,input,solver,optimizer);
        const auto sourceStart=starts.at(fit.bestStart);const auto initial=initialize(fit.best.parameters);
        // Obtain the current fitted state without exposing the reconstructed past
        // as if it had been a causal estimate at those earlier times.
        auto historyOptions=solver;historyOptions.initialTime=times.front();historyOptions.finalTime=times[k];
        historyOptions.outputStep=times[k]-times.front();
        const auto history=integrateKnownInput(system,initial,historyOptions,input);
        const auto terminal=history.state(0,history.sampleCount-1);
        std::vector<double> current(terminal,terminal+dimension);
        if(end==k){
            result.time.push_back(times[k]);result.informationTime.push_back(times[k]);
            result.states.insert(result.states.end(),current.begin(),current.end());
        }else{
            auto options=solver;options.initialTime=times[k];options.finalTime=times[end];
            const auto intervals=(end-k)*o.predictionSubsteps;
            options.outputStep=std::nextafter((times[end]-times[k])/double(intervals),std::numeric_limits<double>::infinity());
            const auto segment=integrateKnownInput(system,current,options,input);
            if(segment.sampleCount!=intervals+1)throw std::runtime_error("Sequential forecast clock mismatch");
            const auto count=j+1==o.updateIndices.size()?segment.sampleCount:segment.sampleCount-1;
            for(std::size_t i=0;i<count;++i){
                result.time.push_back(segment.time[i]);result.informationTime.push_back(times[k]);
                const auto state=segment.state(0,i);result.states.insert(result.states.end(),state,state+dimension);
            }
        }
        auto warm=sourceStart;warm.initial=fit.best.parameters;
        starts=o.retainOriginalStarts?original:std::vector<ShootingStart>{};starts.insert(starts.begin(),std::move(warm));
        // The initial multi-start limit is 256; adding a warm start retains the
        // best-informed candidate and drops only the final original if needed.
        if(starts.size()>256)starts.resize(256);
        result.updates.push_back({k,times[k],std::move(fit)});
    }
    return result;
}
}

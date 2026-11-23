#include <dynsyslab/core/ShootingIdentification.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
ShootingIdentificationResult identifyByShooting(ODESystem& system,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<double>& variances,const std::vector<ShootingStart>& starts,
    const InitialStateMap& initialize,const InputAwareReadout& read,
    const KnownScalarInput& input,IntegrationOptions solver,
    const numerics::NonlinearLeastSquaresOptions& optimizer){
    const auto nm=variances.size();
    if(!initialize||!read||!input||times.size()<2||times.size()>1000000||!nm||
       nm>10000000/times.size()||measurements.size()!=times.size()*nm||
       starts.empty()||starts.size()>256)throw std::invalid_argument("Invalid shooting dimensions/callbacks");
    const double dt=(times.back()-times.front())/double(times.size()-1);
    if(!std::isfinite(dt)||dt<=0)throw std::invalid_argument("Positive measurement step required");
    for(std::size_t j=0;j<times.size();++j){
        const double expected=times[0]+j*dt;
        if(!std::isfinite(times[j])||(j&&times[j]<=times[j-1])||
           std::abs(times[j]-expected)>1e-10*std::max({1.,std::abs(expected),dt}))
            throw std::invalid_argument("Shooting requires a finite uniform measurement clock");
    }
    for(double v:variances)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Positive measurement variances required");
    for(double v:measurements)if(!std::isfinite(v))throw std::invalid_argument("Finite measurements required");
    for(const auto& s:starts){
        if(s.initial.empty()||s.initial.size()!=s.lower.size()||s.initial.size()!=s.upper.size())
            throw std::invalid_argument("Invalid shooting start dimensions");
        for(std::size_t j=0;j<s.initial.size();++j)
            if(!std::isfinite(s.initial[j])||!std::isfinite(s.lower[j])||!std::isfinite(s.upper[j])||
               s.lower[j]>=s.upper[j]||s.initial[j]<s.lower[j]||s.initial[j]>s.upper[j])
                throw std::invalid_argument("Invalid shooting start bounds");
    }
    // Integrator uses ceil(span/step); round the normalized step upward by one
    // ULP to prevent a duplicate final sample from decimal-clock roundoff.
    solver.initialTime=times.front();solver.finalTime=times.back();solver.outputStep=std::nextafter(dt,std::numeric_limits<double>::infinity());
    solver.truncateOnFailure=false;
    ShootingIdentificationResult result;result.best.squaredResidual=std::numeric_limits<double>::infinity();
    auto objective=[&](const std::vector<double>& p){
        ++result.trajectoryEvaluations;
        const auto initial=initialize(p);
        if(initial.size()!=system.descriptor().stateDimension)throw std::invalid_argument("Shooting initial-state map dimension mismatch");
        const auto solution=integrateKnownInput(system,initial,solver,input);
        if(solution.sampleCount!=times.size())throw std::runtime_error("Shooting output clock mismatch");
        std::vector<double> residual;residual.reserve(measurements.size());
        for(std::size_t i=0;i<times.size();++i){
            if(std::abs(solution.time[i]-times[i])>1e-10*std::max(1.,std::abs(times[i])))throw std::runtime_error("Shooting sample times differ");
            const double u=input.measurement(times[i]);
            const auto y=read(times[i],solution.state(0,i),u);
            if(y.size()!=nm)throw std::invalid_argument("Shooting readout dimension mismatch");
            for(std::size_t c=0;c<nm;++c){
                const double r=(y[c]-measurements[i*nm+c])/std::sqrt(variances[c]);
                if(!std::isfinite(r))throw std::domain_error("Nonfinite shooting residual");residual.push_back(r);
            }
        }return residual;
    };
    for(const auto& start:starts){ShootingAttempt attempt;
        try{
            const auto r=objective(start.initial);double cost=0;for(double v:r)cost+=v*v;
            if(!std::isfinite(cost))throw std::overflow_error("Shooting objective overflow");
            attempt.fit.parameters=start.initial;attempt.fit.squaredResidual=cost;
            attempt.fit.reason="initial_candidate";attempt.valid=true;
            try{
                const auto fit=numerics::boundedNonlinearLeastSquares(objective,start.initial,start.lower,start.upper,optimizer);
                if(fit.squaredResidual<=cost)attempt.fit=fit;
            }catch(const std::exception& e){attempt.failure=e.what();attempt.fit.reason="optimizer_failed_initial_retained";}
        }catch(const std::exception& e){attempt.failure=e.what();}
        if(attempt.valid&&attempt.fit.squaredResidual<result.best.squaredResidual){result.best=attempt.fit;result.bestStart=result.attempts.size();}
        result.attempts.push_back(std::move(attempt));
    }
    if(!std::isfinite(result.best.squaredResidual)){
        throw std::runtime_error("No valid shooting candidate: "+result.attempts.front().failure);
    }
    return result;
}
}

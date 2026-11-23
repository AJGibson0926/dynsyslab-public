#include <dynsyslab/core/PeriodicEstimation.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
numerics::NonlinearLeastSquaresResult identifyPeriodicOrbit(const math::PeriodicAtlas& atlas,
    const std::vector<double>& times,const std::vector<double>& measurements,
    const std::vector<std::vector<double>>& starts,const PeriodicFilterOptions& options,
    const numerics::NonlinearLeastSquaresOptions& optimization){
    const auto nm=options.measurementChannels.size(),nc=atlas.knots()[0].channels.size();
    if(times.empty()||nm==0||measurements.size()!=times.size()*nm||options.measurementVariances.size()!=nm||starts.empty())throw std::invalid_argument("Invalid periodic identification arrays");
    for(std::size_t i=0;i<times.size();++i)if(!std::isfinite(times[i])||(i&&times[i]<=times[i-1]))throw std::invalid_argument("Identification times must increase");
    for(auto c:options.measurementChannels)if(c>=nc)throw std::invalid_argument("Identification channel out of range");
    for(double v:options.measurementVariances)if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Positive finite measurement variances required");
    for(double v:measurements)if(!std::isfinite(v))throw std::invalid_argument("Finite measurements required");
    if(!std::isfinite(options.frozenFrequency)||options.frozenFrequency<0)throw std::invalid_argument("Nonnegative frozen frequency required");
    const math::PeriodicAtlasTable table(atlas,options.tableSamples);
    const auto objective=[&](const std::vector<double>& p){
        const double omega=options.frozenFrequency>0?options.frozenFrequency:table.frequency(p[0]);
        std::vector<double> residual;residual.reserve(measurements.size());
        for(std::size_t i=0;i<times.size();++i){const auto y=table.evaluate(p[0],p[1]+omega*(times[i]-times.front()));
            for(std::size_t j=0;j<nm;++j){const auto value=y[options.measurementChannels[j]];
                if(std::abs(value.imag())>1e-9*(1+std::abs(value.real())))throw std::invalid_argument("Periodic identification requires real readouts");
                residual.push_back((value.real()-measurements[i*nm+j])/std::sqrt(options.measurementVariances[j]));}
        }return residual;
    };
    numerics::NonlinearLeastSquaresResult best;best.squaredResidual=std::numeric_limits<double>::infinity();
    const double pi=std::acos(-1.);
    for(const auto& initial:starts){if(initial.size()!=2)throw std::invalid_argument("Periodic start needs [parameter,phase]");
        auto result=numerics::boundedNonlinearLeastSquares(objective,initial,
            {atlas.knots().front().parameter,initial[1]-pi},{atlas.knots().back().parameter,initial[1]+pi},optimization);
        if(result.squaredResidual<best.squaredResidual)best=std::move(result);
    }
    best.parameters[1]=std::remainder(best.parameters[1],2*pi);return best;
}
}

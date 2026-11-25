#include "MeasuredJetControl.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace dynsyslab::workflow::detail {
namespace c = dynsyslab::core;

c::ScalarControlLaw measuredJetControlLaw(const c::KoopmanLQRController& prototype,
    c::Observable measurement,c::DerivativeSpec derivative,double sampleTime,
    std::size_t count,std::size_t warmupSamples,double warmupInput) {
    if(!count||!(sampleTime>0)||!std::isfinite(warmupInput)||measurement.outputDimension()!=1)
        throw std::invalid_argument("online_measurement needs a scalar signal, a sample time and a finite warmup input");
    if(!c::isCausalDerivativeMethod(derivative.method))
        throw std::invalid_argument("online_measurement needs a causal derivative method");
    struct Stream { c::KoopmanLQRController controller; c::CausalDifferentiator differentiator; std::size_t samples{}; };
    auto streams=std::make_shared<std::vector<Stream>>();
    streams->reserve(count);
    for(std::size_t j=0;j<count;++j)streams->push_back({prototype,c::CausalDifferentiator(derivative,sampleTime,1.0)});
    return [streams,measurement=std::move(measurement),warmupSamples,warmupInput](
        std::size_t stream,double time,const double* state) {
        auto& s=streams->at(stream);
        const double y=measurement.evaluateReal(state).front();  // the only use of the state
        const double rate=s.differentiator.step(y);
        if(++s.samples<=warmupSamples)return warmupInput;
        const std::array<double,2> coordinates{y,rate};
        return s.controller.valueFromCoordinates(time,coordinates.data()).front();
    };
}
}

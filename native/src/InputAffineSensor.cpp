#include <dynsyslab/core/InputAffineSensor.hpp>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core {
InputAffineSensor::InputAffineSensor(std::vector<Observable> base,std::vector<double> gain,
    std::optional<Observable> scale,std::size_t component)
    :base_(std::move(base)),gain_(std::move(gain)),scale_(std::move(scale)),component_(component){
    if(base_.empty())throw std::invalid_argument("Sensor needs base observables");
    for(const auto& b:base_){if(b.stateDimension()!=base_[0].stateDimension())throw std::invalid_argument("Sensor state dimension mismatch");width_+=b.outputDimension();}
    if(!gain_.empty()&&gain_.size()!=width_)throw std::invalid_argument("Sensor feedthrough needs one coefficient per raw output");
    for(double d:gain_)if(!std::isfinite(d))throw std::invalid_argument("Finite sensor feedthrough required");
    if(scale_&&(gain_.empty()||scale_->stateDimension()!=base_[0].stateDimension()||component_>=scale_->outputDimension()))
        throw std::invalid_argument("Invalid sensor feedthrough scale");
    if(!scale_&&component_)throw std::invalid_argument("Sensor scale component requires a scale observable");
}
std::vector<double> InputAffineSensor::evaluate(const double* state,double input) const {
    if(!state||!std::isfinite(input))throw std::invalid_argument("Finite input and state pointer required");
    std::vector<double> y;y.reserve(width_);
    for(const auto& b:base_){const auto v=b.evaluateReal(state);y.insert(y.end(),v.begin(),v.end());}
    if(!gain_.empty()){
        const double scale=scale_?scale_->evaluateReal(state).at(component_):1.;
        for(std::size_t k=0;k<width_;++k)y[k]+=gain_[k]*input*scale;
    }
    for(double v:y)if(!std::isfinite(v))throw std::domain_error("Nonfinite simulated sensor output");
    return y;
}
}

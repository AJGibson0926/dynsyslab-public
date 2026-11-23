#pragma once
#include <dynsyslab/core/InputAffineSensor.hpp>
#include <dynsyslab/workflow/Json.hpp>
#include <cmath>
#include <memory>
#include <stdexcept>
namespace dynsyslab::workflow {
inline std::shared_ptr<const core::InputAffineSensor> bindInputSensor(
    const json::Value& spec,const std::vector<core::Observable>& sensors,const core::ObservableCatalog& catalog){
    const auto& v=spec.at("sensor_feedthrough");std::vector<double> gain;std::optional<core::Observable> scale;std::size_t component=0;
    if(!v.isNull()){
        if(!std::holds_alternative<json::Value::Array>(v.at("coefficients").storage()))
            throw std::invalid_argument("sensor_feedthrough.coefficients array required");
        for(const auto& d:v.at("coefficients").array()){
            if(!std::holds_alternative<double>(d.storage())||!std::isfinite(d.number()))throw std::invalid_argument("Finite sensor coefficients required");gain.push_back(d.number());}
        if(gain.empty())throw std::invalid_argument("Explicit sensor feedthrough cannot be empty");
        if(!v.at("scale_observable").isNull())scale=catalog.at(v.at("scale_observable").string());
        const double c=v.at("scale_component").isNull()?0:v.at("scale_component").number(-1);
        if(!std::isfinite(c)||c<0||c!=std::floor(c)||c>1000000)throw std::invalid_argument("Invalid sensor scale component");component=std::size_t(c);
    }
    return std::make_shared<const core::InputAffineSensor>(sensors,std::move(gain),std::move(scale),component);
}
}

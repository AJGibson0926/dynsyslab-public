#include <dynsyslab/workflow/ObservableQuery.hpp>
#include <dynsyslab/core/System.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace dynsyslab::workflow {
json::Value evaluateObservableQuery(const json::Value& request) {
    if(request.at("schema").string()!="dynsyslab.observable-query/v1")
        throw std::invalid_argument("unsupported observable query schema");
    const auto& spec=request.at("system");
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier=spec.at("name").string(spec.at("identifier").string());
    configuration.preset=spec.at("preset").string();
    for(const auto& item:spec.at("parameters").object())configuration.parameters[item.first]=item.second.number();
    auto system=dynsyslab::core::createSystem(configuration);
    const auto catalog=system->observableCatalog();
    if (request.at("operation").string()=="list") {
        json::Value::Array entries;
        for (const auto& name:catalog.names()) {
            const auto& obs=catalog.at(name);
            json::Value::Object entry{{"name",name},{"units",obs.descriptor().units},
                {"coordinates",obs.descriptor().stateCoordinateSystem},{"outputs",double(obs.outputDimension())},
                {"real",!obs.descriptor().complexValued}};
            if (const auto* link=catalog.timeDerivative(name))
                entry["time_derivative"]=json::Value::Object{{"observable",link->target},
                    {"field",link->zeroInputOnly?"zero-input":"all-inputs"}};
            entries.emplace_back(std::move(entry));
        }
        return json::Value::Object{{"schema","dynsyslab.observable-catalog/v1"},{"observables",entries}};
    }
    if(request.at("operation").string()=="evaluate_batch") {
        if(request.at("flow_order").number()!=0||!request.at("partial_axes").array().empty())
            throw std::invalid_argument("Batch readout evaluates values only; use a derived observable for derivatives");
        if(!std::holds_alternative<json::Value::Array>(request.at("states").storage())||request.at("states").array().size()>10000000)
            throw std::invalid_argument("Bounded flat states array required");
        std::vector<double> states;
        for(const auto& v:request.at("states").array()){
            if(!std::holds_alternative<double>(v.storage())||!std::isfinite(v.number()))throw std::invalid_argument("Finite batch states required");states.push_back(v.number());}
        const auto dimension=system->descriptor().stateDimension;
        if(states.empty()||states.size()%dimension)throw std::invalid_argument("Batch state dimensions differ");
        std::vector<core::Observable> observables;
        for(const auto& n:request.at("observables").array())observables.push_back(catalog.at(n.string()));
        if(observables.empty())throw std::invalid_argument("Batch readout requires named observables");
        const auto combined=core::Observable::stack("batch","Batch readout",observables);
        const auto samples=states.size()/dimension;
        if(combined.outputDimension()>10000000/samples)throw std::invalid_argument("Batch readout output resource limit");
        json::Value::Array values;for(double v:combined.evaluateRealBatch(states,samples)){
            if(!std::isfinite(v))throw std::domain_error("Nonfinite batch readout");values.emplace_back(v);}
        return json::Value::Object{{"schema","dynsyslab.observable-result/v1"},{"values",values},
            {"samples",double(samples)},{"output_dimension",double(combined.outputDimension())},
            {"state_dimension",double(dimension)},{"value_type","real"},{"layout","sample-major"}};
    }
    std::vector<double> x;
    for(const auto& item:request.at("state").array())x.push_back(item.number());
    if(x.size()!=system->descriptor().stateDimension)
        throw std::invalid_argument("observable query state dimension mismatch");
    for(const auto value:x)if(!std::isfinite(value))throw std::invalid_argument("query state must be finite");
    const auto name=request.at("observable").string();
    auto observable=request.at("expression").isNull()
        ?(catalog.coordinatePresets().count(name)?catalog.coordinatePreset(name):catalog.at(name))
        :dynsyslab::core::Observable::fromExpressionDefinition(request.at("expression"));
    if(observable.stateDimension()!=x.size())throw std::invalid_argument("expression and query dimensions differ");
    for(const auto& axis:request.at("partial_axes").array()) {
        const double value=axis.number(-1);
        if(value<0||value>=x.size()||value!=std::floor(value))throw std::invalid_argument("invalid partial derivative axis");
        observable=observable.partialDerivative(std::size_t(value));
    }
    const double order=request.at("flow_order").number();
    if(order<0||order>8||order!=std::floor(order))throw std::invalid_argument("flow_order must be an integer from 0 to 8");
    const bool includeInputs=request.at("include_input_terms").boolean(false);
    bool varyingInput=false;
    if(order) {
        const auto field=system->driftExpressions();
        if(field.empty())throw std::invalid_argument("system has no analytic zero-input drift expressions");
        if(includeInputs) {
            const auto parseInput=[](const json::Value& value) {
                std::vector<double> result;
                if(std::holds_alternative<double>(value.storage()))result.push_back(value.number());
                else for(const auto& item:value.array())result.push_back(item.number());
                return result;
            };
            const auto inputs=parseInput(request.at("input"));
            const auto rates=parseInput(request.at("input_time_derivative"));
            varyingInput=std::any_of(rates.begin(),rates.end(),[](double v){return v!=0;});
            observable=observable.inputFlowDerivative(field,system->inputExpressions(),inputs,rates,std::size_t(order));
        } else observable=observable.flowDerivative(field,std::size_t(order));
    }
    json::Value::Array values,jacobian;
    const bool complex=observable.descriptor().complexValued;
    auto append=[&](json::Value::Array& output,std::complex<double> value) {
        if(complex)output.emplace_back(json::Value::Array{value.real(),value.imag()});
        else output.emplace_back(value.real());
    };
    if(complex)for(const auto value:observable.evaluate(x.data()))append(values,value);
    else for(const auto value:observable.evaluateReal(x.data()))values.emplace_back(value);
    for(const auto value:observable.jacobian(x.data()))append(jacobian,value);
    json::Value::Object result{{"schema","dynsyslab.observable-result/v1"},
        {"observable",observable.descriptor().name},{"values",std::move(values)},
        {"jacobian",std::move(jacobian)},{"state_dimension",double(x.size())},
        {"output_dimension",double(observable.outputDimension())},
        {"value_type",complex?"complex":"real"},
        {"gradient_backend",observable.hasAnalyticJacobian()?"analytic":"finite-difference"},
        {"flow_field",includeInputs&&order?"supplied-input":"autonomous-zero-input"},
        {"input_terms_included",includeInputs&&order!=0},
        {"input_time_assumption",includeInputs&&order?(varyingInput?"supplied-first-derivative":"held-constant"):"ignored"},
        {"flow_order",order}};
    if(!observable.expressions().empty())result["expression"]=observable.expressionDefinition();
    return result;
}
}

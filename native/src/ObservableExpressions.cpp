#include <dynsyslab/core/Observable.hpp>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace dynsyslab::core {
Observable Observable::fromExpressions(ObservableDescriptor descriptor,
    std::vector<AnalyticExpression> expressions) {
    if(expressions.empty()||descriptor.complexValued)
        throw std::invalid_argument("analytic observables require nonempty real expressions");
    const auto dimension=expressions.front().stateDimension();
    std::vector<AnalyticExpression> partials;
    for(const auto& expression:expressions) {
        if(expression.stateDimension()!=dimension)throw std::invalid_argument("observable expression dimensions differ");
        for(std::size_t axis=0;axis<dimension;++axis)partials.push_back(expression.derivative(axis));
    }
    Observable result(std::move(descriptor),dimension,expressions.size(),
        [expressions](const double* x) {
            Values values; values.reserve(expressions.size());
            for(const auto& expression:expressions)values.emplace_back(expression.evaluate(x));
            return values;
        },[partials](const double* x) {
            Values values; values.reserve(partials.size());
            for(const auto& expression:partials)values.emplace_back(expression.evaluate(x));
            return values;
        });
    result.expressions_=std::move(expressions);
    return result;
}
Observable Observable::partialDerivative(std::size_t axis) const {
    if(expressions_.empty())throw std::logic_error("observable has no differentiable expression definition");
    std::vector<AnalyticExpression> partials;
    for(const auto& item:expressions_)partials.push_back(item.derivative(axis));
    auto descriptor=descriptor_; descriptor.name+="_dx"+std::to_string(axis);
    descriptor.label="Partial derivative of "+descriptor_.label;
    descriptor.outputNames.clear(); descriptor.units.clear();
    return fromExpressions(std::move(descriptor),std::move(partials));
}
Observable Observable::flowDerivative(const std::vector<AnalyticExpression>& field,std::size_t order) const {
    if(expressions_.empty())throw std::logic_error("observable has no differentiable expression definition");
    std::vector<AnalyticExpression> derivatives;
    for(const auto& item:expressions_)derivatives.push_back(item.flowDerivative(field,order));
    auto descriptor=descriptor_; descriptor.name+="_dt"+std::to_string(order);
    descriptor.label="Flow derivative of "+descriptor_.label;
    descriptor.outputNames.clear(); descriptor.units.clear();
    return fromExpressions(std::move(descriptor),std::move(derivatives));
}
Observable Observable::inputFlowDerivative(
    const std::vector<AnalyticExpression>& drift,
    const std::vector<AnalyticExpression>& inputFields,
    const std::vector<double>& inputs,const std::vector<double>& inputRates,std::size_t order) const {
    if(inputs.empty()||drift.size()!=stateDimension_||inputFields.size()!=stateDimension_*inputs.size()
        ||(!inputRates.empty()&&inputRates.size()!=inputs.size()))
        throw std::invalid_argument("analytic input field and input dimensions differ");
    for(const auto value:inputs)if(!std::isfinite(value))throw std::invalid_argument("input must be finite");
    for(const auto value:inputRates)if(!std::isfinite(value))throw std::invalid_argument("input rate must be finite");
    if(order>2&&std::any_of(inputRates.begin(),inputRates.end(),[](double v){return v!=0;}))
        throw std::invalid_argument("nonconstant input derivatives above order two require higher input derivatives");
    auto field=drift;
    for(std::size_t axis=0;axis<stateDimension_;++axis)
        for(std::size_t input=0;input<inputs.size();++input)
            field[axis]=field[axis]+inputFields[axis*inputs.size()+input]*inputs[input];
    auto result=flowDerivative(field,order);
    if(order==2&&!inputRates.empty()) {
        auto terms=result.expressions();
        for(std::size_t output=0;output<terms.size();++output)
            for(std::size_t axis=0;axis<stateDimension_;++axis)
                for(std::size_t input=0;input<inputs.size();++input)
                    terms[output]=terms[output]+expressions_[output].derivative(axis)
                        *inputFields[axis*inputs.size()+input]*inputRates[input];
        result=fromExpressions(result.descriptor(),std::move(terms));
    }
    return result;
}

json::Value Observable::expressionDefinition() const {
    if(expressions_.empty())throw std::logic_error("observable has no serializable expression definition");
    json::Value::Array expressions,names;
    for(const auto& item:expressions_)expressions.push_back(item.toJson());
    for(const auto& item:descriptor_.outputNames)names.emplace_back(item);
    return json::Value::Object{{"schema","dynsyslab.observable-expression/v1"},
        {"name",descriptor_.name},{"label",descriptor_.label},{"units",descriptor_.units},
        {"coordinates",descriptor_.stateCoordinateSystem},{"output_names",std::move(names)},
        {"expressions",std::move(expressions)}};
}
Observable Observable::fromExpressionDefinition(const json::Value& value) {
    if(value.at("schema").string()!="dynsyslab.observable-expression/v1")
        throw std::invalid_argument("unsupported observable expression encoding");
    ObservableDescriptor descriptor{value.at("name").string(),value.at("label").string(),
        value.at("units").string(),value.at("coordinates").string(),{},false};
    for(const auto& name:value.at("output_names").array())descriptor.outputNames.push_back(name.string());
    std::vector<AnalyticExpression> expressions;
    for(const auto& item:value.at("expressions").array())expressions.push_back(AnalyticExpression::fromJson(item));
    return fromExpressions(std::move(descriptor),std::move(expressions));
}
}  // namespace dynsyslab::core

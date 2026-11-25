#include <dynsyslab/core/FunctionExpansion.hpp>
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
    // The formulas as a FunctionExpansion: one expression term per output, identity coefficients.
    std::vector<DictionaryTerm> terms;
    for(auto& expression:expressions) {
        DictionaryTerm term; term.kind="expression";
        term.expression=std::make_shared<const AnalyticExpression>(std::move(expression));
        terms.push_back(std::move(term));
    }
    const auto outputs=terms.size();
    std::vector<std::complex<double>> identity(outputs*outputs);
    for(std::size_t r=0;r<outputs;++r)identity[r*outputs+r]=1.0;
    result.representation_=std::make_shared<const FunctionExpansion>(
        std::make_shared<const ExpansionDictionary>(dimension,std::move(terms)),std::move(identity),
        FunctionDomain::physicalState(dimension),
        result.descriptor_.outputNames.size()==outputs?result.descriptor_.outputNames:std::vector<std::string>{});
    return result;
}
Observable Observable::partialDerivative(std::size_t axis) const {
    const auto own=expressions();
    if(own.empty())throw std::logic_error("observable has no differentiable expression definition");
    std::vector<AnalyticExpression> partials;
    for(const auto& item:own)partials.push_back(item.derivative(axis));
    auto descriptor=descriptor_; descriptor.name+="_dx"+std::to_string(axis);
    descriptor.label="Partial derivative of "+descriptor_.label;
    descriptor.outputNames.clear(); descriptor.units.clear(); descriptor.decompositions.clear();
    return fromExpressions(std::move(descriptor),std::move(partials));
}
Observable Observable::flowDerivative(const std::vector<AnalyticExpression>& field,std::size_t order) const {
    const auto own=expressions();
    if(own.empty())throw std::logic_error("observable has no differentiable expression definition");
    std::vector<AnalyticExpression> derivatives;
    for(const auto& item:own)derivatives.push_back(item.flowDerivative(field,order));
    auto descriptor=descriptor_; descriptor.name+="_dt"+std::to_string(order);
    descriptor.label="Flow derivative of "+descriptor_.label;
    descriptor.outputNames.clear(); descriptor.units.clear(); descriptor.decompositions.clear();
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
        const auto own=expressions();
        for(std::size_t output=0;output<terms.size();++output)
            for(std::size_t axis=0;axis<stateDimension_;++axis)
                for(std::size_t input=0;input<inputs.size();++input)
                    terms[output]=terms[output]+own[output].derivative(axis)
                        *inputFields[axis*inputs.size()+input]*inputRates[input];
        result=fromExpressions(result.descriptor(),std::move(terms));
    }
    return result;
}

json::Value Observable::expressionDefinition() const {
    const auto own=this->expressions();
    if(own.empty())throw std::logic_error("observable has no serializable expression definition");
    json::Value::Array expressions,names,references;
    for(const auto& ref:descriptor_.decompositions)references.emplace_back(json::Value::Object{
        {"artifact",ref.artifact},{"model_identity",ref.modelIdentity},{"component",ref.component},{"kind",ref.kind},{"index",double(ref.index)}});
    for(const auto& item:own)expressions.push_back(item.toJson());
    for(const auto& item:descriptor_.outputNames)names.emplace_back(item);
    json::Value::Object definition{{"schema","dynsyslab.observable-expression/v1"},
        {"name",descriptor_.name},{"label",descriptor_.label},{"units",descriptor_.units},
        {"coordinates",descriptor_.stateCoordinateSystem},{"output_names",std::move(names)},
        {"expressions",std::move(expressions)}};
    if(!references.empty())definition["decompositions"]=std::move(references);
    return definition;
}
Observable Observable::fromExpressionDefinition(const json::Value& value) {
    if(value.at("schema").string()!="dynsyslab.observable-expression/v1")
        throw std::invalid_argument("unsupported observable expression encoding");
    ObservableDescriptor descriptor{value.at("name").string(),value.at("label").string(),
        value.at("units").string(),value.at("coordinates").string(),{},false};
    for(const auto& name:value.at("output_names").array())descriptor.outputNames.push_back(name.string());
    for(const auto& item:value.at("decompositions").array()) {
        const double index=item.at("index").number(-1);
        if(!std::isfinite(index)||index<0||index!=std::floor(index)||index>9007199254740991.)throw std::invalid_argument("Invalid observable decomposition index");
        descriptor.decompositions.push_back({item.at("artifact").string(),item.at("model_identity").string(),
            item.at("component").string(),item.at("kind").string("reconstruction"),std::size_t(index)});
    }
    std::vector<AnalyticExpression> expressions;
    for(const auto& item:value.at("expressions").array())expressions.push_back(AnalyticExpression::fromJson(item));
    return fromExpressions(std::move(descriptor),std::move(expressions));
}
}  // namespace dynsyslab::core

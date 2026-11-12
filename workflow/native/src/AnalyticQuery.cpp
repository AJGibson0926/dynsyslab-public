#include <cmath>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>
#include <stdexcept>
namespace dynsyslab::workflow {
namespace {
json::Value vectorJson(const std::vector<double> &values) {
    json::Value::Array result;
    for (double x : values) {
        if (!std::isfinite(x))
            result.emplace_back();
        else
            result.emplace_back(x);
    }
    return result;
}
double argument(const json::Value &r, const char *key) {
    const auto &x = r.at(key);
    if (!std::holds_alternative<double>(x.storage()) || !std::isfinite(x.number()))
        throw std::invalid_argument(std::string(key) + " must be a finite number");
    return x.number();
}
std::vector<double> finiteVector(const json::Value &value, const char *name) {
    if (!std::holds_alternative<json::Value::Array>(value.storage()))
        throw std::invalid_argument(std::string(name) + " must be an array of finite numbers");
    std::vector<double> result;
    for (const auto &element : value.array()) {
        if (!std::holds_alternative<double>(element.storage()) || !std::isfinite(element.number()))
            throw std::invalid_argument(std::string(name) + " must contain finite numbers");
        result.push_back(element.number());
    }
    return result;
}
} // namespace
json::Value analyticDescriptorJson(const dynsyslab::core::AnalyticSolutionDescriptor &d) {
    json::Value::Array references;
    for (const auto &url : d.references)
        references.emplace_back(url);
    return json::Value::Object{{"id", d.identifier},         {"label", d.label},
                               {"formula_latex", d.formula}, {"conditions", d.conditions},
                               {"evaluation", d.evaluation}, {"version", d.version},
                               {"references", references}};
}
json::Value evaluateAnalyticQuery(const json::Value &request) {
    if (request.at("schema").string() != "dynsyslab.analytic-query/v1")
        throw std::invalid_argument("Unsupported analytic query schema");
    const auto &spec = request.at("system");
    dynsyslab::core::SystemConfiguration config;
    config.identifier = spec.at("name").string(spec.at("identifier").string());
    config.preset = spec.at("preset").string();
    for (const auto &p : spec.at("parameters").object())
        config.parameters[p.first] = argument(spec.at("parameters"), p.first.c_str());
    const auto &driver = spec.at("driver");
    config.driver = {driver.at("amplitude").number(), driver.at("frequency_hz").number(),
                     driver.at("phase_radians").number()};
    auto definition = resolveSystemDefinition(config);
    std::vector<dynsyslab::core::AnalyticSolutionFamily> families;
    auto knowledge = definition.descriptor.research;
    if (definition.descriptor.family == "ode") {
        auto system = dynsyslab::core::createSystem(config);
        families = system->analyticSolutions();
        knowledge = system->dynamicalKnowledge();
    } else if (definition.descriptor.analyticSolutions) {
        families = definition.descriptor.analyticSolutions(definition.configuration);
    }
    json::Value::Array catalog;
    for (const auto &family : families)
        catalog.push_back(analyticDescriptorJson(family.descriptor));
    json::Value::Object result{
        {"schema", "dynsyslab.analytic-result/v1"},
        {"system", definition.configuration.identifier},
        {"system_family", definition.descriptor.family},
        {"families", catalog},
        {"coordinate_system", definition.descriptor.preferredCoordinates},
        {"equations_fingerprint", dynsyslab::core::equationsFingerprint(definition.configuration)}};
    result["dynamical_knowledge"] = dynamicalKnowledgeJson(knowledge);
    if (request.at("initial_state").isNull())
        return result;
    auto initial = finiteVector(request.at("initial_state"), "initial_state");
    auto time = finiteVector(request.at("time"), "time");
    if (initial.size() != definition.descriptor.stateDimension)
        throw std::invalid_argument("Analytic IC dimension mismatch");
    if (time.empty() || time.size() > 1000000)
        throw std::invalid_argument("Analytic query requires 1..1000000 times");
    auto comparison = dynsyslab::core::compareAnalyticSolution(
        families, initial,
        request.at("initial_time").isNull() ? 0 : argument(request, "initial_time"), time,
        std::vector<double>(time.size() * initial.size()), request.at("family").string("auto"));
    result["status"] = comparison.status;
    result["reason"] = comparison.reason;
    result["solution"] = analyticDescriptorJson(comparison.descriptor);
    result["time"] = vectorJson(time);
    result["values"] = vectorJson(comparison.values);
    result["valid"] = vectorJson(comparison.valid);
    result["state_dimension"] = double(initial.size());
    return result;
}
json::Value evaluateSpecialFunctionQuery(const json::Value &request) {
    if (request.at("schema").string() != "dynsyslab.special-function-query/v1")
        throw std::invalid_argument("Unsupported special-function query schema");
    namespace s = dynsyslab::core::special;
    json::Value::Object result{{"schema", "dynsyslab.special-function-result/v1"}};
    const auto name = request.at("function").string();
    result["function"] = name;
    if (name == "jacobi") {
        auto j = s::jacobi(argument(request, "u"), argument(request, "m"));
        result["sn"] = j.sn;
        result["cn"] = j.cn;
        result["dn"] = j.dn;
    } else if (name == "elliptic_k")
        result["value"] = s::ellipticK(argument(request, "m"));
    else if (name == "elliptic_f")
        result["value"] = s::ellipticF(argument(request, "phi"), argument(request, "m"));
    else if (name == "jacobi_amplitude")
        result["value"] = s::jacobiAmplitude(argument(request, "u"), argument(request, "m"));
    else if (name == "theta") {
        double index = argument(request, "index");
        if (index < 1 || index > 4 || std::floor(index) != index)
            throw std::invalid_argument("Theta index must be 1..4");
        result["value"] =
            s::theta(unsigned(index), argument(request, "z"), argument(request, "nome"));
    } else if (name == "weierstrass") {
        auto w = s::weierstrass(argument(request, "z"), argument(request, "g2"),
                                argument(request, "g3"));
        result["value"] = w.value;
        result["derivative"] = w.derivative;
        result["real_period"] =
            std::isfinite(w.realPeriod) ? json::Value(w.realPeriod) : json::Value();
    } else if (name == "lauricella_fd") {
        result["value"] = s::lauricellaFD(argument(request, "a"),
            finiteVector(request.at("b"), "b"), argument(request, "c"),
            finiteVector(request.at("x"), "x"));
    } else if (name == "linear_delay") {
        result["value"] = dynsyslab::core::linearDelayConstantHistory(
            argument(request, "beta"), argument(request, "gamma"), argument(request, "delay"),
            argument(request, "history"), argument(request, "elapsed"));
    } else
        throw std::invalid_argument("Unknown special function");
    return result;
}
} // namespace dynsyslab::workflow

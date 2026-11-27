#include <cmath>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <dynsyslab/core/PeriodicAnalysis.hpp>
#include <dynsyslab/core/CoefficientShift.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/workflow/AnalyticQuery.hpp>
#include <dynsyslab/workflow/SystemKnowledge.hpp>
#include <dynsyslab/workflow/SystemResolution.hpp>
#include <stdexcept>
#include "TorusQuery.hpp"
namespace dynsyslab::workflow {
json::Value evaluateAnalyticHistoryQuery(const json::Value&,const dynsyslab::core::ODESystem&);
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
    if (request.at("operation").string() == "phase_cover") return evaluatePhaseCoverQuery(request);  // PhaseCoverQuery.cpp
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
    json::Value::Object parameters,groups;
    for(const auto& [name,value]:definition.configuration.parameters)parameters[name]=value;
    for(const auto& [name,value]:definition.configuration.dimensionlessGroups)groups[name]=value;
    result["resolved_parameters"]=parameters;result["dimensionless_groups"]=groups;
    if (request.at("include_eigenfunctions").boolean()) {
        if (definition.descriptor.family != "ode")
            throw std::invalid_argument("Known polynomial eigenfunctions require an ODE system");
        const auto system = dynsyslab::core::createSystem(config);
        result["known_eigenfunctions"] = dynsyslab::core::knownEigenfunctionsJson(
            dynsyslab::core::knownPolynomialEigenfunctions(*system));
    }
    if(request.at("operation").string()=="histories") {
        if(definition.descriptor.family!="ode")throw std::invalid_argument("Analytic histories currently require an ODE system");
        const auto system=dynsyslab::core::createSystem(config);
        const auto histories=evaluateAnalyticHistoryQuery(request,*system);
        for(const auto&[key,value]:histories.object())result[key]=value;
        return result;
    }
    if (request.at("initial_state").isNull())
        return result;
    auto initial = finiteVector(request.at("initial_state"), "initial_state");
    if(request.at("operation").string()=="torus") {
        if(definition.descriptor.family!="ode")throw std::invalid_argument("This torus adapter requires an ODE system");
        const auto system=dynsyslab::core::createSystem(config);
        const auto analysis=evaluateTorusAnalysis(request,*system);
        for(const auto& [key,value]:analysis.object())result[key]=value;
        return result;
    }
    const bool shiftCheck=request.at("operation").string()=="coefficient_shift";
    if(request.at("operation").string()=="periodic"||shiftCheck) {
        if(definition.descriptor.family!="ode")throw std::invalid_argument("Periodic observable sampling currently requires an ODE");
        const auto system=dynsyslab::core::createSystem(config);
        if(initial.size()!=definition.descriptor.stateDimension)throw std::invalid_argument("Periodic IC dimension mismatch");
        const double t0=request.at("initial_time").isNull()?0:argument(request,"initial_time");
        dynsyslab::core::AnalyticTrajectory trajectory;
        const auto selection=request.at("family").string("auto");
        std::string reason;
        dynsyslab::core::AnalyticSolutionFamily selectedFamily;
        for(const auto& family:families){
            if(selection!="auto"&&family.descriptor.identifier!=selection)continue;
            trajectory=family.instantiate(initial,t0);
            if(trajectory){result["solution"]=analyticDescriptorJson(family.descriptor);selectedFamily=family;break;}
            reason+=trajectory.reason+"; ";
        }
        if(!trajectory)throw std::invalid_argument("No applicable periodic analytic evaluator: "+reason);
        const double count=argument(request,"samples");
        if(count<3||count>1000000||count!=std::floor(count))throw std::invalid_argument("Periodic samples must be an integer from 3 to 1000000");
        std::vector<dynsyslab::core::Observable> observables;
        const auto catalog=system->observableCatalog();json::Value::Array names;
        for(const auto& name:request.at("observables").array()){
            const auto& o=catalog.at(name.string());
            if(o.outputDimension()!=1)throw std::invalid_argument("Periodic query requires scalar named channels");
            observables.push_back(o);names.push_back(o.descriptor().name);
        }
        if(shiftCheck) {
            const double maximum=argument(request,"maximum_harmonic");
            if(maximum<0||2*maximum>=count||maximum!=std::floor(maximum))
                throw std::invalid_argument("Need integer maximum_harmonic below Nyquist");
            const auto check=dynsyslab::core::checkPeriodicCoefficientShift(selectedFamily,initial,observables,
                std::size_t(count),int(maximum),finiteVector(request.at("shifts"),"shifts"),t0,
                request.at("period").isNull()?0:argument(request,"period"));
            json::Value::Array rows;
            for(const auto& row:check.shifts)rows.emplace_back(json::Value::Object{
                {"shift",row.shift},{"relative_frequency_error",row.relativeFrequencyError},
                {"initial_state_error",row.initialStateError},{"closure_error",row.closureError},
                {"maximum_coefficient_error",vectorJson(row.maximumCoefficientError)},
                {"coefficient_l2_error",vectorJson(row.coefficientL2Error)},
                {"reference_l2_norm",vectorJson(row.referenceL2Norm)}});
            result["status"]="complete";result["period"]=check.period;
            result["angular_frequency"]=check.angularFrequency;result["closure_error"]=check.closureError;
            result["observables"]=names;result["checks"]=rows;result["certified"]=false;
            result["scope"]="Numerical finite-Fourier time-shift check; requires autonomous periodic flow; no tail or global eigenfunction proof";
            return result;
        }
        const auto sample=dynsyslab::core::samplePeriodicTrajectory(trajectory,initial.size(),observables,std::size_t(count),t0,
            request.at("period").isNull()?0:argument(request,"period"));
        result["period"]=sample.period;result["angular_frequency"]=sample.angularFrequency;
        result["closure_error"]=sample.closureError;result["time"]=vectorJson(sample.times);
        result["values"]=vectorJson(sample.states);result["state_dimension"]=double(sample.stateDimension);
        result["observables"]=names;json::Value::Array measurements,series;
        for(auto c:sample.observables)measurements.emplace_back(json::Value::Array{c.real(),c.imag()});
        result["measurements"]=measurements;
        if(!request.at("maximum_harmonic").isNull()){
            const double maximum=argument(request,"maximum_harmonic");
            if(maximum<0||2*maximum>=count||maximum!=std::floor(maximum))throw std::invalid_argument("Need integer 0 <= maximum_harmonic < samples/2");
            for(std::size_t channel=0;channel<sample.channels;++channel){
                std::vector<std::complex<double>> data;for(std::size_t i=0;i<count;++i)data.push_back(sample.observables[i*sample.channels+channel]);
                const auto s=dynsyslab::math::periodicFourierSeries(data,sample.angularFrequency,int(maximum));
                json::Value::Array coefficients;for(const auto& [n,c]:s.coefficients)coefficients.emplace_back(json::Value::Array{double(n),c.real(),c.imag()});
                series.emplace_back(json::Value::Object{{"angular_frequency",s.angularFrequency},{"coefficients",coefficients}});
            }
            result["fourier_series"]=series;
        }
        result["status"]="complete";result["sample_convention"]="uniform-model-time-no-endpoint";
        return result;
    }
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
    result["period"] = comparison.period>0?json::Value(comparison.period):json::Value();
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
    } else if (name == "quartic_orbit") {
        const auto c=finiteVector(request.at("coefficients"),"coefficients");
        if(c.size()!=5)throw std::invalid_argument("Quartic orbit needs five ascending coefficients");
        const s::QuarticOrbit orbit({c[0],c[1],c[2],c[3],c[4]},argument(request,"initial"),argument(request,"velocity"));
        const auto time=finiteVector(request.at("elapsed"),"elapsed");
        if(time.size()>1000000)throw std::invalid_argument("Quartic sample count exceeds limit");
        json::Value::Array values;for(double t:time){const auto x=orbit.evaluate(t);values.emplace_back(json::Value::Array{x[0],x[1]});}
        result["values"]=values;result["parameter_period"]=std::isfinite(orbit.realPeriod())?json::Value(orbit.realPeriod()):json::Value();
        result["clock"]="elapsed is the quartic equation parameter, not necessarily the physical system time";
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

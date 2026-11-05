#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
json::Value numbers(const std::vector<double>& v) {
    json::Value::Array a; for (double x:v) a.emplace_back(x); return a;
}
json::Value complexNumbers(const std::vector<std::complex<double>>& v) {
    json::Value::Array a;
    for (auto x:v) a.emplace_back(json::Value::Array{x.real(),x.imag()});
    return a;
}
void different(std::vector<std::string>& out, bool differs, const std::string& field) {
    if(differs)out.push_back(field);
}
}
std::vector<std::string> systemCompatibilityDifferences(
    const SystemConfiguration& expected, const SystemConfiguration& actual) {
    const auto e=createSystem(expected)->configuration();
    const auto a=createSystem(actual)->configuration();
    std::vector<std::string> out;
    different(out,e.identifier!=a.identifier,"system.identifier");
    std::set<std::string> keys;
    for(const auto& p:e.parameters)keys.insert(p.first);
    for(const auto& p:a.parameters)keys.insert(p.first);
    for(const auto& k:keys)
        different(out,!e.parameters.count(k)||!a.parameters.count(k)
            || e.parameters.at(k)!=a.parameters.at(k),"parameters."+k);
    // Retain the archived equation assertion; resolving current code cannot erase it.
    const auto ev=expected.equationsVersion.empty()?e.equationsVersion:expected.equationsVersion;
    const auto av=actual.equationsVersion.empty()?a.equationsVersion:actual.equationsVersion;
    different(out,ev!=av,"equations_version");
    different(out,e.coordinateScaling.secondsPerModelTime!=a.coordinateScaling.secondsPerModelTime,
        "seconds_per_model_time");
    different(out,e.coordinateScaling.physicalInputPerModelInput!=a.coordinateScaling.physicalInputPerModelInput,
        "physical_input_per_model_input");
    // Explicit archived groups also detect a changed nondimensionalization.
    for(const auto& g:expected.dimensionlessGroups)
        different(out,!a.dimensionlessGroups.count(g.first)||a.dimensionlessGroups.at(g.first)!=g.second,
            "dimensionless_groups."+g.first);
    return out;
}
void requireCompatibleSystems(const SystemConfiguration& e,const SystemConfiguration& a) {
    const auto differences=systemCompatibilityDifferences(e,a);
    if(differences.empty())return;
    std::string message="Scientific system mismatch:";
    for(const auto& d:differences)message+=" "+d;
    throw std::invalid_argument(message);
}
std::string koopmanModelIdentity(const KoopmanModel& model) {
    auto coordinates=eigenfunctionCoordinatesJson(model.coordinates());
    // The frozen provider's content identity, rather than its filesystem location,
    // belongs in the model identity. Keep the path only for unversioned providers.
    if(!model.coordinates().providerIdentity.empty())coordinates.object().erase("provider_artifact");
    const auto& s=model.spectrum();
    json::Value::Object m{{"schema","dynsyslab.model-identity/v1"},
        {"coordinates",coordinates},{"inputs",double(model.inputDimension())},
        {"domain",s.domain()==EigenvalueDomain::discreteTime?"discrete":"continuous"},
        {"dt",s.sampleInterval()},{"spectrum",complexNumbers(s.values())}};
    json::Value::Array branches;for(auto b:s.logarithmBranches())branches.emplace_back(double(b));
    m["branches"]=branches;
    if(s.hasContinuousRates())m["rates"]=complexNumbers(s.continuousRates());
    if(model.hasPhysicalSystem()) {
        const auto& system=model.systemConfiguration();
        json::Value::Object parameters;for(const auto& p:system.parameters)parameters[p.first]=p.second;
        m["system"]=system.identifier;m["parameters"]=parameters;m["equations"]=system.equationsVersion;
        m["time_scale"]=system.coordinateScaling.secondsPerModelTime;
        m["input_scale"]=system.coordinateScaling.physicalInputPerModelInput;
    }
    if(const auto* scattered=model.scatteredEigenfunctions()) {
        m["states"]=numbers(scattered->trainingStates);
        m["values"]=complexNumbers(scattered->trainingValues);
        m["neighbors"]=double(scattered->neighbors);
        m["gradient_step"]=scattered->finiteDifferenceRelativeStep;
    } else {
        const auto& p=model.polynomialEigenfunctions();const auto& d=p.dictionary();
        json::Value::Array exponents;for(auto n:d.exponents)exponents.emplace_back(double(n));
        m["exponents"]=exponents;m["center"]=numbers(d.center);m["scale"]=numbers(d.scale);
        m["coefficients"]=complexNumbers(p.coefficients());
    }
    return "sha256:"+sha256Text(json::dump(m));
}
void requireContinuousControlModel(const KoopmanModel& model,const std::string& representation) {
    if(!model.spectrum().hasContinuousRates())
        throw std::invalid_argument("Continuous controller requires finite continuous rates");
    if(!model.hasControlMap(representation)) {
        if(model.discreteInputMaps().count(representation))
            throw std::invalid_argument("Discrete input map cannot be used by a continuous controller; no implicit conversion is permitted");
        throw std::invalid_argument("Continuous control map is unavailable");
    }
}
} // namespace dynsyslab::core

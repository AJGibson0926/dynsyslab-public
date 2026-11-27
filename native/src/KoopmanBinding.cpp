#include "KoopmanBinding.hpp"
#include "KoopmanArchiveMetadata.hpp"
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
SystemConfiguration readSystem(const json::Value& metadata);
bool coordinateOnly(const NpzRecord& a) {
    // A K-ROM or eigenfunctions object of data with no producing system is a function of its
    // coordinates only.
    return a.metadataJson().empty() || readSystem(json::parse(a.metadataJson())).identifier.empty();
}
SystemConfiguration readSystem(const json::Value& metadata) {
    const auto& encoded=metadata.at("system_configuration");
    SystemConfiguration s;
    s.identifier=encoded.at("identifier").string();
    const bool legacy=s.identifier.empty();
    const auto& source=legacy?metadata:encoded;
    if(legacy)s.identifier=metadata.at("system").string();
    s.preset=legacy?source.at("system_preset").string(source.at("preset").string()):source.at("preset").string();
    s.equationsVersion=source.at("equations_version").string();
    const auto& parameters=legacy
        ? (source.at("resolved_system_parameters").object().empty()
            ?source.at("system_parameters"):source.at("resolved_system_parameters"))
        :source.at("parameters");
    for(const auto& p:parameters.object())s.parameters[p.first]=p.second.number();
    for(const auto& p:source.at("parameter_units").object())s.parameterUnits[p.first]=p.second.string();
    for(const auto& p:source.at("dimensionless_groups").object())s.dimensionlessGroups[p.first]=p.second.number();
    for(const auto& p:source.at("dimensionless_group_units").object())s.dimensionlessGroupUnits[p.first]=p.second.string();
    for(const auto& p:source.at("initial_state").array())s.initialState.push_back(p.number());
    if(s.identifier=="rayleigh-plesset")s.identifier="rpe";
    return s;
}
bool meaningful(const std::string& s) { return s.find_first_not_of(" \t\r\n")!=std::string::npos; }
}
SystemConfiguration systemFromMetadata(const NpzRecord& archive) {
    if(coordinateOnly(archive))return {};
    auto s=readSystem(archive.metadataJson().empty()?json::Value{}:json::parse(archive.metadataJson()));
    if(s.identifier.empty())throw std::invalid_argument(
        "KROM is not self-contained; use the explicit-system loader for a legacy archive");
    return s;
}
ArchiveBinding resolveArchiveBinding(const NpzRecord& a,
    const SystemConfiguration& requested,const ModelBindingOptions& options) {
    const auto metadata=a.metadataJson().empty()?json::Value{}:json::parse(a.metadataJson());
    auto archived=readSystem(metadata);
    json::Value::Object evidence{{"schema","dynsyslab.model-binding/v1"}};
    if(!metadata.at("binding_evidence").isNull())evidence["prior_binding"]=metadata.at("binding_evidence");
    if(coordinateOnly(a)) {
        // Do not invent a physical model for measured coordinates at deployment.
        evidence["status"]="coordinate-only";
        return {{},json::dump(evidence)};
    }
    const bool known=!archived.identifier.empty()&&!archived.parameters.empty();
    if(!known) {
        // An archive that names no physics (e.g. eigenfunctions computed from measured channels of a
        // solve) is bound to the system the caller names, and the binding says so.
        if(requested.identifier.empty())throw std::invalid_argument("this archive records no physical system; name one explicitly");
        // Even partial metadata must not be contradicted without an explicit transfer.
        if(!archived.identifier.empty()) {
            auto partial=requested;partial.identifier=archived.identifier;
            partial.equationsVersion=archived.equationsVersion;
            for(const auto& p:archived.parameters)partial.parameters[p.first]=p.second;
            const auto differences=systemCompatibilityDifferences(partial,requested);
            if(!differences.empty()&&!meaningful(options.systemTransferReason))
                requireCompatibleSystems(partial,requested);
        }
        evidence["status"]="requested-system";
        evidence["transfer_reason"]=options.systemTransferReason;
        evidence["requested_system"]=koopmanSystemJson(requested);
        return {requested,json::dump(evidence)};
    }
    const auto target=requested.identifier.empty()?archived:requested;
    const auto differences=systemCompatibilityDifferences(archived,target);
    if(!differences.empty()&&!meaningful(options.systemTransferReason))
        requireCompatibleSystems(archived,target);
    evidence["status"]=differences.empty()?"verified":"explicit-system-transfer";
    evidence["reason"]=differences.empty()?std::string{}:options.systemTransferReason;
    evidence["training_system"]=koopmanSystemJson(archived);
    evidence["requested_system"]=koopmanSystemJson(target);
    json::Value::Array fields;for(const auto& field:differences)fields.emplace_back(field);
    evidence["differing_fields"]=fields;
    // Transfer never rewrites the trained model's physics or its exact Jacobian.
    return {archived,json::dump(evidence)};
}
}

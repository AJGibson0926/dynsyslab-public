#include <dynsyslab/workflow/KromValidationPlan.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include "KromDatasetSupport.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace core=dynsyslab::core;
namespace json=core::json;
json::Value planAutonomousKoopmanValidation(const core::KoopmanModel& model,
    const std::filesystem::path& modelPath,const std::filesystem::path& root,std::size_t maximum) {
    if(!maximum)throw std::invalid_argument("Validation lookup requires a positive dataset limit");
    if(!model.hasPhysicalSystem())return json::Value::Object{
        {"status","coordinate-provider-required"},{"simulation_executed",false},
        {"message","Coordinate-only KROM: provide compatible unforced measurements and the frozen coordinate provider. No physical system is assumed."}};
    const auto system=core::createSystem(model.systemConfiguration());
    std::set<std::filesystem::path> paths;
    std::set<std::string> preferred,training;
    for(const auto& ref:model.datasetReferences()) {
        auto path=std::filesystem::u8path(ref.path);
        if(path.is_relative()) {
            if(std::filesystem::is_regular_file(root/path))path=root/path;
            else if(std::filesystem::is_regular_file(root.parent_path()/path))path=root.parent_path()/path;
            else path=modelPath.parent_path()/path;
        }
        if(std::filesystem::is_regular_file(path)) {
            path=std::filesystem::absolute(path).lexically_normal();
            if(ref.role=="training")training.insert(path.generic_u8string());
            else {paths.insert(path);preferred.insert(path.generic_u8string());}
        }
    }
    for(const auto& item:core::loadArchiveIndex(root))
        if(std::find(item.arrays.begin(),item.arrays.end(),"time")!=item.arrays.end()
            && (std::find(item.arrays.begin(),item.arrays.end(),"values")!=item.arrays.end()
            || std::find(item.arrays.begin(),item.arrays.end(),"states")!=item.arrays.end())
            && std::filesystem::is_regular_file(item.path))
            paths.insert(std::filesystem::absolute(item.path).lexically_normal());
    // The read-only index may be absent or stale after an external import.
    // Enumerate solution pointers, then inflate metadata only, not state arrays.
    if(std::filesystem::is_directory(root/"solutions"))
        for(const auto& item:std::filesystem::recursive_directory_iterator(root/"solutions"))
            if(item.is_regular_file()&&item.path().extension()==".npz")
                paths.insert(std::filesystem::absolute(item.path()).lexically_normal());
    struct Candidate {std::filesystem::path path;json::Value meta;bool preferred;};
    std::vector<Candidate> matching;json::Value::Array rejected;
    for(const auto& path:paths)try {
        if(training.count(path.generic_u8string()))continue;
        const auto index=core::readNpzIndex(path);
        if(!index.contains("time")||(!index.contains("values")&&!index.contains("states")))continue;
        const auto meta=json::parse(index.metadataJson);
        postprocess_detail::SourceData source;source.path=path;source.metadataJson=index.metadataJson;
        // Unknown forcing is not silently accepted, even if input arrays are absent.
        if(postprocess_detail::datasetDynamics(source)!="unforced")continue;
        const auto cfg=postprocess_detail::systemConfiguration(source);
        const auto resolved=core::createSystem(cfg)->configuration();
        if(resolved.identifier!=system->configuration().identifier
            ||resolved.parameters!=system->configuration().parameters
            ||(!cfg.equationsVersion.empty()&&cfg.equationsVersion!=system->configuration().equationsVersion))continue;
        const double step=meta.at("time_configuration").at("output_step").number();
        const double dt=model.spectrum().sampleInterval();
        const auto coordinates=meta.at("archive_coordinate_system").string(meta.at("integration_form").string());
        if(step>0 && dt>0 && coordinates==system->descriptor().preferredCoordinates) {
            const double ratio=dt/step;
            if(ratio<1-1e-8||std::abs(ratio-std::round(ratio))>1e-7*std::max(1.,ratio))continue;
        }
        matching.push_back({path,meta,preferred.count(path.generic_u8string())!=0});
    } catch(const std::exception& e) {
        rejected.emplace_back(json::Value::Object{{"path",path.generic_u8string()},{"reason",e.what()}});
    }
    std::sort(matching.begin(),matching.end(),[](const auto& a,const auto& b){
        if(a.preferred!=b.preferred)return a.preferred>b.preferred;
        const auto ad=a.meta.at("created_utc").string(),bd=b.meta.at("created_utc").string();
        return ad==bd?a.path.generic_u8string()<b.path.generic_u8string():ad>bd;});
    json::Value::Array candidates,sources;std::vector<core::KoopmanDatasetReference> refs;
    for(const auto& item:matching) {
        candidates.emplace_back(json::Value::Object{{"path",item.path.generic_u8string()},
            {"dynamics","unforced"},{"initial_state",item.meta.at("initial_state")},
            {"time_configuration",item.meta.at("time_configuration")}});
        if(sources.size()<maximum) {
            sources.emplace_back(item.path.generic_u8string());
            refs.push_back({item.path.generic_u8string(),"validation","unforced",{}});
        }
    }
    json::Value::Object result{{"schema","dynsyslab.krom-validation-plan/v1"},
        {"source_krom",modelPath.empty()?"":std::filesystem::absolute(modelPath).generic_u8string()},
        {"status",matching.empty()?"simulation-offer":"existing-data-available"},
        {"matched_count",double(matching.size())},{"inspected_count",double(paths.size())},
        {"candidates",std::move(candidates)},{"unreadable_or_unsupported",std::move(rejected)},
        {"dataset_references",core::koopmanDatasetReferencesJson(refs)},
        {"simulation_executed",false},{"requires_execution_verification",true}};
    if(!sources.empty()&&!modelPath.empty())result["validation_job"]=json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","Autonomous eigenfunction validation"},{"sources",std::move(sources)},
        {"methods",json::Value::Array{json::Value::Object{{"name","koopman_validation"},
            {"settings",json::Value::Object{{"krom",std::filesystem::absolute(modelPath).generic_u8string()},
                {"dataset_role","validation"}}}}}},{"archive",json::Value::Object{{"enabled",true}}}};
    else if(sources.empty()) {
        json::Value::Object parameters;
        for(const auto& [name,value]:system->configuration().parameters)parameters[name]=value;
        result["simulation_offer"]=json::Value::Object{
            {"message","No compatible unforced dataset was found. Offer to simulate from the equations of motion after choosing representative ICs and a validation duration; do not run automatically."},
            {"system",json::Value::Object{{"name",system->descriptor().identifier},
                {"preset",system->configuration().preset},{"parameters",std::move(parameters)}}},
            {"controller_kind","none"},{"output_step",model.spectrum().sampleInterval()},
            {"required_choices",json::Value::Array{"initial_conditions","duration","approval_to_simulate"}}};
    }
    if(model.requiresExternalCoordinates())result["coordinate_warning"]=
        "Reuse candidates are raw trajectories; apply the frozen history coordinate provider before eigenfunction validation.";
    return result;
}

json::Value planAutonomousKoopmanValidation(const std::filesystem::path& path,
    const std::filesystem::path& root,std::size_t maximum) {
    return planAutonomousKoopmanValidation(core::loadKoopmanModel(path),path,root,maximum);
}
}

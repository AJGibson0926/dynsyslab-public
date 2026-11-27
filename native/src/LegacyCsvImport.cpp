#include <dynsyslab/core/LegacyCsv.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Hierarchy.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::string timestamp() {
    const auto now=std::chrono::system_clock::now();const auto raw=std::chrono::system_clock::to_time_t(now);std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc,&raw);
#else
    gmtime_r(&raw,&utc);
#endif
    std::ostringstream out;out<<std::put_time(&utc,"%Y%m%dT%H%M%S")<<'Z';return out.str();
}

std::string safe(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char item){return static_cast<char>(std::tolower(item));});
    for(char& item:value)if(!std::isalnum(static_cast<unsigned char>(item))&&item!='-'&&item!='_')item='-';
    while(value.find("--")!=std::string::npos)value.replace(value.find("--"),2,"-");return value.empty()?"legacy-system":value;
}

json::Value numericObject(const std::map<std::string,double>& values) {
    json::Value::Object result;for(const auto& item:values)result[item.first]=item.second;return json::Value(std::move(result));
}

json::Value stringObject(const std::map<std::string,std::string>& values) {
    json::Value::Object result;for(const auto& item:values)result[item.first]=item.second;return json::Value(std::move(result));
}

json::Value stringArray(const std::vector<std::string>& values) {
    json::Value::Array result;for(const auto& value:values)result.emplace_back(value);return json::Value(std::move(result));
}

std::string groupKey(const LegacyCsvSystem& system,const LegacyCsvImportOptions& options) {
    std::ostringstream out;out<<(options.mappedSystem.identifier.empty()?system.dynamicsName:options.mappedSystem.identifier);
    for(const auto& field:system.trailingFields)out<<'|'<<field;return out.str();
}

struct Identity {
    std::string identifier,label,preset,coordinateSystem,equationsVersion,fingerprint;
    std::vector<std::string> stateNames;
    std::map<std::string,double> parameters;
    std::map<std::string,std::string> parameterUnits;
    std::map<std::string,double> groups;
    std::map<std::string,std::string> groupUnits;
};

Identity identity(const LegacyCsvSystem& legacy,const LegacyCsvImportOptions& options,
                  std::size_t dimension) {
    Identity result;result.coordinateSystem=options.coordinateSystem;
    if(options.mappedSystem.identifier.empty()){result.identifier="legacy-"+safe(legacy.dynamicsName);result.label=legacy.dynamicsName+" (legacy import)";result.preset="unresolved legacy parameters";result.equationsVersion="legacy-external-unresolved";result.fingerprint=stableIdentifier(groupKey(legacy,options));
        for(std::size_t axis=0;axis<dimension;++axis)result.stateNames.push_back("state_"+std::to_string(axis));return result;}
    auto system=createSystem(options.mappedSystem);if(system->descriptor().stateDimension!=dimension)throw std::invalid_argument("mapped native system dimension does not match legacy trajectory state dimension");
    result.identifier=system->descriptor().identifier;result.label=system->descriptor().label;result.preset=system->configuration().preset;result.coordinateSystem=system->descriptor().preferredCoordinates;
    if(options.coordinateSystem!="unknown"&&options.coordinateSystem!=result.coordinateSystem)throw std::invalid_argument("selected import coordinates differ from mapped native system coordinates");
    result.equationsVersion=system->configuration().equationsVersion;result.fingerprint=equationsFingerprint(system->configuration());result.stateNames=system->descriptor().stateNames;result.parameters=system->configuration().parameters;result.parameterUnits=system->configuration().parameterUnits;result.groups=system->configuration().dimensionlessGroups;result.groupUnits=system->configuration().dimensionlessGroupUnits;return result;
}

void atomicJson(const std::filesystem::path& path,const json::Value& value) {
    std::filesystem::create_directories(path.parent_path());const auto temporary=path.parent_path()/(path.filename().wstring()+L".partial");
    {std::ofstream stream(temporary,std::ios::binary|std::ios::trunc);if(!stream)throw std::runtime_error("unable to create legacy import manifest");stream<<json::dump(value)<<'\n';}
    std::error_code ignored;std::filesystem::rename(temporary,path,ignored);if(ignored){std::filesystem::remove(temporary);throw std::runtime_error("unable to publish legacy import manifest");}
}

std::vector<double> scaled(const std::vector<double>& values,double scale) {
    std::vector<double> result=values;for(double& value:result)value*=scale;return result;
}

json::Value importMetadata(const LegacyCsvBatch& batch,const LegacyCsvSystem& source,
                           const LegacyCsvTrajectory& trajectory,const Identity& model,
                           const LegacyCsvImportOptions& options,const std::string& label,
                           const std::string& batchId,const std::string& runId,
                           const std::string& clusterId,const std::string& trajectoryId) {
    json::Value::Object legacy;legacy["source_path"]=batch.source.generic_u8string();legacy["source_sha256"]=batch.sourceSha256;legacy["system_index"]=static_cast<double>(source.index);legacy["cluster_index"]=static_cast<double>(trajectory.clusterIndex);legacy["trajectory_index"]=static_cast<double>(trajectory.trajectoryIndex);legacy["dynamics_name"]=source.dynamicsName;legacy["declared_step"]=source.declaredStep;legacy["declared_transitions"]=static_cast<double>(source.declaredTransitions);legacy["trailing_fields"]=stringArray(source.trailingFields);legacy["solution_source_lines"]=json::parse("["+std::to_string(trajectory.firstSourceLine)+","+std::to_string(trajectory.lastSourceLine)+"]");
    json::Value::Object result;result["schema"]="dynsyslab.solution/v1";result["label"]=label;result["system"]=model.identifier;result["system_label"]=model.label;result["system_preset"]=model.preset;result["resolved_system_parameters"]=numericObject(model.parameters);result["parameter_units"]=stringObject(model.parameterUnits);result["dimensionless_groups"]=numericObject(model.groups);result["dimensionless_group_units"]=stringObject(model.groupUnits);result["equations_version"]=model.equationsVersion;result["equations_fingerprint"]=model.fingerprint;result["integration_form"]=model.coordinateSystem;result["forcing_type"]=source.controlValues.empty()?"free":"forced";result["data_role"]=source.controlValues.empty()?"imported_free":"control_training";result["legacy_import"]=json::Value(std::move(legacy));result["legacy_forcing"]=source.forcing;result["legacy_subforcing"]=source.subforcing;result["control_alignment"]=source.controlValues.empty()?"none":"stored samples 2:end are interval inputs";result["control_input_scale"]=options.inputScale;result["execution_backend"]=json::parse(R"({"name":"native-cpp-import","compiled":true,"python_required":false})");result["batch_id"]=batchId;result["system_run_id"]=runId;result["cluster_id"]=clusterId;result["trajectory_id"]=trajectoryId;return json::Value(std::move(result));
}

}  // namespace

json::Value LegacyCsvImportResult::summary() const {
    json::Value::Array paths;for(const auto& path:archives)paths.emplace_back(path.generic_u8string());json::Value::Object result;
    result["schema"]="dynsyslab.legacy-csv-import-result/v1";result["source"]=source.generic_u8string();result["source_sha256"]=sourceSha256;result["system_count"]=static_cast<double>(systemCount);result["trajectory_count"]=static_cast<double>(trajectoryCount);result["solution_archives"]=json::Value(std::move(paths));result["hierarchy_manifest"]=hierarchyManifest.generic_u8string();result["import_manifest"]=importManifest.generic_u8string();result["python_required"]=false;return json::Value(std::move(result));
}

LegacyCsvImportResult importLegacyDynBatchCsv(const std::filesystem::path& source,
    const std::filesystem::path& archiveRoot,const LegacyCsvImportOptions& options) {
    if(!std::isfinite(options.inputScale))throw std::invalid_argument("legacy import input scale must be finite");
    const auto batch=readLegacyDynBatchCsv(source,options.strict);const std::string stamp=timestamp();const std::string batchId=stableIdentifier(options.label+"|"+batch.sourceSha256+"|batch");
    BatchRecord hierarchy;hierarchy.id=batchId;hierarchy.label=options.label;hierarchy.createdUtc=stamp;hierarchy.dataRole="legacy_import";hierarchy.tags={"imported","legacy-csv"};hierarchy.properties={{"source_path",batch.source.generic_u8string()},{"source_sha256",batch.sourceSha256},{"python_required","false"},{"mapping",options.mappedSystem.identifier.empty()?"preserved legacy identity":"explicit native system"}};
    LegacyCsvImportResult result;result.source=batch.source;result.sourceSha256=batch.sourceSha256;result.trajectoryCount=batch.trajectoryCount();std::map<std::string,std::size_t> groups;std::map<std::pair<std::size_t,std::string>,std::size_t> clusters;
    for(const auto& legacy:batch.systems){if(legacy.trajectories.empty())continue;const std::size_t dimension=legacy.trajectories.front().stateDimension;const auto model=identity(legacy,options,dimension);const std::string key=groupKey(legacy,options);std::size_t runIndex{};
        const auto found=groups.find(key);if(found==groups.end()){runIndex=hierarchy.systems.size();groups[key]=runIndex;SystemRunRecord run;run.id=stableIdentifier(batchId+"|"+key);run.label=model.label;run.systemIdentifier=model.identifier;run.systemLabel=model.label;run.systemFamily="ode";run.preset=model.preset;run.presetLabel=model.preset;run.parameterSetKey=model.fingerprint;run.coordinateSystem=model.coordinateSystem;run.resolvedParameters=model.parameters;run.parameterUnits=model.parameterUnits;run.dimensionlessGroups=model.groups;run.dimensionlessGroupUnits=model.groupUnits;run.equationsVersion=model.equationsVersion;run.equationsFingerprint=model.fingerprint;run.stateNames=model.stateNames;hierarchy.systems.push_back(std::move(run));}else runIndex=found->second;
        auto& run=hierarchy.systems[runIndex];
        for(const auto& trajectory:legacy.trajectories){if(trajectory.stateDimension!=dimension)throw std::invalid_argument("legacy system contains inconsistent state dimensions");const std::string clusterKey=std::to_string(legacy.index)+"|"+std::to_string(trajectory.clusterIndex);const auto lookup=std::make_pair(runIndex,clusterKey);std::size_t clusterIndex;
            const auto clusterFound=clusters.find(lookup);if(clusterFound==clusters.end()){clusterIndex=run.clusters.size();clusters[lookup]=clusterIndex;ClusterRecord cluster;cluster.id=stableIdentifier(run.id+"|cluster|"+clusterKey);cluster.label="legacy-system-"+std::to_string(legacy.index)+"-cluster-"+std::to_string(trajectory.clusterIndex);run.clusters.push_back(std::move(cluster));}else clusterIndex=clusterFound->second;
            const std::string trajectoryId=stableIdentifier(run.clusters[clusterIndex].id+"|trajectory|"+std::to_string(trajectory.trajectoryIndex));const std::string itemLabel=options.label+"-s"+std::to_string(legacy.index)+"-c"+std::to_string(trajectory.clusterIndex)+"-t"+std::to_string(trajectory.trajectoryIndex);
            NpzRecord archive;archive.arrays["time"]=NpyArray::fromDoubles({trajectory.samples},trajectory.time);archive.arrays["values"]=NpyArray::fromDoubles({trajectory.samples,dimension},trajectory.states);
            if(!legacy.controlValues.empty()){if(legacy.controlValues.size()!=trajectory.samples)throw std::invalid_argument("legacy control and trajectory samples are not aligned");const auto controls=scaled(legacy.controlValues,options.inputScale);archive.arrays["control_values"]=NpyArray::fromDoubles({trajectory.samples,1},controls);archive.arrays["control_interval_values"]=NpyArray::fromDoubles({trajectory.samples-1,1},std::vector<double>(controls.begin()+1,controls.end()));}
            archive.arrays["metadata_json"]=NpyArray::fromString(json::dump(importMetadata(batch,legacy,trajectory,model,options,itemLabel,batchId,run.id,run.clusters[clusterIndex].id,trajectoryId)));
            const auto path=archiveRoot/L"solutions"/std::filesystem::u8path(model.identifier)/std::filesystem::u8path(safe(itemLabel)+"__"+stamp+".npz");writeNpzAtomic(path,archive);result.archives.push_back(path);
            TrajectoryRecord record;record.id=trajectoryId;record.label=itemLabel;record.initialState.assign(trajectory.states.begin(),trajectory.states.begin()+static_cast<std::ptrdiff_t>(dimension));record.solutionArchive=path;record.tags={"imported",legacy.controlValues.empty()?"free":"forced",legacy.controlValues.empty()?"imported_free":"control_training"};record.timeSeries.push_back({"state","time","values"});if(!legacy.controlValues.empty())record.timeSeries.push_back({"control","time","control_values"});run.clusters[clusterIndex].trajectories.push_back(std::move(record));}
    }
    result.systemCount=hierarchy.systems.size();result.hierarchyManifest=archiveRoot/L"hierarchies"/L"imports"/std::filesystem::u8path(safe(options.label)+"__"+stamp+".hierarchy.json");writeHierarchyAtomic(result.hierarchyManifest,hierarchy);
    result.importManifest=archiveRoot/L"imports"/L"legacy-csv"/std::filesystem::u8path(safe(options.label)+"__"+stamp+".import.json");json::Value::Object manifest;manifest["schema"]="dynsyslab.legacy-csv-import/v1";manifest["source_summary"]=batch.summary();manifest["result"]=result.summary();manifest["input_scale"]=options.inputScale;manifest["coordinate_system"]=options.coordinateSystem;manifest["mapped_system"]=options.mappedSystem.identifier;manifest["mapped_preset"]=options.mappedSystem.preset;manifest["mapped_parameters"]=numericObject(options.mappedSystem.parameters);atomicJson(result.importManifest,json::Value(std::move(manifest)));return result;
}

}  // namespace dynsyslab::core

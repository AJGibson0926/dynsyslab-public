#include "KromControlSearchSupport.hpp"
#include <dynsyslab/core/OutputPath.hpp>

#include <dynsyslab/core/Archive.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow::search_detail {

std::vector<double> numbers(const json::Value& value){std::vector<double> result;
    for(const auto& item:value.array())result.push_back(item.number());return result;}
json::Value values(const std::vector<double>& source){json::Value::Array result;
    for(const double value:source)result.emplace_back(value);return json::Value(std::move(result));}
json::Value indices(const std::vector<std::size_t>& source){json::Value::Array result;
    for(const auto value:source)result.emplace_back(static_cast<double>(value));
    return json::Value(std::move(result));}

json::Value candidateJson(const dynsyslab::core::KromCandidateScore& source){
    json::Value::Object result;result["mode_indices"]=indices(source.modeIndices);
    result["mode_count"]=static_cast<double>(source.modeIndices.size());
    result["selection_score"]=source.selectionScore;
    result["lifted_target_displacement"]=source.liftedTargetDisplacement;
    result["worst_normalized_controllability"]=source.worstNormalizedControllability;
    result["minimum_exact_input_authority"]=source.minimumExactInputAuthority;
    result["frequency_spread"]=source.frequencySpread;
    result["minimum_controllability_rank"]=static_cast<double>(source.minimumControllabilityRank);
    result["split_dimension"]=static_cast<double>(source.splitDimension);
    return json::Value(std::move(result));}

json::Value extractionJob(const json::Value& search,const std::filesystem::path& parent,
    const std::vector<std::size_t>& modes,std::size_t attempt){json::Value::Object settings;
    settings["selected_modes"]=indices(modes);json::Value::Object method;
    method["name"]="krom_extract";method["enabled"]=true;method["settings"]=json::Value(settings);
    json::Value::Object job;job["schema"]="dynsyslab.compute-job/v1";job["kind"]="postprocess";
    job["label"]=search.at("label").string("krom-control-search")+"-attempt-"+std::to_string(attempt);
    job["sources"]=json::Value::Array{json::Value(parent.generic_u8string())};
    job["methods"]=json::Value::Array{json::Value(std::move(method))};
    json::Value::Object archive;archive["enabled"]=true;job["archive"]=json::Value(std::move(archive));
    return json::Value(std::move(job));}

json::Value regressionJob(const json::Value& search,const std::filesystem::path& controlled,
    const std::filesystem::path& krom,std::size_t attempt){json::Value settings=search.at("regression");
    settings.object()["krom"]=krom.generic_u8string();json::Value::Array additional;
    const auto& sources=search.at("controlled_sources").array();
    for(std::size_t i=1;i<sources.size();++i)additional.emplace_back(sources[i].string());
    settings.object()["additional_sources"]=json::Value(std::move(additional));
    json::Value::Object method;method["name"]="control_regression";method["enabled"]=true;
    method["settings"]=std::move(settings);json::Value::Object job;
    job["schema"]="dynsyslab.compute-job/v1";job["kind"]="postprocess";
    job["label"]=search.at("label").string("krom-control-search")+"-regression-"+std::to_string(attempt);
    job["sources"]=json::Value::Array{json::Value(controlled.generic_u8string())};
    job["methods"]=json::Value::Array{json::Value(std::move(method))};
    json::Value::Object archive;archive["enabled"]=true;job["archive"]=json::Value(std::move(archive));
    return json::Value(std::move(job));}

json::Value solveJob(const json::Value& search,const dynsyslab::core::KoopmanModel& parent,
    const std::filesystem::path& krom,const std::string& representation,std::size_t attempt){
    const auto& configuration=parent.systemConfiguration();json::Value::Object parameters;
    for(const auto& item:configuration.parameters)parameters[item.first]=item.second;
    json::Value::Object system;system["name"]=configuration.identifier;
    system["preset"]=configuration.preset;system["parameters"]=json::Value(std::move(parameters));
    json::Value::Object initial;initial["mode"]="explicit";initial["count"]=1.0;
    initial["clusters"]=1.0;initial["values"]=search.at("initial_state");
    json::Value::Object target;target["kind"]="constant";target["state"]=search.at("target_state");
    json::Value settings=search.at("controller");settings.object()["krom"]=krom.generic_u8string();
    settings.object()["control_representation"]=representation;
    settings.object()["active_time_window"]=search.at("time").at("active_window");
    json::Value::Object controller;controller["kind"]="klqr";
    controller["settings"]=std::move(settings);controller["target"]=json::Value(std::move(target));
    json::Value::Object campaign;campaign["enabled"]=false;
    json::Value::Object archive;archive["enabled"]=true;archive["data_role"]="control_attempt";
    json::Value::Object preflight;preflight["enabled"]=false;
    json::Value::Object job;job["schema"]="dynsyslab.compute-job/v1";job["kind"]="solve";
    job["label"]=search.at("label").string("krom-control-search")+"-attempt-"
        +std::to_string(attempt)+"-"+representation;job["system"]=json::Value(std::move(system));
    job["initial_conditions"]=json::Value(std::move(initial));job["time"]=search.at("time");
    job["solver"]=search.at("solver");job["controller"]=json::Value(std::move(controller));
    job["campaign"]=json::Value(std::move(campaign));job["archive"]=json::Value(std::move(archive));
    job["preflight"]=json::Value(std::move(preflight));return json::Value(std::move(job));}

CaseMetrics metrics(const std::filesystem::path& path,double release,
    const std::vector<double>& target,const std::vector<double>& finalState,double settlingStart){
    const auto archive=dynsyslab::core::readNpz(path);const auto time=archive.at("time").doubles();
    const auto state=archive.at("values").doubles();const auto& shape=archive.at("values").shape;
    if(shape.size()!=2||shape[0]!=time.size()||shape[1]!=target.size()
        ||finalState.size()!=target.size())throw std::runtime_error("control result shape mismatch");
    const std::size_t dimension=shape[1];const auto found=std::lower_bound(time.begin(),time.end(),release);
    const std::size_t releaseIndex=std::min<std::size_t>(found-time.begin(),time.size()-1);
    auto distance=[&](std::size_t row,const std::vector<double>& reference){double sum=0.0;
        for(std::size_t column=0;column<dimension;++column){const double delta=
            state[row*dimension+column]-reference[column];sum+=delta*delta;}return std::sqrt(sum);};
    double squared=0.0;std::size_t count=0;
    for(std::size_t row=0;row<time.size();++row)if(time[row]>=settlingStart&&time[row]<=release){
        const double value=distance(row,target);squared+=value*value;++count;}
    double peak=0.0;if(archive.contains("control_values"))for(const double value:
        archive.at("control_values").doubles())peak=std::max(peak,std::abs(value));
    return {distance(releaseIndex,target),distance(time.size()-1,finalState),peak,
        count?std::sqrt(squared/static_cast<double>(count)):std::numeric_limits<double>::quiet_NaN()};}

json::Value metricsJson(const CaseMetrics& source){json::Value::Object result;
    result["distance_to_target_at_release"]=source.releaseError;
    result["final_distance_to_reference"]=source.finalError;
    result["control_peak_abs"]=source.controlPeak;
    result["controlled_state_rmse_after_settling"]=source.activeRmse;
    return json::Value(std::move(result));}

void writeJsonAtomic(const std::filesystem::path& path,const json::Value& value){
    (void)dynsyslab::core::prepareOutputFile(path, "campaign ledger");const auto temporary=path.wstring()+L".partial";
    {std::ofstream stream(temporary,std::ios::binary|std::ios::trunc);
        if(!stream)throw std::runtime_error("unable to create campaign ledger");
        const auto text=json::dump(value);stream.write(text.data(),static_cast<std::streamsize>(text.size()));}
#ifdef _WIN32
    if(!MoveFileExW(temporary.c_str(),path.c_str(),MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish campaign ledger");
#else
    std::filesystem::rename(temporary,path);
#endif
}

}  // namespace dynsyslab::workflow::search_detail

#include <dynsyslab/workflow/KromControlSearch.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KromScreening.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>

#include "KromControlSearchSupport.hpp"
#include "WorkflowSupport.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

std::size_t positiveInteger(const json::Value& value,std::size_t fallback,
    std::size_t maximum,const char* name){const double number=value.number(
        static_cast<double>(fallback));
    if(number<1.0||number>static_cast<double>(maximum)||std::floor(number)!=number)
        throw std::invalid_argument(std::string(name)+" must be an integer in [1, "
            +std::to_string(maximum)+"]");return static_cast<std::size_t>(number);}

json::Value::Array paths(const json::Value& source){json::Value::Array result;
    for(const auto& item:source.array()){const auto path=std::filesystem::absolute(
        std::filesystem::u8path(item.string()));if(!std::filesystem::is_regular_file(path))
        throw std::invalid_argument("campaign source does not exist: "+path.u8string());
        result.emplace_back(path.generic_u8string());}return result;}

bool succeeded(const json::Value& record,double releaseTolerance,double finalTolerance){
    return record.at("status").string()=="complete"
        &&record.at("metrics").at("distance_to_target_at_release").number(
            std::numeric_limits<double>::infinity())<=releaseTolerance
        &&record.at("metrics").at("final_distance_to_reference").number(
            std::numeric_limits<double>::infinity())<=finalTolerance;}

void report(const std::filesystem::path& path,const dynsyslab::core::ExecutionContext* execution,
    double fraction,std::size_t complete,std::size_t planned,const std::string& message){
    detail::atomicProgress(path,fraction,complete,planned,message);
    if(execution)execution->report({"krom-control-search",message,fraction,complete,planned});}

dynsyslab::core::ExecutionContext stageExecution(
    const dynsyslab::core::ExecutionContext* parent,const std::filesystem::path& progressPath,
    std::size_t complete,std::size_t planned,const std::string& label,
    std::chrono::steady_clock::time_point deadline=std::chrono::steady_clock::time_point::max()){
    dynsyslab::core::ExecutionContext child;if(parent)child.limits=parent->limits;
    child.cancellation=[parent,deadline]{return (parent&&parent->cancelled())
        ||std::chrono::steady_clock::now()>=deadline;};
    child.progress=[parent,progressPath,complete,planned,label](
        const dynsyslab::core::ExecutionProgress& event){
        const double local=std::clamp(event.fraction,0.0,1.0);
        const double global=(static_cast<double>(complete)+local)/planned;
        const auto completed=complete+static_cast<std::size_t>(local>=1.0);
        const std::string message=label+": "+event.message;
        detail::atomicProgress(progressPath,global,completed,planned,message);
        if(parent)parent->report({"krom-control-search/"+event.stage,message,
            global,completed,planned});};
    return child;}

json::Value runCase(const json::Value& job,const dynsyslab::core::KoopmanModel& parent,
    const std::filesystem::path& krom,const std::string& representation,std::size_t attempt,
    const std::filesystem::path& archiveRoot,const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution,std::size_t complete,std::size_t planned){
    json::Value::Object record;record["representation"]=representation;
    const double timeout=job.at("case_timeout_seconds").number(300.0);
    const auto deadline=std::chrono::steady_clock::now()+std::chrono::duration_cast<
        std::chrono::steady_clock::duration>(std::chrono::duration<double>(timeout));
    auto child=stageExecution(execution,progressPath,complete,planned,
        "Attempt "+std::to_string(attempt)+" "+representation+" KLQR",deadline);
    try{const auto run=executeJob(search_detail::solveJob(job,parent,krom,
            representation,attempt),archiveRoot,{},&child);
        if(run.archives.empty())throw std::runtime_error("control solve produced no archive");
        const auto target=search_detail::numbers(job.at("target_state"));
        const auto finalState=search_detail::numbers(job.at("final_reference_state"));
        const auto active=search_detail::numbers(job.at("time").at("active_window"));
        const auto metrics=search_detail::metrics(run.archives.front(),active[1],target,
            finalState,job.at("success").at("settling_start").number(active[0]));
        record["status"]="complete";record["archive"]=run.archives.front().generic_u8string();
        record["metrics"]=search_detail::metricsJson(metrics);
    }catch(const dynsyslab::core::ExecutionCancelled&){if(execution&&execution->cancelled())throw;
        record["status"]="timeout";record["message"]="case exceeded cooperative wall-time limit";
    }catch(const std::exception& error){record["status"]="failed";record["message"]=error.what();}
    return json::Value(std::move(record));}

json::Value ledger(const json::Value& job,const std::filesystem::path& parent,
    const json::Value::Array& controlled,const json::Value::Array& attempts,
    const std::string& status,const std::string& stopReason){json::Value::Object result;
    result["schema"]="dynsyslab.krom-control-search/v1";result["label"]=job.at("label").string();
    result["status"]=status;result["stop_reason"]=stopReason;
    result["parent_krom"]=parent.generic_u8string();result["controlled_sources"]=json::Value(controlled);
    result["policy"]=job.at("screening");result["success_policy"]=job.at("success");
    result["attempt_count"]=static_cast<double>(attempts.size());result["attempts"]=json::Value(attempts);
    return json::Value(std::move(result));}

std::filesystem::path archiveCampaign(const json::Value& job,
    const std::filesystem::path& destination,const std::filesystem::path& parent,
    const json::Value::Array& controlled,const json::Value::Array& attempts,
    const std::string& stopReason){const std::size_t count=attempts.size();
    const std::string cases[] {"exact","constant","bilinear"};
    std::vector<double> attemptNumber,modeCount,selection,modeOffsets{0.0},modeIndices;
    std::vector<double> release(count*3,std::numeric_limits<double>::quiet_NaN());
    std::vector<double> finalError=release,controlPeak=release,success(count*3,0.0);
    std::vector<double> constantMap(count,std::numeric_limits<double>::quiet_NaN());
    std::vector<double> bilinearMap=constantMap;
    const double releaseTolerance=job.at("success").at("release_error_tolerance").number();
    const double finalTolerance=job.at("success").at("final_error_tolerance").number();
    json::Value::Array resultArchives;
    for(std::size_t row=0;row<count;++row){const auto& attempt=attempts[row];
        attemptNumber.push_back(static_cast<double>(row+1));
        const auto indices=attempt.at("candidate").at("mode_indices").array();
        modeCount.push_back(static_cast<double>(indices.size()));for(const auto& value:indices)
            modeIndices.push_back(value.number());modeOffsets.push_back(static_cast<double>(modeIndices.size()));
        selection.push_back(attempt.at("candidate").at("selection_score").number());
        constantMap[row]=attempt.at("constant_error_mean").number(
            std::numeric_limits<double>::quiet_NaN());
        bilinearMap[row]=attempt.at("bilinear_error_mean").number(
            std::numeric_limits<double>::quiet_NaN());
        for(std::size_t column=0;column<3;++column){const auto& item=attempt.at("cases").at(cases[column]);
            if(item.at("status").string()=="complete"){const auto& metrics=item.at("metrics");
                release[row*3+column]=metrics.at("distance_to_target_at_release").number();
                finalError[row*3+column]=metrics.at("final_distance_to_reference").number();
                controlPeak[row*3+column]=metrics.at("control_peak_abs").number();
                success[row*3+column]=succeeded(item,releaseTolerance,finalTolerance)?1.0:0.0;
                resultArchives.emplace_back(item.at("archive").string());}}
        const auto regression=attempt.at("regressed_krom").string();
        if(!regression.empty())resultArchives.emplace_back(regression);}
    json::Value::Object metadata;metadata["format"]="dynsyslab.derived/v1";
    metadata["method"]="krom-control-search";metadata["label"]=job.at("label").string();
    metadata["system"]=dynsyslab::core::loadKoopmanModel(parent)
        .systemConfiguration().identifier;metadata["source_archive"]=parent.generic_u8string();
    metadata["source_krom"]=parent.generic_u8string();metadata["source_archives"]=json::Value(controlled);
    metadata["result_archives"]=json::Value(std::move(resultArchives));
    metadata["stop_reason"]=stopReason;metadata["attempt_count"]=static_cast<double>(count);
    metadata["control_cases"]=json::Value::Array{json::Value("exact"),json::Value("constant"),json::Value("bilinear")};
    metadata["initial_state"]=job.at("initial_state");metadata["target_state"]=job.at("target_state");
    metadata["final_reference_state"]=job.at("final_reference_state");
    metadata["success_policy"]=job.at("success");metadata["screening_policy"]=job.at("screening");
    metadata["attempts"]=json::Value(attempts);json::Value::Object backend;
    backend["name"]="native-cpp";backend["compiled"]=true;backend["python_required"]=false;
    metadata["execution_backend"]=json::Value(std::move(backend));
    dynsyslab::core::NpzRecord archive;archive.arrays["metadata_json"]=
        dynsyslab::core::NpyArray::fromString(json::dump(json::Value(std::move(metadata))));
    archive.arrays["attempt_number"]=dynsyslab::core::NpyArray::fromDoubles({count},attemptNumber);
    archive.arrays["candidate_mode_count"]=dynsyslab::core::NpyArray::fromDoubles({count},modeCount);
    archive.arrays["candidate_selection_score"]=dynsyslab::core::NpyArray::fromDoubles({count},selection);
    archive.arrays["candidate_mode_offsets"]=dynsyslab::core::NpyArray::fromDoubles({count+1},modeOffsets);
    archive.arrays["candidate_mode_indices"]=dynsyslab::core::NpyArray::fromDoubles({modeIndices.size()},modeIndices);
    archive.arrays["release_error"]=dynsyslab::core::NpyArray::fromDoubles({count,3},release);
    archive.arrays["final_error"]=dynsyslab::core::NpyArray::fromDoubles({count,3},finalError);
    archive.arrays["control_peak"]=dynsyslab::core::NpyArray::fromDoubles({count,3},controlPeak);
    archive.arrays["case_success"]=dynsyslab::core::NpyArray::fromDoubles({count,3},success);
    archive.arrays["constant_error_mean"]=dynsyslab::core::NpyArray::fromDoubles({count},constantMap);
    archive.arrays["bilinear_error_mean"]=dynsyslab::core::NpyArray::fromDoubles({count},bilinearMap);
    dynsyslab::core::writeNpzAtomic(destination,archive);return destination;}

}  // namespace

json::Value planKromControlSearch(const json::Value& job){const auto parent=
    std::filesystem::absolute(std::filesystem::u8path(job.at("parent_krom").string()));
    if(!std::filesystem::is_regular_file(parent))throw std::invalid_argument("parent KROM does not exist");
    const auto model=dynsyslab::core::loadKoopmanModel(parent);const auto initial=
        search_detail::numbers(job.at("initial_state"));const auto target=
        search_detail::numbers(job.at("target_state"));const auto finalState=
        search_detail::numbers(job.at("final_reference_state"));
    if(initial.size()!=model.stateDimension()||target.size()!=model.stateDimension()
        ||finalState.size()!=model.stateDimension())throw std::invalid_argument(
            "search initial, target, and final-reference states must match the KROM system");
    const auto controlled=paths(job.at("controlled_sources"));if(controlled.empty())
        throw std::invalid_argument("KROM control search requires controlled regression data");
    const auto active=search_detail::numbers(job.at("time").at("active_window"));
    const double start=job.at("time").at("start").number(),stop=job.at("time").at("stop").number();
    const double step=job.at("time").at("output_step").number();
    if(active.size()!=2||!(start<=active[0]&&active[0]<active[1]&&active[1]<stop)||!(step>0.0))
        throw std::invalid_argument("search time requires start <= active start < active stop < stop");
    const auto attempts=positiveInteger(job.at("screening").at("max_attempts"),25,25,"max_attempts");
    const auto modeCount=positiveInteger(job.at("screening").at("maximum_mode_count"),4,10,"maximum_mode_count");
    const double release=job.at("success").at("release_error_tolerance").number(0.02);
    const double finalTolerance=job.at("success").at("final_error_tolerance").number(0.1);
    if(!(release>0.0)||!(finalTolerance>0.0))throw std::invalid_argument("success tolerances must be positive");
    json::Value::Object result;result["schema"]="dynsyslab.workflow-plan/v1";
    result["status"]="planned";result["kind"]="krom_control_search";
    result["system"]=model.systemConfiguration().identifier;
    result["parent_mode_count"]=static_cast<double>(model.liftedDimension());
    result["maximum_mode_count"]=static_cast<double>(std::min(modeCount,model.liftedDimension()));
    result["maximum_attempts"]=static_cast<double>(attempts);
    result["control_cases_per_attempt"]=3.0;result["controlled_source_count"]=
        static_cast<double>(controlled.size());result["python_required"]=false;
    result["archive_policy"]="durable ledger plus immutable campaign/attempt archives";
    return json::Value(std::move(result));}

RunResult executeKromControlSearch(const json::Value& job,
    const std::filesystem::path& archiveRoot,const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution){(void)planKromControlSearch(job);
    const auto parentPath=std::filesystem::absolute(std::filesystem::u8path(job.at("parent_krom").string()));
    const auto parent=dynsyslab::core::loadKoopmanModel(parentPath);const auto initial=
        search_detail::numbers(job.at("initial_state"));const auto target=
        search_detail::numbers(job.at("target_state"));const auto controlled=paths(job.at("controlled_sources"));
    const std::size_t maxAttempts=positiveInteger(job.at("screening").at("max_attempts"),25,25,"max_attempts");
    dynsyslab::core::KromScreeningOptions options;
    options.maximumModeCount=positiveInteger(job.at("screening").at("maximum_mode_count"),4,10,"maximum_mode_count");
    options.maximumVirtualCandidates=positiveInteger(job.at("screening").at("maximum_virtual_candidates"),5000,100000,"maximum_virtual_candidates");
    options.retainedCandidates=std::max<std::size_t>(maxAttempts,positiveInteger(
        job.at("screening").at("retained_candidates"),100,10000,"retained_candidates"));
    options.seedModeCount=positiveInteger(job.at("screening").at("seed_mode_count"),12,1000,"seed_mode_count");
    const auto candidates=dynsyslab::core::rankKromControlCandidates(parent,initial,target,options);
    const auto directory=archiveRoot/L"derived"/L"krom-control-search"/
        dynsyslab::core::wide(detail::safeName(job.at("label").string(),"krom-control-search")
            +"__"+detail::timestamp());std::filesystem::create_directories(directory);
    const auto ledgerPath=directory/L"campaign-ledger.json";json::Value::Array attempts;
    const double releaseTolerance=job.at("success").at("release_error_tolerance").number(0.02);
    const double finalTolerance=job.at("success").at("final_error_tolerance").number(0.1);
    const bool stopOnSuccess=job.at("screening").at("stop_on_bilinear_success").boolean(true);
    const std::size_t planned=std::min(maxAttempts,candidates.size())*5,attemptLimit=std::min(maxAttempts,candidates.size());
    if(!planned)throw std::runtime_error("KROM screening produced no control candidates");
    std::size_t complete=0;std::string stopReason="hard-attempt-limit-reached";
    report(progressPath,execution,0.0,0,planned,"Screening KROM candidates");
    for(std::size_t attempt=0;attempt<attemptLimit;++attempt){const auto& candidate=candidates[attempt];
        json::Value::Object record;record["attempt_number"]=static_cast<double>(attempt+1);
        record["candidate"]=search_detail::candidateJson(candidate);json::Value::Object cases;
        std::filesystem::path extracted,regressed;std::string regressionFailure;
        try{auto child=stageExecution(execution,progressPath,complete,planned,
                "Attempt "+std::to_string(attempt+1)+" extraction");
            const auto result=executePostprocess(search_detail::extractionJob(job,parentPath,
                candidate.modeIndices,attempt+1),archiveRoot,{},&child);
            extracted=result.archives.at(0);}catch(const std::exception& error){regressionFailure=error.what();}
        report(progressPath,execution,static_cast<double>(++complete)/planned,complete,planned,"Extracted candidate KROM");
        if(!extracted.empty())try{auto child=stageExecution(execution,progressPath,complete,planned,
                "Attempt "+std::to_string(attempt+1)+" regression");
            const auto result=executePostprocess(search_detail::regressionJob(
                job,std::filesystem::u8path(controlled.front().string()),extracted,attempt+1),
                archiveRoot,{},&child);regressed=result.archives.at(0);
            const auto metadata=json::parse(dynsyslab::core::readNpz(regressed).metadataJson());
            record["constant_error_mean"]=metadata.at("constant_error_mean");
            record["bilinear_error_mean"]=metadata.at("bilinear_error_mean");
        }catch(const std::exception& error){regressionFailure=error.what();}
        record["extracted_krom"]=extracted.generic_u8string();record["regressed_krom"]=regressed.generic_u8string();
        if(!regressionFailure.empty())record["regression_failure"]=regressionFailure;
        report(progressPath,execution,static_cast<double>(++complete)/planned,complete,planned,"Regressed candidate control maps");
        for(const std::string representation:{"exact","constant","bilinear"}){
            if(extracted.empty()||(representation!="exact"&&regressed.empty())){json::Value::Object failed;
                failed["status"]="failed";failed["representation"]=representation;
                failed["message"]=regressionFailure.empty()?"candidate extraction failed":regressionFailure;
                cases[representation]=json::Value(std::move(failed));}
            else cases[representation]=runCase(job,parent,representation=="exact"?extracted:regressed,
                representation,attempt+1,archiveRoot,progressPath,execution,complete,planned);
            report(progressPath,execution,static_cast<double>(++complete)/planned,complete,planned,
                "Evaluated "+representation+" KLQR control");}
        record["cases"]=json::Value(cases);attempts.emplace_back(json::Value(std::move(record)));
        search_detail::writeJsonAtomic(ledgerPath,ledger(job,parentPath,controlled,attempts,"running","in-progress"));
        if(stopOnSuccess&&succeeded(cases.at("bilinear"),releaseTolerance,finalTolerance)){
            stopReason="bilinear-success";break;}
        if(attempt+1==candidates.size())stopReason="candidate-list-exhausted";}
    const auto campaignPath=archiveCampaign(job,directory/L"campaign.npz",parentPath,
        controlled,attempts,stopReason);search_detail::writeJsonAtomic(ledgerPath,
            ledger(job,parentPath,controlled,attempts,"complete",stopReason));
    report(progressPath,execution,1.0,planned,planned,"KROM control search archived");
    std::wstring summary=L"Native KROM control search completed\r\n\r\nAttempts: "+
        std::to_wstring(attempts.size())+L"\r\nStop reason: "+dynsyslab::core::wide(stopReason)
        +L"\r\nCampaign: "+campaignPath.wstring();
    return {"krom_control_search",std::move(summary),{campaignPath},ledgerPath};}

}  // namespace dynsyslab::workflow

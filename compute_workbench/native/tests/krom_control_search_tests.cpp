#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
std::string escape(const std::string& value){std::string result;for(char item:value){
    if(item=='\\'||item=='"')result.push_back('\\');result.push_back(item);}return result;}

std::filesystem::path makeKrom(const std::filesystem::path& root){
    dynsyslab::core::KoopmanSINDyResult fit;fit.continuousEigenvalues={{-1.0,0.0}};
    fit.coefficients={{1.0,0.0}};fit.exponents={0,1};fit.center={0,0};fit.scale={1,1};
    fit.stateDimension=2;fit.featureCount=1;dynsyslab::core::SystemConfiguration system;
    system.identifier="duffing";system.preset="double_well_reference";
    dynsyslab::core::KoopmanModel model(std::move(fit),system);
    model.setReconstruction("state",{{{0,0},{1,0}},2,1,0,0,0});const auto path=root/"velocity-krom.npz";
    dynsyslab::core::savePolynomialKoopmanModel(path,model);return path;}

dynsyslab::core::json::Value searchJob(const std::filesystem::path& krom,
    const std::filesystem::path& controlled){return dynsyslab::core::json::parse(
    R"({"schema":"dynsyslab.compute-job/v1","kind":"krom_control_search","label":"native-search-test","parent_krom":")"
    +escape(krom.u8string())+R"(","controlled_sources":[")"+escape(controlled.u8string())
    +R"("],"initial_state":[0,0],"target_state":[0,0],"final_reference_state":[0,0],"time":{"start":0,"active_window":[0,0.01],"stop":0.02,"output_step":0.01,"coordinate_system":"automatic"},"solver":{"method":"RK45","rtol":1e-8,"atol":1e-10,"max_step":0.005,"control_evaluation":"continuous_rhs"},"controller":{"Q":1,"R":1},"regression":{"constant_method":"multistep-output","bilinear_method":"controlled-residual","reconstruction":"state","ridge":1e-10,"conjugate_complete":false,"map_sample_stride":1,"grid_columns":5,"grid_rows":5},"screening":{"max_attempts":1,"maximum_mode_count":1,"maximum_virtual_candidates":10,"retained_candidates":1,"seed_mode_count":1,"stop_on_bilinear_success":false},"success":{"release_error_tolerance":0.01,"final_error_tolerance":0.01,"settling_start":0},"case_timeout_seconds":10,"archive":{"enabled":true}})");}

void run(){const auto root=std::filesystem::temp_directory_path()/"dynsyslab-krom-search-test";
    std::error_code ignored;std::filesystem::remove_all(root,ignored);std::filesystem::create_directories(root);
    const auto controlledJob=dynsyslab::core::json::parse(R"({
      "schema":"dynsyslab.compute-job/v1","kind":"solve","label":"search-regression-source",
      "system":{"name":"duffing","preset":"double_well_reference","parameters":{}},
      "time":{"start":0,"stop":0.04,"output_step":0.01,"coordinate_system":"automatic"},
      "initial_conditions":{"mode":"explicit","count":1,"clusters":1,"values":[[0,0]],"ranges":[]},
      "controller":{"kind":"waveform","settings":{"amplitude":0.2,"frequency":1},"target":{"kind":"constant","state":[0,0]}},
      "solver":{"method":"RK4","step_size":0.001},"archive":{"enabled":true}})");
    const auto controlled=dynsyslab::workflow::executeJob(controlledJob,root).archives.at(0);
    auto job=searchJob(makeKrom(root),controlled);const auto plan=dynsyslab::workflow::planJob(job);
    require(plan.at("maximum_attempts").number()==1,"bounded campaign plan");
    auto invalid=job;invalid.at("screening").at("max_attempts")=dynsyslab::core::json::Value(26.0);
    bool rejected=false;try{(void)dynsyslab::workflow::planJob(invalid);}
    catch(const std::invalid_argument&){rejected=true;}require(rejected,"hard 25-attempt cap");
    std::vector<double> progress;dynsyslab::core::ExecutionContext execution;
    execution.progress=[&](const dynsyslab::core::ExecutionProgress& event){progress.push_back(event.fraction);};
    const auto result=dynsyslab::workflow::executeJob(job,root,root/"progress.json",&execution);
    require(result.kind=="krom_control_search"&&result.archives.size()==1,"campaign result");
    const auto archive=dynsyslab::core::readNpz(result.archives.front());
    require(archive.at("release_error").shape==std::vector<std::size_t>({1,3})
        &&archive.contains("candidate_mode_indices"),"campaign archive contract");
    require(std::filesystem::is_regular_file(result.manifest),"durable campaign ledger");
    require(!progress.empty()&&std::is_sorted(progress.begin(),progress.end())
        &&progress.back()==1.0,"campaign progress is monotonic and complete");
    std::filesystem::remove_all(root,ignored);}

}  // namespace

int main(){try{run();std::cout<<"dynsyslab_krom_control_search_tests: PASS\n";return 0;}
    catch(const std::exception& error){std::cerr<<"dynsyslab_krom_control_search_tests: FAIL: "
        <<error.what()<<'\n';return 1;}}

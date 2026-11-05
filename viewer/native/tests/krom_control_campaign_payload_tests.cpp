#include "dynsyslab_viewer/NativePayload.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(){const auto path=std::filesystem::temp_directory_path()/
    "dynsyslab-krom-control-campaign-payload.npz";try{dynsyslab::core::NpzRecord archive;
    archive.arrays["attempt_number"]=dynsyslab::core::NpyArray::fromDoubles({2},{1,2});
    archive.arrays["candidate_mode_count"]=dynsyslab::core::NpyArray::fromDoubles({2},{1,2});
    archive.arrays["candidate_selection_score"]=dynsyslab::core::NpyArray::fromDoubles({2},{3,2});
    archive.arrays["release_error"]=dynsyslab::core::NpyArray::fromDoubles({2,3},{.1,.2,.3,.01,.02,.03});
    archive.arrays["final_error"]=dynsyslab::core::NpyArray::fromDoubles({2,3},{.2,.3,.4,.02,.03,.04});
    archive.arrays["control_peak"]=dynsyslab::core::NpyArray::fromDoubles({2,3},{1,2,3,2,3,4});
    archive.arrays["case_success"]=dynsyslab::core::NpyArray::fromDoubles({2,3},{0,0,0,1,1,1});
    archive.arrays["constant_error_mean"]=dynsyslab::core::NpyArray::fromDoubles({2},{.4,.2});
    archive.arrays["bilinear_error_mean"]=dynsyslab::core::NpyArray::fromDoubles({2},{.3,.1});
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"Native candidate campaign","method":"krom-control-search","success_policy":{"release_error_tolerance":0.02,"final_error_tolerance":0.05}})");
    dynsyslab::core::writeNpzAtomic(path,archive);const auto payload=dsl::buildNativePayload(
        path,"krom_control_campaign");if(payload.figures.size()!=6||payload.figures[0].traces.size()!=4
        ||!payload.warnings.empty())throw std::runtime_error("campaign dashboard contract");
    std::filesystem::remove(path);std::cout<<"dynsyslab_krom_control_campaign_payload_tests: PASS\n";return 0;
    }catch(const std::exception& error){std::error_code ignored;std::filesystem::remove(path,ignored);
        std::cerr<<"dynsyslab_krom_control_campaign_payload_tests: FAIL: "<<error.what()<<'\n';return 1;}}

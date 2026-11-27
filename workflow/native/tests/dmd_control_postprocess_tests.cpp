#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/EdmdControlArchive.hpp>
#include <dynsyslab/workflow/Plan.hpp>
#include <dynsyslab/workflow/Run.hpp>
#include <dynsyslab/workflow_c_api.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
namespace json=core::json;
void require(bool ok,const char* msg) { if(!ok) throw std::runtime_error(msg); }
template<class F> void rejects(F f) {
    bool failed=false; try { f(); } catch(const std::exception&) { failed=true; }
    require(failed,"invalid workflow input accepted");
}
json::Value job(const std::filesystem::path& source,const std::string& method="dmdc") {
    return json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},{"label","dmdc-test"},
        {"sources",json::Value::Array{source.generic_u8string()}},
        {"methods",json::Value::Array{json::Value::Object{
            {"name",method},{"settings",json::Value::Object{{"snapshot_array","measured_latents"}}}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}
core::NpzRecord sourceRecord(bool interval) {
    constexpr std::size_t count=120;
    std::vector<double> time(count),values(count),control(count);
    values[0]=.4;
    for(std::size_t k=0;k<count;++k) {
        time[k]=.01*k; control[k]=std::sin(.3*k)+.4*std::cos(.8*k);
        if(k+1<count) values[k+1]=.85*values[k]+.2*control[k];
    }
    // A trap: the terminal point input must never be used in identification.
    control.back()=1e20;
    if(interval) control.pop_back();
    core::NpzRecord out;
    out.arrays["time"]=core::NpyArray::fromDoubles({count},time);
    out.arrays["measured_latents"]=core::NpyArray::fromDoubles({count,1},values);
    out.arrays[interval?"control_interval_values":"control_values"]=
        core::NpyArray::fromDoubles({control.size(),1},control);
    // No system or raw state present: proves the path is measurement-only.
    out.arrays["metadata_json"]=core::NpyArray::fromString("{\"measurement\":\"echo-history latents\"}");
    return out;
}
} // namespace
int main() {
    std::filesystem::path root;
    try {
        root=std::filesystem::temp_directory_path()/("dynsyslab-dmdc-"+std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        for(const std::string method:{"dmdc","edmdc"}) for(const bool interval:{false,true}) {
            const auto source=root/(method+(interval?"-interval.npz":"-sample.npz"));
            auto record=sourceRecord(interval);
            if(interval) record.arrays["control_values"]=core::NpyArray::fromDoubles({120,1},std::vector<double>(120,999));
            core::writeNpzAtomicNew(source,record);
            const auto before=std::filesystem::file_size(source);
            const auto request=job(source,method);
            require(dynsyslab::workflow::planJob(request).at("stage_count").number()==1,"DMDc not planned");
            const auto result=dynsyslab::workflow::executeJob(request,root/"archives",{});
            const auto archive=core::readNpz(result.archives.at(0));
            const auto meta=json::parse(archive.metadataJson());
            require(meta.at("method").string()==method&&meta.at("full_row_rank").boolean(),"bad metadata");
            const std::size_t features=method=="edmdc"?5:1;
            require(meta.at(method=="edmdc"?"snapshot_one_step_relative_error":"one_step_relative_error").number(1)<1e-10,"wrong input alignment");
            if(method=="dmdc") {
                require(std::abs(archive.at("operator").doubles()[0]-.85)<1e-12,"wrong A");
                require(std::abs(std::abs(archive.at("control_matrix").doubles()[0])-.2)<1e-12,"wrong B");
            } else {
                require(meta.at("polynomial_degree").number()==4&&meta.at("include_constant").boolean()
                    &&meta.at("normalize").boolean(),"EDMDc defaults differ from SINDy");
                require(archive.at("polynomial_exponents").doubles()==std::vector<double>({0,1,2,3,4}),"dictionary missing");
                require(archive.at("state_center").shape==std::vector<std::size_t>({1})
                    &&archive.at("state_scale").doubles()[0]>0,"normalization not archived");
                require(!archive.contains("lifted_values"),"large dictionary training data duplicated");
                const auto loaded=core::loadEdmdControlModel(result.archives.at(0));
                const auto restoredValues=core::evaluateEdmdControlEigenfunctions(loaded,
                    record.at("measured_latents").doubles(),120);
                const auto savedValues=archive.at("eigenfunction_values").complexDoubles();
                require(restoredValues.size()==savedValues.size(),"EDMDc reload shape mismatch");
                for(std::size_t i=0;i<savedValues.size();++i)
                    require(std::abs(restoredValues[i]-savedValues[i])<1e-12,"EDMDc reload changed polynomial evaluation");
            }
            require(archive.at("eigenfunction_coefficients").shape==std::vector<std::size_t>({features,features}),"coefficients absent");
            require(archive.at("eigenfunction_values").shape==std::vector<std::size_t>({120,features}),"psi samples absent");
            require(meta.at("source_archive").string()==source.generic_u8string(),"missing source pointer");
            require(!archive.contains("values")&&!archive.contains("control_values"),"raw data duplicated");
            require(std::filesystem::file_size(source)==before,"source mutated");
            const auto again=dynsyslab::workflow::executeJob(request,root/"archives",{});
            require(result.archives[0]!=again.archives[0],"archive was overwritten");
            dsl_workflow_result_handle* handle=nullptr;
            const auto encoded=json::dump(request),destination=(root/"abi-archives").u8string();
            const auto status=dsl_workflow_run_json(encoded.c_str(),destination.c_str(),nullptr,nullptr,&handle);
            if(status!=DSL_STATUS_OK) {
                if(handle) dsl_workflow_result_destroy(handle);
                throw std::runtime_error("DMDc failed through the public C ABI");
            }
            std::size_t size=0;
            dsl_workflow_result_json(handle,nullptr,0,&size);
            std::vector<char> response(size);
            const auto copied=dsl_workflow_result_json(handle,response.data(),size,nullptr);
            dsl_workflow_result_destroy(handle);
            require(copied==DSL_STATUS_OK&&std::string(response.data()).find("koopman-"+method)!=std::string::npos,
                "DMDc C ABI omitted its derived model");
        }
        for(const auto badSetting:{"polynomial_degree","dictionary"}) {
            const auto source=root/"edmdc-sample.npz";
            auto request=job(source,"edmdc");
            request.at("methods").array()[0].at("settings").at(badSetting)=
                std::string(badSetting)=="dictionary"?json::Value("rbf"):json::Value(0.0);
            rejects([&]{(void)dynsyslab::workflow::executeJob(request,root/"archives",{});});
        }
        {   // time_window from each record's own metadata: a bank fitted on its active intervals
            auto record=sourceRecord(false);
            record.arrays["metadata_json"]=core::NpyArray::fromString(R"({"measurement":"echo-history latents","drive_start":0.195,"drive_stop":0.605})");
            const auto path=root/"windowed.npz"; core::writeNpzAtomicNew(path,record);
            auto request=job(path);
            request.at("methods").array()[0].at("settings").at("time_window")=
                json::Value::Object{{"from_metadata",json::Value::Array{"drive_start","drive_stop"}}};
            const auto meta=json::parse(core::readNpz(dynsyslab::workflow::executeJob(request,root/"archives",{}).archives.at(0)).metadataJson());
            require(meta.at("pairs").number()==40,"time_window from metadata selects the active interval");
        }
        {   // pool_sources fit records together; validation_sources score frozen one-step predictions
            const auto a=root/"pool-a.npz",b=root/"pool-b.npz",v=root/"pool-v.npz";
            for(const auto& path:{a,b,v}) core::writeNpzAtomicNew(path,sourceRecord(false));
            for(const std::string method:{"dmdc","edmdc"}) {
                auto request=job(a,method);
                auto& s=request.at("methods").array()[0].at("settings");
                s.at("pool_sources")=json::Value::Array{b.generic_u8string()};
                s.at("validation_sources")=json::Value::Array{v.generic_u8string()};
                const auto meta=json::parse(core::readNpz(dynsyslab::workflow::executeJob(request,root/"archives",{}).archives.at(0)).metadataJson());
                require(meta.at("pairs").number()==2*119,"pooled records add their pairs");
                require(meta.at("validation").at("one_step_relative_l2").number()<1e-8,"held-out one-step score of an exact model");
            }
        }
        for(int mode=0;mode<5;++mode) {
            auto record=sourceRecord(false);
            if(mode==0) record.arrays.erase("control_values");
            if(mode==1) record.arrays["control_values"]=core::NpyArray::fromDoubles({119,1},std::vector<double>(119,1));
            if(mode==2) record.arrays["control_values"]=core::NpyArray::fromDoubles({120,1},std::vector<double>(120));
            if(mode==3) { auto t=record.at("time").doubles(); t[20]+=.005; record.arrays["time"]=core::NpyArray::fromDoubles({120},t); }
            const auto path=root/("bad-"+std::to_string(mode)+".npz");
            core::writeNpzAtomicNew(path,record); auto request=job(path);
            if(mode==4) request.at("methods").array()[0].at("settings").at("rank")=-1.0;
            rejects([&]{(void)dynsyslab::workflow::executeJob(request,root/"archives",{});});
        }
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_dmd_control_postprocess_tests: PASS\n"; return 0;
    } catch(const std::exception& e) {
        // Test-owned, unique temporary directory only.
        if(!root.empty()) { std::error_code ignored; std::filesystem::remove_all(root,ignored); }
        std::cerr<<"dynsyslab_dmd_control_postprocess_tests: FAIL: "<<e.what()<<'\n'; return 1;
    }
}

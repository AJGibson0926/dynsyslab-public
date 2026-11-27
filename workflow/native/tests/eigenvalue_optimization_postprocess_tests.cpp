#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;namespace json=core::json;

void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}

void trajectory(const std::filesystem::path& path,double amplitude){
    constexpr std::size_t samples=24;std::vector<double> values;
    for(std::size_t sample=0;sample<samples;++sample){
        values.push_back(amplitude*std::pow(.82,static_cast<double>(sample)));
        values.push_back(0.0);
    }
    core::NpzRecord archive;
    std::vector<double> time;
    for(std::size_t sample=0;sample<samples;++sample)time.push_back(.1*sample);
    archive.arrays["time"]=core::NpyArray::fromDoubles({samples},time);
    archive.arrays["values"]=core::NpyArray::fromDoubles({samples,2},values);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(
        R"({"schema":"dynsyslab.trajectory/v1","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
    core::writeNpzAtomic(path,archive);
}

json::Value job(const std::filesystem::path& first,
                const std::filesystem::path& second,
                const std::filesystem::path& spectrum){
    return json::Value::Object{{"schema","dynsyslab.compute-job/v1"},
        {"kind","postprocess"},{"label","optimized spectrum"},
        {"sources",json::Value::Array{first.generic_u8string()}},
        {"methods",json::Value::Array{json::Value::Object{
            {"name","eigenvalue_optimization"},{"enabled",true},
            {"settings",json::Value::Object{
                {"initial_spectrum",spectrum.generic_u8string()},
                {"observable","state_0"},
                {"additional_sources",json::Value::Array{second.generic_u8string()}},
                {"mode_groups",json::Value::Array{}},{"lower_bound",-1.05},
                {"upper_bound",1.05},{"max_iterations",200.0},{"memory",8.0},
                {"tolerance",1.0e-11},{"preserve_conjugates",true}}}},
            json::Value::Object{{"name","boundary_eigenfunctions"},{"enabled",true},
                {"settings",json::Value::Object{{"spectrum",""},
                    {"observable","state_0"},
                    {"additional_sources",json::Value::Array{second.generic_u8string()}},
                    {"mode_groups",json::Value::Array{}},
                    {"interpolation_neighbors",4.0},{"gradient_relative_step",1.0e-6},
                    {"svd_tolerance",1.0e-12}}}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}

}  // namespace

int main(){std::filesystem::path root;
    try{root=std::filesystem::temp_directory_path()/("dynsyslab-eigenvalue-workflow-"
        +std::to_string(std::chrono::high_resolution_clock::now()
            .time_since_epoch().count()));
        std::filesystem::create_directories(root);const auto first=root/"first.npz";
        const auto second=root/"second.npz";trajectory(first,1.0);trajectory(second,2.0);
        const auto spectrum=root/"seed.npz";core::NpzRecord seed;
        seed.arrays["discrete_eigenvalues"]=core::NpyArray::fromComplex({1},{{.65,0.}});
        seed.arrays["metadata_json"]=core::NpyArray::fromString(
            R"({"schema":"dynsyslab.derived/v1","method":"hankel-dmd","delta_t":0.1})");
        core::writeNpzAtomic(spectrum,seed);
        const auto result=dynsyslab::workflow::executeJob(
            job(first,second,spectrum),root/"archives");
        std::vector<std::filesystem::path> kroms;  // the stage publishes one krom
        for(const auto& path:result.archives)
            if(core::dataObjectKind(core::readNpz(path).metadataJson())==core::object_kind::krom)kroms.push_back(path);
        require(result.archives.size()==2&&kroms.size()==1,"ordered spectral pipeline did not publish its stages");
        const auto archive=core::readNpz(result.archives.front());
        const auto optimized=archive.at("discrete_eigenvalues").complexDoubles();
        require(optimized.size()==1&&std::abs(optimized[0]-.82)<1.0e-6,
            "workflow did not recover the expected spectrum");
        require(archive.at("optimization_metrics").shape==
            std::vector<std::size_t>({1,5}),"optimization diagnostics are incomplete");
        const auto metadata=json::parse(archive.metadataJson());
        require(metadata.at("stage_kind").string()=="eigenvalue-extraction"
            &&!metadata.at("eigenfunctions_constructed").boolean(true),
            "spectrum-only stage boundary was not retained");
        const auto krom=core::loadKoopmanModel(kroms.front());
        require(krom.liftedDimension()==1
            &&krom.eigenfunctionRepresentation()=="scattered-interpolation",
            "optimized spectrum did not chain into boundary eigenfunctions");
        require(std::abs(std::exp(krom.continuousEigenvalues()[0]*.1)-.82)<1.0e-6,
            "downstream eigenfunction stage changed the optimized spectrum");
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_eigenvalue_optimization_postprocess_tests: PASS\n";
        return 0;
    }catch(const std::exception& error){if(!root.empty()){std::error_code ignored;
        std::filesystem::remove_all(root,ignored);}std::cerr
        <<"dynsyslab_eigenvalue_optimization_postprocess_tests: FAIL: "
        <<error.what()<<'\n';return 1;}
}

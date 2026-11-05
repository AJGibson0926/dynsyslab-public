#include <dynsyslab/core/Archive.hpp>
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
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

void trajectory(const std::filesystem::path& path,double steady,double transient){
    constexpr std::size_t samples=5;std::vector<double> values;
    for(std::size_t sample=0;sample<samples;++sample){
        values.push_back(steady+transient*std::pow(0.5,static_cast<double>(sample)));
        values.push_back(0.1*static_cast<double>(sample));}
    core::NpzRecord archive;
    archive.arrays["time"]=core::NpyArray::fromDoubles({samples},{0,.1,.2,.3,.4});
    archive.arrays["values"]=core::NpyArray::fromDoubles({samples,2},values);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(
        R"({"format":"dynsyslab.trajectory/v1","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
    core::writeNpzAtomic(path,archive);
}

json::Value job(const std::filesystem::path& first,const std::filesystem::path& second,
                const std::filesystem::path& spectrum){
    return json::Value::Object{{"schema","dynsyslab.compute-job/v1"},
        {"kind","postprocess"},{"label","boundary postprocess"},
        {"sources",json::Value::Array{first.generic_u8string()}},
        {"methods",json::Value::Array{json::Value::Object{
            {"name","boundary_eigenfunctions"},{"enabled",true},
            {"settings",json::Value::Object{
                {"spectrum",spectrum.generic_u8string()},{"observable","state_0"},
                {"additional_sources",json::Value::Array{second.generic_u8string()}},
                {"mode_groups",json::Value::Array{}},
                {"interpolation_neighbors",4.0},{"gradient_relative_step",1.0e-6},
                {"svd_tolerance",1.0e-12}}}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}
}

int main(){std::filesystem::path root;
    try{root=std::filesystem::temp_directory_path()/(
        "dynsyslab-boundary-workflow-"+std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);const auto first=root/"first.npz";
        const auto second=root/"second.npz";trajectory(first,1.0,2.0);trajectory(second,2.0,.5);
        const auto spectrumPath=root/"spectrum.npz";core::NpzRecord spectrum;
        spectrum.arrays["continuous_eigenvalues"]=core::NpyArray::fromComplex(
            {2},{{0.,0.},{std::log(.5)/.1,0.}});
        spectrum.arrays["metadata_json"]=core::NpyArray::fromString(
            R"({"format":"dynsyslab.derived/v1","method":"trajectory-eigenvalue-optimization","delta_t":0.1})");
        core::writeNpzAtomic(spectrumPath,spectrum);
        const auto result=dynsyslab::workflow::executeJob(
            job(first,second,spectrumPath),root/"archives");
        require(result.archives.size()==1,"boundary stage did not publish one KROM");
        const auto model=core::loadKoopmanModel(result.archives.front());
        const auto archive=core::readNpz(result.archives.front());
        const auto observable=archive.at("reconstruction_names").string();
        require(model.liftedDimension()==2&&model.hasReconstruction(observable)
            &&model.eigenfunctionRepresentation()=="scattered-interpolation",
            "boundary stage did not publish a usable scattered KROM");
        require(archive.at("training_states").shape==std::vector<std::size_t>({10,2})
            &&archive.at("fit_relative_residual_by_output").doubles()[0]<1.0e-11,
            "boundary stage lost its ensemble evidence");
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_boundary_postprocess_tests: PASS\n";return 0;
    }catch(const std::exception& error){if(!root.empty()){std::error_code ignored;
        std::filesystem::remove_all(root,ignored);}std::cerr
        <<"dynsyslab_boundary_postprocess_tests: FAIL: "<<error.what()<<'\n';return 1;}}

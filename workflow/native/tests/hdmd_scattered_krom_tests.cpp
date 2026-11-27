#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/DataObject.hpp>
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
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}

json::Value job(const std::filesystem::path& source) {
    return json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","H-DMD scattered KROM"},
        {"sources",json::Value::Array{source.generic_u8string()}},
        {"methods",json::Value::Array{json::Value::Object{
            {"name","hdmd"},{"enabled",true},{"settings",json::Value::Object{
                {"observable","state_0"},{"matrix_size",6.0},
                {"sampling_window",16.0},{"sample_stride",1.0},{"rank",2.0},
                {"svd_tolerance",1.0e-12},
                {"mode_sort","terminal-amplitude-energy"}}}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root=std::filesystem::temp_directory_path()/(
            "dynsyslab-hdmd-scattered-"+std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto sourcePath=root/"trajectory.npz";
        constexpr std::size_t samples=30;std::vector<double> time(samples),states;
        states.reserve(samples*2);
        for(std::size_t index=0;index<samples;++index){
            const double t=0.05*static_cast<double>(index);time[index]=t;
            states.push_back(std::cos(1.7*t)+0.2*std::cos(3.1*t));
            states.push_back(-1.7*std::sin(1.7*t)-0.62*std::sin(3.1*t));
        }
        core::NpzRecord source;
        source.arrays["time"]=core::NpyArray::fromDoubles({samples},time);
        source.arrays["values"]=core::NpyArray::fromDoubles({samples,2},states);
        source.arrays["metadata_json"]=core::NpyArray::fromString(
            R"({"schema":"dynsyslab.trajectory/v1","label":"HDMD source","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
        core::writeNpzAtomic(sourcePath,source);

        const auto run=dynsyslab::workflow::executeJob(job(sourcePath),root/"archives");
        // H-DMD publishes one data object: its K-ROM (spectrum, eigenfunctions, readouts C).
        require(run.archives.size()==1,"H-DMD did not publish its krom");
        const auto archive=core::readNpz(run.archives.front());
        require(core::dataObjectKind(archive.metadataJson())=="krom","a krom declares its kind");
        require(core::isKrom(archive)
            &&!archive.contains("training_states")
            &&archive.contains("eigenfunction_values")&&archive.contains("window_map")
            &&!archive.contains("modes")&&archive.contains("reconstruction_names"),
            "the krom holds eigenfunctions and readouts, and no copied state");
        require(core::json::parse(archive.metadataJson()).at("data_object").at("produced_by").at("method").string()=="hdmd",
                "the producing method is a provenance tag");
        const auto observable=archive.at("reconstruction_names").string().substr(0,archive.at("reconstruction_names").string().find('\n'));
        const auto model=core::loadKoopmanModel(run.archives.front());
        {   // One eigenfunction on its own evaluates a new window without H-DMD.
            const auto psi=model.eigenfunction(0).as("time-delay");
            const auto layout=core::hankelWindowLayoutFromJson(core::json::parse(archive.metadataJson())
                .at("krom").at("representations").array()[1]);
            const core::TimeSeries series(time,states,2);
            const auto window=layout.window_at(series,static_cast<std::size_t>(archive.at("training_sample_indices").doubles().back()));
            const auto value=psi.evaluate(window.data());
            const auto stored=archive.at("eigenfunction_values").complexDoubles();
            const std::size_t last=archive.at("eigenfunction_values").shape[0]-1,modes=archive.at("eigenfunction_values").shape[1];
            require(std::abs(value[0]-stored[last*modes])<1e-8*(1+std::abs(stored[last*modes])),"a stand-alone eigenfunction reproduces its value");
        }
        require(model.liftedDimension()==2
            &&model.eigenfunctionRepresentation()=="scattered-interpolation"&&!model.hasControlMap("constant")
            &&model.hasReconstruction(observable),
            "H-DMD archive did not load as a control-ready scattered KROM");
        // The state domain resolves through the input record: evaluate at a training state.
        const auto label=static_cast<std::size_t>(archive.at("training_sample_indices").doubles().front());
        const auto lifted=model.evaluate(states.data()+2*label);
        require(lifted.size()==2&&std::isfinite(lifted[0].real()),
            "interpolated H-DMD eigenfunction cannot be evaluated");
        auto jointJob=job(sourcePath);
        auto& settings=jointJob.object()["methods"].array()[0].object()["settings"].object();
        settings["observables"]=json::Value::Array{"state_0","state_1"};
        settings["rank"]=4.0;settings["window_alignment"]="end";
        settings["reconstruction_observables"]=json::Value::Array{"state"};
        const auto jointRun=dynsyslab::workflow::executeJob(jointJob,root/"joint");
        const auto joint=core::loadKoopmanModel(jointRun.archives.front());
        require(joint.hasReconstruction("position")&&joint.hasReconstruction("velocity")
            &&joint.hasReconstruction("state"),"joint KROM lost named reconstruction maps");
        const auto jointArchive=core::readNpz(jointRun.archives.front());
        require(std::abs(jointArchive.at("training_time").doubles().front()-0.75)<1e-12,
            "end-aligned joint KROM time is wrong");
        auto reconstructionJob=job(jointRun.archives.front());
        reconstructionJob.object()["methods"]=json::Value::Array{json::Value::Object{
            {"name","observable_reconstruction"},{"enabled",true},{"settings",json::Value::Object{
                {"krom",jointRun.archives.front().generic_u8string()},
                {"observable","measured_state"},{"observable_array","training_states"}}}}};
        const auto reconstructed=dynsyslab::workflow::executeJob(reconstructionJob,root/"reconstructed");
        require(core::loadKoopmanModel(reconstructed.archives.front()).hasReconstruction("measured_state"),
            "post-fit arbitrary measured observable was not retained");
        {   // The readout step publishes a child krom: the new C only, its parent by reference.
            const auto child=core::readNpz(reconstructed.archives.front());
            const auto parent=core::json::parse(child.metadataJson()).at("krom").at("parent");
            require(!parent.isNull()&&!child.contains("eigenfunction_values")&&!child.contains("window_map")
                &&child.at("reconstruction_names").string()=="measured_state","a readout krom stores only its readout");
            require(core::loadKoopmanModel(reconstructed.archives.front()).hasReconstruction("position"),
                "the child keeps its parent's readouts");
        }
        const auto matricesPath=root/"matrices.npz";core::NpzRecord matrices;
        matrices.arrays["A"]=core::NpyArray::fromDoubles({2,2},{-1,0,0,-2});
        matrices.arrays["B"]=core::NpyArray::fromDoubles({2,1},{1,1});
        matrices.arrays["C"]=core::NpyArray::fromDoubles({1,2},{1,1});
        core::writeNpzAtomic(matricesPath,matrices);
        auto analysisJob=job(matricesPath);
        analysisJob.object()["methods"]=json::Value::Array{json::Value::Object{
            {"name","linear_analysis"},{"enabled",true},{"settings",json::Value::Object{
                {"analysis_method","continuous-infinite-gramian"}}}}};
        const auto analyzed=dynsyslab::workflow::executeJob(analysisJob,root/"analyzed");
        const auto analyzedRecord=core::readNpz(analyzed.archives.front());
        require(analyzedRecord.contains("observability_gramian")&&analyzedRecord.contains("controllability_gramian"),
            "matrix diagnostics were not published");
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_hdmd_scattered_krom_tests: PASS\n";return 0;
    } catch(const std::exception& error) {
        if(!root.empty()){std::error_code ignored;std::filesystem::remove_all(root,ignored);}
        std::cerr<<"dynsyslab_hdmd_scattered_krom_tests: FAIL: "
                 <<error.what()<<'\n';return 1;
    }
}

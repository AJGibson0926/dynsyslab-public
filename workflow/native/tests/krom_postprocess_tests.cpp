#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
namespace core=dynsyslab::core;
namespace json=core::json;

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

core::KoopmanModel parentModel() {
    auto fit=core::polynomialEigenfunctions(2,{1,0,0,1},
        {{1.0,0.0},{0.0,0.0},{0.0,0.0},{1.0,0.0}},
        core::Spectrum::fromContinuousRates({{-1.0,0.0},{-2.0,0.0}}),{0.0,0.0},{1.0,1.0});
    core::SystemConfiguration system;
    system.identifier="duffing"; system.preset="double_well_reference";
    core::KoopmanModel model(std::move(fit),std::move(system));
    model.setControlMap("constant",core::ConstantControlMap{
        {{1.0,0.0},{2.0,0.0}},2,1,"parent regression"});
    model.setControlMap("bilinear",core::BilinearControlMap{
        {{3.0,0.0},{4.0,0.0},{5.0,0.0},{6.0,0.0}},
        {{7.0,0.0},{8.0,0.0}},{},2,1,"parent regression"});
    model.setReconstruction("state",{{
        {10.0,0.0},{11.0,0.0},{12.0,0.0},{13.0,0.0}},2,2,0.0,0.0,0.0});
    return model;
}

json::Value job(const std::filesystem::path& source) {
    return json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","independent-minimal-krom"},
        {"sources",json::Value::Array{source.generic_u8string()}},
        {"methods",json::Value::Array{json::Value::Object{
            {"name","krom_extract"},{"enabled",true},
            {"settings",json::Value::Object{{"selected_modes",json::Value::Array{1.0}}}}}}},
        {"archive",json::Value::Object{{"enabled",true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root=std::filesystem::temp_directory_path()/(
            "dynsyslab-krom-postprocess-"+std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto source=root/"parent.npz";
        std::map<std::string,core::NpyArray> extra;
        extra["training_time"]=core::NpyArray::fromDoubles({2},{0.0,0.1});
        extra["training_states"]=core::NpyArray::fromDoubles(
            {2,2},{1.0,0.0,0.9,-0.1});
        extra["discrete_eigenvalues"]=core::NpyArray::fromComplex(
            {2},{{std::exp(-0.1),0.0},{std::exp(-0.2),0.0}});
        extra["eigenvalues"]=extra["discrete_eigenvalues"];
        extra["source_mode_indices"]=core::NpyArray::fromDoubles({2},{431.0,433.0});
        auto parent=parentModel();
        parent.setSamples({{0.0,0.1},{0,1},{{1.0,0.0},{2.0,0.0},{0.9,0.0},{1.8,0.0}},2});
        extra.erase("training_time");  // the samples' own
        core::saveKoopmanModel(source,parent,
            R"({"schema":"dynsyslab.derived/v1","label":"Parent KROM","delta_t":0.1})",
            std::move(extra));
        const auto result=dynsyslab::workflow::executeJob(job(source),root/"archives",{});
        require(result.archives.size()==1,"KROM extraction did not publish one archive");
        const auto archive=core::readNpz(result.archives.front());
        require(archive.at("source_mode_indices").doubles()==std::vector<double>{433.0},
            "ancestral source mode identity was not archived");
        const auto reference=core::json::parse(archive.metadataJson()).at("krom").at("parent");
        require(!archive.contains("training_states") && !archive.contains("eigenfunction_values")
            && !archive.contains("discrete_eigenvalues") && !archive.contains("expansion_coefficients")
            && !archive.contains("control_constant_matrix") && reference.at("modes").array().size()==1
            && std::filesystem::equivalent(std::filesystem::u8path(reference.at("path").string()),source),
            "extracted KROM must reference its parent, not copy it");
        const auto metadata=json::parse(archive.metadataJson());
        require(metadata.at("independent_archive").boolean()
            && metadata.at("system_configuration").at("identifier").string()=="duffing",
            "extracted KROM is not self-contained");
        const auto model=core::loadKoopmanModel(result.archives.front());
        require(model.liftedDimension()==1 && model.hasControlMap("constant")
            && model.hasControlMap("bilinear") && model.hasReconstruction("state"),
            "extracted KROM lost inherited objects");
        require(std::get<core::ConstantControlMap>(model.storedControlMap("constant"))
                    .matrix.front()==std::complex<double>{2.0,0.0},
            "constant map did not reduce by output row");
        const auto& bilinear=std::get<core::BilinearControlMap>(
            model.storedControlMap("bilinear"));
        require(bilinear.operators.front()==std::complex<double>{6.0,0.0}
            && bilinear.offsets.front()==std::complex<double>{8.0,0.0},
            "bilinear map did not reduce to the principal block");
        require(model.reconstruction("state").matrix==std::vector<std::complex<double>>{
            {11.0,0.0},{13.0,0.0}},"reconstruction map did not reduce by mode");
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_krom_postprocess_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) { std::error_code ignored; std::filesystem::remove_all(root,ignored); }
        std::cerr<<"dynsyslab_krom_postprocess_tests: FAIL: "<<error.what()<<'\n'; return 1;
    }
}

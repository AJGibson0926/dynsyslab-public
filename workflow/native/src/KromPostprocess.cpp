#include "PostprocessInternal.hpp"
#include "KromEvidence.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/KromValidationPlan.hpp>

#include <filesystem>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

std::vector<std::size_t> modeIndices(const json::Value& settings) {
    std::vector<std::size_t> result;
    for (const auto& value : settings.at("selected_modes").array()) {
        const double number=value.number(-1.0);
        const auto index=static_cast<std::size_t>(number);
        if (number<0.0 || number!=static_cast<double>(index))
            throw std::invalid_argument("selected_modes must contain nonnegative integers");
        result.push_back(index);
    }
    if (result.empty())
        throw std::invalid_argument("KROM extraction requires selected_modes");
    return result;
}

json::Value extractionMetadata(const std::filesystem::path& source,
    const std::string& label,const std::vector<std::size_t>& indices,
    const dynsyslab::core::KoopmanModel& model,const json::Value& parent) {
    json::Value::Array selected;
    for (const auto index : indices) selected.emplace_back(static_cast<double>(index));
    json::Value::Array maps;
    if(model.hasControlMap("exact"))maps.emplace_back("exact");
    for (const auto& map : model.controlMaps()) maps.emplace_back(map.first);
    json::Value::Array reconstructions;
    for (const auto& item : model.reconstructions()) reconstructions.emplace_back(item.first);
    json::Value::Object backend;
    backend["name"]="native-cpp"; backend["compiled"]=true;
    backend["python_required"]=false;
    json::Value::Object metadata;
    metadata["format"]="dynsyslab.derived/v1";
    metadata["label"]=label;
    metadata["method"]="krom-extract";
    metadata["source_archive"]=source.generic_u8string();
    metadata["source_krom"]=source.generic_u8string();
    metadata["selected_modes"]=json::Value(std::move(selected));
    metadata["control_maps"]=json::Value(std::move(maps));
    metadata["reconstruction_observables"]=json::Value(std::move(reconstructions));
    metadata["independent_archive"]=true;
    const double deltaTime=parent.at("delta_t").number();
    if (deltaTime>0.0) metadata["delta_t"]=deltaTime;
    metadata["eigenvalue_domain"]=model.spectrum().domain()==dynsyslab::core::EigenvalueDomain::discreteTime
        ?"discrete-time":"continuous-time";
    metadata["execution_backend"]=json::Value(std::move(backend));
    return json::Value(std::move(metadata));
}

}  // namespace

std::filesystem::path runKromExtraction(
    const std::filesystem::path& source,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    if (!std::filesystem::is_regular_file(source))
        throw std::invalid_argument("KROM extraction source does not exist");
    const auto indices=modeIndices(settings);
    const auto parentArchive=dynsyslab::core::readNpz(source);
    const auto parentMetadata=parentArchive.metadataJson().empty()
        ? json::Value{}:json::parse(parentArchive.metadataJson());
    const auto parent=dynsyslab::core::loadKoopmanModel(source);
    auto reduced=parent.extract(indices);
    const auto directory=root/L"derived"/L"koopman-krom-extract";
    std::filesystem::create_directories(directory);
    const auto path=directory/dynsyslab::core::wide(safeName(
        source.stem().u8string()+"__krom-extract__"+timestamp())+".npz");
    const auto inherited=inheritedSourceModes(parentArchive,indices,reduced.liftedDimension());
    std::vector<double> sourceIndices;
    for (const auto index : inherited) sourceIndices.push_back(static_cast<double>(index));
    std::map<std::string,dynsyslab::core::NpyArray> extra;
    extra["source_mode_indices"]=dynsyslab::core::NpyArray::fromDoubles(
        {sourceIndices.size()},sourceIndices);
    inheritKromEvidence(parentArchive,indices,extra);
    auto meta=extractionMetadata(source,label,indices,reduced,parentMetadata);
    const auto origin=parentMetadata.at("extraction_algorithm").string(parentMetadata.at("method").string());
    meta.object()["extraction_algorithm"]=origin;
    if(origin=="edmdc") {
        const auto plan=planAutonomousKoopmanValidation(reduced,path,root);
        for(const auto& ref:dynsyslab::core::koopmanDatasetReferencesFromMetadata(plan))reduced.addDatasetReference(ref);
        meta.object()["autonomous_validation_plan"]=plan;
    }
    dynsyslab::core::saveKoopmanModel(path,reduced,json::dump(meta),std::move(extra));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail

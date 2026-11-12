#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/ArchiveIndex.hpp>

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

std::vector<double> metrics(const dynsyslab::core::KoopmanErrorStatistics& value) {
    return {static_cast<double>(value.samplePairs),value.absoluteRms,
        value.maximumAbsolute,value.relativeL2,value.medianPointwiseRelative,
        value.p95PointwiseRelative};
}

void appendMetrics(std::vector<double>& target,
                   const dynsyslab::core::KoopmanErrorStatistics& value) {
    const auto row=metrics(value);target.insert(target.end(),row.begin(),row.end());
}

std::pair<std::vector<double>,std::vector<double>> trainingBounds(
    const dynsyslab::core::NpzRecord& archive,std::size_t dimension) {
    if(archive.contains("training_minimum")&&archive.contains("training_maximum"))
        return {archive.at("training_minimum").doubles(),archive.at("training_maximum").doubles()};
    if(!archive.contains("training_states"))return {};
    const auto& states=archive.at("training_states");
    if(states.shape.size()!=2||states.shape[1]!=dimension||!states.shape[0])return {};
    const auto values=states.doubles();
    std::vector<double> minimum(dimension,std::numeric_limits<double>::infinity());
    std::vector<double> maximum(dimension,-std::numeric_limits<double>::infinity());
    for(std::size_t row=0;row<states.shape[0];++row)
        for(std::size_t axis=0;axis<dimension;++axis){
            const double value=values[row*dimension+axis];
            minimum[axis]=std::min(minimum[axis],value);
            maximum[axis]=std::max(maximum[axis],value);
        }
    return {std::move(minimum),std::move(maximum)};
}

}  // namespace

std::filesystem::path runKromValidation(
    const SourceData& inputSource,const json::Value& settings,
    const std::filesystem::path& root,const std::string& label) {
    const auto kromPath=std::filesystem::u8path(settings.at("krom").string());
    if(kromPath.empty()||!std::filesystem::is_regular_file(kromPath))
        throw std::invalid_argument("Koopman validation requires an existing KROM path");
    const std::string role=settings.at("dataset_role").string("test");
    if(role!="training"&&role!="test"&&role!="validation")
        throw std::invalid_argument("Koopman validation dataset_role must be training, validation, or test");
    std::vector<double> horizons;
    for(const auto& value:settings.at("horizons").array())horizons.push_back(value.number());
    const auto kromArchive=dynsyslab::core::readNpz(kromPath);
    auto model=dynsyslab::core::loadKoopmanModel(kromPath);
    requireMatchingUnforcedDataset(model,inputSource);
    const auto source=validationSamples(inputSource,model,settings);
    model.addDatasetReference({std::filesystem::absolute(source.path).lexically_normal().generic_u8string(),
        role,"unforced",{}});
    const auto bounds=trainingBounds(kromArchive,source.dimension);
    const auto validation=dynsyslab::core::validateEigenfunctions(
        model,source.time,source.values,source.samples,source.dimension,
        dynsyslab::core::EigenfunctionValidationInput::physicalState,horizons,
        bounds.first,bounds.second);

    std::vector<std::complex<double>> eigenvalues;
    std::vector<double> oneStep,anchored,generator,horizonValues,changeRelative,changeDefined,rateAvailable;
    for(const auto& mode:validation.modes){
        eigenvalues.push_back(mode.eigenvalue);appendMetrics(oneStep,mode.oneStep);
        appendMetrics(anchored,mode.anchoredLongTerm);
        appendMetrics(generator,mode.generator);
        changeRelative.push_back(mode.oneStepChangeRelativeL2);
        changeDefined.push_back(mode.oneStepChangeDefined?1.:0.);
        rateAvailable.push_back(mode.continuousRateAvailable?1.:0.);
    }
    if(!validation.modes.empty())
        for(const auto& item:validation.modes.front().horizons)
            horizonValues.push_back(item.horizon);
    std::vector<double> horizonRelative,horizonRms;
    for(const auto& mode:validation.modes)for(const auto& item:mode.horizons){
        horizonRelative.push_back(item.statistics.relativeL2);
        horizonRms.push_back(item.statistics.absoluteRms);
    }
    const std::size_t modes=validation.modes.size(),count=horizonValues.size();
    dynsyslab::core::NpzRecord archive;
    archive.arrays["continuous_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {modes},eigenvalues);
    archive.arrays["continuous_rate_available"]=dynsyslab::core::NpyArray::fromDoubles({modes},rateAvailable);
    archive.arrays["generator_available"]=dynsyslab::core::NpyArray::fromDoubles({modes},rateAvailable);
    archive.arrays["one_step_metrics"]=dynsyslab::core::NpyArray::fromDoubles(
        {modes,6},oneStep);
    archive.arrays["anchored_long_term_metrics"]=
        dynsyslab::core::NpyArray::fromDoubles({modes,6},anchored);
    archive.arrays["generator_metrics"]=dynsyslab::core::NpyArray::fromDoubles(
        {modes,6},generator);
    archive.arrays["one_step_change_relative_l2"]=dynsyslab::core::NpyArray::fromDoubles({modes},changeRelative);
    archive.arrays["one_step_change_defined"]=dynsyslab::core::NpyArray::fromDoubles({modes},changeDefined);
    archive.arrays["horizons"]=dynsyslab::core::NpyArray::fromDoubles(
        {count},horizonValues);
    archive.arrays["horizon_relative_l2"]=dynsyslab::core::NpyArray::fromDoubles(
        {modes,count},horizonRelative);
    archive.arrays["horizon_absolute_rms"]=dynsyslab::core::NpyArray::fromDoubles(
        {modes,count},horizonRms);
    auto validationMeta=json::parse(metadata(
        source,"koopman-validation",label,
        "\"source_krom\":\""+escaped(kromPath.generic_u8string())
        +"\",\"dataset_role\":\""+role+"\",\"mode_count\":"
        +std::to_string(modes)+",\"sample_count\":"
        +std::to_string(validation.samples)+",\"training_domain_available\":"
        +(validation.trainingDomainAvailable?"true":"false")
        +",\"outside_training_samples\":"
        +std::to_string(validation.outsideTrainingSamples)
        +",\"outside_training_fraction\":"
        +std::to_string(validation.outsideTrainingFraction)
        +",\"maximum_relative_excursion\":"
        +std::to_string(validation.maximumRelativeExcursion)
        +",\"metric_columns\":[\"sample_pairs\",\"absolute_rms\","
         "\"maximum_absolute\",\"relative_l2\","
         "\"median_pointwise_relative\",\"p95_pointwise_relative\"]"));
    validationMeta.object()["validation_dynamics"]="unforced";
    validationMeta.object()["validation_delta_t"]=uniformDeltaTime(source);
    validationMeta.object()["validation_start"]=source.time.front();
    validationMeta.object()["validation_stop"]=source.time.back();
    validationMeta.object()["source_samples"]=double(inputSource.samples);
    validationMeta.object()["sample_selection"]="existing samples only; no interpolation, source mutation, or refitting";
    validationMeta.object()["generator_method"]="finite differences of sampled eigenfunction; diagnostic for the stored continuous logarithm branch";
    validationMeta.object()["dataset_references"]=dynsyslab::core::koopmanDatasetReferencesJson(model.datasetReferences());
    // Completed models are immutable. Publish a coefficient-identical revision
    // carrying the validation pointer rather than rewriting the parent archive.
    if(settings.at("retain_validation_pointer").boolean(true)) {
        auto linkedMeta=json::parse(kromArchive.metadataJson());
        linkedMeta.object()["extraction_algorithm"]=linkedMeta.at("extraction_algorithm").string(linkedMeta.at("method").string());
        linkedMeta.object()["method"]="krom-validation-link";
        linkedMeta.object()["source_krom"]=std::filesystem::absolute(kromPath).generic_u8string();
        linkedMeta.object()["label"]=label+" / validated KROM references";
        linkedMeta.object()["dataset_references"]=validationMeta.at("dataset_references");
        std::map<std::string,dynsyslab::core::NpyArray> evidence;
        if(!bounds.first.empty()) {
            evidence["training_minimum"]=dynsyslab::core::NpyArray::fromDoubles({bounds.first.size()},bounds.first);
            evidence["training_maximum"]=dynsyslab::core::NpyArray::fromDoubles({bounds.second.size()},bounds.second);
        }
        if(kromArchive.contains("source_mode_indices"))evidence["source_mode_indices"]=kromArchive.at("source_mode_indices");
        const auto directory=root/"derived"/"koopman-validation";
        std::filesystem::create_directories(directory);
        const auto stem=safeName(source.path.stem().u8string()+"__krom-link__"+timestamp());
        std::filesystem::path linked;
        for(std::size_t suffix=0;suffix<10000;++suffix) {
            auto candidate=directory/(stem+"__"+std::to_string(suffix)+".npz");
            if(!std::filesystem::exists(candidate)){linked=std::move(candidate);break;}
        }
        if(linked.empty())throw std::runtime_error("Cannot allocate validation-linked KROM path");
        dynsyslab::core::saveKoopmanModel(linked,model,json::dump(linkedMeta),std::move(evidence));
        dynsyslab::core::indexArchiveFiles(root,{linked});
        validationMeta.object()["validated_krom"]=std::filesystem::absolute(linked).generic_u8string();
    }
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(json::dump(validationMeta));
    return publish(root,source,"koopman-validation",label,std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail

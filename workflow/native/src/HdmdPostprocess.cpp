#include "HankelAnalysis.hpp"
#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <set>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dynsyslab::workflow::postprocess_detail {
namespace {
std::size_t countSetting(const json::Value& settings,const char* name,double fallback,std::size_t minimum=0) {
    const double value=settings.at(name).number(fallback);
    if(!std::isfinite(value)||value<double(minimum)||std::floor(value)!=value
        ||value>=double((std::numeric_limits<std::size_t>::max)()))
        throw std::invalid_argument(std::string(name)+" must be a valid integer count");
    return static_cast<std::size_t>(value);
}
}
std::filesystem::path runHdmd(const SourceData& source, const json::Value& settings,
                              const std::filesystem::path& root, const std::string& label) {
    dynsyslab::core::HankelDMDOptions options;
    auto input=prepareHankel(source,settings,true);
    options.matrixSize=input.matrix.columns;
    options.samplingWindow=countSetting(settings,"sampling_window",input.matrix.window.delays,1);
    options.snapshotStride=input.matrix.window.snapshotStride;
    options.sampleStride=input.matrix.window.sampleStride;
    options.alignToWindowEnd=input.alignEnd;
    options.channelScales=input.scales;
    options.rank=countSetting(settings,"rank",0);
    options.singularValueTolerance=settings.at("svd_tolerance").number(1.0e-10);
    options.terminalAmplitudeEnergySort=settings.at("mode_sort").string("terminal-amplitude-energy")!="none";
    dynsyslab::core::NpzRecord archive;
    archiveHankelLayout(archive,input,settings.at("save_matrix").boolean());
    const auto result=dynsyslab::core::fitHankelDMD(std::move(input.matrix),input.names,uniformDeltaTime(source),options);
    archive.arrays["eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.discreteEigenvalues);
    archive.arrays["discrete_eigenvalues"] = archive.arrays["eigenvalues"];
    archive.arrays["continuous_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.continuousEigenvalues);
    archive.arrays["modes"] = dynsyslab::core::NpyArray::fromComplex({result.featureDimension, result.rank}, result.modes);
    archive.arrays["eigenfunction_values"] = dynsyslab::core::NpyArray::fromComplex(
        {result.samples,result.rank},result.coordinates);
    std::vector<double> alignedTime,alignedStates;
    alignedTime.reserve(result.samples);alignedStates.reserve(result.samples*source.dimension);
    for(std::size_t sample=0;sample<result.samples;++sample){
        const std::size_t sourceSample=result.coordinateSampleOffsets.at(sample);
        alignedTime.push_back(source.time.at(sourceSample));
        alignedStates.insert(alignedStates.end(),
            source.values.begin()+static_cast<std::ptrdiff_t>(sourceSample*source.dimension),
            source.values.begin()+static_cast<std::ptrdiff_t>((sourceSample+1)*source.dimension));
    }
    archive.arrays["training_time"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.samples},alignedTime);
    archive.arrays["training_states"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.samples,source.dimension},alignedStates);
    const auto channels=result.observableNames.size();
    archive.arrays["reconstruction_matrix"]=dynsyslab::core::NpyArray::fromComplex(
        {channels,result.rank},result.reconstructionMatrix);
    archive.arrays["koopman_modes"]=archive.arrays["reconstruction_matrix"];
    std::string reconstructionNames;
    std::set<std::string> keys;
    for(std::size_t channel=0;channel<channels;++channel) {
        const auto& name=result.observableNames[channel];
        const auto key=safeArrayKey(name);
        if(!keys.insert(key).second)throw std::invalid_argument("observable archive keys collide");
        if(!reconstructionNames.empty())reconstructionNames+="\n";
        reconstructionNames+=name;
        archive.arrays["reconstruction__"+key]=dynsyslab::core::NpyArray::fromComplex(
            {1,result.rank},std::vector<std::complex<double>>(
                result.reconstructionMatrix.begin()+channel*result.rank,
                result.reconstructionMatrix.begin()+(channel+1)*result.rank));
    }
    if(!settings.at("reconstruction_observables").array().empty()) {
        auto system=dynsyslab::core::createSystem(systemConfiguration(source));
        const auto catalog=system->observableCatalog();
        for(const auto& item:settings.at("reconstruction_observables").array()) {
            const auto name=item.string();const auto key=safeArrayKey(name);
            if(!keys.insert(key).second)throw std::invalid_argument("duplicate reconstruction name");
            const auto& observable=catalog.at(name);
            const auto real=observable.evaluateRealBatch(alignedStates,result.samples);
            const std::vector<std::complex<double>> values(real.begin(),real.end());
            const auto fitted=dynsyslab::core::fitObservableReconstruction(result.coordinates,
                result.samples,result.rank,values,observable.outputDimension(),
                settings.at("reconstruction_ridge").number());
            archive.arrays["reconstruction__"+key]=dynsyslab::core::NpyArray::fromComplex(
                {fitted.outputs,result.rank},fitted.matrix);
            archive.arrays["reconstruction_error__"+key]=dynsyslab::core::NpyArray::fromDoubles(
                {3},{fitted.residualRms,fitted.residualMaximum,fitted.relativeResidual});
            reconstructionNames+="\n"+name;
        }
    }
    archive.arrays["reconstruction_names"]=dynsyslab::core::NpyArray::fromString(reconstructionNames);
    archive.arrays["hankel_singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.singularValues.size()},result.singularValues);
    archive.arrays["terminal_modal_amplitudes"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.modalAmplitudes);
    std::vector<double> indices; for (auto value : result.sourceModeIndices) indices.push_back(static_cast<double>(value));
    archive.arrays["energy_sort_source_indices"] = dynsyslab::core::NpyArray::fromDoubles({result.rank}, indices);
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(metadata(source, "hankel-dmd", label,
        "\"eigenvalue_domain\":\"discrete-time\",\"full_spectrum_retained\":"
        + std::string(result.rank==result.singularValues.size()?"true":"false")
        + ",\"mode_sort\":\"" + std::string(options.terminalAmplitudeEnergySort ? "terminal-amplitude-energy" : "none")
        + "\",\"matrix_size\":" + std::to_string(options.matrixSize)
        + ",\"observable\":\"" + escaped(result.observableName) + "\""
        + ",\"eigenfunction_representation\":\"scattered-interpolation\""
        + ",\"gradient_backend\":\"finite-difference\",\"interpolation_neighbors\":8"
        + ",\"gradient_relative_step\":0.000001"
        + ",\"eigenfunction_provenance\":\"H-DMD trajectory coordinates with state-space interpolation\""
        + ",\"hankel_window\":"+json::dump(input.provenance)
        + ",\"delta_t\":" + json::dump(json::Value(result.deltaTime))
        + ",\"sampling_window\":" + std::to_string(options.samplingWindow)
        + ",\"matrix_rows\":" + std::to_string(result.matrixRows)
        + ",\"matrix_columns\":" + std::to_string(result.matrixColumns)
        + ",\"hankel_contract\":\"joint-hankel/v1\",\"row_layout\":\"delay-major\""
        + ",\"observables\":"+jsonStringArray(result.observableNames)
        + ",\"channel_scales\":"+jsonArray(options.channelScales)
        + ",\"snapshot_stride\":"+std::to_string(options.snapshotStride)
        + ",\"reference_delay\":"+std::to_string(result.referenceDelay)
        + ",\"delay_count\":"+std::to_string(result.delayCount)
        + ",\"sample_stride\":"+std::to_string(options.sampleStride)
        + ",\"rank\":" + std::to_string(result.rank)));
    return publish(root, source, "koopman-hankel-dmd", label, std::move(archive));
}


}

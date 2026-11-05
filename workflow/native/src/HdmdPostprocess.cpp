#include "PostprocessInternal.hpp"
#include <dynsyslab/core/Spectral.hpp>
namespace dynsyslab::workflow::postprocess_detail {
std::filesystem::path runHdmd(const SourceData& source, const json::Value& settings,
                              const std::filesystem::path& root, const std::string& label) {
    dynsyslab::core::HankelDMDOptions options;
    options.matrixSize = static_cast<std::size_t>(settings.at("matrix_size").number(2000));
    options.samplingWindow = static_cast<std::size_t>(settings.at("sampling_window").number(20000));
    options.sampleStride = static_cast<std::size_t>(settings.at("sample_stride").number(1));
    options.observableIndex = static_cast<std::size_t>(settings.at("observable").number());
    options.rank = static_cast<std::size_t>(settings.at("rank").number());
    options.singularValueTolerance = settings.at("svd_tolerance").number(1.0e-10);
    options.terminalAmplitudeEnergySort = settings.at("mode_sort").string(
        "terminal-amplitude-energy") != "none";
    const std::string observableName = settings.at("observable").string();
    dynsyslab::core::HankelDMDResult result;
    if (observableName.empty()) {
        result = dynsyslab::core::hankelDynamicModeDecomposition(
            source.values, source.samples, source.dimension, uniformDeltaTime(source), options);
    } else {
        auto system = dynsyslab::core::createSystem(systemConfiguration(source));
        result = dynsyslab::core::hankelDynamicModeDecomposition(
            source.values, source.samples, source.dimension, uniformDeltaTime(source),
            system->observableCatalog().at(observableName), options);
    }
    dynsyslab::core::NpzRecord archive;
    archive.arrays["eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.discreteEigenvalues);
    archive.arrays["discrete_eigenvalues"] = archive.arrays["eigenvalues"];
    archive.arrays["continuous_eigenvalues"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.continuousEigenvalues);
    archive.arrays["modes"] = dynsyslab::core::NpyArray::fromComplex({result.featureDimension, result.rank}, result.modes);
    archive.arrays["eigenfunction_values"] = dynsyslab::core::NpyArray::fromComplex(
        {result.samples,result.rank},result.coordinates);
    std::vector<double> alignedTime,alignedStates;
    alignedTime.reserve(result.samples);alignedStates.reserve(result.samples*source.dimension);
    for(std::size_t sample=0;sample<result.samples;++sample){
        const std::size_t sourceSample=sample*options.sampleStride;
        alignedTime.push_back(source.time.at(sourceSample));
        alignedStates.insert(alignedStates.end(),
            source.values.begin()+static_cast<std::ptrdiff_t>(sourceSample*source.dimension),
            source.values.begin()+static_cast<std::ptrdiff_t>((sourceSample+1)*source.dimension));
    }
    archive.arrays["training_time"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.samples},alignedTime);
    archive.arrays["training_states"]=dynsyslab::core::NpyArray::fromDoubles(
        {result.samples,source.dimension},alignedStates);
    archive.arrays["reconstruction_matrix"] = dynsyslab::core::NpyArray::fromComplex(
        {1,result.rank},result.reconstructionMatrix);
    archive.arrays["koopman_modes"] = archive.arrays["reconstruction_matrix"];
    archive.arrays["reconstruction_names"]=dynsyslab::core::NpyArray::fromString(
        result.observableName);
    archive.arrays["reconstruction__"+safeArrayKey(result.observableName)]=
        archive.arrays["reconstruction_matrix"];
    archive.arrays["terminal_modal_amplitudes"] = dynsyslab::core::NpyArray::fromComplex({result.rank}, result.modalAmplitudes);
    std::vector<double> indices; for (auto value : result.sourceModeIndices) indices.push_back(static_cast<double>(value));
    archive.arrays["energy_sort_source_indices"] = dynsyslab::core::NpyArray::fromDoubles({result.rank}, indices);
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(metadata(source, "hankel-dmd", label,
        "\"eigenvalue_domain\":\"discrete-time\",\"full_spectrum_retained\":true,"
        "\"mode_sort\":\"" + std::string(options.terminalAmplitudeEnergySort ? "terminal-amplitude-energy" : "none")
        + "\",\"matrix_size\":" + std::to_string(options.matrixSize)
        + ",\"observable\":\"" + escaped(result.observableName) + "\""
        + ",\"eigenfunction_representation\":\"scattered-interpolation\""
        + ",\"gradient_backend\":\"finite-difference\",\"interpolation_neighbors\":8"
        + ",\"gradient_relative_step\":0.000001"
        + ",\"eigenfunction_provenance\":\"H-DMD trajectory coordinates with state-space interpolation\""
        + ",\"delta_t\":" + std::to_string(result.deltaTime)
        + ",\"sampling_window\":" + std::to_string(options.samplingWindow)
        + ",\"matrix_rows\":" + std::to_string(result.matrixRows)
        + ",\"matrix_columns\":" + std::to_string(result.matrixColumns)
        + ",\"rank\":" + std::to_string(result.rank)));
    return publish(root, source, "koopman-hankel-dmd", label, std::move(archive));
}


}

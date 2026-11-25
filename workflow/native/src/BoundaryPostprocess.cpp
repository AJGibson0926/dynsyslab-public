#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"

#include <dynsyslab/core/BoundaryEigenfunctions.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>

#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

bool sameSystem(const dynsyslab::core::SystemConfiguration& left,
                const dynsyslab::core::SystemConfiguration& right) {
    return left.identifier==right.identifier&&left.parameters==right.parameters;
}

std::vector<std::complex<double>> spectrum(
    const dynsyslab::core::NpzRecord& archive,const json::Value& metadata,
    double& deltaTime) {
    deltaTime=metadata.at("delta_t").number();
    if(archive.contains("continuous_eigenvalues"))
        return archive.at("continuous_eigenvalues").complexDoubles();
    const auto& discrete=archive.contains("discrete_eigenvalues")
        ?archive.at("discrete_eigenvalues"):archive.at("eigenvalues");
    if(!(deltaTime>0.0))throw std::invalid_argument(
        "discrete eigenvalue spectrum requires positive delta_t");
    std::vector<std::complex<double>> result;
    for(const auto value:discrete.complexDoubles())result.push_back(std::log(value)/deltaTime);
    return result;
}

std::vector<std::vector<std::size_t>> modeGroups(const json::Value& settings) {
    std::vector<std::vector<std::size_t>> result;
    for(const auto& group:settings.at("mode_groups").array()){
        std::vector<std::size_t> row;
        for(const auto& value:group.array())row.push_back(
            static_cast<std::size_t>(value.number()));
        result.push_back(std::move(row));
    }
    return result;
}

}  // namespace

std::filesystem::path runBoundaryEigenfunctions(
    const SourceData& primary,const json::Value& settings,
    const std::filesystem::path& precedingSpectrum,
    const std::filesystem::path& root,const std::string& label) {
    auto spectrumPath=std::filesystem::u8path(settings.at("spectrum").string());
    if(spectrumPath.empty())spectrumPath=precedingSpectrum;
    if(spectrumPath.empty()||!std::filesystem::is_regular_file(spectrumPath))
        throw std::invalid_argument(
            "boundary eigenfunction fitting requires a spectrum archive");
    const auto spectralArchive=dynsyslab::core::readNpz(spectrumPath);
    const auto spectralMetadata=spectralArchive.metadataJson().empty()?json::Value{}:
        json::parse(spectralArchive.metadataJson());
    double deltaTime=0.0;const auto eigenvalues=spectrum(
        spectralArchive,spectralMetadata,deltaTime);
    const double trajectoryStep=uniformDeltaTime(primary);
    if(!(deltaTime>0.0))deltaTime=trajectoryStep;
    if(std::abs(deltaTime-trajectoryStep)>1.0e-9*std::max(1.0,trajectoryStep))
        throw std::invalid_argument(
            "boundary spectrum and trajectory sample intervals differ");

    std::vector<SourceData> sources;sources.push_back(primary);
    for(const auto& path:settings.at("additional_sources").array())
        sources.push_back(loadSource(std::filesystem::u8path(path.string())));
    const auto configuration=systemConfiguration(primary);
    dynsyslab::core::KoopmanTrajectoryEnsemble ensemble;
    ensemble.time=primary.time;ensemble.samples=primary.samples;
    ensemble.trajectories=sources.size();ensemble.stateDimension=primary.dimension;
    json::Value::Array sourcePaths;
    for(const auto& source:sources){
        if(source.samples!=ensemble.samples||source.dimension!=ensemble.stateDimension
            ||source.time!=ensemble.time
            ||!sameSystem(configuration,systemConfiguration(source)))
            throw std::invalid_argument(
                "boundary fitting sources must share system and uniform sample grid");
        ensemble.states.insert(ensemble.states.end(),source.values.begin(),source.values.end());
        sourcePaths.emplace_back(source.path.generic_u8string());
    }
    auto system=dynsyslab::core::createSystem(configuration);
    const std::string observableName=settings.at("observable").string("state_0");
    const auto observableCatalog=system->observableCatalog();
    const auto& observable=observableCatalog.at(observableName);
    dynsyslab::core::BoundaryEigenfunctionOptions options;
    options.modeGroups=modeGroups(settings);
    options.interpolationNeighbors=static_cast<std::size_t>(
        settings.at("interpolation_neighbors").number(8.0));
    options.gradientRelativeStep=settings.at("gradient_relative_step").number(1.0e-6);
    options.singularValueTolerance=settings.at("svd_tolerance").number(1.0e-12);
    auto fit=dynsyslab::core::fitBoundaryEigenfunctions(
        eigenvalues,deltaTime,ensemble,observable,options);
    dynsyslab::core::KoopmanModel model(std::move(fit.eigenfunctions),
        dynsyslab::core::Spectrum::fromContinuousRates(fit.continuousEigenvalues),configuration);
    for(const auto& source:sources)model.addDatasetReference({
        std::filesystem::absolute(source.path).lexically_normal().generic_u8string(),
        "training",datasetDynamics(source),{}});
    model.setReconstruction(observable.descriptor().name,
        std::move(fit.observableReconstruction));

    std::vector<double> repeatedTime;repeatedTime.reserve(
        ensemble.samples*ensemble.trajectories);
    for(std::size_t trajectory=0;trajectory<ensemble.trajectories;++trajectory)
        repeatedTime.insert(repeatedTime.end(),ensemble.time.begin(),ensemble.time.end());
    std::vector<std::complex<double>> discrete;
    for(const auto value:eigenvalues)discrete.push_back(std::exp(value*deltaTime));
    std::map<std::string,dynsyslab::core::NpyArray> extra;
    extra["training_time"]=dynsyslab::core::NpyArray::fromDoubles(
        {repeatedTime.size()},repeatedTime);
    extra["discrete_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {discrete.size()},discrete);extra["eigenvalues"]=extra["discrete_eigenvalues"];
    extra["fit_relative_residual_by_output"]=
        dynsyslab::core::NpyArray::fromDoubles(
            {fit.relativeResidualByOutput.size()},fit.relativeResidualByOutput);
    const auto directory=objectFolder(root, "koopman-boundary-eigenfunctions");
    std::filesystem::create_directories(directory);
    const auto output=directory/dynsyslab::core::wide(safeName(
        primary.path.stem().u8string()+"__boundary-eigenfunctions__"+timestamp())+".npz");
    const std::string details="\"source_spectrum\":\""
        +escaped(spectrumPath.generic_u8string())+"\",\"source_archives\":"
        +json::dump(json::Value(std::move(sourcePaths)))
        +",\"observable\":\""+escaped(observable.descriptor().name)
        +"\",\"trajectory_count\":"+std::to_string(ensemble.trajectories)
        +",\"delta_t\":"+std::to_string(deltaTime)
        +",\"eigenvalue_algorithm\":\""
        +escaped(spectralMetadata.at("method").string("provided"))
        +"\",\"eigenfunction_representation\":\"scattered-interpolation\""
        +",\"gradient_backend\":\"finite-difference\"";
    dynsyslab::core::saveKoopmanModel(output,model,
        metadata(primary,"fitted-boundary-eigenfunctions",label,details),std::move(extra));
    return output;
}

}  // namespace dynsyslab::workflow::postprocess_detail

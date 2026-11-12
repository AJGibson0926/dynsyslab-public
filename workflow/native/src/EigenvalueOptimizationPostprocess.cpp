#include "PostprocessInternal.hpp"

#include <dynsyslab/core/EigenvalueOptimization.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

bool sameSystem(const dynsyslab::core::SystemConfiguration& left,
                const dynsyslab::core::SystemConfiguration& right){
    return left.identifier==right.identifier&&left.parameters==right.parameters;
}

std::vector<std::vector<std::size_t>> groups(const json::Value& settings,
                                             std::size_t modes,
                                             std::size_t outputs){
    std::vector<std::vector<std::size_t>> result;
    for(const auto& group:settings.at("mode_groups").array()){
        std::vector<std::size_t> row;
        for(const auto& value:group.array()){
            const auto index=static_cast<std::size_t>(value.number());
            if(index>=modes)throw std::invalid_argument(
                "eigenvalue optimization mode group exceeds seed spectrum");
            row.push_back(index);
        }
        result.push_back(std::move(row));
    }
    if(result.empty()&&outputs==1){result.emplace_back();
        for(std::size_t index=0;index<modes;++index)result[0].push_back(index);}
    if(result.size()!=outputs)throw std::invalid_argument(
        "eigenvalue optimization needs one mode group per observable output");
    return result;
}

std::vector<std::complex<double>> discreteSpectrum(
    const dynsyslab::core::NpzRecord& archive,double dt){
    if(archive.contains("discrete_eigenvalues"))
        return archive.at("discrete_eigenvalues").complexDoubles();
    if(archive.contains("eigenvalues"))
        return archive.at("eigenvalues").complexDoubles();
    if(!archive.contains("continuous_eigenvalues"))throw std::invalid_argument(
        "seed archive contains no Koopman spectrum");
    std::vector<std::complex<double>> result;
    for(const auto value:archive.at("continuous_eigenvalues").complexDoubles())
        result.push_back(std::exp(value*dt));
    return result;
}

}  // namespace

std::filesystem::path runEigenvalueOptimization(
    const SourceData& primary,const json::Value& settings,
    const std::filesystem::path& precedingSpectrum,
    const std::filesystem::path& root,const std::string& label){
    auto seedPath=std::filesystem::u8path(settings.at("initial_spectrum").string());
    if(seedPath.empty())seedPath=precedingSpectrum;
    if(seedPath.empty()||!std::filesystem::is_regular_file(seedPath))
        throw std::invalid_argument(
            "eigenvalue optimization requires an initial spectrum archive");
    const double dt=uniformDeltaTime(primary);
    const auto seedArchive=dynsyslab::core::readNpz(seedPath);
    const auto initial=discreteSpectrum(seedArchive,dt);
    if(initial.empty())throw std::invalid_argument("initial spectrum is empty");

    std::vector<SourceData> sources{primary};
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
                "optimization sources must share system and uniform sample grid");
        ensemble.states.insert(ensemble.states.end(),source.values.begin(),source.values.end());
        sourcePaths.emplace_back(source.path.generic_u8string());
    }
    auto system=dynsyslab::core::createSystem(configuration);
    const auto catalog=system->observableCatalog();
    const auto& observable=catalog.at(settings.at("observable").string("state_0"));
    const auto modeGroups=groups(settings,initial.size(),observable.outputDimension());
    dynsyslab::core::EigenvalueOptimizationOptions options;
    for(const auto& group:modeGroups){options.initialDiscreteGroups.emplace_back();
        for(const auto index:group)options.initialDiscreteGroups.back().push_back(initial[index]);}
    options.lowerBound=settings.at("lower_bound").number(-1.05);
    options.upperBound=settings.at("upper_bound").number(1.05);
    options.maximumIterations=static_cast<std::size_t>(
        settings.at("max_iterations").number(500));
    options.memory=static_cast<std::size_t>(settings.at("memory").number(8));
    options.tolerance=settings.at("tolerance").number(1.0e-10);
    options.preserveConjugatePairs=settings.at("preserve_conjugates").boolean(true);
    const auto fit=dynsyslab::core::optimizeTrajectoryEigenvalues(
        ensemble,observable,options);

    std::vector<std::complex<double>> selectedInitial;
    std::vector<double> sourceModeIndices;
    std::vector<double> metrics,offsets{0.0};std::string messages;
    for(std::size_t output=0;output<modeGroups.size();++output){
        for(const auto index:modeGroups[output]){
            selectedInitial.push_back(initial[index]);
            sourceModeIndices.push_back(static_cast<double>(index));
        }
        offsets.push_back(static_cast<double>(selectedInitial.size()));
        const auto& diagnostic=fit.diagnostics[output];
        metrics.insert(metrics.end(),{diagnostic.normalizedSquaredResidual,
            diagnostic.relativeResidual,static_cast<double>(diagnostic.iterations),
            static_cast<double>(diagnostic.functionEvaluations),
            diagnostic.converged?1.0:0.0});
        if(!messages.empty())messages.push_back('\n');messages+=diagnostic.message;
    }
    dynsyslab::core::NpzRecord archive;
    archive.arrays["initial_discrete_eigenvalues"]=
        dynsyslab::core::NpyArray::fromComplex({selectedInitial.size()},selectedInitial);
    archive.arrays["eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {fit.discreteEigenvalues.size()},fit.discreteEigenvalues);
    archive.arrays["discrete_eigenvalues"]=archive.arrays["eigenvalues"];
    archive.arrays["continuous_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
        {fit.continuousEigenvalues.size()},fit.continuousEigenvalues);
    archive.arrays["optimization_metrics"]=dynsyslab::core::NpyArray::fromDoubles(
        {fit.diagnostics.size(),5},metrics);
    archive.arrays["mode_group_offsets"]=dynsyslab::core::NpyArray::fromDoubles(
        {offsets.size()},offsets);
    archive.arrays["source_mode_indices"]=dynsyslab::core::NpyArray::fromDoubles(
        {sourceModeIndices.size()},sourceModeIndices);
    archive.arrays["optimization_messages"]=dynsyslab::core::NpyArray::fromString(messages);
    json::Value::Array encodedGroups;
    for(const auto& group:modeGroups){json::Value::Array row;
        for(const auto index:group)row.emplace_back(static_cast<double>(index));
        encodedGroups.emplace_back(json::Value(std::move(row)));}
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(metadata(
        primary,"trajectory-eigenvalue-optimization",label,
        "\"stage_kind\":\"eigenvalue-extraction\",\"eigenfunctions_constructed\":false,"
        "\"source_spectrum\":\""+escaped(seedPath.generic_u8string())+"\","
        "\"source_archives\":"+json::dump(json::Value(std::move(sourcePaths)))+","
        "\"observable\":\""+escaped(observable.descriptor().name)+"\","
        "\"trajectory_count\":"+std::to_string(ensemble.trajectories)+","
        "\"mode_groups\":"+json::dump(json::Value(std::move(encodedGroups)))+","
        "\"delta_t\":"+std::to_string(dt)+",\"eigenvalue_domain\":\"discrete-time\","
        "\"generator_conversion\":\"lambda=log(mu)/delta_t\""));
    return publish(root,primary,"koopman-eigenvalue-optimization",label,std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail

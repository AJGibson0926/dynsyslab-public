#include <dynsyslab/workflow/EchoHistoryKoopman.hpp>

#include "PostprocessInternal.hpp"
#include "ObservationSource.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/HistoryKoopman.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::workflow {
namespace {

using postprocess_detail::SourceData;

struct ObservationRun {
    SourceData source;
    std::vector<double> time;
    std::vector<double> observation;
};

std::size_t positiveInteger(const json::Value& settings, const char* name,
                            std::size_t fallback) {
    const double value = settings.at(name).number(static_cast<double>(fallback));
    if (!(value >= 1.0) || !std::isfinite(value) || std::floor(value) != value)
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(value);
}

void requireUnforced(const SourceData& source, bool training) {
    if (source.metadataJson.empty())
        throw std::invalid_argument("echo-history data require provenance metadata");
    const auto metadata = json::parse(source.metadataJson);
    if (metadata.at("forcing_type").string() != "free"
        || metadata.at("controller_kind").string() != "none")
        throw std::invalid_argument(
            "echo-history representation data must be unforced and uncontrolled");
    if (training && metadata.at("data_role").string() != "free_training")
        throw std::invalid_argument(
            "echo-history training sources must have data_role=free_training");
}

ObservationRun loadObservation(const std::filesystem::path& path,
                               const json::Value& settings,
                               bool training) {
    auto source = postprocess_detail::loadObservationSource(path,
        settings.at("observation_array").string("echo"),
        settings.at("observable").string("far_field_acoustic_echo"));
    requireUnforced(source, training);
    if(source.dimension!=1)throw std::invalid_argument("Echo-history observations must be scalar");
    auto time=source.time,observation=source.values;
    return {std::move(source),std::move(time),std::move(observation)};
}

std::vector<double> decimate(const std::vector<double>& values,
                             std::size_t stride) {
    std::vector<double> result;
    result.reserve((values.size() + stride - 1) / stride);
    for (std::size_t index = 0; index < values.size(); index += stride)
        result.push_back(values[index]);
    return result;
}

struct CoordinateRun {
    std::vector<double> time;
    std::vector<double> values;
};

CoordinateRun fittedSuffix(const ObservationRun& run,
                           const dynsyslab::core::HistoryProjection& projection,
                           std::size_t decimation,
                           std::size_t delays,
                           std::size_t retainedRank,
                           double start,
                           double end) {
    const auto sampledTime = decimate(run.time, decimation);
    if (sampledTime.size() != projection.windows + delays - 1)
        throw std::runtime_error("history coordinate time alignment changed");
    CoordinateRun result;
    for (std::size_t window = 0; window < projection.windows; ++window) {
        const double time = sampledTime[window + delays - 1];
        if (time < start || time > end) continue;
        result.time.push_back(time);
        const auto begin = projection.coordinates.begin()
            + static_cast<std::ptrdiff_t>(window * projection.rank);
        result.values.insert(result.values.end(), begin,
            begin + static_cast<std::ptrdiff_t>(retainedRank));
    }
    if (result.time.size() < 3)
        throw std::invalid_argument("echo-history fit interval contains too few samples");
    return result;
}


std::vector<double> validationRow(
    const dynsyslab::core::KoopmanValidationResult& result) {
    const auto& statistics = result.modes.at(0).oneStep;
    return {statistics.absoluteRms, statistics.maximumAbsolute,
        statistics.relativeL2, statistics.medianPointwiseRelative,
        statistics.p95PointwiseRelative};
}

json::Value::Array pathsJson(const std::vector<std::filesystem::path>& paths) {
    json::Value::Array result;
    for (const auto& path : paths) result.emplace_back(path.generic_u8string());
    return result;
}

}  // namespace

std::filesystem::path runEchoHistoryKoopman(
    const std::vector<std::filesystem::path>& trainingSources,
    const json::Value& settings,
    const std::filesystem::path& archiveRoot,
    const std::string& label) {
    if (trainingSources.empty())
        throw std::invalid_argument("echo-history Koopman fitting needs training sources");
    const std::size_t decimation = positiveInteger(settings, "decimation", 1);
    const std::size_t delays = positiveInteger(settings, "delays", 100);
    const std::size_t hankelRank = positiveInteger(settings, "hankel_rank", 4);
    const std::size_t retainedRank = positiveInteger(settings, "coordinate_rank", 2);
    if (retainedRank > hankelRank)
        throw std::invalid_argument("coordinate_rank cannot exceed hankel_rank");
    const double fitStart = settings.at("fit_start").number(
        -std::numeric_limits<double>::infinity());
    const double fitEnd = settings.at("fit_end").number(
        std::numeric_limits<double>::infinity());
    if (!(fitEnd > fitStart))
        throw std::invalid_argument("echo-history fit_end must exceed fit_start");

    std::vector<ObservationRun> training;
    std::vector<std::vector<double>> sampledObservations;
    for (const auto& path : trainingSources) {
        training.push_back(loadObservation(path, settings, true));
        sampledObservations.push_back(decimate(training.back().observation, decimation));
    }
    const double sourceStep = postprocess_detail::uniformDeltaTime(training.front().source);
    for (const auto& run : training)
    {
        postprocess_detail::requireSameObservationSystem(training.front().source,run.source);
        if (std::abs(postprocess_detail::uniformDeltaTime(run.source) - sourceStep)
            > 1.0e-12 * std::max(1.0, std::abs(sourceStep)))
            throw std::invalid_argument(
                "all echo-history trajectories must share one output step");
    }

    dynsyslab::core::PooledHankelOptions hankelOptions;
    hankelOptions.delays = delays;
    hankelOptions.rank = hankelRank;
    hankelOptions.oversampling = positiveInteger(settings, "hankel_oversampling", 4);
    hankelOptions.subspaceIterations = positiveInteger(
        settings, "hankel_subspace_iterations", 8);
    hankelOptions.seed = positiveInteger(settings, "seed", 20261002);
    const auto basis = dynsyslab::core::fitPooledHankelBasis(
        sampledObservations, hankelOptions);
    std::vector<std::size_t> retained;for(std::size_t j=0;j<retainedRank;++j)retained.push_back(j);
    const dynsyslab::core::FrozenHistoryTransform transform(basis,sourceStep*double(decimation),retained,
        settings.at("observable").string("far_field_acoustic_echo"),settings.at("time_units").string("model-time"),
        settings.at("observation_units").string("unspecified"));

    std::vector<CoordinateRun> fitRuns;
    std::vector<std::vector<double>> fitValues;
    for (std::size_t index = 0; index < training.size(); ++index) {
        const auto projection = dynsyslab::core::projectHistory(
            sampledObservations[index], basis);
        fitRuns.push_back(fittedSuffix(training[index], projection,
            decimation, delays, retainedRank, fitStart, fitEnd));
        fitValues.push_back(fitRuns.back().values);
    }
    dynsyslab::core::PooledPolynomialEDMDOptions edmdOptions;
    edmdOptions.degree = positiveInteger(settings, "polynomial_degree", 5);
    edmdOptions.pairStride = positiveInteger(settings, "pair_stride", 1);
    edmdOptions.sampleTime = sourceStep * static_cast<double>(decimation);
    edmdOptions.relativeRidge = settings.at("relative_ridge").number(1.0e-10);
    const auto edmd = dynsyslab::core::fitPooledPolynomialEDMD(
        fitValues, retainedRank, edmdOptions);
    dynsyslab::core::EDMDModeSearchOptions search;
    search.minimumDiscreteMagnitude = settings.at("minimum_discrete_magnitude").number(0.90);
    search.maximumDiscreteMagnitude = settings.at("maximum_discrete_magnitude").number(1.005);
    search.minimumContinuousReal = settings.at("minimum_continuous_real").number(-1.0);
    search.maximumContinuousReal = settings.at("maximum_continuous_real").number(0.1);
    search.minimumContinuousImaginary = settings.at("minimum_continuous_imaginary").number(4.0);
    search.maximumContinuousImaginary = settings.at("maximum_continuous_imaginary").number(7.0);
    const auto candidates = dynsyslab::core::selectPolynomialEDMDModes(
        edmd, fitValues, search);
    if (candidates.empty())
        throw std::runtime_error("no admissible echo-history Koopman mode was found");
    const auto eigenfunction = dynsyslab::core::selectedPolynomialEigenfunction(
        edmd, candidates.front());
    auto bank=dynsyslab::core::polynomialEigenfunctionsFromLegacy(eigenfunction);
    auto coordinateContract=bank.coordinates();coordinateContract.physicalStateDimension.reset();
    coordinateContract.providerIdentity=transform.identity();
    // Preserve the same coefficients/rates, while retaining their discrete origin.
    bank=dynsyslab::core::PolynomialEigenfunctionSet(bank.dictionary(),bank.coefficients(),
        dynsyslab::core::EigenfunctionSpectrum::fromStoredRates({candidates.front().continuousEigenvalue},edmd.pairTime,
            dynsyslab::core::EigenvalueDomain::discreteTime),coordinateContract);
    dynsyslab::core::EigenfunctionFitDiagnostics diagnostics;
    diagnostics.algorithm="pooled-polynomial-edmd";diagnostics.sampleCount=eigenfunction.samples;
    const auto physical=postprocess_detail::observationSystem(training.front().source);
    auto model=physical?dynsyslab::core::KoopmanModel(bank,*physical,diagnostics):
        dynsyslab::core::KoopmanModel(bank,positiveInteger(settings,"input_dimension",1),diagnostics);

    std::vector<double> trainingCoordinates, trainingTime, offsets;
    std::vector<double> trainingValidation;
    for (const auto& run : fitRuns) {
        offsets.push_back(static_cast<double>(trainingTime.size()));
        trainingCoordinates.insert(trainingCoordinates.end(),
            run.values.begin(), run.values.end());
        trainingTime.insert(trainingTime.end(), run.time.begin(), run.time.end());
        const auto validation = dynsyslab::core::validateKoopmanSemigroup(
            model, run.time, run.values, run.time.size(), retainedRank,
            {0.01, 0.1, 0.5, 1.0, 2.0, 4.0});
        const auto row = validationRow(validation);
        trainingValidation.insert(trainingValidation.end(), row.begin(), row.end());
    }

    std::vector<std::filesystem::path> validationPaths;
    std::vector<double> heldoutValidation;
    for (const auto& value : settings.at("validation_sources").array()) {
        const auto path = std::filesystem::u8path(value.string());
        validationPaths.push_back(path);
        const auto run = loadObservation(path, settings, false);
        postprocess_detail::requireSameObservationSystem(training.front().source,run.source);
        if (std::abs(postprocess_detail::uniformDeltaTime(run.source) - sourceStep)
            > 1.0e-12 * std::max(1.0, std::abs(sourceStep)))
            throw std::invalid_argument("held-out echo history uses a different output step");
        const auto observation = decimate(run.observation, decimation);
        const auto projection = dynsyslab::core::projectHistory(observation, basis);
        const auto coordinates = fittedSuffix(run, projection, decimation,
            delays, retainedRank, fitStart, fitEnd);
        const auto validation = dynsyslab::core::validateKoopmanSemigroup(
            model, coordinates.time, coordinates.values, coordinates.time.size(),
            retainedRank, {0.01, 0.1, 0.5, 1.0, 2.0, 4.0});
        const auto row = validationRow(validation);
        heldoutValidation.insert(heldoutValidation.end(), row.begin(), row.end());
    }

    dynsyslab::core::NpzRecord archive;
    archive.arrays["eigenvalues"]=dynsyslab::core::NpyArray::fromComplex({1},model.spectrum().discreteMultipliers());
    archive.arrays["history_decimation"] = dynsyslab::core::NpyArray::fromDoubles(
        {1}, {static_cast<double>(decimation)});
    archive.arrays["history_sample_time"] = dynsyslab::core::NpyArray::fromDoubles(
        {1}, {edmdOptions.sampleTime});
    archive.arrays["edmd_operator"] = dynsyslab::core::NpyArray::fromDoubles(
        {eigenfunction.featureCount, eigenfunction.featureCount},
        edmd.operatorMatrix);
    archive.arrays["edmd_feature_scale"] = dynsyslab::core::NpyArray::fromDoubles(
        {eigenfunction.featureCount}, edmd.featureScale);
    archive.arrays["training_states"] = dynsyslab::core::NpyArray::fromDoubles(
        {trainingTime.size(), retainedRank}, trainingCoordinates);
    archive.arrays["training_time"] = dynsyslab::core::NpyArray::fromDoubles(
        {trainingTime.size()}, trainingTime);
    archive.arrays["training_trajectory_offsets"] =
        dynsyslab::core::NpyArray::fromDoubles({offsets.size()}, offsets);
    archive.arrays["training_validation_metrics"] =
        dynsyslab::core::NpyArray::fromDoubles({fitRuns.size(), 5}, trainingValidation);
    if (!heldoutValidation.empty())
        archive.arrays["heldout_validation_metrics"] =
            dynsyslab::core::NpyArray::fromDoubles({validationPaths.size(), 5},
                heldoutValidation);

    auto metadata = json::parse(postprocess_detail::metadata(
        training.front().source, "echo-history-koopman", label,
        "\"eigenvalue_domain\":\"discrete-time\"" )).object();
    metadata["coordinate_provider"] = "causal-history";
    metadata["eigenfunction_representation"] = "analytic-polynomial";
    metadata["delta_t"] = edmd.pairTime;
    metadata["sample_time"] = edmd.pairTime;
    metadata["training_archives"] = pathsJson(trainingSources);
    metadata["validation_archives"] = pathsJson(validationPaths);
    std::vector<dynsyslab::core::KoopmanDatasetReference> datasetRefs;
    for(const auto& path:trainingSources) dynsyslab::core::mergeKoopmanDatasetReference(datasetRefs,
        {std::filesystem::absolute(path).lexically_normal().generic_u8string(),"training","unforced",{}});
    for(const auto& path:validationPaths) dynsyslab::core::mergeKoopmanDatasetReference(datasetRefs,
        {std::filesystem::absolute(path).lexically_normal().generic_u8string(),"validation","unforced",{}});
    metadata["dataset_references"]=dynsyslab::core::koopmanDatasetReferencesJson(datasetRefs);
    metadata["data_role"] = "reduced_order_model";
    metadata["history_observation"] = settings.at("observable").string(
        "far_field_acoustic_echo");
    metadata["history_delays"] = static_cast<double>(delays);
    metadata["history_decimation"] = static_cast<double>(decimation);
    metadata["history_coordinate_rank"] = static_cast<double>(retainedRank);
    metadata["fit_start"] = std::isfinite(fitStart)?json::Value(fitStart):json::Value{};
    metadata["fit_end"] = std::isfinite(fitEnd)?json::Value(fitEnd):json::Value{};
    metadata["polynomial_degree"] = static_cast<double>(edmdOptions.degree);
    metadata["pair_stride"] = static_cast<double>(edmdOptions.pairStride);
    metadata["native_history_pipeline"] = true;
    metadata["selected_mean_pair_residual"]=candidates.front().meanPairResidual;
    metadata["selected_minimum_phase_r_squared"]=candidates.front().minimumPhaseRSquared;
    dynsyslab::core::writeFrozenHistoryTransform(archive,transform);
    archive=dynsyslab::core::makeKoopmanArchive(model,json::dump(metadata),std::move(archive.arrays));
    return postprocess_detail::publish(archiveRoot, training.front().source,
        "echo-history-koopman", label, std::move(archive));
}

}  // namespace dynsyslab::workflow

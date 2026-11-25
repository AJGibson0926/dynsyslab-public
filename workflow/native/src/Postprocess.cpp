#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/EstimationPostprocess.hpp>
#include <dynsyslab/workflow/DecompositionPostprocess.hpp>
#include <dynsyslab/workflow/EchoHistoryKoopman.hpp>
#include <dynsyslab/workflow/HistoryPostprocess.hpp>
#include "HankelShorthand.hpp"
#include "CoordinateChangeKinds.hpp"
#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "HankelAnalysis.hpp"
#include "ObservationSource.hpp"
#include "KromEvidence.hpp"

#include <dynsyslab/core/Text.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <optional>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow {
namespace postprocess_detail {

std::string jsonArray(const std::vector<double>& values) {
    std::ostringstream out; out << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) out << ','; out << std::setprecision(17) << values[index];
    }
    out << ']'; return out.str();
}

std::string jsonStringArray(const std::vector<std::string>& values) {
    std::ostringstream out; out << '[';
    for (std::size_t index=0;index<values.size();++index) {
        if (index) out << ','; out << '"' << escaped(values[index]) << '"';
    }
    out << ']'; return out.str();
}

std::string systemIdentifier(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char item) {
        return static_cast<char>(std::tolower(item));
    });
    if (value.find("rayleigh") != std::string::npos || value == "rpe") return "rpe";
    if (value.find("duffing") != std::string::npos) return "duffing";
    if (value.find("marmottant") != std::string::npos) return "marmottant";
    if (value.find("reddy") != std::string::npos) return "reddy_szeri";
    for (const auto& descriptor : dynsyslab::core::systemCatalog())
        if (value == descriptor.identifier) return value;
    return {};
}

dynsyslab::core::SystemConfiguration systemConfiguration(const SourceData& source) {
    if (source.metadataJson.empty()) throw std::invalid_argument(
        "source archive lacks the system identity required by this operation");
    const auto metadata = json::parse(source.metadataJson);
    const auto& encoded=metadata.at("system_configuration");
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier=systemIdentifier(encoded.at("identifier").string(metadata.at("system").string()));
    configuration.preset=encoded.at("preset").string(metadata.at("system_preset").string(metadata.at("preset").string()));
    const auto copyNumbers=[&](const json::Value& values) {
        for(const auto& item:values.object())
            if(std::holds_alternative<double>(item.second.storage()))configuration.parameters[item.first]=item.second.number();
    };
    copyNumbers(metadata.at("system_parameters").at("preset"));
    copyNumbers(metadata.at("system_parameters"));
    copyNumbers(metadata.at("resolved_system_parameters"));
    copyNumbers(encoded.at("parameters"));
    const auto& driver=encoded.at("driver").isNull()?metadata.at("driver"):encoded.at("driver");
    configuration.driver={driver.at("amplitude").number(),driver.at("frequency_hz").number(),driver.at("phase_radians").number()};
    configuration.equationsVersion=encoded.at("equations_version").string(metadata.at("equations_version").string());
    if(configuration.identifier.empty())throw std::invalid_argument("source system is not native");
    return configuration;
}

std::vector<double> exactVelocities(const SourceData& source, bool includeInputTerms) {
    auto system=dynsyslab::core::createSystem(systemConfiguration(source));
    if(system->descriptor().stateDimension!=source.dimension)
        throw std::invalid_argument("source state dimension differs from model");
    std::vector<double> result(source.values.size());
    if(!source.inputs.empty()&&source.inputDimension!=1)
        throw std::invalid_argument("system-RHS training currently requires scalar control input");
    for(std::size_t sample=0;sample<source.samples;++sample) {
        const auto x=source.values.data()+sample*source.dimension;
        const auto f=result.data()+sample*source.dimension;
        if(!includeInputTerms)system->rhsWithInput(source.time[sample],x,0.,f);
        else if(source.inputs.empty())system->rhs(source.time[sample],x,f);
        else system->rhsWithInput(source.time[sample],x,
            source.inputs[std::min(sample,source.inputs.size()-1)],f);
    }
    return result;
}

std::filesystem::path runControlRegression(
    const SourceData& source, const json::Value& settings,
    const std::filesystem::path& root, const std::string& label) {
    const std::string kromPath = settings.at("krom").string();
    if (kromPath.empty()) throw std::invalid_argument("control regression requires a KROM path");
    std::vector<std::size_t> selected;
    for (const auto value : settings.at("selected_modes").array())
        selected.push_back(static_cast<std::size_t>(value.number()));
    auto model = dynsyslab::core::loadKoopmanModel(   // readouts may be separate reconstructions objects
            dynsyslab::core::wide(kromPath), systemConfiguration(source), selected, composedBinding(settings));
    std::vector<SourceData> sources{source};
    const auto referenceConfiguration=systemConfiguration(source);
    const auto referenceFingerprint=dynsyslab::core::equationsFingerprint(
        dynsyslab::core::createSystem(referenceConfiguration)->configuration());
    for(const auto& value:settings.at("additional_sources").array()){
        const auto path=std::filesystem::absolute(
            std::filesystem::u8path(value.string()));
        if(std::any_of(sources.begin(),sources.end(),[&](const SourceData& item){
            return std::filesystem::absolute(item.path)==path;}))continue;
        auto additional=loadSource(path);
        const auto configuration=systemConfiguration(additional);
        const auto fingerprint=dynsyslab::core::equationsFingerprint(
            dynsyslab::core::createSystem(configuration)->configuration());
        if(additional.dimension!=source.dimension||fingerprint!=referenceFingerprint)
            throw std::invalid_argument(
                "all control-regression trajectories must share identical dynamics");
        sources.push_back(std::move(additional));
    }
    std::vector<dynsyslab::core::ControlledData> data;
    for(const auto& item:sources)data.push_back({item.time,item.values,item.inputs,
        item.samples,item.dimension,item.inputDimension});
    const double ridge = settings.at("ridge").number(0.0);
    const bool conjugate = settings.at("conjugate_complete").boolean(false);
    const std::string constantMethod=settings.at("constant_method").string(
        "multistep-output");
    const std::string bilinearMethod=settings.at("bilinear_method").string(
        "controlled-residual");
    const double strideValue=settings.at("map_sample_stride").number(1);
    if(!std::isfinite(strideValue)||strideValue<1.0||std::floor(strideValue)!=strideValue)
        throw std::invalid_argument("map_sample_stride must be a positive integer");
    const std::size_t stride=static_cast<std::size_t>(strideValue);
    std::vector<double> mapStates;
    for(const auto& item:sources)for(std::size_t sample=0;sample<item.samples;sample+=stride)
        mapStates.insert(mapStates.end(),item.values.begin()+sample*item.dimension,
                         item.values.begin()+(sample+1)*item.dimension);
    if(mapStates.size()<3*source.dimension)
        throw std::invalid_argument("control-map regression needs at least three sampled states");
    const auto constant=[&]{
        if(constantMethod=="exact-map")return dynsyslab::core::fitConstantExactControlMap(
            model,mapStates,mapStates.size()/source.dimension);
        if(std::any_of(sources.begin(),sources.end(),[](const SourceData& item){
            return item.inputs.empty();}))throw std::invalid_argument(
            "selected constant regression requires an archived control history");
        if(constantMethod=="controlled-residual")
            return dynsyslab::core::fitConstantControlResidual(model,data,ridge);
        if(constantMethod!="multistep-output")throw std::invalid_argument(
            "constant_method must be multistep-output, controlled-residual, or exact-map");
        std::string reconstruction=settings.at("reconstruction").string("state");
        if(reconstruction.empty()&&!model.reconstructions().empty())
            reconstruction=model.reconstructions().begin()->first;
        if(!model.hasReconstruction(reconstruction))throw std::invalid_argument(
            "multistep output regression requires the selected KROM reconstruction C");
        auto system=dynsyslab::core::createSystem(model.systemConfiguration());
        return dynsyslab::core::fitConstantControlMultistepOutput(model,data,
            system->observableCatalog().at(reconstruction),
            model.reconstruction(reconstruction),ridge);
    }();
    const auto bilinear=[&]{
        if(bilinearMethod=="exact-map")return dynsyslab::core::fitBilinearExactControlMap(
            model,mapStates,mapStates.size()/source.dimension,ridge,conjugate);
        if(bilinearMethod!="controlled-residual")throw std::invalid_argument(
            "bilinear_method must be controlled-residual or exact-map");
        if(std::any_of(sources.begin(),sources.end(),[](const SourceData& item){
            return item.inputs.empty();}))throw std::invalid_argument(
            "controlled-residual bilinear regression requires an archived control history");
        return dynsyslab::core::fitBilinearControlResidual(model,data,ridge,conjugate);
    }();
    model.setControlMap("constant", constant.representation);
    model.setControlMap("bilinear", bilinear.representation);
    double x0 = source.values[0], x1 = x0, y0 = source.values[1], y1 = y0;
    for(const auto& item:sources)for(std::size_t sample=0;sample<item.samples;++sample){
        x0=std::min(x0,item.values[sample*item.dimension]);
        x1=std::max(x1,item.values[sample*item.dimension]);
        y0=std::min(y0,item.values[sample*item.dimension+1]);
        y1=std::max(y1,item.values[sample*item.dimension+1]);
    }
    const auto bounds = settings.at("bounds").array();
    if (bounds.size() == 4) {
        x0 = bounds[0].number(); x1 = bounds[1].number();
        y0 = bounds[2].number(); y1 = bounds[3].number();
    } else {
        const double xPad = std::max(1e-9, 0.05 * (x1-x0));
        const double yPad = std::max(1e-9, 0.05 * (y1-y0));
        x0 -= xPad; x1 += xPad; y0 -= yPad; y1 += yPad;
    }
    const auto columns = static_cast<std::size_t>(settings.at("grid_columns").number(101));
    const auto rows = static_cast<std::size_t>(settings.at("grid_rows").number(101));
    const auto constantError = dynsyslab::core::evaluateControlMapError2D(
        model, "constant", x0, x1, y0, y1, columns, rows);
    const auto bilinearError = dynsyslab::core::evaluateControlMapError2D(
        model, "bilinear", x0, x1, y0, y1, columns, rows);
    std::map<std::string, dynsyslab::core::NpyArray> extra;
    const auto sourceKromArchive=dynsyslab::core::readNpz(
        dynsyslab::core::wide(kromPath));
    inheritKromEvidence(sourceKromArchive,selected,extra);
    const auto sourceModes=inheritedSourceModes(
        sourceKromArchive,selected,model.liftedDimension());
    std::vector<double> encodedSourceModes;
    for (const auto index:sourceModes)
        encodedSourceModes.push_back(static_cast<double>(index));
    extra["source_mode_indices"]=dynsyslab::core::NpyArray::fromDoubles(
        {encodedSourceModes.size()},encodedSourceModes);
    extra["control_error_axis0"] = dynsyslab::core::NpyArray::fromDoubles({columns}, constantError.axis0);
    extra["control_error_axis1"] = dynsyslab::core::NpyArray::fromDoubles({rows}, constantError.axis1);
    extra["constant_absolute_error"] = dynsyslab::core::NpyArray::fromDoubles({rows,columns}, constantError.absoluteError);
    extra["constant_relative_error"] = dynsyslab::core::NpyArray::fromDoubles({rows,columns}, constantError.relativeError);
    extra["bilinear_absolute_error"] = dynsyslab::core::NpyArray::fromDoubles({rows,columns}, bilinearError.absoluteError);
    extra["bilinear_relative_error"] = dynsyslab::core::NpyArray::fromDoubles({rows,columns}, bilinearError.relativeError);
    const auto directory = objectFolder(root, "koopman-control-regression");
    const auto path = directory / dynsyslab::core::wide(safeName(source.path.stem().u8string()
        + "__control-regression__" + timestamp()) + ".npz");
    const std::string details = "\"source_krom\":\"" + escaped(kromPath)
        + "\",\"constant_regression_method\":\""+escaped(constantMethod)
        + "\",\"bilinear_regression_method\":\""+escaped(bilinearMethod)
        + "\",\"reconstruction_observable\":\""
        + escaped(settings.at("reconstruction").string("state"))
        + "\",\"controlled_trajectory_count\":"+std::to_string(sources.size())
        + ",\"constant_relative_residual\":" + std::to_string(constant.diagnostics.relativeResidual)
        + ",\"bilinear_relative_residual\":" + std::to_string(bilinear.diagnostics.relativeResidual)
        + ",\"constant_error_mean\":" + std::to_string(constantError.relativeMean)
        + ",\"bilinear_error_mean\":" + std::to_string(bilinearError.relativeMean)
        + ",\"control_map_time_domain\":\"continuous-time\"";
    dynsyslab::core::savePolynomialKoopmanModel(path, model,
        metadata(source, "koopman-control-regression", label, details), std::move(extra));
    return path;
}

namespace {
// Methods that read all sources together (one archive per stage, not one per source).
bool pooledMethod(const std::string& name) {
    return name == "echo_history_koopman" || name == "history_basis" || name == "history_control_regression"
        || name == "input_map" || name == "bilinear_certificate" || name == "control_target" || name == "derivative_benchmark"
        || name == "mode_selection" || name == "krom_recenter" || name == "krom_readout"
        || name == "reference_signal" || name == "control_score" || name == "figure_dataset" || name == "io_delay_model" || name == "animation";
}
// SINDy fits on the rows it is given. "sample_stride" and "time_window" select
// those rows explicitly so the fit can use the same clock as the H-DMD stage
// (whose Hankel sample_stride does not otherwise propagate).
SourceData sindyRows(const SourceData& full, const json::Value& settings) {
    const double strideValue = settings.at("sample_stride").number(1);
    if (!std::isfinite(strideValue) || strideValue < 1 || std::floor(strideValue) != strideValue)
        throw std::invalid_argument("SINDy sample_stride must be a positive integer");
    const auto stride = static_cast<std::size_t>(strideValue);
    const auto& window = settings.at("time_window").array();
    if (!window.empty() && window.size() != 2)
        throw std::invalid_argument("SINDy time_window must be [start, stop]");
    if (stride == 1 && window.empty()) return full;
    // Measured derivative arrays are matched to the kept rows by time (coordinate SINDy).
    const double start = window.empty() ? -INFINITY : window[0].number();
    const double stop = window.empty() ? INFINITY : window[1].number();
    SourceData rows = full; rows.time.clear(); rows.values.clear(); rows.inputs.clear();
    std::size_t kept = 0;
    for (std::size_t sample = 0; sample < full.samples; ++sample) {
        if (full.time[sample] < start || full.time[sample] > stop) continue;
        if (kept++ % stride) continue;
        rows.time.push_back(full.time[sample]);
        rows.values.insert(rows.values.end(), full.values.begin() + sample * full.dimension,
                           full.values.begin() + (sample + 1) * full.dimension);
        if (!full.inputs.empty())
            rows.inputs.insert(rows.inputs.end(), full.inputs.begin() + sample * full.inputDimension,
                               full.inputs.begin() + (sample + 1) * full.inputDimension);
    }
    rows.samples = rows.time.size();
    if (rows.samples < 3) throw std::invalid_argument("SINDy row selection leaves fewer than three samples");
    return rows;
}
}  // namespace

std::filesystem::path runSindy(const SourceData& givenSource, const json::Value& givenSettings,
                               const std::filesystem::path& spectrumPath,
                               const std::filesystem::path& root, const std::string& label) {
    const auto input = analyticCoordinateInput(givenSource, givenSettings);
    const auto& fullSource = input.source;
    const auto& settings = input.settings;
    const auto selectedRows = sindyRows(fullSource, settings);
    const auto& source = selectedRows;
    const auto requestedSpectrum=settings.at("spectrum_archive").string();
    const auto activeSpectrum=requestedSpectrum.empty()?spectrumPath:std::filesystem::u8path(requestedSpectrum);
    if (activeSpectrum.empty()) throw std::invalid_argument(
        "SINDy requires a preceding native spectrum in the ordered pipeline");
    const auto spectrum = dynsyslab::core::readNpz(activeSpectrum);
    std::vector<std::complex<double>> continuous;
    if (spectrum.contains("continuous_eigenvalues"))
        continuous = spectrum.at("continuous_eigenvalues").complexDoubles();
    else {
        const auto discrete = spectrum.at("eigenvalues").complexDoubles();
        const double dt = uniformDeltaTime(source);
        for (const auto value : discrete) continuous.push_back(std::log(value) / dt);
    }
    dynsyslab::core::KoopmanSINDyOptions options;
    options.polynomialDegree = static_cast<std::size_t>(settings.at("polynomial_degree").number(4));
    options.threshold = settings.at("threshold").number(0.05);
    options.maximumIterations = static_cast<std::size_t>(settings.at("max_iterations").number(20));
    options.normalizeStates = settings.at("normalize").boolean(true);
    options.includeConstant = settings.at("include_constant").boolean(true);
    options.reuseConjugateFits = settings.at("reuse_conjugate_fits").boolean(true);
    options.materializeTrainingValues = settings.at("compute_training_values").boolean(false);
    const auto nullspace = settings.at("nullspace_solver").string("gram");
    if (nullspace != "gram" && nullspace != "svd")
        throw std::invalid_argument("SINDy nullspace_solver must be gram or svd");
    options.gramNullspace = nullspace == "gram";
    const auto sourceMetadata=source.metadataJson.empty()?json::Value{}:json::parse(source.metadataJson);
    if(sourceMetadata.at("value_space").string()=="history-coordinates" || settings.at("input_space").string()=="coordinates")
        return runCoordinateSindy(source,settings,activeSpectrum,continuous,options,root,label);
    const bool includeInputTerms=settings.at("include_input_terms").boolean(true);
    const auto velocities = settings.at("gamma_method").string("system_rhs") == "finite_difference"
        ? std::vector<double>{} : exactVelocities(source,includeInputTerms);
    auto system=dynsyslab::core::createSystem(systemConfiguration(source));
    const auto catalog=system->observableCatalog();
    std::vector<std::string> coordinateNames;
    for (const auto& value : settings.at("coordinate_observables").array())
        coordinateNames.push_back(value.string());
    const std::string coordinateName=settings.at("coordinate_observable").string();
    dynsyslab::core::Observable coordinates = !settings.at("coordinate_expression").isNull()
        ? dynsyslab::core::Observable::fromExpressionDefinition(settings.at("coordinate_expression"))
        : !coordinateNames.empty()
        ? catalog.coordinates(coordinateNames,"requested_coordinates")
        : (!coordinateName.empty() && catalog.coordinatePresets().count(coordinateName)
            ? catalog.coordinatePreset(coordinateName)
            : catalog.at(coordinateName.empty() ? "state" : coordinateName));
    const auto result = dynsyslab::core::fitKoopmanEigenfunctionsInCoordinates(
        source.values, source.samples, source.dimension, source.time, continuous,
        coordinates, velocities, options);
    std::vector<std::complex<double>> discrete;
    const double dt = uniformDeltaTime(source);
    for (const auto value : continuous) discrete.push_back(std::exp(value * dt));
    const bool archiveDerivatives=settings.at("archive_coordinate_derivatives").boolean(false);
    if (archiveDerivatives && (velocities.empty() || !coordinates.hasAnalyticJacobian()))
        throw std::invalid_argument("coordinate derivative evidence requires system_rhs and an analytic coordinate Jacobian");
    // SINDy publishes data objects only: its eigenfunctions (spectrum + one polynomial expansion,
    // each eigenfunction's coefficients its own row) and the reconstructions written over them. The
    // training series is referenced by its origin, never copied.
    namespace c = dynsyslab::core;
    c::NpzRecord archive, readouts;
    c::DataOrigin origin{"sindy", settings, {dataReference(fullSource.path, "time_series")}};
    origin.inputs.push_back(dataReference(activeSpectrum, "spectrum"));
    auto fit = c::polynomialEigenfunctionsFromLegacy(result);
    c::EigenfunctionSamples samples;
    samples.modes = continuous.size();
    if (!result.trainingValues.empty()) {
        samples.values = result.trainingValues;
        for (std::size_t k = 0; k < source.samples; ++k) {  // rows of the input series, matched by time
            const auto found = std::lower_bound(fullSource.time.begin(), fullSource.time.end(), source.time[k] - 1e-12);
            samples.sampleIndices.push_back(std::size_t(found - fullSource.time.begin()));
            samples.time.push_back(source.time[k]);
        }
    }
    const c::EigenfunctionSet set(c::Spectrum::fromDiscrete(discrete, continuous, dt, origin), samples, origin,
        {std::make_shared<c::SharedRepresentation>(std::make_shared<const c::FunctionExpansion>(fit.functions))});
    auto objects = c::writeEigenfunctionSet(archive, set, timeSeriesRecord(fullSource).channels);
    archive.arrays["sample_interval"] = c::NpyArray::fromDoubles({1}, {dt});
    archive.arrays["training_residual_rms"] = c::NpyArray::fromDoubles({continuous.size()}, result.residualRms);
    archive.arrays["training_residual_maximum"] = c::NpyArray::fromDoubles({continuous.size()}, result.residualMaximum);
    if (archiveDerivatives) {  // requested evidence: the coordinates' exact rates and Jacobians at the samples
        archive.arrays["training_coordinate_velocities"]=c::NpyArray::fromDoubles(
            {source.samples,result.stateDimension},
            coordinates.pushForwardRealBatch(source.values,velocities,source.samples,
                c::ObservableJacobianBackend::analytic));
        std::vector<double> gradients;
        gradients.reserve(source.samples*result.stateDimension*source.dimension);
        for (std::size_t row=0;row<source.samples;++row)
            for (const auto value:coordinates.jacobian(source.values.data()+row*source.dimension,
                    c::ObservableJacobianBackend::analytic))
                gradients.push_back(value.real());
        archive.arrays["training_coordinate_jacobians"]=c::NpyArray::fromDoubles(
            {source.samples,result.stateDimension,source.dimension},gradients);
    }
    dynsyslab::core::KoopmanModel fittedModel(result,systemConfiguration(source));
    std::vector<std::string> reconstructionNames;
    for (const auto& value : settings.at("reconstruction_observables").array())
        reconstructionNames.push_back(value.string());
    // "reconstruction":"none" skips the all-mode least-squares C fit (samples x
    // retained modes), which dominates SINDy runtime for a full H-DMD spectrum
    // and is unused by state-Jacobian/coordinate-target control.
    const bool skipReconstruction = settings.at("reconstruction").string() == "none";
    if (skipReconstruction && !reconstructionNames.empty())
        throw std::invalid_argument("reconstruction \"none\" conflicts with reconstruction_observables");
    if (reconstructionNames.empty() && !skipReconstruction) reconstructionNames.push_back("state");
    std::string encodedReconstructions;
    for (const auto& name : reconstructionNames) {
        const auto projection=dynsyslab::core::fitObservableReconstruction(
            fittedModel,source.values,source.samples,catalog.at(name),
            settings.at("reconstruction_ridge").number(0.0));
        readouts.arrays["reconstruction__"+safeArrayKey(name)]=
            dynsyslab::core::NpyArray::fromComplex(
                {projection.outputs,projection.modes},projection.matrix);
        if (!encodedReconstructions.empty()) encodedReconstructions.push_back('\n');
        encodedReconstructions+=name;
    }
    if (!encodedReconstructions.empty())
        readouts.arrays["reconstruction_names"]=dynsyslab::core::NpyArray::fromString(encodedReconstructions);
    const std::string gamma = velocities.empty() ? "finite_difference" : "system_rhs";
    const std::vector<double> multipliers(result.stateDimension, 1.0);
    auto sindyMeta = json::parse(metadata(source,
        "koopman-with-sindy", label, "\"eigenvalue_domain\":\"discrete-time\","
        "\"generator_eigenvalue_domain\":\"continuous-time\",\"generator_conversion\":"
        "\"lambda=log(mu)/delta_t\",\"delta_t\":" + std::to_string(dt)
        + ",\"polynomial_degree\":" + std::to_string(options.polynomialDegree)
        + ",\"feature_count\":" + std::to_string(result.featureCount)
        + ",\"threshold\":" + std::to_string(options.threshold)
        + ",\"sample_stride\":" + std::to_string(static_cast<long long>(settings.at("sample_stride").number(1)))
        + ",\"training_samples\":" + std::to_string(source.samples)
        + ",\"source_samples\":" + std::to_string(fullSource.samples)
        + ",\"gamma_method\":\"" + gamma + "\",\"state_preprocessing\":{\"center\":"
        + jsonArray(result.center) + ",\"scale\":" + jsonArray(result.scale)
        + ",\"multiplier\":" + jsonArray(multipliers)
        + "},\"eigenfunction_representation\":\"analytic-polynomial\","
        "\"coordinate_observables\":" + jsonStringArray(result.coordinateObservableNames) + ","
        "\"reconstruction_observables\":" + jsonStringArray(reconstructionNames) + ","
        "\"source_spectrum_archive\":\"" + escaped(activeSpectrum.generic_u8string()) + "\","
        "\"include_input_terms\":" + (includeInputTerms?"true":"false") + ","
        "\"coordinate_derivative_evidence\":\"" + (archiveDerivatives?"analytic-push-forward/v1":"none") + "\","
        "\"coordinate_gradient_backend\":\"" + (coordinates.hasAnalyticJacobian()?"analytic":"finite-difference") + "\","
        "\"analytic_gradient_available\":" + (coordinates.hasAnalyticJacobian()?"true":"false") + ",\"direct_eigenvalue_fits\":"
        + std::to_string(result.directFits)));
    sindyMeta.object()["eigenfunction_objects"] = std::move(objects);
    sindyMeta.object()["object_properties"] = json::Value::Object{
        {"per_eigenfunction", json::Value::Array{"training_residual_rms", "training_residual_maximum"}}};
    archive.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(json::dump(sindyMeta));
    const auto eigenfunctions = publish(root, source, "sindy-eigenfunctions", label, std::move(archive));
    if (readouts.contains("reconstruction_names")) {
        auto readoutMeta = json::parse(metadata(source, "sindy-reconstructions", label, ""));
        readoutMeta.object()["reconstructions"] = c::reconstructionsBlock(eigenfunctions, continuous.size());
        readouts.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(readoutMeta));
        registerAdditionalObject(publish(root, source, "reconstructions", label, std::move(readouts)));
    }
    return eigenfunctions;
}

}  // namespace postprocess_detail

std::wstring PostprocessResult::summary() const {
    std::wostringstream out; out << L"Native post-processing completed\r\n\r\nStages: ";
    for (std::size_t index = 0; index < methods.size(); ++index) {
        if (index) out << L", "; out << dynsyslab::core::wide(methods[index]);
    }
    out << L"\r\n\r\nDerived archives:";
    for (const auto& path : archives) out << L"\r\n" << path.wstring();
    return out.str();
}

namespace {
// A step's objects, primary first, each with its data-object envelope.
void addStepObjects(PostprocessResult& result, const std::filesystem::path& primary, const std::string& method) {
    for (const auto& path : postprocess_detail::collectStepObjects(primary, method)) {
        result.archives.push_back(path);
        result.methods.push_back(method);
    }
}
}  // namespace

PostprocessResult executePostprocess(const json::Value& submitted,
    const std::filesystem::path& archiveRoot, const std::filesystem::path& progressPath,
    const dynsyslab::core::ExecutionContext* execution) {
    using namespace postprocess_detail;
    // Observables named in a Hankel method become an explicit coordinate change.
    const auto job = expandHankelCoordinateShorthand(submitted);
    if (execution) execution->report({"postprocess", "Loading source archives", 0.0, 0, 1});
    std::vector<std::filesystem::path> sources;
    for (const auto& value : job.at("sources").array())
        sources.emplace_back(dynsyslab::core::wide(value.string()));
    std::vector<const json::Value*> stages;
    for (const auto& stage : job.at("methods").array()) if (stage.at("enabled").boolean(true)) stages.push_back(&stage);
    if (sources.empty() || stages.empty()) throw std::invalid_argument("post-processing needs sources and enabled stages");
    const auto pooledStageCount = static_cast<std::size_t>(std::count_if(
        stages.begin(), stages.end(), [](const json::Value* stage) {
            const auto name=stage->at("name").string();
            return pooledMethod(name);
        }));
    const std::size_t ordinaryStageCount = stages.size() - pooledStageCount;
    PostprocessResult result; const std::size_t planned =
        sources.size() * ordinaryStageCount + pooledStageCount;
    std::size_t complete = 0; atomicProgress(progressPath, 0.0, 0, planned, "Starting native post-processing");
    for (const auto* stage : stages) {
        const auto name=stage->at("name").string();
        if (!pooledMethod(name)) continue;
        for (const auto& path : sources) requireTrainingDataRole(name, path, stage->at("settings"));
        const std::string label = job.at("label").string("native-postprocess")
            + " / " + name;
        const StageScope scope(name, stage->at("settings"));
        const auto output = name=="history_basis"?runHistoryBasis(sources,stage->at("settings"),archiveRoot,label):
            name=="history_control_regression"?runHistoryControlRegression(sources,stage->at("settings"),archiveRoot,label):
            name=="input_map"?runInputMap(sources,stage->at("settings"),archiveRoot,label):
            name=="bilinear_certificate"?runBilinearCertificate(sources,stage->at("settings"),archiveRoot,label):
            name=="control_target"?runControlTarget(sources,stage->at("settings"),archiveRoot,label):
            name=="derivative_benchmark"?runDerivativeBenchmark(sources,stage->at("settings"),archiveRoot,label):
            name=="mode_selection"?runModeSelection(sources,stage->at("settings"),archiveRoot,label):
            name=="krom_recenter"?runKromRecenter(sources,stage->at("settings"),archiveRoot,label):
            name=="krom_readout"?runKromReadout(sources,stage->at("settings"),archiveRoot,label):
            name=="reference_signal"?runReferenceSignal(sources,stage->at("settings"),archiveRoot,label):
            name=="control_score"?runControlScore(sources,stage->at("settings"),archiveRoot,label):
            name=="figure_dataset"?runFigureDataset(sources,stage->at("settings"),archiveRoot,label):
            name=="animation"?runAnimation(sources,stage->at("settings"),archiveRoot,label):
            name=="io_delay_model"?runIoDelayModel(sources,stage->at("settings"),archiveRoot,label):
            runEchoHistoryKoopman(sources, stage->at("settings"), archiveRoot, label);
        addStepObjects(result, output, name);
        ++complete;
        atomicProgress(progressPath, static_cast<double>(complete) / planned,
            complete, planned, "Completed "+name);
        if (execution) execution->report({"postprocess",
            "Completed "+name, static_cast<double>(complete) / planned,
            complete, planned});
    }
    for (const auto& sourcePath : sources) {
      std::optional<SourceData> loadedSource;
      std::filesystem::path precedingHankel;
      std::filesystem::path precedingSpectrum;
      std::filesystem::path precedingCoordinates;
      std::filesystem::path precedingHankelInput; // the series the stored matrix was built from
      for (const auto* stage : stages) {
        if (execution) execution->checkpoint();
        const std::string method = stage->at("name").string();
        const auto& originalSettings=stage->at("settings");
        std::optional<json::Value> inherited;
        if (originalSettings.at("use_previous_hankel").boolean()) {
            if (method!="hdmd" && method!="hankel_dmd" && method!="hankel_svd" && method!="hankel_pod"
                && method!="koopman_spectral" && method!="resdmd" && method!="mpedmd"
                && method!="hermitian_dmd" && method!="tlsdmd" && method!="multdmd" && method!="stochastic_resdmd" && method!="specrkhs")
                throw std::invalid_argument("use_previous_hankel requires a Hankel decomposition");
            if (precedingHankel.empty() || !originalSettings.at("hankel_archive").isNull())
                throw std::invalid_argument("select one preceding Hankel matrix for this trajectory");
            inherited=originalSettings;
            inherited->at("hankel_archive")=precedingHankel.generic_u8string();
        }
        const auto& settings=inherited?*inherited:originalSettings;
        if (pooledMethod(method)) continue;
        {   // Eigenfunctions/spectra from unforced data only; B_psi regressions from forced data only.
            const bool previous = settings.at("use_previous_output").boolean(false) && !precedingCoordinates.empty();
            requireTrainingDataRole(method, previous ? precedingCoordinates : sourcePath, settings);
        }
        const std::string label = job.at("label").string("native-postprocess") + " / " + method;
        const StageScope scope(method, settings);
        std::filesystem::path output;
        if (method == "linear_analysis" || method == "linear_systems" || method == "controllability" || method == "observable_reconstruction" || method == "krom_extract" || method == "dmdc" || method == "edmdc" || method=="history_projection" || method=="eigenfunction_expansion" || method=="eigenfunction_errors" || method=="domain_field"
            || method=="klqe" || method=="observability" || method=="empirical_observability") {
            if(method=="history_projection")output=runHistoryProjection(sourcePath,settings,archiveRoot,label);
            else if(method=="eigenfunction_expansion")output=runEigenfunctionExpansion(sourcePath,settings,archiveRoot,label);
            else if(method=="domain_field")output=runDomainField(sourcePath,settings,archiveRoot,label);
            else if(method=="eigenfunction_errors")output=runEigenfunctionErrors(sourcePath,settings,archiveRoot,label);
            else if(method=="linear_analysis"||method=="controllability")output=runLinearAnalysis(sourcePath,settings,archiveRoot,label);
            else if(method=="linear_systems")output=runLinearSystems(sourcePath,settings,archiveRoot,label);
            else if(method=="empirical_observability")output=runEmpiricalObservability(sourcePath,settings,archiveRoot,label);
            else if(method=="observable_reconstruction")output=runObservableReconstruction(sourcePath,settings,archiveRoot,label);
            else if(method=="klqe"||method=="observability")output=runEstimationPostprocess(method,sourcePath,settings,archiveRoot,label);
            else output = method != "krom_extract"
                ? runDmdControl(sourcePath, settings, archiveRoot, label, method == "edmdc")
                : runKromExtraction(sourcePath, settings, archiveRoot, label);
            addStepObjects(result, output, method); ++complete;
            atomicProgress(progressPath, static_cast<double>(complete) / planned, complete, planned,
                           "Completed " + method);
            if (execution) execution->report({"postprocess", "Completed " + method,
                static_cast<double>(complete) / planned, complete, planned});
            continue;
        }
        if(method=="sindy_dynamics" || method=="coordinate_transform" || method=="dynamics_validation") {
            const bool previous=settings.at("use_previous_output").boolean(false);
            if(previous&&precedingCoordinates.empty())
                throw std::invalid_argument("No preceding coordinate dataset to consume");
            const auto path=previous?precedingCoordinates:sourcePath;
            const auto input=method=="coordinate_transform"?loadCoordinateInput(path,settings):loadSource(path);
            output=method=="coordinate_transform"?runCoordinateTransform(input,settings,archiveRoot,label):
                method=="dynamics_validation"?runDynamicsValidation(input,settings,archiveRoot,label):
                runSindyDynamics(input,settings,archiveRoot,label);
            if(method=="coordinate_transform")precedingCoordinates=output;
            addStepObjects(result, output, method);++complete;
            atomicProgress(progressPath,double(complete)/planned,complete,planned,"Completed "+method);
            if(execution)execution->report({"postprocess","Completed "+method,double(complete)/planned,complete,planned});
            continue;
        }
        // A method may read the preceding coordinate change instead of the source.
        const bool previousOutput = settings.at("use_previous_output").boolean(false);
        if (previousOutput && precedingCoordinates.empty())
            throw std::invalid_argument(method + " use_previous_output needs a preceding coordinate_transform");
        std::optional<SourceData> coordinateInput;
        if (previousOutput) coordinateInput = loadSource(precedingCoordinates);
        else if (inherited && !precedingHankelInput.empty() && precedingHankelInput != sourcePath)
            coordinateInput = loadSource(precedingHankelInput);
        const bool readsCoordinates = coordinateInput.has_value();
        if (!readsCoordinates && !loadedSource) loadedSource = loadSource(sourcePath);
        const auto& source = readsCoordinates ? *coordinateInput : *loadedSource;
        if (method=="hankel_matrix") {
            output=runHankelMatrix(source,settings,archiveRoot,label);
            precedingHankel=output;
            precedingHankelInput=previousOutput?precedingCoordinates:sourcePath;
        }
        else if(method=="koopman_spectral"||method=="resdmd"||method=="mpedmd"||method=="hermitian_dmd"||method=="tlsdmd"||method=="multdmd"||method=="stochastic_resdmd"||method=="specrkhs")
            output=runKoopmanSpectral(source,settings,archiveRoot,label,method);
        else if (method == "dmd") output = runDmd(source, settings, archiveRoot, label);
        else if (method == "hdmd" || method == "hankel_dmd") output = runHdmd(source, settings, archiveRoot, label);
        else if (method == "sindy") output = runSindy(source, settings, precedingSpectrum, archiveRoot, label);
        else if (method == "control_regression")
            output = runControlRegression(source, settings, archiveRoot, label);
        else if (method == "koopman_validation")
            output = runKromValidation(source, settings, archiveRoot, label);
        else if (method == "boundary_eigenfunctions")
            output = runBoundaryEigenfunctions(
                source,settings,precedingSpectrum,archiveRoot,label);
        else if (method == "eigenvalue_optimization")
            output = runEigenvalueOptimization(
                source,settings,precedingSpectrum,archiveRoot,label);
        else if (method == "edmd" || method == "havok")
            output = runDecompositionPostprocess(
                method, source.path, settings, archiveRoot, label);
        else if (method == "hankel_pod" || method == "hankel_svd")
            output = runHankelSvd(method, source, settings, archiveRoot, label);
        else if (method == "svd" || method == "pca" || method == "pod"
                 || method == "bagged_pca" || method == "robust_pca")
            output = runPcaPostprocess(method, source, settings, archiveRoot, label);
        else {
            const std::string extension = stage->at("extension").string();
            if (extension.empty())
                throw std::runtime_error(method
                    + " has not yet been ported to the native post-processing core");
            output = runExtensionPostprocess(
                source, settings, archiveRoot, label, method, extension);
        }
        addStepObjects(result, output, method);
        if (method == "dmd" || method == "hdmd" || method == "hankel_dmd"
            || method == "eigenvalue_optimization") precedingSpectrum = output;
        ++complete;
        atomicProgress(progressPath, static_cast<double>(complete) / planned, complete, planned,
                       "Completed " + method);
        if (execution) execution->report({"postprocess", "Completed " + method,
            static_cast<double>(complete) / planned, complete, planned});
      }
    }
    return result;
}

}  // namespace dynsyslab::workflow

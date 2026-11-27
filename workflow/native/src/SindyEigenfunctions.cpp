// The "sindy" post-processing method: polynomial Koopman eigenfunctions fitted to the generator
// equation grad(psi)·f = lambda psi over a monomial dictionary (SINDy: Brunton, Proctor & Kutz, PNAS
// 113:3932, 2016; Koopman eigenfunctions as its null vectors: Kaiser, Kutz & Brunton, Mach. Learn.:
// Sci. Technol. 2:035023, 2021). What it is given decides what it publishes (see
// docs/architecture/EIGENFUNCTION_OBJECTS.md, "Adding representations with SINDy"):
//   eigenvalues (settings), a spectrum, or training data alone -> new eigenfunctions objects;
//   eigenfunctions objects or a K-ROM -> the same eigenfunctions with a calibrated polynomial
//   representation added (an extension object; the input's own representations are kept by reference).
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/SpectrumArchive.hpp>
#include "CoordinateChangeKinds.hpp"
#include "ObjectPublishing.hpp"
#include "SampledDomain.hpp"
#include "SindyTraining.hpp"
#include "DictionaryTerms.hpp"
#include "TimeSeriesSource.hpp"

#include <dynsyslab/core/DataObject.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/SINDy.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
namespace {
using Complex = std::complex<double>;

// What SINDy was given, in order of precedence: "eigenvalues" (a list, or "from_data"), "krom",
// "eigenfunctions", "spectrum_archive", the preceding method's object; nothing at all means from data.
struct Given {
    enum class Kind { eigenvalues, spectrum, krom, data } kind = Kind::data;
    std::filesystem::path path;            // the object read
    std::vector<Complex> rates;            // eigenvalues/spectrum
    std::vector<std::size_t> modes;        // eigenfunctions/krom: which of the object's modes
};

std::vector<std::size_t> indexList(const json::Value& values, const char* what) {
    std::vector<std::size_t> out;
    for (const auto& value : values.array()) {
        const double number = value.number(-1.0);
        if (!(number >= 0) || std::floor(number) != number) throw std::invalid_argument(std::string(what) + " must be nonnegative integers");
        out.push_back(std::size_t(number));
    }
    return out;
}

Given resolveGiven(const json::Value& settings, const std::filesystem::path& preceding, double dt) {
    Given given;
    const auto& eigenvalues = settings.at("eigenvalues");
    if (std::holds_alternative<std::string>(eigenvalues.storage())) {
        if (eigenvalues.string() != "from_data") throw std::invalid_argument("SINDy eigenvalues must be a list or \"from_data\"");
        return given;
    }
    if (!eigenvalues.isNull()) {
        given.kind = Given::Kind::eigenvalues;
        for (const auto& pair : eigenvalues.array()) {
            const auto& parts = pair.array();
            given.rates.emplace_back(parts.at(0).number(), parts.size() > 1 ? parts[1].number() : 0.0);
        }
        if (given.rates.empty()) throw std::invalid_argument("SINDy eigenvalues is empty");
        return given;
    }
    const auto named = [&](const char* key) { return settings.at(key).string(); };
    const std::string path = !named("krom").empty() ? named("krom") : !named("eigenfunctions").empty() ? named("eigenfunctions")
                           : !named("spectrum_archive").empty() ? named("spectrum_archive") : preceding.u8string();
    if (path.empty()) return given;  // from data
    given.path = std::filesystem::u8path(path);
    const auto archive = c::readNpz(given.path);
    if (c::isKrom(archive)) { given.kind = Given::Kind::krom; return given; }  // its spectrum and eigenfunctions
    given.kind = Given::Kind::spectrum;  // any archive with eigenvalues
    given.rates = c::readSpectrumArrays(archive).continuousRates();
    (void)dt;
    return given;
}

const char* kindName(Given::Kind kind) {
    switch (kind) {
        case Given::Kind::eigenvalues: return "eigenvalues";
        case Given::Kind::spectrum: return "spectrum";
        case Given::Kind::krom: return "krom";
        default: return "data";
    }
}

std::vector<std::size_t> selectModes(const json::Value& settings, std::size_t count, std::vector<std::size_t> given) {
    if (!settings.at("modes").isNull()) given = indexList(settings.at("modes"), "SINDy modes");
    if (given.empty()) { given.resize(count); for (std::size_t j = 0; j < count; ++j) given[j] = j; }
    for (const auto mode : given) if (mode >= count) throw std::invalid_argument("SINDy mode index out of range");
    return given;
}

// p_j at the given dynamical states (rows x modes), through the fit's coordinates.
std::vector<Complex> polynomialValues(const c::EigenfunctionFitResult& result, const c::Observable& coordinates,
                                      const std::vector<double>& states, std::size_t rows) {
    const auto q = coordinates.evaluateRealBatch(states, rows);
    const std::size_t m = result.spectrum.size(), d = coordinates.outputDimension();
    std::vector<Complex> out;
    out.reserve(rows * m);
    for (std::size_t k = 0; k < rows; ++k) {
        const auto values = result.functions.evaluate(q.data() + k * d);
        out.insert(out.end(), values.begin(), values.end());
    }
    return out;
}

// Multiplies each mode's coefficients by its factor. The generator residual and its maximum are linear
// in psi, so they scale by |factor| and keep describing the published eigenfunction.
void scaleModes(c::EigenfunctionFitResult& result, const std::vector<Complex>& factors) {
    auto coefficients = result.functions.coefficients();
    const std::size_t features = result.functions.dictionary().size();
    for (std::size_t j = 0; j < factors.size(); ++j) {
        for (std::size_t f = 0; f < features; ++f) coefficients[j * features + f] *= factors[j];
        if (auto residual = result.spectrum.residual(j)) {
            residual->value *= std::abs(factors[j]);
            result.spectrum.setResidual(j, std::move(residual));
        }
        if (j < result.diagnostics.residualMaximum.size()) result.diagnostics.residualMaximum[j] *= std::abs(factors[j]);
    }
    result.functions = c::FunctionExpansion(result.functions.dictionaryHandle(), std::move(coefficients),
                                            result.functions.coordinates(), result.functions.outputNames());
}

// Calibration against reference values psi (rows x modes): c_j = p_j^H psi_j / p_j^H p_j, and the
// agreement ||c_j p_j - psi_j|| / ||psi_j|| after it. Without references, or for a mode whose
// polynomial disagrees with its references beyond `tolerance` (the references then describe another
// function, so their factor is arbitrary): unit RMS over the normalization rows q and a real positive
// value at the first of them (a convention; eigenfunctions are defined up to a factor). `apply` false
// keeps every factor at 1 (an anchored fit is already on its references).
struct Calibration { std::vector<Complex> scale; std::vector<double> mismatch; std::vector<bool> referenced; };
Calibration calibrate(c::EigenfunctionFitResult& result, const std::vector<Complex>& p, const std::vector<Complex>* psi,
                      std::size_t rows, const std::vector<Complex>& q, std::size_t normalizationRows,
                      double tolerance, bool apply) {
    const std::size_t m = result.spectrum.size();
    Calibration out{std::vector<Complex>(m, 1.0), std::vector<double>(m, 0.0), std::vector<bool>(m, false)};
    for (std::size_t j = 0; j < m; ++j) {
        Complex cross{}; double pp = 0, refNorm = 0;
        if (psi)
            for (std::size_t k = 0; k < rows; ++k) {
                const auto pk = p[k * m + j];
                pp += std::norm(pk);
                cross += std::conj(pk) * (*psi)[k * m + j]; refNorm += std::norm((*psi)[k * m + j]);
            }
        Complex factor = 1.0;
        if (psi && pp > 0) {
            factor = cross / pp;
            if (refNorm > 0) {
                double error = 0;
                for (std::size_t k = 0; k < rows; ++k) error += std::norm(factor * p[k * m + j] - (*psi)[k * m + j]);
                out.mismatch[j] = std::sqrt(error / refNorm);
            }
            out.referenced[j] = out.mismatch[j] <= tolerance;
        }
        if (!out.referenced[j]) {
            double qq = 0;
            for (std::size_t k = 0; k < normalizationRows; ++k) qq += std::norm(q[k * m + j]);
            factor = qq > 0 ? std::polar(1.0 / std::sqrt(qq / double(normalizationRows)), -std::arg(q[j])) : Complex(1.0);
        }
        out.scale[j] = apply ? factor : Complex(1.0);
    }
    scaleModes(result, out.scale);
    return out;
}

std::vector<double> realParts(const std::vector<Complex>& values) { std::vector<double> out; for (auto v : values) out.push_back(v.real()); return out; }

}  // namespace

std::filesystem::path runSindy(const SourceData& givenSource, const json::Value& givenSettings,
                               const std::filesystem::path& spectrumPath,
                               const std::filesystem::path& root, const std::string& label) {
    const auto input = analyticCoordinateInput(givenSource, givenSettings);
    const auto& fullSource = input.source;
    const auto& settings = input.settings;
    const auto selectedRows = sindyRows(fullSource, settings);
    const auto& source = selectedRows;
    const double dt = uniformDeltaTime(source);
    c::KoopmanSINDyOptions options;
    options.polynomialDegree = static_cast<std::size_t>(settings.at("polynomial_degree").number(4));
    options.threshold = settings.at("threshold").number(0.05);
    options.maximumIterations = static_cast<std::size_t>(settings.at("max_iterations").number(20));
    options.normalizeStates = settings.at("normalize").boolean(true);
    options.includeConstant = settings.at("include_constant").boolean(true);
    options.reuseConjugateFits = settings.at("reuse_conjugate_fits").boolean(true);
    const bool trainingValues = settings.at("compute_training_values").boolean(false);
    const auto nullspace = settings.at("nullspace_solver").string("gram");
    if (nullspace != "gram" && nullspace != "svd") throw std::invalid_argument("SINDy nullspace_solver must be gram or svd");
    options.gramNullspace = nullspace == "gram";
    // Options of the fit itself (core/SINDy.hpp KoopmanSINDyOptions).
    if (!settings.at("rescale").isNull()) {
        if (std::holds_alternative<double>(settings.at("rescale").storage())) options.rescale = {settings.at("rescale").number()};
        else for (const auto& value : settings.at("rescale").array()) options.rescale.push_back(value.number());
    }
    options.sparsification = settings.at("sparsification_policy").string("iterative-relative");
    options.directionalDerivative = settings.at("directional_derivative").string("analytic");
    options.directionalStep = settings.at("directional_step").number(1e-6);
    auto given = resolveGiven(settings, spectrumPath, dt);
    const auto sourceMetadata = source.metadataJson.empty() ? json::Value{} : json::parse(source.metadataJson);
    if (sourceMetadata.at("value_space").string() == "history-coordinates" || settings.at("input_space").string() == "coordinates") {
        if (given.kind == Given::Kind::krom)  // coordinate SINDy reads a spectrum
            given.rates = c::loadKoopmanModel(given.path).spectrum().continuousRates();
        if (given.rates.empty()) throw std::invalid_argument("coordinate SINDy needs eigenvalues, a spectrum or eigenfunctions");
        return runCoordinateSindy(source, settings, given.path, given.rates, options, root, label);
    }
    const bool includeInputTerms = settings.at("include_input_terms").boolean(true);
    const auto gammaMethod = settings.at("gamma_method").string("system_rhs");
    if (gammaMethod != "system_rhs" && gammaMethod != "finite_difference")
        throw std::invalid_argument("SINDy gamma_method must be system_rhs or finite_difference");
    // One trajectory, several pooled with their boundaries kept, or a state-space domain (SindyTraining.hpp).
    const auto training = sindyTraining(source, settings, gammaMethod == "system_rhs", includeInputTerms);
    const auto& velocities = training.velocities;
    if (!training.oneTrajectory()) options.segments = training.segments;
    auto system = c::createSystem(systemConfiguration(source));
    const auto catalog = system->observableCatalog();
    std::vector<std::string> coordinateNames;
    for (const auto& value : settings.at("coordinate_observables").array()) coordinateNames.push_back(value.string());
    const std::string coordinateName = settings.at("coordinate_observable").string();
    c::Observable coordinates = !settings.at("coordinate_expression").isNull()
        ? c::Observable::fromExpressionDefinition(settings.at("coordinate_expression"))
        : !coordinateNames.empty() ? catalog.coordinates(coordinateNames, "requested_coordinates")
        : (!coordinateName.empty() && catalog.coordinatePresets().count(coordinateName)
            ? catalog.coordinatePreset(coordinateName) : catalog.at(coordinateName.empty() ? "state" : coordinateName));
    // A dictionary other than the polynomial library, on the fit's coordinates.
    options.dictionary = jobDictionary(settings.at("dictionary"), coordinates.outputDimension(), coordinates.descriptor().outputNames);

    // The eigenvalues, and for eigenfunctions or a K-ROM the eigenfunctions themselves.
    json::Value::Object estimateRecord;
    std::optional<c::KoopmanModel> base;  // a krom: the SINDy krom is its child (its modes, one more representation)
    std::vector<std::size_t> modes;
    SampleStates anchors;
    if (given.kind == Given::Kind::data) {
        // Kaiser, Kutz & Brunton (2021), Algorithm 1 (SINDy.hpp): eigenvalues from the library generator.
        c::KoopmanSINDyEigenvalueOptions estimate;
        estimate.validationTolerance = settings.at("validation_tolerance").number(0.1);
        estimate.refinementIterations = std::size_t(settings.at("refinement_iterations").number(50));
        if (!training.oneTrajectory()) throw std::invalid_argument("SINDy eigenvalues from data need one trajectory");
        if (options.dictionary) throw std::invalid_argument("SINDy eigenvalues from data use the polynomial library; give eigenvalues for a dictionary fit");
        const auto found = c::estimateKoopmanEigenvaluesInCoordinates(training.states, training.samples, training.dimension,
            training.time, coordinates, velocities, options, estimate);
        if (found.eigenvalues.empty())
            throw std::invalid_argument("SINDy found no eigenvalue that passes validation on this data ("
                + std::to_string(found.candidates) + " candidates); raise validation_tolerance or the polynomial degree");
        given.rates = found.eigenvalues;
        json::Value::Array initial, errors;
        for (const auto value : found.initial) initial.emplace_back(json::Value::Array{value.real(), value.imag()});
        for (const auto value : found.validationError) errors.emplace_back(value);
        estimateRecord = {{"algorithm", "Kaiser, Kutz & Brunton 2021, Algorithm 1 (library generator, null vectors, Rayleigh update)"},
                          {"candidates", double(found.candidates)}, {"rejected", double(found.rejected)},
                          {"initial", initial}, {"validation_error", errors}, {"validation_tolerance", estimate.validationTolerance}};
    } else if (given.kind == Given::Kind::krom) {
        const auto all = c::loadKoopmanModel(given.path);
        modes = selectModes(settings, all.liftedDimension(), given.modes);
        base.emplace(all.extract(modes));
        given.rates = base->spectrum().continuousRates();
        anchors = trainingStates(*base);  // the dynamical states at the eigenfunctions' samples
        if (anchors.dimension != training.dimension)
            throw std::invalid_argument("the eigenfunctions' samples are not states of this training series");
        for (const auto& name : base->representationNames())
            if (name == "polynomial" && !settings.at("refit").boolean(false))
                throw std::invalid_argument("these eigenfunctions already have a polynomial representation; give \"refit\": true to fit another object");
    }
    const auto& continuous = given.rates;

    // Fit: the null space of the generator equation (default; Kaiser, Kutz & Brunton 2021), or, on
    // request, anchored to the eigenfunctions' own values. The anchored fit trades generator residual
    // for agreement with the references, so where the dictionary cannot represent them (H-DMD values
    // on a late window, say) its polynomial is no longer an eigenfunction; it is opt-in.
    const bool reference = base.has_value();
    const auto fitMode = settings.at("fit").string("nullspace");
    if (fitMode != "anchored" && fitMode != "nullspace") throw std::invalid_argument("SINDy fit must be anchored or nullspace");
    if (fitMode == "anchored" && !reference) throw std::invalid_argument("an anchored SINDy fit needs eigenfunctions to anchor to");
    const std::size_t anchorRows = reference ? base->samples().samples() : 0;
    if (fitMode == "anchored") {
        options.anchorStates = &anchors.values;
        options.anchorValues = &base->samples().values;
        options.anchorWeight = settings.at("anchor_weight").number(1.0);
    }
    auto result = c::fitKoopmanEigenfunctionsInCoordinates(training.states, training.samples, training.dimension, training.time,
                                                           continuous, coordinates, velocities, options);
    // Calibrate the scale and phase to the reference values (anchored fits are already on them; their
    // scale is still the least-squares one), or normalize when there are none or they disagree.
    const double tolerance = settings.at("agreement_tolerance").number(0.1);
    Calibration calibration;
    {
        const std::size_t stride = std::max<std::size_t>(1, training.samples / 2000);
        std::vector<double> rows;
        for (std::size_t k = 0; k < training.samples; k += stride)
            rows.insert(rows.end(), training.states.begin() + k * training.dimension, training.states.begin() + (k + 1) * training.dimension);
        const std::size_t count = rows.size() / training.dimension;
        const auto q = polynomialValues(result, coordinates, rows, count);
        if (reference) {
            const auto p = polynomialValues(result, coordinates, anchors.values, anchorRows);
            calibration = calibrate(result, p, &base->samples().values, anchorRows, q, count, tolerance, fitMode == "nullspace");
        } else {
            calibration = calibrate(result, {}, nullptr, 0, q, count, tolerance, true);
        }
    }
    json::Value::Array flagged, calibrationKinds;
    for (std::size_t j = 0; j < calibration.mismatch.size(); ++j) {
        if (reference && calibration.mismatch[j] > tolerance) flagged.emplace_back(double(j));
        calibrationKinds.emplace_back(calibration.referenced[j] ? "reference" : "normalized");
    }

    c::NpzRecord archive;
    c::DataOrigin origin{"sindy", settings, {dataReference(fullSource.path, "time_series")}};
    if (!given.path.empty()) origin.inputs.push_back(dataReference(given.path, reference ? "krom" : "spectrum"));
    std::shared_ptr<const c::FunctionRepresentation> polynomial = std::make_shared<const c::FunctionExpansion>(result.functions);
    const bool makeDefault = settings.at("make_default").boolean(true);
    std::optional<c::KoopmanModel> krom;
    if (reference) {  // the same eigenfunctions, one more representation; the base's by reference (its parent)
        krom.emplace(std::move(*base));
        const auto names = krom->representationNames();
        if (std::find(names.begin(), names.end(), c::representationName(*polynomial)) != names.end()) {
            krom->replaceRepresentation(polynomial);
            if (makeDefault) krom->useRepresentation(c::representationName(*polynomial));
        } else {
            krom->addRepresentation(polynomial, makeDefault);
        }
        krom->setOrigin(origin);
        // The pairs' residuals describe the evaluated eigenfunction: with the polynomial evaluated they
        // are this fit's generator residuals; with make_default false the parent's stay.
        if (makeDefault)
            for (std::size_t j = 0; j < krom->liftedDimension(); ++j) krom->setResidual(j, result.spectrum.residual(j));
    } else {
        std::vector<Complex> discrete;
        for (const auto value : continuous) discrete.push_back(std::exp(value * dt));
        c::EigenfunctionSamples samples;
        samples.modes = continuous.size();
        if (trainingValues && training.oneTrajectory()) {  // the calibrated eigenfunctions at the training rows
            samples.values = polynomialValues(result, coordinates, source.values(), source.samples());
            for (std::size_t k = 0; k < source.samples(); ++k) {  // rows of the input series, matched by time
                const auto found = std::lower_bound(fullSource.time().begin(), fullSource.time().end(), source.time()[k] - 1e-12);
                samples.sampleIndices.push_back(std::size_t(found - fullSource.time().begin()));
                samples.time.push_back(source.time()[k]);
            }
        }
        auto spectrum = c::Spectrum::fromDiscrete(discrete, continuous, dt, origin);
        spectrum.adoptEvidence(result.spectrum);  // the pairs' generator residuals
        krom.emplace(std::move(spectrum), samples, origin,
                     std::vector<std::shared_ptr<const c::FunctionRepresentation>>{polynomial});
        krom->bindSystem(systemConfiguration(source));
        archive.arrays["sample_interval"] = c::NpyArray::fromDoubles({1}, {dt});
    }
    const std::size_t m = continuous.size();
    archive.arrays["training_residual_maximum"] = c::NpyArray::fromDoubles({m}, result.diagnostics.residualMaximum);
    archive.arrays["calibration_scale"] = c::NpyArray::fromComplex({m}, calibration.scale);
    if (reference) archive.arrays["agreement_mismatch"] = c::NpyArray::fromDoubles({m}, calibration.mismatch);
    const bool archiveDerivatives = settings.at("archive_coordinate_derivatives").boolean(false);
    if (archiveDerivatives && (velocities.empty() || !coordinates.hasAnalyticJacobian()))
        throw std::invalid_argument("coordinate derivative evidence requires system_rhs and an analytic coordinate Jacobian");
    if (archiveDerivatives) {  // requested evidence: the coordinates' exact rates and Jacobians at the samples
        archive.arrays["training_coordinate_velocities"] = c::NpyArray::fromDoubles({training.samples, coordinates.outputDimension()},
            coordinates.pushForwardRealBatch(training.states, velocities, training.samples, c::ObservableJacobianBackend::analytic));
        std::vector<double> gradients;
        gradients.reserve(training.samples * coordinates.outputDimension() * training.dimension);
        for (std::size_t row = 0; row < training.samples; ++row)
            for (const auto value : coordinates.jacobian(training.states.data() + row * training.dimension,
                                                         c::ObservableJacobianBackend::analytic))
                gradients.push_back(value.real());
        archive.arrays["training_coordinate_jacobians"] = c::NpyArray::fromDoubles(
            {training.samples, coordinates.outputDimension(), training.dimension}, gradients);
    }

    c::KoopmanModel fittedModel(result, systemConfiguration(source));
    std::vector<std::string> reconstructionNames;
    for (const auto& value : settings.at("reconstruction_observables").array()) reconstructionNames.push_back(value.string());
    // "reconstruction":"none" skips the all-mode least-squares C fit, which dominates SINDy runtime for
    // a full H-DMD spectrum and is unused by state-Jacobian/coordinate-target control.
    const bool skipReconstruction = settings.at("reconstruction").string() == "none";
    if (skipReconstruction && !reconstructionNames.empty())
        throw std::invalid_argument("reconstruction \"none\" conflicts with reconstruction_observables");
    if (reconstructionNames.empty() && !skipReconstruction) reconstructionNames.push_back("state");
    // Readout transfer: given a krom, its readouts C apply to these same eigenfunctions once their scale
    // and phase are calibrated to its values (eigenfunctions are defined up to a factor; the anchored
    // calibration fixes it). The transfer is accepted when every mode agrees with its reference within
    // "readout_transfer" tolerance (default 0.2); then the inherited C are kept, else they are refitted.
    json::Value transferRecord;
    const auto& transferSetting = settings.at("readout_transfer");
    if (reference && (transferSetting.boolean(false) || std::holds_alternative<json::Value::Object>(transferSetting.storage()))) {
        const double tolerance = transferSetting.at("tolerance").number(0.2);
        double worst = 0;
        for (const double value : calibration.mismatch) worst = std::max(worst, value);
        const bool passed = std::isfinite(worst) && worst <= tolerance;
        json::Value::Array transferred, refitted;
        for (const auto& [name, readout] : krom->reconstructions()) {
            (void)readout;
            if (passed) transferred.emplace_back(name);
            else if (std::find(reconstructionNames.begin(), reconstructionNames.end(), name) == reconstructionNames.end()
                     && catalog.contains(name)) { reconstructionNames.push_back(name); refitted.emplace_back(name); }
        }
        if (passed)  // the transferred readouts are the parent's; fit only those it lacks
            reconstructionNames.erase(std::remove_if(reconstructionNames.begin(), reconstructionNames.end(),
                [&](const std::string& name) { return krom->hasReconstruction(name); }), reconstructionNames.end());
        transferRecord = json::Value::Object{{"method", "eigenfunction-scale-aligned C transfer"}, {"passed", passed},
            {"tolerance", tolerance}, {"maximum_relative_alignment_residual", worst},
            {"transferred", transferred}, {"refitted", refitted}};
    }
    for (const auto& name : reconstructionNames)  // C of each named observable, fitted on the polynomial
        krom->setReconstruction(name, c::fitObservableReconstruction(fittedModel, training.states, training.samples, catalog.at(name),
                                                                     settings.at("reconstruction_ridge").number(0.0)));
    const std::string gamma = velocities.empty() ? "finite_difference" : "system_rhs";
    std::vector<double> multipliers(coordinates.outputDimension(), 1.0);  // rescale after normalization
    for (std::size_t i = 0; i < multipliers.size() && !options.rescale.empty(); ++i)
        multipliers[i] = options.rescale[options.rescale.size() == 1 ? 0 : i];
    auto sindyMeta = json::parse(metadata(source,
        "koopman-with-sindy", label, "\"eigenvalue_domain\":\"discrete-time\","
        "\"generator_eigenvalue_domain\":\"continuous-time\",\"generator_conversion\":"
        "\"lambda=log(mu)/delta_t\",\"delta_t\":" + std::to_string(dt)
        + ",\"polynomial_degree\":" + std::to_string(options.polynomialDegree)
        + ",\"feature_count\":" + std::to_string(result.functions.dictionary().size())
        + ",\"threshold\":" + std::to_string(options.threshold)
        + ",\"sample_stride\":" + std::to_string(static_cast<long long>(settings.at("sample_stride").number(1)))
        + ",\"training_samples\":" + std::to_string(training.samples)
        + ",\"source_samples\":" + std::to_string(fullSource.samples())
        + ",\"gamma_method\":\"" + gamma + "\",\"state_preprocessing\":{\"center\":"
        + jsonArray(result.functions.monomial() ? result.functions.monomialTable().center : std::vector<double>{})
        + ",\"scale\":" + jsonArray(result.functions.monomial() ? result.functions.monomialTable().scale : std::vector<double>{})
        + ",\"multiplier\":" + jsonArray(multipliers)
        + "},\"eigenfunction_representation\":\"analytic-polynomial\","
        "\"coordinate_observables\":" + jsonStringArray(result.functions.coordinates().names) + ","
        "\"reconstruction_observables\":" + jsonStringArray(reconstructionNames) + ","
        "\"source_spectrum_archive\":\"" + escaped(given.path.generic_u8string()) + "\","
        "\"include_input_terms\":" + (includeInputTerms ? "true" : "false") + ","
        "\"coordinate_derivative_evidence\":\"" + (archiveDerivatives ? "analytic-push-forward/v1" : "none") + "\","
        "\"coordinate_gradient_backend\":\"" + (coordinates.hasAnalyticJacobian() ? "analytic" : "finite-difference") + "\","
        "\"analytic_gradient_available\":" + (coordinates.hasAnalyticJacobian() ? "true" : "false") + ",\"direct_eigenvalue_fits\":"
        + std::to_string(result.diagnostics.directFitCount)));
    sindyMeta.object()["sindy_input"] = json::Value::Object{{"kind", kindName(given.kind)}, {"fit", fitMode}};
    sindyMeta.object()["sindy_training"] = training.describe();
    if (!transferRecord.isNull()) sindyMeta.object()["readout_transfer"] = transferRecord;
    sindyMeta.object()["sindy_options"] = json::Value::Object{{"sparsification_policy", options.sparsification},
        {"directional_derivative", options.directionalDerivative}, {"rescale", json::parse(jsonArray(options.rescale))},
        {"dictionary", options.dictionary ? options.dictionary->describe() : json::Value("polynomial")}};
    if (!estimateRecord.empty()) sindyMeta.object()["eigenvalue_estimate"] = json::Value(estimateRecord);
    if (reference)
        sindyMeta.object()["agreement"] = json::Value::Object{{"reference", "eigenfunction_values (the eigenfunctions' sample values)"},
            {"tolerance", tolerance}, {"flagged_modes", flagged}, {"calibration", calibrationKinds},
            {"measure", "||c p_j - psi_j|| / ||psi_j|| at the eigenfunctions' samples"}};
    json::Value::Array perMode{"training_residual_maximum", "calibration_scale"};  // the residual RMS is the pairs' own
    if (reference) perMode.emplace_back("agreement_mismatch");
    sindyMeta.object()["object_properties"] = json::Value::Object{{"per_eigenfunction", perMode}};
    const auto published = publishKrom(root, source, *krom, "sindy-eigenfunctions", json::dump(sindyMeta),
                                       std::move(archive.arrays), fullSource.record.channels);
    if (!flagged.empty())
        std::fprintf(stderr, "SINDy: %zu eigenfunction(s) disagree with their reference values beyond %.3g "
                     "(agreement.flagged_modes)\n", flagged.size(), tolerance);
    return published;
}

}  // namespace dynsyslab::workflow::postprocess_detail

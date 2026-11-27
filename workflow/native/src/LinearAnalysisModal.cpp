// linear_analysis methods on a K-ROM source:
//   "modal-gramian": the closed-form observability Gramian of the diagonal modal dynamics over
//     `horizon` (core analyzeModalObservability), its eigenvalues, directions and singular values
//     (finite-horizon observability Gramian: Kalman, Ho & Narendra 1963, "Controllability of
//     linear dynamical systems", Contrib. Differ. Equ. 1:189-213).
//   "klqe-assessment": a finite-window linear-Gaussian benchmark of a mode subset for echo-only
//     estimation (core assessEstimationSubset): state and measurement reconstruction MSE,
//     posterior state-noise MSE and autonomous-propagation MSE on a trajectory, with the score.
//     Keys: trajectory (an archive with time and values), measurement, output, output_scales,
//     measurement_standard_deviation, horizon, sample_time, dynamics_weight,
//     measurement_mismatch_weight, selected_modes.
#include "LinearAnalysisModal.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/EstimationSelection.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/core/System.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

std::vector<std::size_t> indices(const json::Value& values) {
    std::vector<std::size_t> result;
    for (const auto& v : values.array()) result.push_back(static_cast<std::size_t>(v.number()));
    return result;
}

dynsyslab::core::NpyArray complexArray(const std::vector<std::complex<double>>& v, std::size_t a, std::size_t b) {
    std::vector<double> flat;
    for (const auto& x : v) { flat.push_back(x.real()); flat.push_back(x.imag()); }
    return dynsyslab::core::NpyArray::fromDoubles({a, b, 2}, flat);
}

std::filesystem::path finish(const std::filesystem::path& root, const SourceData& source, const std::string& label,
                             dynsyslab::core::NpzRecord out, json::Value::Object details) {
    auto serialized = json::dump(json::Value(std::move(details)));
    out.arrays["metadata_json"] = dynsyslab::core::NpyArray::fromString(
        metadata(source, "linear-analysis", label, serialized.substr(1, serialized.size() - 2)));
    return publishObject(root, source, dynsyslab::core::object_kind::report, "linear-analysis", std::move(out));
}

}  // namespace

std::filesystem::path runModalGramian(const std::filesystem::path& sourcePath, const SourceData& source,
                                      const json::Value& settings, const std::filesystem::path& root,
                                      const std::string& label) {
    using namespace dynsyslab::core;
    const auto model = loadKoopmanModel(sourcePath, indices(settings.at("selected_modes")));
    const auto measurement = settings.at("measurement").string();
    if (measurement.empty()) throw std::invalid_argument("modal-gramian needs a measurement");
    ObservabilityOptions options;
    options.horizon = settings.at("horizon").number(1.0);
    options.relativeTolerance = settings.at("relative_tolerance").number(options.relativeTolerance);
    std::vector<double> scales;
    for (const auto& v : settings.at("coordinate_scales").array()) scales.push_back(v.number());
    const auto r = analyzeModalObservability(model, measurement, options, scales);
    NpzRecord out;
    out.arrays["modal_gramian"] = complexArray(r.gramian, r.modes, r.modes);
    out.arrays["modal_gramian_eigenvalues"] = NpyArray::fromDoubles({r.gramianEigenvalues.size()}, r.gramianEigenvalues);
    out.arrays["modal_gramian_directions"] = complexArray(r.directions, r.modes, r.directions.size() / std::max<std::size_t>(1, r.modes));
    out.arrays["modal_singular_values"] = NpyArray::fromDoubles({r.singularValues.size()}, r.singularValues);
    return finish(root, source, label, std::move(out),
        {{"schema", "dynsyslab.linear-analysis/v1"}, {"analysis_method", "modal-gramian"},
         {"source_krom", sourcePath.generic_u8string()}, {"measurement", measurement},
         {"horizon", options.horizon}, {"rank", double(r.rank)}, {"rank_tolerance", r.rankTolerance}});
}

std::filesystem::path runKlqeAssessment(const std::filesystem::path& sourcePath, const SourceData& source,
                                        const json::Value& settings, const std::filesystem::path& root,
                                        const std::string& label) {
    using namespace dynsyslab::core;
    const auto model = loadKoopmanModel(sourcePath, indices(settings.at("selected_modes")));
    const auto trajectory = settings.at("trajectory").string();
    if (trajectory.empty()) throw std::invalid_argument("klqe-assessment needs a trajectory archive");
    const auto run = readNpz(std::filesystem::u8path(trajectory));
    const auto time = run.at("time").doubles();
    const auto values = run.at("values").doubles();
    const std::size_t n = model.stateDimension(), samples = time.size();
    if (values.size() != samples * n) throw std::invalid_argument("klqe-assessment trajectory does not match the K-ROM state");
    EstimationSelectionOptions options;
    options.measurement = settings.at("measurement").string(options.measurement);
    options.output = settings.at("output").string(options.output);
    for (const auto& v : settings.at("output_scales").array()) options.outputScales.push_back(v.number());
    options.measurementStandardDeviation = settings.at("measurement_standard_deviation").number(options.measurementStandardDeviation);
    options.horizon = settings.at("horizon").number(options.horizon);
    options.sampleTime = settings.at("sample_time").number(samples > 1 ? time[1] - time[0] : options.sampleTime);
    options.dynamicsWeight = settings.at("dynamics_weight").number(options.dynamicsWeight);
    options.measurementMismatchWeight = settings.at("measurement_mismatch_weight").number(options.measurementMismatchWeight);
    auto system = createSystem(model.systemConfiguration());
    const auto catalog = system->observableCatalog();
    const auto& measured = catalog.at(settings.at("measurement_observable").string(options.measurement));
    const auto& output = catalog.at(settings.at("output_observable").string(options.output));
    std::vector<std::complex<double>> psi;
    std::vector<double> outputs, measurements;
    for (std::size_t k = 0; k < samples; ++k) {
        const double* x = values.data() + k * n;
        const auto p = model.evaluate(x);
        psi.insert(psi.end(), p.begin(), p.end());
        const auto y = measured.evaluateReal(x), o = output.evaluateReal(x);
        measurements.insert(measurements.end(), y.begin(), y.end());
        outputs.insert(outputs.end(), o.begin(), o.end());
    }
    const auto a = assessEstimationSubset(model, time, psi, outputs, measurements, options);
    NpzRecord out;
    out.arrays["mode_scales"] = NpyArray::fromDoubles({a.modeScales.size()}, a.modeScales);
    out.arrays["observability_singular_values"] = NpyArray::fromDoubles({a.observability.singularValues.size()},
                                                                        a.observability.singularValues);
    json::Value::Object details{{"schema", "dynsyslab.linear-analysis/v1"}, {"analysis_method", "klqe-assessment"},
         {"source_krom", sourcePath.generic_u8string()}, {"trajectory", trajectory},
         {"state_reconstruction_mse", a.stateReconstructionMse},
         {"measurement_reconstruction_mse", a.measurementReconstructionMse},
         {"posterior_state_noise_mse", a.posteriorStateNoiseMse},
         {"autonomous_propagation_mse", a.autonomousPropagationMse}, {"score", a.score}};
    // "role": "validation" when the trajectory is held out: the scores are then also under
    // "validation", where experiment select stages read held-out metrics (EXPERIMENTS.md).
    if (settings.at("role").string() == "validation")
        details["validation"] = json::Value::Object{{"state_reconstruction_mse", a.stateReconstructionMse},
            {"measurement_reconstruction_mse", a.measurementReconstructionMse}, {"posterior_state_noise_mse", a.posteriorStateNoiseMse},
            {"autonomous_propagation_mse", a.autonomousPropagationMse}, {"score", a.score}};
    return finish(root, source, label, std::move(out), std::move(details));
}

}  // namespace dynsyslab::workflow::postprocess_detail

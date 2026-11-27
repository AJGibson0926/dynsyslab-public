// input_map fit "symmetric_pulse": B_psi samples from paired +-u first responses
// (core/SymmetricPulseControlMap.hpp). Each source is one trajectory solved with the
// symmetric_pulse prescribed input, whose metadata.input_design carries pair_id and polarity;
// its first two samples give the anchor and the response after one held sample. The samples
// fit representation "constant" (their mean) or "bilinear" (least squares in psi(anchor)).
#include "PostprocessInternal.hpp"
#include "InputMapPulses.hpp"

#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/SymmetricPulseControlMap.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {

json::Value::Object symmetricPulseRegression(dynsyslab::core::KoopmanModel& model,
                                             const std::vector<std::filesystem::path>& sources,
                                             const json::Value& settings, const std::string& representation,
                                             const std::string& mapName) {
    std::vector<dynsyslab::core::SymmetricPulseRun> runs;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        dynsyslab::core::requireCompatibleSystems(model.systemConfiguration(), systemConfiguration(source));
        const auto meta = source.metadataJson.empty() ? json::Value{} : json::parse(source.metadataJson);
        const auto& design = meta.at("input_design");
        // pulse_branches: unpaired branches (the zero branch, unmatched amplitudes) are skipped and
        // a pair's pulse must last exactly the first sample interval.
        const bool branches = design.at("kind").string() == "pulse_branches";
        if (branches && design.at("pair_id").isNull()) continue;
        if (branches && meta.at("control_evaluation").string() != "sampled")
            throw std::invalid_argument("pulse_branches pairs need solver.control_evaluation sampled (held pulses): " + path.generic_u8string());
        if (branches && source.samples() >= 2
            && std::abs(design.at("width").number() - (source.time()[1] - source.time()[0])) > 1e-9 * design.at("width").number())
            throw std::invalid_argument("symmetric_pulse fit of pulse_branches needs pulses one sample long: " + path.generic_u8string());
        if (design.at("pair_id").isNull() || design.at("polarity").isNull())
            throw std::invalid_argument("symmetric_pulse sources need input_design pair_id and polarity: "
                                        + path.generic_u8string());
        if (source.samples() < 2 || source.inputs().empty())
            throw std::invalid_argument("symmetric_pulse sources need two samples and their input: " + path.generic_u8string());
        const std::size_t n = source.dimension();
        dynsyslab::core::SymmetricPulseRun run;
        run.pairId = static_cast<std::size_t>(design.at("pair_id").number());
        run.polarity = static_cast<int>(design.at("polarity").number());
        run.anchor.assign(source.values().begin(), source.values().begin() + n);
        run.response.assign(source.values().begin() + n, source.values().begin() + 2 * n);
        run.input = source.inputs().front();
        run.sampleTime = source.time()[1] - source.time()[0];
        runs.push_back(std::move(run));
    }
    const auto samples = dynsyslab::core::estimateSymmetricPulseControlMap(model, runs);
    const auto fitted = representation == "constant"
        ? dynsyslab::core::fitConstantControlMapSamples(model, samples)
        : dynsyslab::core::fitBilinearControlMapSamples(model, samples, settings.at("ridge").number(0.0));
    model.setControlMap(mapName, fitted.map);
    const auto norms = samples.responseNorms;
    return {{"method", json::Value("symmetric-opposite-pulse-control-map-estimator/v1")},
            {"pair_count", json::Value(double(samples.pairs))}, {"analytic_Bpsi_used", json::Value(false)},
            {"relative_residual", json::Value(fitted.relativeResidual)},
            {"complex_rmse", json::Value(fitted.complexRmse)},
            {"lifted_response_norm_minimum", json::Value(*std::min_element(norms.begin(), norms.end()))},
            {"lifted_response_norm_maximum", json::Value(*std::max_element(norms.begin(), norms.end()))}};
}

}  // namespace dynsyslab::workflow::postprocess_detail

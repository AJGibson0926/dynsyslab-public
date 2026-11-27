// krom_control_search extensions (job keys under "screening"):
//   "schedule": [[mode indices], ...], "schedule_mode": "prepend" | "only"
//       an explicit candidate schedule (historical anchors, repeats) before or instead of the
//       ranked order;
//   "closure_weight": w   reorder the ranked candidates by selection_score * (1 + w q), with q the
//       structural bilinear-closure quality 1 - (relative residual of the best b0 + sum b_i psi_i
//       fit to grad(psi) . g on the controlled states): how far the exact input map stays in
//       span(psi, 1) (Goswami & Paley 2017, "Global bilinearization and controllability of
//       control-affine nonlinear systems: a Koopman spectral approach", IEEE CDC 6107-6112);
//   "map_variants": [{"name", "family": "bilinear"|"constant", "conjugate_complete",
//       "ridge", "sampling": {"path": "controlled"|"prior_exact", "arc_points": N,
//       "tube_radius": r, "tube_count": k, "seed": s}}, ...]
//       extra per-attempt control cases, each an exact-map fit on the sampled states
//       (default [{"name": "bilinear_map"}]: bilinear on the controlled path);
// and top-level "resume_from": a campaign-ledger.json whose attempts are kept and not rerun.
#include "KromControlSearchSupport.hpp"
#include "KromControlSearchVariants.hpp"
#include "ObjectPublishing.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/PathSampling.hpp>
#include <dynsyslab/core/Text.hpp>

#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow::search_detail {
namespace {

std::vector<std::size_t> modes(const json::Value& value) {
    std::vector<std::size_t> result;
    for (const auto& item : value.array()) result.push_back(static_cast<std::size_t>(item.number()));
    std::sort(result.begin(), result.end());
    return result;
}

double closureQuality(const core::KoopmanModel& parent, const std::vector<std::size_t>& indices,
                      const std::vector<double>& states, std::size_t samples) {
    try {
        const auto model = parent.extract(indices);
        const auto fit = core::fitBilinearExactControlMap(model, states, samples, 0.0, false);
        return std::max(0.0, 1.0 - fit.diagnostics.relativeResidual);
    } catch (const std::exception&) { return 0.0; }
}

}  // namespace

std::vector<double> controlledStates(const std::filesystem::path& source, std::size_t& samples) {
    const auto archive = core::readNpz(source);
    const auto& values = archive.at("values");
    samples = values.shape.empty() ? 0 : values.shape[0];
    return values.doubles();
}

std::vector<core::KromCandidateScore> scheduleCandidates(const json::Value& job,
    std::vector<core::KromCandidateScore> ranked, const core::KoopmanModel& parent,
    const std::vector<double>& states, std::size_t samples, json::Value& evidence) {
    const auto& screening = job.at("screening");
    const double weight = screening.at("closure_weight").number(0.0);
    json::Value::Array closure;
    if (weight > 0.0) {
        std::vector<std::pair<double, std::size_t>> order;
        for (std::size_t k = 0; k < ranked.size(); ++k) {
            const double q = closureQuality(parent, ranked[k].modeIndices, states, samples);
            closure.emplace_back(q);
            order.emplace_back(-ranked[k].selectionScore * (1.0 + weight * q), k);
        }
        std::stable_sort(order.begin(), order.end());
        std::vector<core::KromCandidateScore> reordered;
        for (const auto& [score, k] : order) reordered.push_back(ranked[k]);
        ranked = std::move(reordered);
    }
    std::vector<core::KromCandidateScore> result;
    std::set<std::vector<std::size_t>> seen;
    for (const auto& entry : screening.at("schedule").array()) {
        core::KromCandidateScore candidate;
        candidate.modeIndices = modes(entry);
        if (candidate.modeIndices.empty()) throw std::invalid_argument("screening.schedule entries need mode indices");
        for (auto j : candidate.modeIndices)
            if (j >= parent.liftedDimension()) throw std::invalid_argument("screening.schedule names a mode the parent lacks");
        result.push_back(candidate);  // a scheduled repeat is kept on purpose
        seen.insert(candidate.modeIndices);
    }
    const auto mode = screening.at("schedule_mode").string("prepend");
    if (mode != "prepend" && mode != "only") throw std::invalid_argument("screening.schedule_mode must be prepend or only");
    if (mode == "prepend")
        for (const auto& candidate : ranked) {
            auto key = candidate.modeIndices; std::sort(key.begin(), key.end());
            if (!seen.count(key)) result.push_back(candidate);
        }
    evidence = json::Value::Object{{"closure_weight", weight}, {"closure_quality", json::Value(std::move(closure))},
                                   {"scheduled", double(screening.at("schedule").array().size())},
                                   {"schedule_mode", mode}};
    return result;
}

json::Value::Array resumedAttempts(const json::Value& job) {
    const auto path = job.at("resume_from").string();
    if (path.empty()) return {};
    std::ifstream stream(std::filesystem::u8path(path), std::ios::binary);
    if (!stream) throw std::invalid_argument("resume_from ledger does not exist: " + path);
    std::ostringstream text; text << stream.rdbuf();
    const auto ledger = json::parse(text.str());
    if (ledger.at("schema").string() != "dynsyslab.krom-control-search/v1")
        throw std::invalid_argument("resume_from is not a krom-control-search ledger");
    return ledger.at("attempts").array();
}

bool attempted(const json::Value::Array& attempts, const std::vector<std::size_t>& indices) {
    auto key = indices; std::sort(key.begin(), key.end());
    for (const auto& attempt : attempts)
        if (modes(attempt.at("candidate").at("mode_indices")) == key) return true;
    return false;
}

std::vector<MapVariant> fitMapVariants(const json::Value& job, const std::filesystem::path& extracted,
    const std::vector<double>& controlled, std::size_t controlledSamples, const std::string& priorExact,
    const std::filesystem::path& directory, std::size_t attempt) {
    json::Value::Array specs = job.at("screening").at("map_variants").array();
    if (job.at("screening").at("map_variants").isNull())
        specs.emplace_back(json::Value::Object{{"name", "bilinear_map"}});
    std::vector<MapVariant> result;
    for (const auto& spec : specs) {
        MapVariant variant;
        variant.name = spec.at("name").string();
        if (variant.name.empty() || variant.name == "exact" || variant.name == "constant" || variant.name == "bilinear")
            throw std::invalid_argument("map_variants need a distinct name");
        try {
            auto model = core::loadKoopmanModel(extracted);
            const std::size_t n = model.stateDimension();
            const auto& sampling = spec.at("sampling");
            std::vector<double> path = controlled;
            std::size_t samples = controlledSamples;
            if (sampling.at("path").string("controlled") == "prior_exact") {
                if (priorExact.empty()) throw std::invalid_argument("no successful prior exact path yet");
                path = controlledStates(std::filesystem::u8path(priorExact), samples);
            }
            if (!sampling.at("arc_points").isNull()) {
                const auto points = static_cast<std::size_t>(sampling.at("arc_points").number());
                path = core::resampleByArcLength(path, samples, n, points);
                samples = points;
            }
            if (sampling.at("tube_count").number(0) > 0) {
                const auto per = static_cast<std::size_t>(sampling.at("tube_count").number());
                path = core::tubeStates(path, samples, n, sampling.at("tube_radius").number(0.05), per,
                                        static_cast<std::uint64_t>(sampling.at("seed").number(0)));
                samples = path.size() / n;
            }
            const auto family = spec.at("family").string("bilinear");
            if (family == "constant") {
                auto fit = core::fitConstantExactControlMap(model, path, samples);
                variant.relativeResidual = fit.diagnostics.relativeResidual;
                model.setControlMap(variant.name, fit.representation);
            } else if (family == "bilinear") {
                auto fit = core::fitBilinearExactControlMap(model, path, samples, spec.at("ridge").number(0.0),
                                                            spec.at("conjugate_complete").boolean(false));
                variant.relativeResidual = fit.diagnostics.relativeResidual;
                model.setControlMap(variant.name, fit.representation);
            } else throw std::invalid_argument("map variant family must be bilinear or constant");
            variant.archive = directory / core::wide("attempt-" + std::to_string(attempt) + "__" + variant.name + ".npz");
            const json::Value meta = json::Value::Object{{"schema", "dynsyslab.derived/v1"}, {"method", "krom-control-map-variant"},
                {"variant", spec}, {"source_krom", extracted.generic_u8string()},
                {"sampled_states", double(samples)}, {"relative_residual", variant.relativeResidual}};
            core::saveKoopmanModel(variant.archive, model, postprocess_detail::withDataObject(
                json::dump(meta), core::object_kind::krom, {extracted}));
        } catch (const std::exception& error) { variant.failure = error.what(); }
        result.push_back(std::move(variant));
    }
    return result;
}

}  // namespace dynsyslab::workflow::search_detail

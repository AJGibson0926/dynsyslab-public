// Offline control-target lifting ("control_target" post-processing method).
//
// A KLQR/KLQE controller tracks a target in its own coordinates. When the target is
// specified as a physical state waveform (e.g. a radius history r(t) with r'(t)), mapping
// it into the controller's coordinates may require a physical model: the far-field echo
// coordinates (e, e') and the causal echo-history coordinates xi are computed from the
// state target through the system's observable expressions. This method performs that
// mapping OFFLINE, once, on the control clock, and writes the coordinate target as an
// archive that the controller reads as an ordinary "trajectory" target. Online control then
// evaluates no physical model to obtain its reference.
//
// Settings:
//   krom          KROM (or input-map) archive whose coordinates the target is lifted into
//                 (defaults to the first source).
//   state_target  state target specification (waveform, constant, trajectory; optional ramp).
//   time          {start, stop, step}; step defaults to the history clock (history KROMs)
//                 or 0.001.
//   measurement_observable  history KROMs: observable sampled along the target (default:
//                 the frozen transform's measurement).
// Output arrays: time, values (coordinate target), state_target_values.
#include "PostprocessInternal.hpp"
#include "HistoryControl.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;

std::filesystem::path runControlTarget(const std::vector<std::filesystem::path>& sources, const json::Value& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    const auto requested = settings.at("krom").string();
    if (requested.empty() && sources.empty()) throw std::invalid_argument("control_target needs a KROM");
    const auto kromPath = requested.empty() ? sources.front() : std::filesystem::u8path(requested);
    const auto model = c::loadKoopmanModel(kromPath);
    const bool history = model.coordinates().kind == c::EigenfunctionCoordinateKind::causalHistory;
    const auto system = c::createSystem(model.systemConfiguration());
    const auto stateDimension = system->initialState().size();
    if (settings.at("state_target").isNull()) throw std::invalid_argument("control_target needs a state_target");
    const auto stateTarget = resolveControlTarget(settings.at("state_target"), stateDimension);

    std::optional<c::FrozenHistoryTransform> transform;
    std::string transformPath;
    c::StateTargetFunction coordinateTarget;
    if (history) {
        transformPath = settings.at("history_transform").string(model.coordinates().providerArtifact.empty()
            ? kromPath.u8string() : model.coordinates().providerArtifact);
        transform = c::loadFrozenHistoryTransform(std::filesystem::u8path(transformPath));
        if (transform->identity() != model.coordinates().providerIdentity)
            throw std::invalid_argument("History transform differs from KROM training coordinates");
        const auto name = settings.at("measurement_observable").string(transform->measurement());
        coordinateTarget = detail::historyCoordinateTarget(stateTarget, *transform,
                                                           system->observableCatalog().at(name));
    } else {
        // Measured-coordinate KROMs carry no state observable; the target names the one that
        // defines their coordinates (e.g. echo_jet2 = (e, e')), evaluated offline on the target.
        const auto named = settings.at("coordinate_observable").string();
        const auto catalog = system->observableCatalog();
        const auto coordinates = !named.empty() ? (catalog.coordinatePresets().count(named) ? catalog.coordinatePreset(named) : catalog.at(named))
                                                : model.coordinateObservable();
        coordinateTarget = [stateTarget, coordinates](double time) {
            const auto x = stateTarget(time); return coordinates.evaluateReal(x.data());
        };
    }

    const auto& time = settings.at("time");
    const double start = time.at("start").number(0.0), stop = time.at("stop").number(NAN);
    const double step = time.at("step").number(history ? transform->sampleInterval() : 0.001);
    if (!std::isfinite(start) || !std::isfinite(stop) || !(stop > start) || !(step > 0))
        throw std::invalid_argument("control_target time needs finite start < stop and a positive step");
    if (history && std::abs(step - transform->sampleInterval()) > 1e-10 * step)
        throw std::invalid_argument("History targets must be sampled on the history transform clock");
    const auto count = static_cast<std::size_t>(std::floor((stop - start) / step + 1e-9)) + 1;
    std::vector<double> t, values, states;
    const auto width = model.coordinateDimension();
    for (std::size_t k = 0; k < count; ++k) {
        const double tk = start + static_cast<double>(k) * step; t.push_back(tk);
        const auto q = coordinateTarget(tk);
        if (q.size() != width) throw std::invalid_argument("Lifted target width differs from the KROM coordinates");
        for (double v : q) { if (!std::isfinite(v)) throw std::invalid_argument("Lifted target is not finite"); values.push_back(v); }
        const auto x = stateTarget(tk); states.insert(states.end(), x.begin(), x.end());
    }

    json::Value::Array names; for (const auto& n : model.coordinates().names) names.emplace_back(n);
    json::Value::Object meta{
        {"format", json::Value("dynsyslab.derived/v1")}, {"method", json::Value("control-target")},
        {"label", json::Value(label)}, {"source_krom", json::Value(kromPath.generic_u8string())},
        {"target_space", json::Value("observable")}, {"coordinate_names", json::Value(std::move(names))},
        {"coordinate_kind", json::Value(history ? "causal-history" : "measured")},
        {"state_target", settings.at("state_target")},
        {"computed", json::Value("offline")},
        {"provenance", json::Value(history
            ? "xi_ref: frozen history projection of the measurement observable evaluated along the state target"
            : "coordinate observables of the KROM evaluated along the state target")},
        {"data_role", json::Value("control_target")}};
    if (!transformPath.empty()) meta["history_transform"] = json::Value(transformPath);
    c::NpzRecord record;
    record.arrays["time"] = c::NpyArray::fromDoubles({count}, t);
    record.arrays["values"] = c::NpyArray::fromDoubles({count, width}, values);
    record.arrays["state_target_values"] = c::NpyArray::fromDoubles({count, stateDimension}, states);
    record.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(json::Value(std::move(meta))));
    const auto directory = root / L"derived" / L"control-target";
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__control-target__" + timestamp()) + ".npz");
    c::writeNpzAtomic(path, record);
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail

#include "HeldInputRecords.hpp"
#include "KromDatasetSupport.hpp"

#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using Complex = std::complex<double>;

enum class Hold { always, ifUnchanged };

Hold holdKind(const SourceData& source, double dt, bool strict) {
    if (source.inputs.size() == (source.samples - 1) * source.inputDimension) return Hold::always;
    if (!strict) return Hold::always;  // legacy history regression: every pair uses u_k
    if (source.metadataJson.empty()) return Hold::ifUnchanged;
    const auto meta = json::parse(source.metadataJson);
    if (meta.at("control_evaluation").string() != "sampled") return Hold::ifUnchanged;
    const double clock = meta.at("control_sample_time").number(0.0);
    const double ratio = clock / dt, whole = std::round(ratio);
    if (!(clock > 0.0) || whole < 1.0 || std::abs(ratio - whole) > 1e-6)
        throw std::invalid_argument("held-input pairs need a control clock that is a whole multiple of the "
                                    "record step: " + source.path.generic_u8string());
    return Hold::always;
}

bool insideWindow(const json::Value& window, const SourceData& s, std::size_t k) {
    if (window.isNull()) return true;
    if (window.array().size() != s.dimension)
        throw std::invalid_argument("state_window needs one [low, high] pair per state component");
    for (std::size_t j = 0; j < s.dimension; ++j) {
        const auto& bounds = window.array().at(j).array();
        const double value = s.values[k * s.dimension + j];
        if (value < bounds.at(0).number() || value > bounds.at(1).number()) return false;
    }
    return true;
}
}  // namespace

HeldInputRecords heldInputRecords(const c::KoopmanModel& model, const c::FrozenHistoryTransform* transform,
    const std::vector<std::filesystem::path>& sources, const json::Value& settings, bool strictHold) {
    const bool history = model.coordinates().computableFromMeasurements();
    if (history && !transform) throw std::invalid_argument("causal-history held-input pairs need the frozen transform");
    const bool external = !history && model.requiresExternalCoordinates();
    const double fitStart = settings.at("fit_start").number(0.0);
    const auto& window = settings.at("state_window");
    if (history && !window.isNull()) throw std::invalid_argument("state_window does not apply to causal-history KROMs");
    HeldInputRecords out;
    if (history) out.sampleTime = transform->sampleInterval();
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        if (source.inputs.empty() || source.inputDimension != 1 || datasetDynamics(source) != "forced")
            throw std::invalid_argument("input_map regression needs forced records with a scalar control history: "
                                        + path.generic_u8string());
        const double dt = uniformDeltaTime(source);
        if (out.sampleTime == 0.0) out.sampleTime = dt;
        if (std::abs(dt - out.sampleTime) > 1e-8 * out.sampleTime)
            throw std::invalid_argument(history ? "forced record must be sampled on the history transform clock"
                                                : "forced records must share one sample time");
        if (external && source.dimension != model.coordinateDimension())
            throw std::invalid_argument("forced record coordinates differ from the KROM coordinates");
        if (!history && !external) c::requireCompatibleSystems(model.systemConfiguration(), systemConfiguration(source));
        const auto hold = holdKind(source, dt, strictHold);

        std::optional<c::Observable> measurement;
        std::optional<c::CausalHistoryProjector> projector;
        if (history) {
            const auto system = c::createSystem(systemConfiguration(source));
            measurement = system->observableCatalog().at(settings.at("measurement_observable").string(transform->measurement()));
            projector.emplace(*transform);
        }
        const auto lift = [&](std::size_t k) -> std::optional<std::vector<Complex>> {
            const double* row = source.values.data() + k * source.dimension;
            if (std::any_of(row, row + source.dimension, [](double v) { return !std::isfinite(v); })) return std::nullopt;
            if (history) {
                const auto q = projector->push(source.time[k], measurement->evaluateReal(row).front());
                if (!q) return std::nullopt;
                return model.evaluateCoordinates(q->data());
            }
            return external ? model.evaluateCoordinates(row) : model.evaluate(row);
        };
        const std::size_t inputRows = source.inputs.size();
        std::optional<std::vector<Complex>> previous;
        std::size_t previousIndex = 0;
        for (std::size_t k = 0; k < source.samples; ++k) {
            auto psi = lift(k);
            if (!psi) { previous.reset(); continue; }
            if (previous && previousIndex + 1 == k && source.time[previousIndex] >= fitStart
                && insideWindow(window, source, previousIndex) && insideWindow(window, source, k)) {
                const double u = previousIndex < inputRows ? source.inputs[previousIndex] : 0.0;
                const bool changed = hold == Hold::ifUnchanged && k < inputRows && source.inputs[k] != u;
                if (changed) ++out.droppedAtInputChange;
                else if (std::isfinite(u)) out.rows.push_back({*previous, *psi, u});
            }
            previous = std::move(psi);
            previousIndex = k;
        }
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail

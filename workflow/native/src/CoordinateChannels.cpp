#include "CoordinateChannels.hpp"
#include "KromDatasetSupport.hpp"
#include <dynsyslab/core/ObservableDerivative.hpp>
#include <dynsyslab/core/KromSpectralObservables.hpp>
#include <dynsyslab/core/System.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;

std::size_t count(const J& v, const char* key, std::size_t fallback) {
    if (v.isNull()) return fallback;
    const double x = v.number(-1.);
    if (!std::isfinite(x) || x < 0 || x > 9007199254740991. || std::floor(x) != x)
        throw std::invalid_argument(std::string("invalid coordinate channel integer: ") + key);
    return static_cast<std::size_t>(x);
}

class Resolver {
public:
    explicit Resolver(const SourceData& s) : source(s) {}
    CoordinateChannel resolve(const J& selector, std::size_t order) {
        if (std::holds_alternative<double>(selector.storage())) {
            const auto k = count(selector, "observable", 0);
            if (k >= source.dimension() || order)
                throw std::invalid_argument("use a named observable for its derivative; check component index");
            CoordinateChannel result{"state_" + std::to_string(k), "state_" + std::to_string(k), "none", {}, {}};
            result.values.reserve(source.samples());
            for (std::size_t i = 0; i < source.samples(); ++i) result.values.push_back(source.values()[i * source.dimension() + k]);
            return result;
        }
        if (std::holds_alternative<J::Object>(selector.storage())) return krom(selector, order);
        if (!system) {
            system = c::createSystem(systemConfiguration(source));
            if (system->descriptor().stateDimension != source.dimension())
                throw std::invalid_argument("Coordinate observable and source state dimensions differ");
            catalog = system->observableCatalog();
        }
        const auto& base = catalog->at(selector.string());
        if (base.outputDimension() != 1 || base.descriptor().complexValued)
            throw std::invalid_argument("Coordinate channels require real scalar observables");
        const bool missingInput = source.inputs().empty() && datasetDynamics(source) == "forced"
            && system->configuration().driver.amplitude == 0.;
        std::optional<double> constantInput;
        if (!source.inputs().empty()) {
            if (source.inputDimension() == 1 && std::all_of(source.inputs().begin(), source.inputs().end(),
                [&](double u) { return u == source.inputs().front(); })) constantInput = source.inputs().front();
        } else if (!missingInput && system->configuration().driver.amplitude == 0.) constantInput = 0.;
        auto plan = c::resolveObservableDerivative(*system, *catalog, base.descriptor().name, order, constantInput);
        CoordinateChannel result{plan.name, base.descriptor().name, plan.method,
                                 plan.pushForward ? "" : plan.map.descriptor().units, {}};
        if (plan.pushForward) {
            if (missingInput) throw std::invalid_argument("Observable derivative needs the recorded forcing or a prescribed driver");
            if (!velocities) velocities = exactVelocities(source, true);
            result.values = plan.map.pushForwardRealBatch(source.values(), *velocities, source.samples(),
                                                         c::ObservableJacobianBackend::analytic);
        } else result.values = plan.map.evaluateRealBatch(source.values(), source.samples());
        return result;
    }

private:
    CoordinateChannel krom(const J& selector, std::size_t order) {
        if (order) throw std::invalid_argument("A KROM decomposition reference has no implicit time derivative");
        c::ObservableDecompositionReference ref;
        ref.artifact = selector.at("artifact").string();
        ref.modelIdentity = selector.at("model_identity").string();
        ref.component = selector.at("component").string();
        ref.kind = selector.at("kind").string("reconstruction");
        ref.index = count(selector.at("index"), "index", 0);
        auto observable = c::existingKromObservable(std::filesystem::u8path(ref.artifact), ref,
                                                    selector.at("name").string("krom-observable"));
        if (observable.stateDimension() != source.dimension())
            throw std::invalid_argument("KROM observable source coordinate dimension differs");
        if (selector.at("coordinate_contract").string() != observable.descriptor().stateCoordinateSystem)
            throw std::invalid_argument("KROM reference requires explicit matching coordinate_contract; dimensions alone do not identify coordinates");
        const auto output = count(selector.at("output"), "output", 0);
        if (output >= observable.outputDimension()) throw std::invalid_argument("KROM observable output index out of range");
        const auto part = selector.at("part").string("real");
        if (part != "real" && part != "imag") throw std::invalid_argument("KROM channel part must be real or imag");
        CoordinateChannel result{observable.descriptor().name, ref.component, "none", "", {}};
        for (std::size_t i = 0; i < source.samples(); ++i) {
            const auto value = observable.evaluate(source.values().data() + i * source.dimension())[output];
            result.values.push_back(part == "real" ? value.real() : value.imag());
        }
        return result;
    }

    const SourceData& source;
    std::unique_ptr<c::ODESystem> system;
    std::optional<c::ObservableCatalog> catalog;
    std::optional<std::vector<double>> velocities;
};

bool isKromReference(const J& v) {
    return std::holds_alternative<J::Object>(v.storage()) && !v.at("artifact").isNull();
}
}  // namespace

std::vector<CoordinateChannel> evaluateCoordinateChannels(const SourceData& source, const J::Array& specs) {
    if (specs.empty()) throw std::invalid_argument("coordinate channels must be a nonempty array");
    Resolver resolver(source);
    std::map<std::string, std::size_t> seen;
    std::vector<CoordinateChannel> out;
    for (const auto& item : specs) {
        const bool object = std::holds_alternative<J::Object>(item.storage()) && !isKromReference(item);
        if (object) for (const auto& entry : item.object())
            if (entry.first != "observable" && entry.first != "derivative_order")
                throw std::invalid_argument("unknown coordinate channel field: " + entry.first);
        const auto& selector = object ? item.at("observable") : item;
        if (std::holds_alternative<std::string>(selector.storage()) && selector.string().empty())
            throw std::invalid_argument("expected a coordinate observable name or index");
        const auto order = object ? count(item.at("derivative_order"), "derivative_order", 0) : 0;
        if (order > 8) throw std::invalid_argument("coordinate derivative_order exceeds 8");
        if (!seen.emplace(json::dump(selector) + "#" + std::to_string(order), out.size()).second) continue;
        out.push_back(resolver.resolve(selector, order));
    }
    return out;
}

std::vector<std::string> timeSeriesColumnNames(const SourceData& source) {
    const auto meta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    const auto& recorded = meta.at("coordinate_names");
    if (std::holds_alternative<J::Array>(recorded.storage()) && recorded.array().size() == source.dimension()) {
        std::vector<std::string> names;
        for (const auto& name : recorded.array()) names.push_back(name.string());
        return names;
    }
    if (meta.at("coordinate_provider").isNull()) {
        try {
            const auto system = c::createSystem(systemConfiguration(source));
            const auto& names = system->descriptor().stateNames;
            if (names.size() == source.dimension()) return names;
        } catch (const std::exception&) {}
    }
    std::vector<std::string> names;
    for (std::size_t k = 0; k < source.dimension(); ++k) names.push_back("state_" + std::to_string(k));
    return names;
}

bool selectorNeedsCoordinateChange(const J& item, const std::vector<std::string>& columns) {
    if (isKromReference(item)) return true;
    const bool object = std::holds_alternative<J::Object>(item.storage());
    if (object && count(item.at("derivative_order"), "derivative_order", 0) > 0) return true;
    const auto& selector = object ? item.at("observable") : item;
    if (isKromReference(selector)) return true;
    if (!std::holds_alternative<std::string>(selector.storage())) return false;
    const auto name = selector.string();
    if (std::find(columns.begin(), columns.end(), name) != columns.end()) return false;
    for (std::size_t k = 0; k < columns.size(); ++k) if (name == "state_" + std::to_string(k)) return false;
    return true;
}

} // namespace dynsyslab::workflow::postprocess_detail

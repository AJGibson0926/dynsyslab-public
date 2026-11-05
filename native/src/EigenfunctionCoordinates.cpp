#include <cmath>
#include <dynsyslab/core/EigenfunctionCoordinates.hpp>
#include <limits>
#include <stdexcept>
namespace dynsyslab::core {
bool EigenfunctionCoordinates::requiresExternalValues() const noexcept {
    return kind == EigenfunctionCoordinateKind::measured ||
           kind == EigenfunctionCoordinateKind::causalHistory;
}
void EigenfunctionCoordinates::validate() const {
    if (!dimension || (!names.empty() && names.size() != dimension) ||
        (physicalStateDimension && !*physicalStateDimension))
        throw std::invalid_argument("Coordinate dimensions/names are inconsistent");
    if (kind == EigenfunctionCoordinateKind::physicalState && physicalStateDimension &&
        *physicalStateDimension != dimension)
        throw std::invalid_argument("Physical-state coordinates must match physical dimension");
    if (requiresExternalValues() && !expressionDefinition.empty())
        throw std::invalid_argument(
            "Measured/history coordinates cannot claim an instantaneous state expression");
    if (kind == EigenfunctionCoordinateKind::analyticObservable && expressionDefinition.empty() &&
        names.empty())
        throw std::invalid_argument("Analytic observables require expressions or registered names");
}
json::Value eigenfunctionCoordinatesJson(const EigenfunctionCoordinates &c) {
    c.validate();
    const char *kind = c.kind == EigenfunctionCoordinateKind::physicalState ? "physical-state"
                       : c.kind == EigenfunctionCoordinateKind::analyticObservable
                           ? "analytic-observable"
                       : c.kind == EigenfunctionCoordinateKind::causalHistory ? "causal-history"
                                                                              : "measured";
    json::Value::Array names;
    for (const auto &n : c.names)
        names.emplace_back(n);
    json::Value::Object out{{"schema", "dynsyslab.eigenfunction-coordinates/v1"},
                            {"kind", kind},
                            {"dimension", double(c.dimension)},
                            {"names", names},
                            {"expression_definition", c.expressionDefinition},
                            {"provider_identity", c.providerIdentity},
                            {"provider_artifact", c.providerArtifact}};
    if (c.physicalStateDimension)
        out["physical_state_dimension"] = double(*c.physicalStateDimension);
    return out;
}
EigenfunctionCoordinates eigenfunctionCoordinatesFromJson(const json::Value &value) {
    if (value.at("schema").string() != "dynsyslab.eigenfunction-coordinates/v1")
        throw std::invalid_argument("Unknown coordinate schema");
    auto dimension = [](const json::Value &n) {
        const auto x = n.number(-1);
        if (!std::isfinite(x) || x < 1 || x != std::floor(x) ||
            x > double(std::numeric_limits<std::size_t>::max() / 2))
            throw std::invalid_argument("Invalid coordinate dimension");
        return std::size_t(x);
    };
    EigenfunctionCoordinates out;
    out.dimension = dimension(value.at("dimension"));
    const auto k = value.at("kind").string();
    if (k == "physical-state")
        out.kind = EigenfunctionCoordinateKind::physicalState;
    else if (k == "analytic-observable")
        out.kind = EigenfunctionCoordinateKind::analyticObservable;
    else if (k == "measured")
        out.kind = EigenfunctionCoordinateKind::measured;
    else if (k == "causal-history")
        out.kind = EigenfunctionCoordinateKind::causalHistory;
    else
        throw std::invalid_argument("Unknown coordinate kind");
    if (!value.at("physical_state_dimension").isNull())
        out.physicalStateDimension = dimension(value.at("physical_state_dimension"));
    for (const auto &n : value.at("names").array())
        out.names.push_back(n.string());
    out.expressionDefinition = value.at("expression_definition").string();
    out.providerIdentity = value.at("provider_identity").string();
    out.providerArtifact = value.at("provider_artifact").string();
    out.validate();
    return out;
}
} // namespace dynsyslab::core

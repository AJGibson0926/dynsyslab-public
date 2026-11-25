#include <dynsyslab/core/FunctionDomain.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {

bool FunctionDomain::isPhysicalState() const noexcept {
    if (!stateMap || stateMap->stateDimension() != dimension || stateMap->outputDimension() != dimension) return false;
    for (std::size_t i = 0; i < dimension; ++i)
        if (stateMap->stateComponentIndex(i) != std::optional<std::size_t>(i)) return false;
    return true;
}

void FunctionDomain::validate() const {
    if (!dimension || (!names.empty() && !stateMapFromCatalog && names.size() != dimension) ||
        (physicalStateDimension && !*physicalStateDimension))
        throw std::invalid_argument("Domain dimensions/names are inconsistent");
    if (stateMap && (stateMap->outputDimension() != dimension ||
                     (physicalStateDimension && stateMap->stateDimension() != *physicalStateDimension)))
        throw std::invalid_argument("The state map's dimensions differ from the domain");
    if (stateMapFromCatalog && names.empty())
        throw std::invalid_argument("A catalog state map needs the observable names");
    if (providerArtifact.size() && providerIdentity.empty())
        throw std::invalid_argument("A measurement-map artifact needs its identity");
    if (window && (stateMap || window->rows() != dimension))
        throw std::invalid_argument("A window domain has no state map and one coordinate per window row");
}

FunctionDomain FunctionDomain::seriesWindow(HankelWindowLayout layout) {
    FunctionDomain out;
    out.dimension = layout.rows();
    out.window = std::make_shared<const HankelWindowLayout>(std::move(layout));
    out.validate();
    return out;
}

std::string FunctionDomain::stateMapDefinition() const {
    if (!stateMap || stateMapFromCatalog || isPhysicalState() || stateMap->expressions().empty()) return {};
    return json::dump(stateMap->expressionDefinition());
}

FunctionDomain FunctionDomain::bound(const ObservableCatalog& catalog) const {
    if (!needsCatalog()) return *this;
    auto out = *this;
    out.stateMap = std::make_shared<const Observable>(catalog.coordinates(names, "domain_coordinates"));
    if (!out.physicalStateDimension) out.physicalStateDimension = out.stateMap->stateDimension();
    out.validate();
    return out;
}

FunctionDomain FunctionDomain::physicalState(std::size_t dimension) {
    std::vector<std::string> labels;
    for (std::size_t i = 0; i < dimension; ++i) labels.push_back("x_" + std::to_string(i));
    FunctionDomain out;
    out.dimension = dimension;
    out.physicalStateDimension = dimension;
    out.stateMap = std::make_shared<const Observable>(Observable::identity(dimension, labels, labels));
    return out;
}

FunctionDomain FunctionDomain::fromStateMap(Observable map) {
    FunctionDomain out;
    out.dimension = map.outputDimension();
    out.physicalStateDimension = map.stateDimension();
    if (map.outputDimension() == 1) out.names = {map.descriptor().name};
    else out.names = map.descriptor().outputNames;
    if (out.names.size() != out.dimension) out.names.clear();
    out.stateMap = std::make_shared<const Observable>(std::move(map));
    return out;
}

FunctionDomain FunctionDomain::fromCatalog(std::vector<std::string> names, std::size_t dimension,
                                           const ObservableCatalog* catalog) {
    FunctionDomain out;
    out.dimension = dimension;
    out.names = std::move(names);
    out.stateMapFromCatalog = true;
    return catalog ? out.bound(*catalog) : out;
}

FunctionDomain FunctionDomain::supplied(std::size_t dimension) {
    FunctionDomain out;
    out.dimension = dimension;
    return out;
}

namespace {
std::size_t count(const json::Value& n) {
    const auto x = n.number(-1);
    if (!std::isfinite(x) || x < 1 || x != std::floor(x) || x > double(std::numeric_limits<std::size_t>::max() / 2))
        throw std::invalid_argument("Invalid domain dimension");
    return std::size_t(x);
}
}  // namespace

json::Value functionDomainJson(const FunctionDomain& d) {
    d.validate();
    json::Value::Array names;
    for (const auto& n : d.names) names.emplace_back(n);
    std::string map = "none", definition;
    if (d.isPhysicalState()) map = "identity";
    else if (d.stateMapFromCatalog) map = "catalog";
    else if (d.stateMap) {
        definition = d.stateMapDefinition();
        if (definition.empty())
            throw std::invalid_argument("This domain's state map has no formula to store; use catalog observables");
        map = "expression";
    }
    json::Value::Object out{{"schema", "dynsyslab.function-domain/v1"}, {"dimension", double(d.dimension)},
                            {"names", names}, {"state_map", map}, {"expression_definition", definition},
                            {"measurement_map", json::Value::Object{{"identity", d.providerIdentity},
                                                                    {"artifact", d.providerArtifact}}}};
    const auto physical = d.stateMap ? std::optional<std::size_t>(d.stateMap->stateDimension()) : d.physicalStateDimension;
    if (physical) out["physical_state_dimension"] = double(*physical);
    if (d.window) out["window"] = hankelWindowLayoutJson(*d.window);
    return out;
}

FunctionDomain functionDomainFromJson(const json::Value& value) {
    const auto schema = value.at("schema").string();
    const auto dimension = count(value.at("dimension"));
    std::optional<std::size_t> physical;
    if (!value.at("physical_state_dimension").isNull()) physical = count(value.at("physical_state_dimension"));
    std::vector<std::string> names;
    for (const auto& n : value.at("names").array()) names.push_back(n.string());
    const auto definition = value.at("expression_definition").string();
    std::string map = "none", identity, artifact;
    if (schema == "dynsyslab.function-domain/v1" || schema == "dynsyslab.eigenfunction-domain/v1") {
        map = value.at("state_map").string("none");
        identity = value.at("measurement_map").at("identity").string();
        artifact = value.at("measurement_map").at("artifact").string();
    } else if (schema == "dynsyslab.eigenfunction-coordinates/v1") {  // former exclusive kinds
        const auto kind = value.at("kind").string();
        if (kind == "physical-state") map = "identity";
        else if (kind == "analytic-observable") map = definition.empty() ? "catalog" : "expression";
        else if (kind != "measured" && kind != "causal-history") throw std::invalid_argument("Unknown coordinate kind");
        identity = value.at("provider_identity").string();
        artifact = value.at("provider_artifact").string();
        if (kind == "causal-history" && identity.empty()) identity = "causal-history";
    } else {
        throw std::invalid_argument("Unknown domain schema");
    }
    FunctionDomain out;
    if (map == "identity") out = FunctionDomain::physicalState(dimension);
    else if (map == "expression") {
        if (definition.empty()) throw std::invalid_argument("An expression state map needs its definition");
        out = FunctionDomain::fromStateMap(Observable::fromExpressionDefinition(json::parse(definition)));
    } else if (map == "catalog") out = FunctionDomain::fromCatalog(names, dimension);
    else if (map == "none") out = FunctionDomain::supplied(dimension);
    else throw std::invalid_argument("Unknown state map " + map);
    if (!names.empty()) out.names = std::move(names);
    if (physical && !out.stateMap) out.physicalStateDimension = physical;
    out.providerIdentity = std::move(identity);
    out.providerArtifact = std::move(artifact);
    if (!value.at("window").isNull())
        out.window = std::make_shared<const HankelWindowLayout>(hankelWindowLayoutFromJson(value.at("window")));
    out.validate();
    return out;
}

}  // namespace dynsyslab::core

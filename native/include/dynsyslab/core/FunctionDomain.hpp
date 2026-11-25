#pragma once
// The domain a function (an observable or an eigenfunction) is evaluated on: what its coordinates q
// are, and the one optional recipe that computes them from the physical state x.
//
//   stateMap         q = q(x), an Observable of the physical state: the identity (q = x), a formula
//                    of x, or catalog observables of the system. Absent when q cannot be computed
//                    from x; the caller then supplies q (measured or history coordinates).
//   stateMapFromCatalog  the map is the system catalog's observables named by `names`. Such a map is
//                    stored by name and bound from the system's ObservableCatalog (bound()) when the
//                    domain is loaded without it.
//   measurement map  q online from measurements: a frozen provider (e.g. a Hankel basis) named by
//                    its identity and artifact. Independent of the state map.
//   window           q is the delay window of a time series at a sample (Hankel.hpp layout): the
//                    series columns read at each delay. A function on such a domain is evaluated
//                    along a series by cutting its window at each sample (Observable::evaluate).
#include <cstddef>
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct FunctionDomain {
    std::size_t dimension{};
    std::optional<std::size_t> physicalStateDimension;
    // Coordinate identities, one per dimension (optional). For a catalog map: the catalog observable
    // names, whose stacked outputs give the coordinates.
    std::vector<std::string> names;
    std::shared_ptr<const Observable> stateMap;
    bool stateMapFromCatalog{};
    std::string providerIdentity;            // measurement map identity, e.g. "sha256:..."
    std::string providerArtifact;            // archive holding that frozen map
    std::shared_ptr<const HankelWindowLayout> window;  // q = the series window at a sample

    // q = x: the state map is the identity on a state of this dimension.
    [[nodiscard]] bool isPhysicalState() const noexcept;
    [[nodiscard]] bool computableFromState() const noexcept { return stateMap || stateMapFromCatalog; }
    [[nodiscard]] bool computableFromMeasurements() const noexcept { return !providerIdentity.empty(); }
    // The caller must supply q (it cannot be computed from a state).
    [[nodiscard]] bool requiresExternalValues() const noexcept { return !computableFromState(); }
    // A catalog map known by name only, waiting for bound().
    [[nodiscard]] bool needsCatalog() const noexcept { return stateMapFromCatalog && !stateMap; }
    // This domain with its catalog map bound (unchanged when it needs no catalog).
    [[nodiscard]] FunctionDomain bound(const ObservableCatalog& catalog) const;
    // The map's formula as expression-definition JSON text when it is a formula of x (not the
    // identity, not catalog observables); empty otherwise.
    [[nodiscard]] std::string stateMapDefinition() const;
    void validate() const;

    [[nodiscard]] static FunctionDomain physicalState(std::size_t dimension);
    // A formula or other observable q(x); its output names become the coordinate names.
    [[nodiscard]] static FunctionDomain fromStateMap(Observable map);
    // The catalog observables `names` (bound later, or now when `catalog` is given).
    [[nodiscard]] static FunctionDomain fromCatalog(std::vector<std::string> names, std::size_t dimension,
                                                    const ObservableCatalog* catalog = nullptr);
    // No state map: q is supplied (measured or history coordinates).
    [[nodiscard]] static FunctionDomain supplied(std::size_t dimension);
    // q is the window `layout` cuts from a series at each sample.
    [[nodiscard]] static FunctionDomain seriesWindow(HankelWindowLayout layout);
};

// dynsyslab.function-domain/v1 (state_map: identity, expression, catalog or none); also reads
// dynsyslab.eigenfunction-domain/v1 and dynsyslab.eigenfunction-coordinates/v1.
[[nodiscard]] json::Value functionDomainJson(const FunctionDomain& domain);
[[nodiscard]] FunctionDomain functionDomainFromJson(const json::Value& value);

}  // namespace dynsyslab::core

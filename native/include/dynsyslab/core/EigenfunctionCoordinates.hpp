#pragma once
#include <cstddef>
#include <dynsyslab/core/Json.hpp>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::core {
enum class EigenfunctionCoordinateKind {
    physicalState,
    analyticObservable,
    measured,
    causalHistory
};
struct EigenfunctionCoordinates {
    EigenfunctionCoordinateKind kind{EigenfunctionCoordinateKind::measured};
    std::size_t dimension{};
    std::optional<std::size_t> physicalStateDimension;
    std::vector<std::string> names;
    std::string expressionDefinition;
    std::string providerIdentity;
    std::string providerArtifact; // frozen transform, not a live ring buffer
    bool requiresExternalValues() const noexcept;
    void validate() const;
};
json::Value eigenfunctionCoordinatesJson(const EigenfunctionCoordinates &coordinates);
EigenfunctionCoordinates eigenfunctionCoordinatesFromJson(const json::Value &value);
} // namespace dynsyslab::core

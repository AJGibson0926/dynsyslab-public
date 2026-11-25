#pragma once
#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/DomainGrid.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dsl {

// Figures for a dynsyslab.grid-field/v1 archive (DomainFieldPayload.cpp).
[[nodiscard]] Payload buildDomainFieldPayload(const std::filesystem::path& source, const dynsyslab::core::NpzRecord& record);

// A 3-D grid field as translucent planes across its third axis at `slices`, in a "three_d" figure.
// Each plane is a trace with encoding "slice" (colorValues: value normalized to [0, 1]) or
// "slice-domain-color" (colorValues: phase, magnitudes: relative magnitude); x, y, z outline the
// plane and zShape gives its nodes (rows = second axis). Invalid points are left transparent.
[[nodiscard]] Figure buildDomainField3D(const dynsyslab::core::GridField& field, const std::vector<std::size_t>& slices,
                                        const std::string& titleSuffix);

}  // namespace dsl

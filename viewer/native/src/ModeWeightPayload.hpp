#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dsl {

// "Mode weights" view: how each reconstructed observable decomposes over the
// Koopman modes, and which modes a reduced model kept.
//
// For observable o and mode j the weight is |c_oj| * rms_t |phi_j(x(t))| over the
// archived training eigenfunction values (|c_oj| alone when none are stored). It
// does not depend on how the eigenfunctions are scaled, and its square is the
// energy that the observable_energy K-ROM selection rule ranks.
//
// Panels are 3-D stems on the complex plane of eigenvalues (continuous lambda when
// stored, otherwise discrete mu with the unit circle): first the total over observables,
// sqrt(sum_o w_oj^2), then one panel per scalar observable. When the archive was
// extracted from a parent model (metadata source_krom + selected_modes) and the
// parent can be read, all parent modes are drawn in blue and the kept modes in
// red on top, at the parent's weights, with the kept share of energy in the title.
// The view is limited to the extent of the modes carrying >= 0.1% of the heaviest
// total weight; the title counts the lighter modes left outside.
// Inline so the archive catalogue can use it without linking the payload builders.
inline bool hasModeWeights(const dynsyslab::core::NpzIndex& archive) {
    if (!archive.contains("discrete_eigenvalues") && !archive.contains("continuous_eigenvalues")) return false;
    if (archive.contains("reconstruction_matrix")) return true;
    for (const auto& name : archive.arrays)
        if (name.rfind("reconstruction__", 0) == 0) return true;
    return false;
}

Payload modeWeightPayload(const std::filesystem::path& source, const std::string& view,
                          const dynsyslab::core::NpzRecord& record,
                          const dynsyslab::core::json::Value& metadata);

}  // namespace dsl

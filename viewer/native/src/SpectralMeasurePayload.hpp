#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <vector>

namespace dsl {

// Stem plots of each Hankel observable's modal weight |c_j phi_j(t0)|, based at
// the eigenvalue of mode j. Because sum_j c_j phi_j(t0) reproduces the observable
// at the first aligned sample, these weights are the archived H-DMD spectral
// measure of that observable; they do not depend on the eigenvector scaling.
// Returns no figures when the archive lacks eigenfunction values or a C map.
std::vector<Figure> spectralMeasureFigures(
    const dynsyslab::core::NpzRecord& record,
    const dynsyslab::core::json::Value& metadata);

}  // namespace dsl

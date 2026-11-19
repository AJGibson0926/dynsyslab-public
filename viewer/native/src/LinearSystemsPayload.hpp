#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// Views of a linear_systems archive (schema dynsyslab.linear-systems/v1):
//   hankel_spectrum    Hankel singular values with the 2*sum(tail) balanced
//                      truncation bound per reduced order; chosen order marked.
//   gramian_spectrum   controllability/observability gramian eigenvalues with
//                      the HSVs, and the two principal gramian ellipses.
//   frequency_response singular values (and SISO phase) of the full and the
//                      reduced model on the archived frequency grid.
//   pole_zero          poles, reduced poles and transmission zeros; feedback
//                      poles of full and reduced (lifted) LQR/LQE when present.
// Colours follow the suite palette: blue = full model, red = reduced model,
// gold = bounds, zeros and estimator quantities.
Payload buildLinearSystemsPayload(const std::filesystem::path& source, const std::string& view,
                                  const dynsyslab::core::NpzRecord& record,
                                  const dynsyslab::core::json::Value& metadata);

}  // namespace dsl

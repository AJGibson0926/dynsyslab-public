#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// Matrix-valued archive views:
//   control_map_error  constant/bilinear control-map error fields over a state grid.
//   observability      log singular spectrum and right singular directions of a K-ROM
//                      observability analysis; with mode participation, observer poles and
//                      the state-estimation error norm when the archive has them.
//   bifurcation        parameter versus post-transient maxima (bifurcation_points method).
//   observable_spectrum_comparison  H-DMD spectra by scalar observable (discrete on the unit
//                      circle, continuous) and the pairwise nearest-spectrum RMS distances.
Payload buildControlErrorPayload(const std::filesystem::path& source, const std::string& view,
                                 const dynsyslab::core::NpzRecord& record,
                                 const dynsyslab::core::json::Value& metadata);
Payload buildObservabilityPayload(const std::filesystem::path& source, const std::string& view,
                                  const dynsyslab::core::NpzRecord& record,
                                  const dynsyslab::core::json::Value& metadata);
Payload buildBifurcationPayload(const std::filesystem::path& source, const std::string& view,
                                const dynsyslab::core::NpzRecord& record,
                                const dynsyslab::core::json::Value& metadata);
Payload buildObservableSpectraPayload(const std::filesystem::path& source, const std::string& view,
                                     const dynsyslab::core::NpzRecord& record,
                                     const dynsyslab::core::json::Value& metadata);

}  // namespace dsl

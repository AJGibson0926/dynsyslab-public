#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// Classical frequency-domain views of a linear_systems archive:
//   bode         magnitude (dB) and unwrapped phase of every plant channel
//                (transfer_function / bode_* arrays).
//   nyquist      L(jw) of the archived loops (loop__<name>__response): control
//                loops broken at the plant input (LQR blue, LQG red, LTR sweep gold)
//                and loops at the plant output (Kalman gold, LQG red), with -1 and
//                the unit disk around it (the LQR / Kalman guarantee |1 + L| >= 1).
//   loop_bode    magnitude and phase of the same loops.
//   sensitivity  |S| = 1 / |1 + L| and |T| = |L / (1 + L)| in dB for each loop.
// Suite palette: blue = state feedback / plant, red = output feedback, gold =
// estimator and recovery sweeps.
bool isClassicalLinearSystemsView(const std::string& view);
Payload buildClassicalLinearSystemsPayload(const std::filesystem::path& source, const std::string& view,
                                           const dynsyslab::core::NpzRecord& record,
                                           const dynsyslab::core::json::Value& metadata);

}  // namespace dsl

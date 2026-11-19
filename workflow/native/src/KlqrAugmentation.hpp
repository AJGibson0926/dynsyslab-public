#pragma once

#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LQR.hpp>

namespace dynsyslab::workflow {
namespace json = dynsyslab::core::json;

// KLQR adapters for the general integral-action (LQI) and PID components.
//
//   "integral_action": {"reconstruction": "state", "components": [0], "penalty": 100,
//                       "integrator_limit": [..], "target_components": [0]}
//   "pid_correction":  {"reconstruction": "state", "components": [0], "Kp": .., "Ki": ..,
//                       "Kd": .., "derivative_filter_time": .., "integrator_limit": [..],
//                       "output_limits": L | [low, high], "conditional_integration": true}
//
// The controlled output is y = y0 + C z: rows `components` of the named KROM
// reconstruction (with its intercept), applied to the lifted state the KLQR acts on (psi(x)
// for full-state feedback, the KLQE estimate z_hat for estimated feedback). It is compared
// with components `target_components` (default: `components`) of the dynamical-state target.
bool klqrIsStateful(const json::Value& settings);
void applyKlqrAugmentations(core::KoopmanLQRController& controller, const core::KoopmanModel& model,
                            const json::Value& settings);
json::Value klqrAugmentationEvidence(const json::Value& settings);

}  // namespace dynsyslab::workflow

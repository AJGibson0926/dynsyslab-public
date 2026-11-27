#pragma once
// Training-domain enforcement for Koopman control (controller.settings.target_domain_policy):
//   "warn" (default): count evaluations whose plant state or state target leaves the K-ROM's
//                     training domain and archive the report (domain_report_* arrays);
//   "error":          refuse the evaluation (the solve fails with the axes and excursion);
//   "ignore":         no check.
// The domain is settings.training_domain {minimum, maximum}, else the K-ROM archive's
// training_minimum/training_maximum or the bounds of its training_states. A K-ROM extrapolates
// outside the data it was fitted on, so leaving the domain is reported, never silent.
#include <dynsyslab/workflow/ControlFactory.hpp>

#include <string>

namespace dynsyslab::workflow::detail {

void applyTargetDomainPolicy(ControllerSetup& setup, const json::Value& settings,
                             const std::string& kromPath, std::size_t stateDimension,
                             dynsyslab::core::StateTargetFunction stateTarget,
                             std::size_t trajectoryCount);

// The K-ROM's axis-aligned training domain (empty bounds when it records none).
[[nodiscard]] std::pair<std::vector<double>, std::vector<double>> kromTrainingDomain(
    const std::string& kromPath, const json::Value& settings, std::size_t stateDimension);

}  // namespace dynsyslab::workflow::detail

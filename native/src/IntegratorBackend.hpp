#pragma once

#include <dynsyslab/core/Integrator.hpp>

namespace dynsyslab::core::detail {

[[nodiscard]] NativeSolution integrateBDF(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t trajectoryCount, const IntegrationOptions& options,
    ScalarControlPolicy controller, ProgressCallback progress,
    const ExecutionContext* execution);

}  // namespace dynsyslab::core::detail

#include "IntegratorBackend.hpp"

#include <stdexcept>

namespace dynsyslab::core::detail {

NativeSolution integrateBDF(ODESystem&, const std::vector<double>&, std::size_t,
                            const IntegrationOptions&, ScalarControlPolicy,
                            ProgressCallback, const ExecutionContext*) {
    throw std::runtime_error("BDF integration is unavailable in this native build");
}

}  // namespace dynsyslab::core::detail

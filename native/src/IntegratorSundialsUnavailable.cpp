#include "IntegratorBackend.hpp"

#include <stdexcept>

namespace dynsyslab::core::detail {

std::unique_ptr<Stepper> makeCvodeStepper(const IntegrationOptions& options, std::size_t) {
    throw std::runtime_error(integrationMethodName(options.method)
        + " integration needs SUNDIALS CVODE, which is not in this native build");
}

}  // namespace dynsyslab::core::detail

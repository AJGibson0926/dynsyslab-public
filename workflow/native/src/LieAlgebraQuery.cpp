#include <dynsyslab/core/LieAlgebra.hpp>
#include <dynsyslab/workflow/LieAlgebraQuery.hpp>
namespace dynsyslab::workflow {
dynsyslab::core::json::Value
evaluateLieAlgebraQuery(const dynsyslab::core::json::Value &request) {
  return dynsyslab::core::lie::analyze(request);
}
} // namespace dynsyslab::workflow

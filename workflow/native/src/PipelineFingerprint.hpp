#pragma once
#include <string>

namespace dynsyslab::workflow::detail {

// Identity of the pipeline that is running: a SHA-256 over the DynSysLab binaries loaded with this
// code (the module holding it and every dynsyslab* library beside it). Any rebuild changes it, so an
// experiment reuses a stage only when the same pipeline computed it: change a method, rerun the same
// experiment, and the changed stages recompute. Computed once per process.
[[nodiscard]] const std::string& pipelineFingerprint();

}  // namespace dynsyslab::workflow::detail

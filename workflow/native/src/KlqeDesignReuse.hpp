#pragma once
// Reload a verified KLQE design (KlqeDesignReuse.cpp).
#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <string>

namespace dynsyslab::workflow {

// Replaces design's gain, covariance and poles with the archive's; returns the evidence.
json::Value reuseKlqeDesign(core::KoopmanLQEDesign& design, const std::string& archivePath,
                            const core::KoopmanModel& model);

}  // namespace dynsyslab::workflow

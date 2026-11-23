#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dynsyslab::workflow {
core::json::Value auditKromTrajectory(const core::KoopmanModel &,
                                      const core::json::Value &);
}

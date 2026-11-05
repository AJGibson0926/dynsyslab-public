#pragma once

#include <dynsyslab/workflow/Json.hpp>

namespace dynsyslab::workflow::detail {

void loadExtensions(const json::Value& job);

}  // namespace dynsyslab::workflow::detail

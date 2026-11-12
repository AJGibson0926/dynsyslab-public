#pragma once
#include <dynsyslab/workflow/PhaseSpaceQuery.hpp>
namespace dynsyslab::workflow::detail {
json::Value phaseSpaceResultJson(const dynsyslab::core::PhaseSpaceResult& result);
}

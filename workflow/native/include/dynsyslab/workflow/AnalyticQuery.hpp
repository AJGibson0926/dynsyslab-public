#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/System.hpp>
namespace dynsyslab::workflow {
namespace json = dynsyslab::core::json;
json::Value analyticDescriptorJson(const dynsyslab::core::AnalyticSolutionDescriptor &descriptor);
json::Value evaluateAnalyticQuery(const json::Value &request);
// Operation "phase_cover": a checked prime-power orbit tower (PhaseCoverQuery.cpp).
json::Value evaluatePhaseCoverQuery(const json::Value &request);
json::Value evaluateSpecialFunctionQuery(const json::Value &request);
} // namespace dynsyslab::workflow

#pragma once
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Json.hpp>
namespace dynsyslab::workflow {
core::json::Value evaluateTorusAnalysis(const core::json::Value&,const core::ODESystem&);
}

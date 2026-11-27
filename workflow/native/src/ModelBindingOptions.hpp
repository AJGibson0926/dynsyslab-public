#pragma once
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/workflow/Json.hpp>
namespace dynsyslab::workflow {
inline dynsyslab::core::ModelBindingOptions modelBindingOptions(const json::Value& settings) {
    dynsyslab::core::ModelBindingOptions options;
    options.systemTransferReason=settings.at("system_transfer_reason").string();
    return options;
}
}

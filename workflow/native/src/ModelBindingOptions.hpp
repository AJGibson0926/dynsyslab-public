#pragma once
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <dynsyslab/workflow/Json.hpp>
namespace dynsyslab::workflow {
inline dynsyslab::core::ModelBindingOptions modelBindingOptions(const json::Value& settings) {
    dynsyslab::core::ModelBindingOptions options;
    options.allowUnverifiedLegacy=false;
    options.unverifiedLegacyReason=settings.at("unverified_legacy_reason").string();
    options.systemTransferReason=settings.at("system_transfer_reason").string();
    return options;
}
}

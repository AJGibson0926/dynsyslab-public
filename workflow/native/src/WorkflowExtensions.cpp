#include "WorkflowExtensions.hpp"

#include <dynsyslab/core/Extension.hpp>
#include <dynsyslab/core/Text.hpp>

#include <stdexcept>
#include <string>

namespace dynsyslab::workflow::detail {

void loadExtensions(const json::Value& job) {
    for (const auto& entry : job.at("extensions").array()) {
        std::string encoded = entry.string();
        if (encoded.empty()) encoded = entry.at("path").string();
        if (encoded.empty())
            throw std::invalid_argument("extension entry requires a library path");
        dynsyslab::core::ExtensionRegistry::instance().load(
            dynsyslab::core::wide(encoded));
    }
}

}  // namespace dynsyslab::workflow::detail

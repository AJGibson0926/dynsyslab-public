#pragma once
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dynsyslab::core {
// Read-only, versioned summary for CLI/clients; no fitting or validation occurs.
json::Value inspectKoopmanModel(const KoopmanModel &model);
} // namespace dynsyslab::core

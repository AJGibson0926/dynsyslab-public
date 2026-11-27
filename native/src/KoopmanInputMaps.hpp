#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::core::detail {
// `ownOnly` writes only the maps the model owns (not those inherited from its parent krom).
void writeNamedInputMaps(NpzRecord& archive, const KoopmanModel& model, bool ownOnly = false);
void readNamedInputMaps(const NpzRecord& archive, KoopmanModel& model);
}  // namespace dynsyslab::core::detail

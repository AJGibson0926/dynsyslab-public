#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
namespace dynsyslab::core::detail {
bool hasPolynomialContract(const NpzRecord &);
KoopmanModel readPolynomialContract(const NpzRecord &, const SystemConfiguration &);
void writePolynomialContract(NpzRecord &, const KoopmanModel &, const std::string& metadataJson);
void readNamedInputMaps(const NpzRecord &, KoopmanModel &);
void writeNamedInputMaps(NpzRecord &, const KoopmanModel &);
} // namespace dynsyslab::core::detail

#pragma once
// Writes the krom data object of an expansion fit (DMD, EDMD, eigenfunction expansions) beside a
// method's own arrays.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <optional>
#include <string>

namespace dynsyslab::workflow {

// The domain of state-snapshot fits: the physical state when the snapshots are a native solve's
// state values of matching width, otherwise data coordinates with no state map.
[[nodiscard]] dynsyslab::core::FunctionDomain snapshotDomain(const std::string& sourceMetadataJson,
    std::size_t dimension, const std::optional<dynsyslab::core::SystemConfiguration>& system);

// Makes `archive` the krom of an expansion fit (core/KoopmanArchive.hpp): the spectrum and one
// expansion whose row j is eigenfunction j's own coefficients, with their values at the samples (rows
// of the input series) when given, stamped as a krom of the running step. Arrays already in
// `archive` stay as the producing method's properties. Returns false (with `reason`, `archive`
// unchanged) when the K-ROM cannot be formed.
bool writeExpansionKrom(dynsyslab::core::NpzRecord& archive, const dynsyslab::core::EigenfunctionFitResult& bank,
    const dynsyslab::core::DataOrigin& origin, const std::vector<double>& time,
    const std::vector<std::complex<double>>& values, const std::vector<dynsyslab::core::ChannelDescriptor>& channels,
    const std::string& metadataJson, std::string& reason);

}  // namespace dynsyslab::workflow

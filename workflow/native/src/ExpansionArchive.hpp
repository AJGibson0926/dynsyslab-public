#pragma once
// Stores a fitted eigenfunction expansion as a loadable K-ROM contract beside a method's own arrays
// (DMD, EDMD, DMDc), so loadKoopmanModel evaluates it in its own basis.
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

// Replaces `archive` with the K-ROM archive of `bank` (keeping every existing array) and returns
// true; on failure (e.g. a zero multiplier) leaves it unchanged, notes why in `reason`, returns false.
bool attachExpansionContract(dynsyslab::core::NpzRecord& archive, const dynsyslab::core::EigenfunctionFitResult& bank,
    const std::optional<dynsyslab::core::SystemConfiguration>& system, const std::string& metadataJson,
    std::string& reason);

// Writes the eigenfunctions data object of an expansion fit (core/DataObject.hpp) into `archive`:
// the spectrum and one expansion whose row j is eigenfunction j's own coefficients, with their values
// at the samples (rows of the input series) when given. Arrays already in `archive` stay as the
// producing method's properties. Returns false (with `reason`) when the set cannot be formed.
bool writeExpansionEigenfunctions(dynsyslab::core::NpzRecord& archive, const dynsyslab::core::EigenfunctionFitResult& bank,
    const dynsyslab::core::DataOrigin& origin, const std::vector<double>& time,
    const std::vector<std::complex<double>>& values, const std::vector<dynsyslab::core::ChannelDescriptor>& channels,
    const std::string& metadataJson, std::string& reason);

}  // namespace dynsyslab::workflow

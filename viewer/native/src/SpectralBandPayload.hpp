#pragma once

#include "dynsyslab_viewer/Model.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <string>

namespace dsl {

// The "spectral_measures" view of a krom: its Spectrum's Lebesgue decomposition of each observable.
// For every observable (and output) with a measure on a band:
//   - a 3-D stem plot on the unit circle: the atoms |c_j|^2 ||psi_j||^2 of the measure (one per
//     eigenpair, at its multiplier mu_j; Spectrum::atoms) and, as curves around the circle, the
//     density of every unit-circle band at e^{i theta} (each band's own grid);
//   - the same against the band parameter (angle or frequency): each band's density, the atoms as
//     stems at arg mu_j, so a reference band and a fitted band of one observable overlay;
//   - the wave packets' amplitudes when a measure carries packets is left to the workflow's
//     packet_validation report (packets are functions of the state, evaluated along records).
// Densities are per unit of the band parameter, atoms are masses: the height axis carries both, as
// the measure's Lebesgue decomposition does (Colbrook, Drysdale & Horning, arXiv:2405.00782, fig. 2).
Payload spectralMeasuresPayload(const std::filesystem::path& source, const std::string& view);

// Whether the archive's spectrum has bands or measures (the catalogue offers the view then).
inline bool hasSpectralMeasures(const dynsyslab::core::NpzIndex& archive) {
    return archive.contains("spectral_measures_json") || archive.contains("spectral_bands_json")
        || archive.contains("continuous_spectrum_json");
}

}  // namespace dsl

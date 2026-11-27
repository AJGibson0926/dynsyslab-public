#pragma once
// The one archive encoding of a Spectrum, used by every object that stores eigenvalues:
//   continuous_eigenvalues  the continuous rates lambda (when every mode has one);
//   discrete_eigenvalues    the sampled multipliers mu (when the spectrum is sampled);
//   sample_interval         dt (when sampled);
//   logarithm_branches      the branch of log(mu) each rate uses (when sampled);
//   eigenvalue_domain       "continuous" or "discrete": the domain the method produced them in.
// and its evidence, when present (Spectrum.hpp):
//   eigenpair_residual      one value per eigenpair (NaN where a pair has none), with
//   eigenpair_residual_json one {definition, scope} or null per eigenpair;
//   pseudospectrum_grid, pseudospectrum_minimum, pseudospectrum_minimizing_coefficients (grid x
//                           features) and pseudospectrum_json {definition, scope, features};
//   continuous_spectrum_json the continuous components (spectral measures).
// Reading rebuilds the spectrum exactly (values, rates, branches and domain), so anything keyed to
// a spectrum (a model's identity, its input maps) survives a round trip.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Spectrum.hpp>

namespace dynsyslab::core {

void writeSpectrumArrays(NpzRecord& archive, const Spectrum& spectrum);
[[nodiscard]] bool hasSpectrumArrays(const NpzRecord& archive);
[[nodiscard]] Spectrum readSpectrumArrays(const NpzRecord& archive, DataOrigin origin = {});
// The evidence alone (a child krom that adds residuals or a pseudospectrum to its parent's spectrum).
void writeSpectralEvidence(NpzRecord& archive, const Spectrum& spectrum);
[[nodiscard]] bool hasSpectralEvidence(const NpzRecord& archive);
// Lays the archive's evidence over `spectrum`; each part it lacks is left as it is.
void readSpectralEvidence(const NpzRecord& archive, Spectrum& spectrum);

}  // namespace dynsyslab::core

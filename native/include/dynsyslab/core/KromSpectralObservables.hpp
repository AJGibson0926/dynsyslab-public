#pragma once
// Observables defined by a krom's spectral fit, referenced by artifact and identity
// (ObservableDecompositionReference), evaluated on the frozen dictionary's coordinates:
//   "reconstruction"  a readout of dictionary features of the compressed operator (component: its name)
//   "eigenfunction"   eigenfunction `index` (a pair of the krom)
//   "wave-packet"     the generalized eigenfunctions of measure `index` of the continuous spectrum, one
//                     output per grid point (Rigged DMD; Colbrook, Drysdale & Horning 2024)
// The reference names the compression's identity (koopmanCompressionIdentity). A krom without a
// compression offers its named readouts ("reconstruction"), under koopmanModelIdentity.
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <filesystem>
#include <string>

namespace dynsyslab::core {
// The identity a reference to this krom's observables names.
[[nodiscard]] std::string kromObservableIdentity(const KoopmanModel& model);
[[nodiscard]] Observable kromSpectralObservable(const KoopmanModel& model, ObservableDecompositionReference reference,
                                                std::string name);
[[nodiscard]] Observable existingKromObservable(const std::filesystem::path& krom, ObservableDecompositionReference reference,
                                                std::string name);
}  // namespace dynsyslab::core

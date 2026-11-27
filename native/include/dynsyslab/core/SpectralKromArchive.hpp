#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SpectralKROM.hpp>
#include <dynsyslab/core/Spectrum.hpp>
namespace dynsyslab::core {
// A spectral fit's output on a Spectrum of its eigenvalues: residuals on the eigenpairs, the residual
// grid as the pseudospectrum, measures as continuous components (Colbrook & Townsend 2024; Colbrook
// 2023). The reverse rebuilds the fit's SpectralData from a krom's spectrum and diagnostics.
void attachSpectralData(Spectrum& spectrum, const SpectralData& data);
SpectralData spectralDataOf(const Spectrum& spectrum, const SpectralDiagnostics& diagnostics);
// The archive of a spectral KROM: its operator contract and spectrum. Inside a krom the contract is
// written without its SpectralData (the krom keeps that once, on its spectrum and diagnostics);
// reading such an archive restores it.
NpzRecord makeSpectralKromArchive(const SpectralKROM&,const std::string& metadataJson="{}");
SpectralKROM readSpectralKromArchive(const NpzRecord&);
SpectralKROM loadSpectralKrom(const std::filesystem::path&);
std::string spectralKromIdentity(const SpectralKROM&);
Observable spectralKromObservable(SpectralKROM,ObservableDecompositionReference,std::string name);
Observable existingKromObservable(const std::filesystem::path&,ObservableDecompositionReference,std::string name);
}

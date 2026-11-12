#pragma once
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/SpectralKROM.hpp>
namespace dynsyslab::core {
NpzRecord makeSpectralKromArchive(const SpectralKROM&,const std::string& metadataJson="{}");
SpectralKROM readSpectralKromArchive(const NpzRecord&);
SpectralKROM loadSpectralKrom(const std::filesystem::path&);
std::string spectralKromIdentity(const SpectralKROM&);
Observable spectralKromObservable(SpectralKROM,ObservableDecompositionReference,std::string name);
Observable existingKromObservable(const std::filesystem::path&,ObservableDecompositionReference,std::string name);
}

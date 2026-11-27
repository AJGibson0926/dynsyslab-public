#pragma once
// A component of the continuous part of a Koopman spectrum: the spectral measure of an observable on
// a support in the complex plane (the unit circle for a measure-preserving map, parametrized by angle),
// as a density on a grid plus any atoms and moments. It has no eigenfunctions. ResDMD/Rigged-DMD
// smoothed measures (Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024; Colbrook, Drysdale &
// Horning, arXiv:2405.00782, 2024) and the measure-preserving EDMD measure (Colbrook, SIAM J. Numer.
// Anal. 61, 2023) are such components; methods that compute none leave the continuous part empty.
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct SpectralMeasure {
    std::string name, domain{"unit-circle"}, method, assumptions;
    std::vector<double> grid, density, atomLocations, atomMasses;
    std::vector<std::complex<double>> moments;
    double totalMass{}, bandwidth{};
    std::size_t order{};
};

}  // namespace dynsyslab::core

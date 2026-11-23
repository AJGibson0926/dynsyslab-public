#pragma once
#include <dynsyslab/math/TorusFourier.hpp>
#include <cstddef>
#include <vector>
namespace dynsyslab::math {
struct SmallDivisorRecord {FourierIndex index;double divisor{},scaledDivisor{};};
// Exhaustive FINITE |k|_1<=maximumOrder search. Never certifies all integer k.
std::vector<SmallDivisorRecord> smallDivisors(const std::vector<double>&frequencies,
    unsigned maximumOrder,double exponent,std::size_t retain=20);
struct NormalFormMonomial {std::size_t component{};std::vector<unsigned> powers;std::complex<double> coefficient;};
struct NormalFormStep {std::vector<NormalFormMonomial> change,retained;double minimumRemovedDivisor{};unsigned degree{};};
// One homogeneous order for a DIAGONAL semisimple linear part. Convention x=y+h(y):
// new nonlinear term = old - (Dh Ay - A h). Lower-order substitutions must already
// have been applied by the caller; this is not a full multi-order normal form.
NormalFormStep homogeneousNormalForm(const std::vector<std::complex<double>>& eigenvalues,
    const std::vector<NormalFormMonomial>&terms,double divisorTolerance);
// Sampled first-order planar Melnikov integral. Missing infinite tails and
// quadrature error are NOT bounded by this function.
double sampledMelnikov(const std::vector<double>&time,const std::vector<double>&drift,
    const std::vector<double>&perturbation);
// A^T J + J A + gamma J for a constant symplectic form J. gamma=0: Hamiltonian
// infinitesimal defect; positive gamma: conformally symplectic damping convention.
double conformalSymplecticDefect(const std::vector<double>&jacobian,
    const std::vector<double>&form,std::size_t dimension,double damping=0);
}

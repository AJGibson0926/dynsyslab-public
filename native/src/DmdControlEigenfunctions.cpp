#include <dynsyslab/core/DmdControl.hpp>
#include <stdexcept>

namespace dynsyslab::core {
EigenfunctionFitResult DmdControlResult::eigenfunctionSet() const {
    if(!hasEigenfunctions||!dimension||discreteEigenvalues.size()!=rank)
        throw std::invalid_argument("DMDc has no exportable linear eigenbasis");
    FunctionDomain q;q.dimension=dimension;
    return {FunctionExpansion({dimension,polynomialExponentTable(dimension,1,false),{},{}},eigenfunctionCoefficients,q),
        Spectrum::fromDiscreteMultipliers(discreteEigenvalues,deltaTime),{}};
}
KoopmanSINDyResult DmdControlResult::expansion() const {
    if(!hasEigenfunctions||!dimension||discreteEigenvalues.size()!=rank)
        throw std::invalid_argument("DMDc has no exportable linear eigenbasis");
    return KoopmanSINDyResult::fromDiscreteCoefficients(dimension,
        polynomialExponentTable(dimension,1,false),eigenfunctionCoefficients,
        discreteEigenvalues,deltaTime);
}
} // namespace dynsyslab::core

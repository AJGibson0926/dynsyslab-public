#include <dynsyslab/core/DmdControl.hpp>
#include <stdexcept>

namespace dynsyslab::core {
PolynomialEigenfunctionSet DmdControlResult::eigenfunctionSet() const {
    if(!hasEigenfunctions||!dimension||discreteEigenvalues.size()!=rank)
        throw std::invalid_argument("DMDc has no exportable linear eigenbasis");
    EigenfunctionCoordinates q;q.dimension=dimension;
    return {{dimension,polynomialExponentTable(dimension,1,false),{},{}},eigenfunctionCoefficients,
        EigenfunctionSpectrum::fromDiscreteMultipliers(discreteEigenvalues,deltaTime),q};
}
KoopmanSINDyResult DmdControlResult::polynomialEigenfunctions() const {
    if(!hasEigenfunctions||!dimension||discreteEigenvalues.size()!=rank)
        throw std::invalid_argument("DMDc has no exportable linear eigenbasis");
    return KoopmanSINDyResult::fromDiscreteCoefficients(dimension,
        polynomialExponentTable(dimension,1,false),eigenfunctionCoefficients,
        discreteEigenvalues,deltaTime);
}
} // namespace dynsyslab::core

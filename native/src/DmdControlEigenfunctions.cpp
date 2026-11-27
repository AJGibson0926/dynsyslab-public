#include <dynsyslab/core/DmdControl.hpp>
#include <stdexcept>

namespace dynsyslab::core {
EigenfunctionFitResult DmdControlResult::eigenfunctionSet() const {
    if(!hasEigenfunctions||!dimension||spectrum.size()!=rank)
        throw std::invalid_argument("DMDc has no exportable linear eigenbasis");
    FunctionDomain q;q.dimension=dimension;
    return {FunctionExpansion({dimension,polynomialExponentTable(dimension,1,false),{},{}},eigenfunctionCoefficients,q),
        spectrum,{}};
}
} // namespace dynsyslab::core

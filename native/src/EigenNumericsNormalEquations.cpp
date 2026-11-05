#include <dynsyslab/numerics/NormalEquations.hpp>

#include <Eigen/Dense>
#include <Eigen/QR>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<double> solveNormalEquations(
    const std::vector<double>& values,std::size_t rows,std::size_t columns,
    const std::vector<double>& rhsValues,std::size_t outputs,double ridge) {
    if(!std::isfinite(ridge)||ridge<0.0)
        throw std::invalid_argument(
            "normal-equation ridge must be finite and nonnegative");
    if(!rows||!columns||!outputs||values.size()!=rows*columns
        ||rhsValues.size()!=rows*outputs)
        throw std::invalid_argument("normal-equation matrix shapes are inconsistent");
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const auto matrix=Eigen::Map<const Matrix>(values.data(),rows,columns);
    const auto rhs=Eigen::Map<const Matrix>(rhsValues.data(),rows,outputs);
    Matrix gram=matrix.transpose()*matrix;
    gram.diagonal().array()+=ridge;
    const Matrix cross=matrix.transpose()*rhs;
    const Matrix result=gram.completeOrthogonalDecomposition().solve(cross);
    return {result.data(),result.data()+result.size()};
}

}  // namespace dynsyslab::numerics

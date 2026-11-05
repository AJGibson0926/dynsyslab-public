#include <dynsyslab/numerics/LeastSquares.hpp>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<double> solveLeastSquares(
    const std::vector<double>& values,std::size_t rows,std::size_t columns,
    const std::vector<double>& rhsValues,std::size_t outputs,
    const LeastSquaresOptions& options) {
    if(!std::isfinite(options.ridge)||options.ridge<0.0
        ||!std::isfinite(options.singularValueTolerance)
        ||options.singularValueTolerance<0.0)
        throw std::invalid_argument(
            "least-squares tolerances must be finite and nonnegative");
    if(!rows||!columns||!outputs||values.size()!=rows*columns
        ||rhsValues.size()!=rows*outputs)
        throw std::invalid_argument("least-squares matrix shapes are inconsistent");
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const auto matrix=Eigen::Map<const Matrix>(values.data(),rows,columns);
    const auto rhs=Eigen::Map<const Matrix>(rhsValues.data(),rows,outputs);
    Matrix result;
    if(options.ridge>0.0) {
        Matrix normal=matrix.transpose()*matrix;
        normal.diagonal().array()+=options.ridge;
        result=normal.ldlt().solve(matrix.transpose()*rhs);
    } else if(options.method==LeastSquaresMethod::divideAndConquerSvd) {
        Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(matrix);
        if(options.singularValueTolerance>0.0)
            svd.setThreshold(options.singularValueTolerance);
        result=svd.solve(rhs);
    } else result=matrix.colPivHouseholderQr().solve(rhs);
    return {result.data(),result.data()+result.size()};
}

}  // namespace dynsyslab::numerics

#include <dynsyslab/numerics/LeastSquares.hpp>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <cmath>
#include <complex>
#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<std::complex<double>> solveLeastSquares(
    const std::vector<std::complex<double>>& values,std::size_t rows,
    std::size_t columns,
    const std::vector<std::complex<double>>& rhsValues,std::size_t outputs,
    const LeastSquaresOptions& options) {
    if(!std::isfinite(options.ridge)||options.ridge<0.0
        ||!std::isfinite(options.singularValueTolerance)
        ||options.singularValueTolerance<0.0)
        throw std::invalid_argument(
            "least-squares tolerances must be finite and nonnegative");
    if(!rows||!columns||!outputs||values.size()!=rows*columns
        ||rhsValues.size()!=rows*outputs)
        throw std::invalid_argument("least-squares matrix shapes are inconsistent");
    using Scalar=std::complex<double>;
    using Matrix=Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const auto matrix=Eigen::Map<const Matrix>(values.data(),rows,columns);
    const auto rhs=Eigen::Map<const Matrix>(rhsValues.data(),rows,outputs);
    Matrix result;
    if(options.ridge>0.0) {
        Matrix normal=matrix.adjoint()*matrix;
        normal.diagonal().array()+=options.ridge;
        result=normal.ldlt().solve(matrix.adjoint()*rhs);
    } else if(options.method==LeastSquaresMethod::divideAndConquerSvd) {
        Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(matrix);
        if(options.singularValueTolerance>0.0)
            svd.setThreshold(options.singularValueTolerance);
        result=svd.solve(rhs);
    } else result=matrix.colPivHouseholderQr().solve(rhs);
    return {result.data(),result.data()+result.size()};
}

}  // namespace dynsyslab::numerics

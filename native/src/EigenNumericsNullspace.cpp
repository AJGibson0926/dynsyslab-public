#include <dynsyslab/numerics/Nullspace.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<std::complex<double>> smallestRightSingularVector(
    const std::vector<std::complex<double>>& values,std::size_t rows,
    std::size_t columns) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument("complex SVD matrix shape does not match its data");
    using Matrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,
        Eigen::Dynamic,Eigen::RowMajor>;
    const auto source=Eigen::Map<const Matrix>(values.data(),rows,columns);
    Eigen::JacobiSVD<Matrix,Eigen::ComputeThinV> svd(source);
    if(!svd.matrixV().cols())
        throw std::runtime_error("complex SVD has no right singular vector");
    const Eigen::VectorXcd vector=svd.matrixV().col(svd.matrixV().cols()-1);
    return {vector.data(),vector.data()+vector.size()};
}

}  // namespace dynsyslab::numerics

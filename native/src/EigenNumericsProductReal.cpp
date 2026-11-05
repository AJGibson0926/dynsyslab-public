#include <dynsyslab/numerics/Product.hpp>

#include <Eigen/Dense>

#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<double> multiply(
    const std::vector<double>& left,std::size_t rows,std::size_t inner,
    const std::vector<double>& right,std::size_t columns) {
    if(!rows||!inner||!columns||left.size()!=rows*inner
        ||right.size()!=inner*columns)
        throw std::invalid_argument("real matrix-product shapes are inconsistent");
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const auto lhs=Eigen::Map<const Matrix>(left.data(),rows,inner);
    const auto rhs=Eigen::Map<const Matrix>(right.data(),inner,columns);
    Matrix result(rows,columns);
    result.noalias()=lhs*rhs;
    return {result.data(),result.data()+result.size()};
}

}  // namespace dynsyslab::numerics

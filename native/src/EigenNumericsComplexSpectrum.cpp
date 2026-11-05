#include <dynsyslab/numerics/Svd.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <stdexcept>

namespace dynsyslab::numerics {

std::vector<double> singularValues(
    const std::vector<std::complex<double>>& values,std::size_t rows,
    std::size_t columns) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument("complex SVD matrix shape does not match its data");
    using Matrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,
        Eigen::Dynamic,Eigen::RowMajor>;
    const auto source=Eigen::Map<const Matrix>(values.data(),rows,columns);
    Eigen::BDCSVD<Matrix> svd(source);
    return {svd.singularValues().data(),
        svd.singularValues().data()+svd.singularValues().size()};
}

}  // namespace dynsyslab::numerics

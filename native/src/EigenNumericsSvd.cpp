#include <dynsyslab/numerics/Svd.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dynsyslab::numerics {
namespace {
using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

Matrix mapped(const std::vector<double>& values,std::size_t rows,
              std::size_t columns,const char* name) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument(std::string(name)+" shape does not match its data");
    return Eigen::Map<const Matrix>(values.data(),
        static_cast<Eigen::Index>(rows),static_cast<Eigen::Index>(columns));
}
}  // namespace

RealSvd singularValueDecomposition(
    const std::vector<double>& values,std::size_t rows,std::size_t columns,
    std::size_t requestedRank,double relativeTolerance) {
    if(!std::isfinite(relativeTolerance)||relativeTolerance<0.0)
        throw std::invalid_argument("SVD tolerance must be finite and nonnegative");
    const auto source=mapped(values,rows,columns,"SVD matrix");
    Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(source);
    const auto& singular=svd.singularValues();
    std::size_t numerical=static_cast<std::size_t>(singular.size());
    if(relativeTolerance>0.0 && singular.size()) {
        const double threshold=relativeTolerance*singular[0]; numerical=0;
        while(numerical<static_cast<std::size_t>(singular.size())
              &&singular[static_cast<Eigen::Index>(numerical)]>threshold)++numerical;
    }
    const std::size_t rank=requestedRank?std::min(requestedRank,numerical):numerical;
    const Matrix left=svd.matrixU().leftCols(static_cast<Eigen::Index>(rank));
    const Matrix right=svd.matrixV().leftCols(
        static_cast<Eigen::Index>(rank)).transpose();
    RealSvd result;result.rows=rows;result.columns=columns;result.rank=rank;
    result.left.assign(left.data(),left.data()+left.size());
    result.rightTranspose.assign(right.data(),right.data()+right.size());
    result.singularValues.assign(singular.data(),
        singular.data()+static_cast<Eigen::Index>(rank));
    return result;
}

}  // namespace dynsyslab::numerics

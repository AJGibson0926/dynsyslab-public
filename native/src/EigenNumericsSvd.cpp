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
thread_local SvdDriver currentDriver = SvdDriver::divideAndConquer;
}  // namespace

void setSvdDriver(SvdDriver driver) { currentDriver = driver; }
SvdDriver svdDriver() { return currentDriver; }
SvdDriver svdDriverFromName(const std::string& name) {
    if (name.empty() || name == "gesdd" || name == "divide_and_conquer") return SvdDriver::divideAndConquer;
    if (name == "gesvd" || name == "jacobi") return SvdDriver::jacobi;
    throw std::invalid_argument("svd_driver must be gesdd, divide_and_conquer, gesvd or jacobi");
}
const char* svdDriverName(SvdDriver driver) { return driver == SvdDriver::jacobi ? "jacobi" : "divide_and_conquer"; }

RealSvd singularValueDecomposition(
    const std::vector<double>& values,std::size_t rows,std::size_t columns,
    std::size_t requestedRank,double relativeTolerance) {
    if(!std::isfinite(relativeTolerance)||relativeTolerance<0.0)
        throw std::invalid_argument("SVD tolerance must be finite and nonnegative");
    const auto source=mapped(values,rows,columns,"SVD matrix");
    Eigen::VectorXd singular;Matrix u,v;
    if(currentDriver==SvdDriver::jacobi){Eigen::JacobiSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(source);
        singular=svd.singularValues();u=svd.matrixU();v=svd.matrixV();}
    else{Eigen::BDCSVD<Matrix,Eigen::ComputeThinU|Eigen::ComputeThinV> svd(source);
        singular=svd.singularValues();u=svd.matrixU();v=svd.matrixV();}
    std::size_t numerical=static_cast<std::size_t>(singular.size());
    if(relativeTolerance>0.0 && singular.size()) {
        const double threshold=relativeTolerance*singular[0]; numerical=0;
        while(numerical<static_cast<std::size_t>(singular.size())
              &&singular[static_cast<Eigen::Index>(numerical)]>threshold)++numerical;
    }
    const std::size_t rank=requestedRank?std::min(requestedRank,numerical):numerical;
    const Matrix left=u.leftCols(static_cast<Eigen::Index>(rank));
    const Matrix right=v.leftCols(
        static_cast<Eigen::Index>(rank)).transpose();
    RealSvd result;result.rows=rows;result.columns=columns;result.rank=rank;
    result.left.assign(left.data(),left.data()+left.size());
    result.rightTranspose.assign(right.data(),right.data()+right.size());
    result.singularValues.assign(singular.data(),
        singular.data()+static_cast<Eigen::Index>(rank));
    return result;
}

ComplexSvd complexSingularValueDecomposition(
    const std::vector<std::complex<double>>& values,std::size_t rows,std::size_t columns) {
    using ComplexMatrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument("complex SVD matrix shape does not match its data");
    const ComplexMatrix source=Eigen::Map<const ComplexMatrix>(values.data(),
        static_cast<Eigen::Index>(rows),static_cast<Eigen::Index>(columns));
    Eigen::VectorXd s;ComplexMatrix v;  // v: columns x rank
    if(currentDriver==SvdDriver::jacobi){Eigen::JacobiSVD<ComplexMatrix,Eigen::ComputeThinV> svd(source);s=svd.singularValues();v=svd.matrixV();}
    else{Eigen::BDCSVD<ComplexMatrix,Eigen::ComputeThinV> svd(source);s=svd.singularValues();v=svd.matrixV();}
    ComplexSvd result;result.rows=rows;result.columns=columns;
    result.rank=static_cast<std::size_t>(s.size());
    result.singularValues.assign(s.data(),s.data()+s.size());
    result.right.assign(v.data(),v.data()+v.size());
    return result;
}

}  // namespace dynsyslab::numerics

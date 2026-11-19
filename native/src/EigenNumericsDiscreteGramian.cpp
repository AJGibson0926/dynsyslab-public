#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include "EigenNumericsLinearSystems.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <stdexcept>
namespace dynsyslab::numerics {
std::vector<double> discreteInfiniteObservabilityGramian(
    const std::vector<double>& transition,const std::vector<double>& output,
    std::size_t states,std::size_t outputs) {
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    if(!states||!outputs||transition.size()!=states*states||output.size()!=outputs*states)
        throw std::invalid_argument("discrete Gramian dimensions mismatch");
    const Matrix A=Eigen::Map<const Matrix>(transition.data(),states,states);
    const Matrix C=Eigen::Map<const Matrix>(output.data(),outputs,states);
    if(!A.allFinite()||!C.allFinite())throw std::invalid_argument("Gramian inputs must be finite");
    Eigen::EigenSolver<Matrix> eigen(A,false);
    if(eigen.info()!=Eigen::Success)throw std::runtime_error("Gramian stability eigensolve failed");
    if(eigen.eigenvalues().cwiseAbs().maxCoeff()>=1.0)
        throw std::invalid_argument("infinite discrete Gramian requires strictly stable transition");
    // Shared Bartels-Stewart Stein solver: A^T W A - W + C^T C = 0.
    const Matrix W=core::linear_detail::lyapunov(A.transpose(),C.transpose()*C,true);
    return {W.data(),W.data()+W.size()};
}
} // namespace dynsyslab::numerics

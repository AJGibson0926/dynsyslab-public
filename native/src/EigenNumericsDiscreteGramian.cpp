#include <dynsyslab/numerics/MatrixFunctions.hpp>
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
    const Matrix Q=C.transpose()*C;
    Matrix W=Q,power=A;
    for(int iteration=0;iteration<100;++iteration) {
        const Matrix contribution=power.transpose()*W*power;
        W+=contribution;power=(power*power).eval();
        if(!W.allFinite()||!power.allFinite())throw std::runtime_error("discrete Gramian overflow");
        const double scale=std::max(1e-300,Q.norm()+W.norm());
        if(contribution.norm()<1e-13*scale) {
            W=(0.5*(W+W.transpose())).eval();
            const Matrix residual=A.transpose()*W*A-W+Q;
            if(residual.norm()>1e-10*scale)throw std::runtime_error("discrete Gramian residual too large");
            return {W.data(),W.data()+W.size()};
        }
    }
    throw std::runtime_error("discrete Gramian doubling did not converge");
}
} // namespace dynsyslab::numerics

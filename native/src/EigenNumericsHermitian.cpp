#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <stdexcept>
namespace dynsyslab::numerics {
HermitianEigensystem hermitianEigensystem(const std::vector<std::complex<double>>& values,
    std::size_t n) {
    using Matrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    if(!n||values.size()!=n*n)throw std::invalid_argument("Hermitian eigensystem dimensions mismatch");
    const Matrix A=Eigen::Map<const Matrix>(values.data(),n,n);
    if(!A.allFinite()||(A-A.adjoint()).norm()>1e-10*std::max(1.0,A.norm()))
        throw std::invalid_argument("matrix must be finite and Hermitian");
    Eigen::SelfAdjointEigenSolver<Matrix> eigen(A);
    if(eigen.info()!=Eigen::Success)throw std::runtime_error("Hermitian eigensolve failed");
    const Matrix V=eigen.eigenvectors();
    return {{eigen.eigenvalues().data(),eigen.eigenvalues().data()+n},
        {V.data(),V.data()+V.size()}};
}
} // namespace dynsyslab::numerics

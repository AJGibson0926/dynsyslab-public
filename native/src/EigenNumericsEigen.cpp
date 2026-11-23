#include <dynsyslab/numerics/Eigensystem.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <stdexcept>

namespace dynsyslab::numerics {

ComplexSchurDecomposition complexSchur(const std::vector<double>& values,std::size_t n) {
    if(!n||values.size()!=n*n)throw std::invalid_argument("Schur matrix shape mismatch");
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using CM=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const Matrix a=Eigen::Map<const Matrix>(values.data(),n,n);
    if(!a.allFinite())throw std::invalid_argument("Finite Schur matrix required");
    // Eigen ComplexSchur: https://libeigen.gitlab.io/eigen/docs-3.4/classEigen_1_1ComplexSchur.html
    Eigen::ComplexSchur<Matrix> solver(a);
    if(solver.info()!=Eigen::Success)throw std::runtime_error("Complex Schur decomposition did not converge");
    const CM u=solver.matrixU(),t=solver.matrixT();
    return {{u.data(),u.data()+u.size()},{t.data(),t.data()+t.size()},n};
}

RealEigensystem eigensystem(const std::vector<double>& values,
                            std::size_t dimension,bool computeEigenvectors) {
    if(!dimension||values.size()!=dimension*dimension)
        throw std::invalid_argument("eigensystem matrix shape does not match its data");
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const auto source=Eigen::Map<const Matrix>(values.data(),
        static_cast<Eigen::Index>(dimension),static_cast<Eigen::Index>(dimension));
    Eigen::ComplexEigenSolver<Matrix> solver(source,computeEigenvectors);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("real eigensystem did not converge");
    RealEigensystem result;result.dimension=dimension;
    result.eigenvalues.assign(solver.eigenvalues().data(),
        solver.eigenvalues().data()+solver.eigenvalues().size());
    if(computeEigenvectors) {
        using CMatrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,
            Eigen::Dynamic,Eigen::RowMajor>;
        const CMatrix vectors=solver.eigenvectors();
        result.rightEigenvectors.assign(vectors.data(),vectors.data()+vectors.size());
    }
    return result;
}

}  // namespace dynsyslab::numerics

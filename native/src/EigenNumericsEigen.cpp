#include <dynsyslab/numerics/Eigensystem.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <stdexcept>

namespace dynsyslab::numerics {

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

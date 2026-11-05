#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::numerics {
namespace {

using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

Matrix mapped(const std::vector<double>& values,std::size_t rows,
              std::size_t columns,const char* name) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument(std::string(name)+" dimensions are inconsistent");
    Matrix result=Eigen::Map<const Matrix>(values.data(),rows,columns);
    if(!result.allFinite())
        throw std::invalid_argument(std::string(name)+" must be finite");
    return result;
}

std::vector<double> values(const Matrix& matrix) {
    return {matrix.data(),matrix.data()+matrix.size()};
}

Matrix symmetric(Matrix matrix,const char* name) {
    matrix=0.5*(matrix+matrix.transpose()).eval();
    if(!matrix.allFinite())
        throw std::runtime_error(std::string(name)+" is not finite");
    return matrix;
}

}  // namespace

SymmetricEigensystem symmetricEigensystem(
    const std::vector<double>& source,std::size_t dimension) {
    const Matrix matrix=mapped(source,dimension,dimension,"symmetric eigensystem");
    Eigen::SelfAdjointEigenSolver<Matrix> solver(matrix);
    if(solver.info()!=Eigen::Success)
        throw std::runtime_error("symmetric eigensystem did not converge");
    const Matrix vectors=solver.eigenvectors();
    SymmetricEigensystem result;
    result.dimension=dimension;
    result.eigenvalues.assign(
        solver.eigenvalues().data(),solver.eigenvalues().data()+dimension);
    result.eigenvectors=values(vectors);
    return result;
}

std::vector<double> matrixExponential(
    const std::vector<double>& source,std::size_t dimension,double scale) {
    if(!std::isfinite(scale))
        throw std::invalid_argument("matrix-exponential scale must be finite");
    return values((mapped(source,dimension,dimension,"matrix exponential")*scale).exp());
}

ZeroOrderHold zeroOrderHold(
    const std::vector<double>& generator,const std::vector<double>& input,
    std::size_t states,std::size_t inputs,double sampleTime) {
    if(!(sampleTime>0.0)||!std::isfinite(sampleTime))
        throw std::invalid_argument("zero-order-hold sample time must be positive");
    const Matrix A=mapped(generator,states,states,"zero-order-hold generator");
    const Matrix B=mapped(input,states,inputs,"zero-order-hold input");
    Matrix augmented=Matrix::Zero(states+inputs,states+inputs);
    augmented.topLeftCorner(states,states)=A;
    augmented.topRightCorner(states,inputs)=B;
    const Matrix propagated=(augmented*sampleTime).exp();
    const Matrix transition=propagated.topLeftCorner(states,states);
    const Matrix response=propagated.topRightCorner(states,inputs);
    return {values(transition),values(response),states,inputs};
}

bool isDetectable(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs) {
    const Matrix A=mapped(generator,states,states,"detectability generator");
    const Matrix C=mapped(output,outputs,states,"detectability output");
    Eigen::ComplexEigenSolver<Matrix> eigen(A,false);
    if(eigen.info()!=Eigen::Success)
        throw std::runtime_error("detectability eigensolve failed");
    const double tolerance=1.0e-10*std::max(1.0,A.norm());
    for(Eigen::Index pole=0;pole<A.rows();++pole) {
        const auto lambda=eigen.eigenvalues()[pole];
        if(lambda.real()<-tolerance)continue;
        Eigen::MatrixXcd pbh(A.rows()+C.rows(),A.rows());
        pbh.topRows(A.rows())=lambda*Eigen::MatrixXcd::Identity(A.rows(),A.rows())
            -A.cast<std::complex<double>>();
        pbh.bottomRows(C.rows())=C.cast<std::complex<double>>();
        Eigen::FullPivLU<Eigen::MatrixXcd> factor(pbh);
        factor.setThreshold(tolerance);
        if(factor.rank()<A.rows())return false;
    }
    return true;
}

std::vector<double> finiteObservabilityGramian(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs,double horizon) {
    if(!(horizon>0.0)||!std::isfinite(horizon))
        throw std::invalid_argument(
            "continuous Gramian horizon must be positive");
    const Matrix A=mapped(generator,states,states,"observability generator");
    const Matrix C=mapped(output,outputs,states,"observability output");
    const Matrix forcing=C.transpose()*C;
    Matrix augmented=Matrix::Zero(2*states,2*states);
    augmented.topLeftCorner(states,states)=-A.transpose();
    augmented.topRightCorner(states,states)=forcing;
    augmented.bottomRightCorner(states,states)=A;
    const Matrix propagated=(augmented*horizon).exp();
    return values(symmetric((A.transpose()*horizon).exp()
        *propagated.topRightCorner(states,states),
        "continuous observability Gramian"));
}

std::vector<double> infiniteObservabilityGramian(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs) {
    const Matrix A=mapped(generator,states,states,"observability generator");
    const Matrix C=mapped(output,outputs,states,"observability output");
    Eigen::ComplexEigenSolver<Matrix> eigen(A,false);
    if(eigen.info()!=Eigen::Success)
        throw std::runtime_error("infinite-Gramian stability eigensolve failed");
    const double scale=std::max(1.0,A.norm());
    if((eigen.eigenvalues().real().array()>=-1.0e-12*scale).any())
        throw std::invalid_argument(
            "infinite observability Gramian requires a strictly stable generator");
    Eigen::ComplexSchur<Matrix> schur(A);
    if(schur.info()!=Eigen::Success)
        throw std::runtime_error("infinite-Gramian Schur decomposition failed");
    const auto U=schur.matrixU();
    const auto T=schur.matrixT();
    const Eigen::MatrixXcd transformed=U.adjoint()
        *(C.transpose()*C).cast<std::complex<double>>()*U;
    Eigen::MatrixXcd solution=Eigen::MatrixXcd::Zero(A.rows(),A.cols());
    for(Eigen::Index row=0;row<A.rows();++row)
        for(Eigen::Index column=0;column<A.cols();++column) {
            std::complex<double> value=-transformed(row,column);
            for(Eigen::Index index=0;index<row;++index)
                value-=std::conj(T(index,row))*solution(index,column);
            for(Eigen::Index index=0;index<column;++index)
                value-=solution(row,index)*T(index,column);
            const auto denominator=std::conj(T(row,row))+T(column,column);
            if(std::abs(denominator)<=std::numeric_limits<double>::epsilon()*scale)
                throw std::runtime_error(
                    "infinite-Gramian Lyapunov solve is singular");
            solution(row,column)=value/denominator;
        }
    const Eigen::MatrixXcd complexResult=U*solution*U.adjoint();
    if(complexResult.imag().norm()
        >1.0e-8*std::max(1.0,complexResult.real().norm()))
        throw std::runtime_error(
            "infinite-Gramian Lyapunov solve was not real");
    const Matrix result=symmetric(
        complexResult.real(),"infinite observability Gramian");
    const Matrix residual=A.transpose()*result+result*A+C.transpose()*C;
    if(residual.norm()>1.0e-8*std::max(1.0,(C.transpose()*C).norm()))
        throw std::runtime_error(
            "infinite-Gramian Lyapunov residual is too large");
    return values(result);
}

}  // namespace dynsyslab::numerics

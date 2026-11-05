#include <dynsyslab/numerics/Control.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::numerics {
namespace {

using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector=Eigen::VectorXd;
using CMatrix=Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>;

Matrix mapped(const std::vector<double>& values,std::size_t rows,
              std::size_t columns,const char* name) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument(std::string(name)+" dimensions are inconsistent");
    Matrix result=Eigen::Map<const Matrix>(values.data(),rows,columns);
    if(!result.allFinite())
        throw std::invalid_argument(std::string(name)+" must be finite");
    return result;
}

Matrix careResidual(const Matrix& A,const Matrix& gram,const Matrix& Q,
                    const Matrix& P) {
    return A.transpose()*P+P*A-P*gram*P+Q;
}

double residualRatio(const Matrix& A,const Matrix& gram,const Matrix& Q,
                     const Matrix& P) {
    const Matrix residual=careResidual(A,gram,Q,P);
    const double scale=std::max({A.norm()*P.norm(),(P*gram*P).norm(),
                                 Q.norm(),1.0});
    return residual.norm()/scale;
}

bool stable(const Matrix& matrix) {
    Eigen::EigenSolver<Matrix> eigenvalues(matrix,false);
    return eigenvalues.info()==Eigen::Success
        && eigenvalues.eigenvalues().real().maxCoeff()<-1.0e-10;
}

void refineRiccati(const Matrix& A,const Matrix& gram,const Matrix& Q,
                   Matrix& P) {
    const std::size_t states=static_cast<std::size_t>(A.rows());
    const std::size_t unknowns=states*states;
    for(std::size_t iteration=0;iteration<10;++iteration) {
        const double current=residualRatio(A,gram,Q,P);
        if(current<=1.0e-11)return;
        const Matrix closed=A-gram*P;
        if(!stable(closed))return;
        Matrix operation=Matrix::Zero(
            static_cast<Eigen::Index>(unknowns),static_cast<Eigen::Index>(unknowns));
        Vector right(static_cast<Eigen::Index>(unknowns));
        const Matrix residual=careResidual(A,gram,Q,P);
        for(std::size_t row=0;row<states;++row)
            for(std::size_t column=0;column<states;++column) {
                const auto equation=static_cast<Eigen::Index>(row*states+column);
                right[equation]=-residual(static_cast<Eigen::Index>(row),
                                          static_cast<Eigen::Index>(column));
                for(std::size_t inner=0;inner<states;++inner) {
                    operation(equation,static_cast<Eigen::Index>(inner*states+column))
                        +=closed(static_cast<Eigen::Index>(inner),
                                 static_cast<Eigen::Index>(row));
                    operation(equation,static_cast<Eigen::Index>(row*states+inner))
                        +=closed(static_cast<Eigen::Index>(inner),
                                 static_cast<Eigen::Index>(column));
                }
            }
        Eigen::FullPivLU<Matrix> factor(operation);
        if(!factor.isInvertible())return;
        const Vector correctionValues=factor.solve(right);
        Matrix correction(static_cast<Eigen::Index>(states),
                          static_cast<Eigen::Index>(states));
        for(std::size_t row=0;row<states;++row)
            for(std::size_t column=0;column<states;++column)
                correction(static_cast<Eigen::Index>(row),
                           static_cast<Eigen::Index>(column))=
                    correctionValues[static_cast<Eigen::Index>(row*states+column)];
        correction=0.5*(correction+correction.transpose()).eval();
        bool accepted=false;
        for(double step=1.0;step>=1.0/128.0;step*=0.5) {
            Matrix candidate=P+step*correction;
            candidate=0.5*(candidate+candidate.transpose()).eval();
            if(candidate.allFinite()&&stable(A-gram*candidate)
                &&residualRatio(A,gram,Q,candidate)<current) {
                P=std::move(candidate);accepted=true;break;
            }
        }
        if(!accepted)return;
    }
}

}  // namespace

ContinuousLqrResult solveContinuousLqr(
    const std::vector<double>& generator,const std::vector<double>& input,
    const std::vector<double>& statePenalty,
    const std::vector<double>& inputPenalty,
    std::size_t states,std::size_t inputs) {
    const Matrix A=mapped(generator,states,states,"continuous LQR A");
    const Matrix B=mapped(input,states,inputs,"continuous LQR B");
    const Matrix Q=mapped(statePenalty,states,states,"continuous LQR Q");
    const Matrix R=mapped(inputPenalty,inputs,inputs,"continuous LQR R");
    Eigen::LDLT<Matrix> rFactor(R);
    if(rFactor.info()!=Eigen::Success||!rFactor.isPositive())
        throw std::invalid_argument("continuous LQR R must be positive definite");
    const Matrix inverseR=rFactor.solve(Matrix::Identity(inputs,inputs));
    const Matrix gram=B*inverseR*B.transpose();
    const double gramNorm=gram.norm(),qNorm=Q.norm();
    const double pScale=(gramNorm>std::numeric_limits<double>::min()
        &&qNorm>std::numeric_limits<double>::min())
        ?std::sqrt(qNorm/gramNorm):1.0;
    Matrix hamiltonian(2*states,2*states);
    hamiltonian.topLeftCorner(states,states)=A;
    hamiltonian.topRightCorner(states,states)=-pScale*gram;
    hamiltonian.bottomLeftCorner(states,states)=-Q/pScale;
    hamiltonian.bottomRightCorner(states,states)=-A.transpose();
    hamiltonian/=std::max(hamiltonian.cwiseAbs().maxCoeff(),1.0);
    Eigen::ComplexEigenSolver<Matrix> decomposition(hamiltonian);
    if(decomposition.info()!=Eigen::Success)
        throw std::runtime_error("continuous LQR Hamiltonian eigensolve failed");
    std::vector<Eigen::Index> stableIndices;
    for(Eigen::Index index=0;index<decomposition.eigenvalues().size();++index)
        if(decomposition.eigenvalues()[index].real()<-1.0e-10)
            stableIndices.push_back(index);
    if(stableIndices.size()!=states)
        throw std::runtime_error(
            "continuous LQR has no unique stabilizing invariant subspace");
    CMatrix upper(states,states),lower(states,states);
    for(std::size_t column=0;column<states;++column) {
        upper.col(static_cast<Eigen::Index>(column))=decomposition.eigenvectors()
            .col(stableIndices[column]).head(static_cast<Eigen::Index>(states));
        lower.col(static_cast<Eigen::Index>(column))=decomposition.eigenvectors()
            .col(stableIndices[column]).tail(static_cast<Eigen::Index>(states));
    }
    Eigen::FullPivLU<CMatrix> upperFactor(upper);
    if(!upperFactor.isInvertible())
        throw std::runtime_error("continuous LQR stable subspace is singular");
    Eigen::JacobiSVD<CMatrix> upperSvd(upper);
    const double largest=upperSvd.singularValues()[0];
    const double smallest=upperSvd.singularValues()[upperSvd.singularValues().size()-1];
    const double reciprocalCondition=smallest/
        std::max(largest,std::numeric_limits<double>::min());
    const CMatrix complexP=lower*upperFactor.inverse();
    const double imaginaryRatio=complexP.imag().norm()/
        std::max(complexP.real().norm(),std::numeric_limits<double>::min());
    if(imaginaryRatio>1.0e-6)
        throw std::runtime_error("continuous LQR Riccati solution is not real");
    Matrix P=pScale*complexP.real();
    P=0.5*(P+P.transpose()).eval();
    refineRiccati(A,gram,Q,P);
    const Matrix K=inverseR*B.transpose()*P;
    const Matrix residual=careResidual(A,gram,Q,P);
    const double scale=std::max({A.norm()*P.norm(),(P*gram*P).norm(),
                                 Q.norm(),1.0});
    Eigen::EigenSolver<Matrix> closedLoop(A-B*K,false);
    if(closedLoop.info()!=Eigen::Success)
        throw std::runtime_error("closed-loop pole solve failed");
    const double tolerance=std::clamp(
        std::numeric_limits<double>::epsilon()/
            std::max(reciprocalCondition*reciprocalCondition,
                     std::numeric_limits<double>::epsilon()),
        1.0e-6,1.0e-2);
    const double relativeResidual=residual.norm()/scale;
    if(closedLoop.eigenvalues().real().maxCoeff()>=-1.0e-10)
        throw std::runtime_error("continuous LQR solution is not stabilizing");
    if(relativeResidual>tolerance)
        throw std::runtime_error("continuous LQR Riccati residual exceeds tolerance");
    ContinuousLqrResult result;
    result.riccati.assign(P.data(),P.data()+P.size());
    result.gain.assign(K.data(),K.data()+K.size());
    result.closedLoopPoles.assign(closedLoop.eigenvalues().data(),
                                  closedLoop.eigenvalues().data()+states);
    result.relativeResidual=relativeResidual;
    return result;
}

}  // namespace dynsyslab::numerics

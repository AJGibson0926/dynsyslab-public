#include <dynsyslab/numerics/Control.hpp>

#include <Eigen/Dense>

#include <cmath>
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

}  // namespace

FiniteHorizonLqrResult solveFiniteHorizonLqr(
    const std::vector<double>& generator,const std::vector<double>& input,
    const std::vector<double>& statePenalty,
    const std::vector<double>& inputPenalty,
    std::size_t states,std::size_t inputs,std::size_t horizon,
    double sampleTime) {
    const Matrix A=mapped(generator,states,states,"finite-horizon LQR A");
    const Matrix B=mapped(input,states,inputs,"finite-horizon LQR B");
    const Matrix Q=mapped(statePenalty,states,states,"finite-horizon LQR Q");
    const Matrix R=mapped(inputPenalty,inputs,inputs,"finite-horizon LQR R");
    if(!horizon||!std::isfinite(sampleTime)||sampleTime<=0.0)
        throw std::invalid_argument(
            "finite-horizon LQR horizon and sample time must be positive");
    const Matrix Ad=Matrix::Identity(states,states)+sampleTime*A;
    const Matrix Bd=sampleTime*B;
    Matrix P=Q,gain(inputs,states);
    for(std::size_t step=0;step<horizon;++step) {
        const Matrix normal=R+Bd.transpose()*P*Bd;
        Eigen::LDLT<Matrix> factor(normal);
        if(factor.info()!=Eigen::Success||!factor.isPositive())
            throw std::runtime_error(
                "finite-horizon LQR normal matrix is not positive definite");
        gain=factor.solve(Bd.transpose()*P*Ad);
        P=Q+Ad.transpose()*P*(Ad-Bd*gain);
    }
    return {{gain.data(),gain.data()+gain.size()},
            {P.data(),P.data()+P.size()}};
}

}  // namespace dynsyslab::numerics

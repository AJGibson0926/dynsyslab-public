#include <dynsyslab/numerics/StreamingQr.hpp>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
// Row-block QR reduction: Demmel, Grigori, Hoemmen & Langou (2012),
// https://doi.org/10.1137/080731992. Retain only R, not the Q tree. This is
// a simple sequential implementation, not the paper's optimized TSQR kernel.
namespace dynsyslab::numerics {
WeightedQrAccumulator::WeightedQrAccumulator(std::size_t columns):columns_(columns){
    if(!columns||columns>50000000/columns)
        throw std::invalid_argument("QR factor exceeds fifty million entries");
}
void WeightedQrAccumulator::append(const std::vector<double>& values,std::size_t rows,
    const std::vector<double>& weights){
    if(!rows||rows>50000000/columns_-rows_||values.size()!=rows*columns_
        ||(!weights.empty()&&weights.size()!=rows)||rows>std::numeric_limits<std::size_t>::max()-samples_)
        throw std::invalid_argument("QR block shape or workspace limit exceeded");
    double addedWeight=0;
    for(std::size_t i=0;i<rows;++i){const double w=weights.empty()?1:weights[i];
        if(!std::isfinite(w)||w<0)throw std::invalid_argument("QR weights must be finite and nonnegative");
        addedWeight+=w;}
    if(!std::isfinite(totalWeight_+addedWeight))throw std::invalid_argument("QR total weight overflow");
    for(double v:values)if(!std::isfinite(v))throw std::invalid_argument("Finite QR entries required");
    // Column-major is used internally for Eigen's Householder kernels.
    Eigen::MatrixXd stacked(rows_+rows,columns_);
    for(std::size_t i=0;i<rows_;++i)for(std::size_t j=0;j<columns_;++j)
        stacked(i,j)=factor_[i*columns_+j];
    for(std::size_t i=0;i<rows;++i){const double scale=std::sqrt(weights.empty()?1:weights[i]);
        for(std::size_t j=0;j<columns_;++j)stacked(rows_+i,j)=scale*values[i*columns_+j];}
    if(!stacked.allFinite())throw std::invalid_argument("Weighted QR entries overflow");
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(stacked);
    const auto kept=std::min(rows_+rows,columns_);
    std::vector<double> result(kept*columns_);
    for(std::size_t i=0;i<kept;++i)for(std::size_t j=i;j<columns_;++j)result[i*columns_+j]=qr.matrixQR()(i,j);
    for(double v:result)if(!std::isfinite(v))throw std::runtime_error("QR factor overflow");
    factor_=std::move(result);rows_=kept;samples_+=rows;totalWeight_+=addedWeight;
}
}

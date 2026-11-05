#include <dynsyslab/core/DmdControl.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include "FlatMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
using R = detail::FlatMatrix<double>;
using C = detail::FlatMatrix<std::complex<double>>;

template<class T> detail::FlatMatrix<T> transpose(const detail::FlatMatrix<T>& a) {
    detail::FlatMatrix<T> out(a.columns,a.rows);
    for(std::size_t i=0;i<a.rows;++i)
        for(std::size_t j=0;j<a.columns;++j) out(j,i)=a(i,j);
    return out;
}
template<class T> detail::FlatMatrix<T> product(
    const detail::FlatMatrix<T>& a,const detail::FlatMatrix<T>& b) {
    return {numerics::multiply(a.values,a.rows,a.columns,b.values,b.columns),a.rows,b.columns};
}
C complex(const R& a) { return {{a.values.begin(),a.values.end()},a.rows,a.columns}; }

void validate(const std::vector<double>& a,std::size_t rows,std::size_t cols) {
    if(!rows||!cols||rows>std::numeric_limits<std::size_t>::max()/cols
        ||a.size()!=rows*cols
        ||!std::all_of(a.begin(),a.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("DMDc requires finite, correctly shaped snapshot/input matrices");
}
std::size_t numericalRank(const numerics::RealSvd& svd,double tolerance) {
    if(svd.singularValues.empty()) return 0;
    const double cutoff=tolerance>0?tolerance:std::numeric_limits<double>::epsilon()
        *static_cast<double>(std::max(svd.rows,svd.columns))*svd.singularValues.front();
    return static_cast<std::size_t>(std::count_if(svd.singularValues.begin(),
        svd.singularValues.end(),[&](double s){return s>cutoff;}));
}
std::size_t retained(std::size_t numerical,std::size_t requested) {
    const auto rank=requested?std::min(numerical,requested):numerical;
    if(!rank) throw std::invalid_argument("DMDc snapshot matrix has numerical rank zero");
    return rank;
}
double relativeError(const R& prediction,const R& reference) {
    double error=0;
    for(std::size_t i=0;i<reference.values.size();++i)
        error+=std::norm(std::complex<double>(prediction.values[i]-reference.values[i]));
    return std::sqrt(error/std::max(detail::squaredNorm(reference),
        std::numeric_limits<double>::min()));
}

void spectrum(DmdControlResult& out,const R& basis,const R& stateAction) {
    const auto r=out.rank;
    const auto eig=numerics::eigensystem(out.operatorMatrix,r);
    out.discreteEigenvalues=eig.eigenvalues;
    const C w(eig.rightEigenvectors,r,r);
    out.modes=product(complex(basis),w).values;
    out.exactModes=product(complex(stateAction),w).values;
    const auto singular=numerics::singularValues(w.values,r,r);
    out.eigenvectorConditionNumber=singular.back()>0
        ?singular.front()/singular.back():std::numeric_limits<double>::max();
    // A defective/nearly defective eigensystem is not a usable modal chart.
    // Keep A/B available, but never advertise a pseudoinverse as an eigenbasis.
    if(out.eigenvectorConditionNumber>1e6) return;
    C identity(r,r);
    for(std::size_t i=0;i<r;++i) identity(i,i)=1;
    numerics::LeastSquaresOptions opts;
    opts.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    const C left(numerics::solveLeastSquares(w.values,r,r,identity.values,r,opts),r,r);
    out.eigenfunctionCoefficients=product(left,complex(transpose(basis))).values;
    out.modalInputMatrix=product(left,complex(R(out.controlMatrix,r,out.inputDimension))).values;
    out.hasEigenfunctions=true;
}
} // namespace

DmdControlResult dynamicModeDecompositionWithControl(
    const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<double>& inputs,std::size_t pairs,std::size_t dimension,
    std::size_t inputDimension,double deltaTime,const DmdControlOptions& options) {
    validate(x,pairs,dimension); validate(y,pairs,dimension);
    validate(inputs,pairs,inputDimension);
    if(pairs<2||!std::isfinite(deltaTime)||deltaTime<=0
        ||!std::isfinite(options.singularValueTolerance)||options.singularValueTolerance<0)
        throw std::invalid_argument("DMDc requires >=2 pairs, positive finite dt and nonnegative SVD tolerance");
    if(std::all_of(inputs.begin(),inputs.end(),[](double v){return v==0;}))
        throw std::invalid_argument("DMDc cannot identify B from unexcited (all-zero) inputs");
    const R xt=transpose(R(x,pairs,dimension)),yt=transpose(R(y,pairs,dimension));
    R omega(dimension+inputDimension,pairs);
    for(std::size_t i=0;i<dimension;++i)
        for(std::size_t k=0;k<pairs;++k) omega(i,k)=xt(i,k);
    for(std::size_t j=0;j<inputDimension;++j)
        for(std::size_t k=0;k<pairs;++k) omega(dimension+j,k)=inputs[k*inputDimension+j];
    const auto aug=numerics::singularValueDecomposition(omega.values,omega.rows,pairs);
    const auto output=numerics::singularValueDecomposition(yt.values,dimension,pairs);
    DmdControlResult out;
    out.dimension=dimension; out.inputDimension=inputDimension;
    out.pairs=pairs; out.deltaTime=deltaTime;
    out.augmentedNumericalRank=numericalRank(aug,options.singularValueTolerance);
    out.augmentedRank=retained(out.augmentedNumericalRank,options.augmentedRank);
    out.rank=retained(numericalRank(output,options.singularValueTolerance),options.rank);
    out.fullRowRank=out.augmentedNumericalRank==omega.rows;
    out.augmentedSingularValues=aug.singularValues;
    out.outputSingularValues=output.singularValues;
    out.retainedConditionNumber=aug.singularValues.front()/aug.singularValues[out.augmentedRank-1];
    const auto p=out.augmentedRank,r=out.rank;
    R basis(dimension,r),vScaled(pairs,p),u1t(p,dimension),u2t(p,inputDimension);
    for(std::size_t i=0;i<dimension;++i)
        for(std::size_t j=0;j<r;++j) basis(i,j)=output.left[i*output.rank+j];
    for(std::size_t j=0;j<p;++j) {
        for(std::size_t k=0;k<pairs;++k)
            vScaled(k,j)=aug.rightTranspose[j*pairs+k]/aug.singularValues[j];
        for(std::size_t i=0;i<dimension;++i) u1t(j,i)=aug.left[i*aug.rank+j];
        for(std::size_t i=0;i<inputDimension;++i) u2t(j,i)=aug.left[(dimension+i)*aug.rank+j];
    }
    const R factor=product(yt,vScaled);
    const R projected=product(transpose(basis),factor);
    const R stateAction=product(factor,product(u1t,basis));
    const R a=product(transpose(basis),stateAction),b=product(projected,u2t);
    out.basis=basis.values; out.operatorMatrix=a.values; out.controlMatrix=b.values;
    R prediction=product(a,product(transpose(basis),xt));
    const R driven=product(b,transpose(R(inputs,pairs,inputDimension)));
    for(std::size_t i=0;i<prediction.values.size();++i) prediction.values[i]+=driven.values[i];
    out.oneStepRelativeError=relativeError(product(basis,prediction),yt);
    out.projectionRelativeError=relativeError(product(basis,product(transpose(basis),yt)),yt);
    spectrum(out,basis,stateAction);
    return out;
}

std::vector<std::complex<double>> evaluateDmdControlEigenfunctions(
    const DmdControlResult& model,const std::vector<double>& snapshots,std::size_t samples) {
    validate(snapshots,samples,model.dimension);
    if(!model.hasEigenfunctions||!model.rank
        ||model.eigenfunctionCoefficients.size()!=model.rank*model.dimension)
        throw std::invalid_argument("DMDc model has no usable linear eigenfunction basis");
    const C coefficients(model.eigenfunctionCoefficients,model.rank,model.dimension);
    return product(complex(R(snapshots,samples,model.dimension)),transpose(coefficients)).values;
}

std::vector<double> predictDmdControl(const DmdControlResult& model,
    const std::vector<double>& snapshots,const std::vector<double>& inputs,std::size_t samples) {
    validate(snapshots,samples,model.dimension); validate(inputs,samples,model.inputDimension);
    validate(model.basis,model.dimension,model.rank);
    validate(model.operatorMatrix,model.rank,model.rank);
    validate(model.controlMatrix,model.rank,model.inputDimension);
    const R basis(model.basis,model.dimension,model.rank);
    auto next=product(product(R(snapshots,samples,model.dimension),basis),
        transpose(R(model.operatorMatrix,model.rank,model.rank)));
    const auto forced=product(R(inputs,samples,model.inputDimension),
        transpose(R(model.controlMatrix,model.rank,model.inputDimension)));
    for(std::size_t i=0;i<next.values.size();++i) next.values[i]+=forced.values[i];
    return product(next,transpose(basis)).values;
}
} // namespace dynsyslab::core

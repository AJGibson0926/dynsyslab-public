#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include "FlatMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace dynsyslab::core {
namespace {

using RMatrix=detail::FlatMatrix<double>;
using CMatrix=detail::FlatMatrix<std::complex<double>>;

RMatrix matrix(const std::vector<double>& values,std::size_t rows,
               std::size_t columns) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument("matrix shape does not match its data");
    return {values,rows,columns};
}

std::size_t retainedRank(const std::vector<double>& singular,std::size_t rows,
                         std::size_t columns,std::size_t requested,
                         double tolerance) {
    if(singular.empty())return 0;
    const double threshold=tolerance>0.0?tolerance:
        std::numeric_limits<double>::epsilon()
            *static_cast<double>(std::max(rows,columns))*singular.front();
    std::size_t numerical=0;
    while(numerical<singular.size()&&singular[numerical]>threshold)++numerical;
    return requested?std::min(requested,numerical):numerical;
}

CMatrix solve(const CMatrix& left,const CMatrix& right) {
    numerics::LeastSquaresOptions options;
    options.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    return {numerics::solveLeastSquares(
        left.values,left.rows,left.columns,right.values,right.columns,options),
        left.columns,right.columns};
}

CMatrix complexValues(const RMatrix& source) {
    CMatrix result(source.rows,source.columns);
    for(std::size_t index=0;index<source.values.size();++index)
        result.values[index]=source.values[index];
    return result;
}

CMatrix transpose(const CMatrix& source) {
    CMatrix result(source.columns,source.rows);
    for(std::size_t row=0;row<source.rows;++row)
        for(std::size_t column=0;column<source.columns;++column)
            result(column,row)=source(row,column);
    return result;
}

RMatrix transpose(const RMatrix& source) {
    RMatrix result(source.columns,source.rows);
    for(std::size_t row=0;row<source.rows;++row)
        for(std::size_t column=0;column<source.columns;++column)
            result(column,row)=source(row,column);
    return result;
}

RMatrix product(const RMatrix& left,const RMatrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

CMatrix product(const CMatrix& left,const CMatrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

DMDResult dmdPair(const RMatrix& x,const RMatrix& y,double deltaTime,
                  const DMDOptions& options,bool coordinates) {
    if(x.rows<2||x.columns<2||x.rows!=y.rows||x.columns!=y.columns)
        throw std::invalid_argument(
            "DMD snapshot matrices require equal dimensions of at least 2x2");
    const auto svd=numerics::singularValueDecomposition(
        x.values,x.rows,x.columns);
    const std::size_t rank=retainedRank(
        svd.singularValues,x.rows,x.columns,options.rank,
        options.singularValueTolerance);
    if(!rank)throw std::runtime_error("DMD snapshot matrix has numerical rank zero");
    RMatrix left(x.rows,rank),rightTranspose(rank,x.columns);
    for(std::size_t row=0;row<x.rows;++row)
        for(std::size_t mode=0;mode<rank;++mode)
            left(row,mode)=svd.left[row*svd.rank+mode];
    for(std::size_t mode=0;mode<rank;++mode)
        for(std::size_t column=0;column<x.columns;++column)
            rightTranspose(mode,column)=
                svd.rightTranspose[mode*svd.columns+column];
    const RMatrix projected=product(transpose(left),y);
    RMatrix rightScaled(x.columns,rank);
    for(std::size_t snapshot=0;snapshot<x.columns;++snapshot)
        for(std::size_t mode=0;mode<rank;++mode)
            rightScaled(snapshot,mode)=rightTranspose(mode,snapshot)
                /svd.singularValues[mode];
    const RMatrix reduced=product(projected,rightScaled);
    const auto eigen=numerics::eigensystem(reduced.values,rank);
    const CMatrix eigenvectors{eigen.rightEigenvectors,rank,rank};
    const RMatrix basis=options.exactModes?product(y,rightScaled):std::move(left);
    const CMatrix modes=product(complexValues(basis),eigenvectors);
    DMDResult result;
    result.featureDimension=x.rows;
    result.samples=x.columns+1;
    result.rank=rank;
    result.deltaTime=deltaTime;
    result.discreteEigenvalues=eigen.eigenvalues;
    result.continuousEigenvalues.reserve(rank);
    for(const auto value:result.discreteEigenvalues)
        result.continuousEigenvalues.push_back(std::log(value)/deltaTime);
    result.modes=modes.values;
    result.operatorMatrix=complexValues(reduced).values;
    if(coordinates)result.coordinates=transpose(solve(modes,complexValues(x))).values;
    return result;
}

void reorder(HankelDMDResult& result,const std::vector<std::size_t>& order) {
    auto vector=[&](auto& source) {
        using Value=typename std::decay_t<decltype(source)>::value_type;
        std::vector<Value> copy(source.size());
        for(std::size_t index=0;index<order.size();++index)
            copy[index]=source[order[index]];
        source=std::move(copy);
    };
    vector(result.discreteEigenvalues);
    vector(result.continuousEigenvalues);
    vector(result.modalAmplitudes);
    const auto modes=result.modes;
    for(std::size_t row=0;row<result.featureDimension;++row)
        for(std::size_t column=0;column<result.rank;++column)
            result.modes[row*result.rank+column]=
                modes[row*result.rank+order[column]];
    result.sourceModeIndices=order;
}

HankelDMDResult hankelSignal(const std::vector<double>& signal,double deltaTime,
                             const HankelDMDOptions& options) {
    if(options.matrixSize<2||options.samplingWindow<3||!options.sampleStride
        ||!(deltaTime>0.0))
        throw std::invalid_argument("invalid Hankel-DMD settings");
    std::vector<double> sampled;
    for(std::size_t sample=0;sample<signal.size();sample+=options.sampleStride)
        sampled.push_back(signal[sample]);
    const std::size_t rows=std::min(options.samplingWindow,
        sampled.size()>options.matrixSize?sampled.size()-options.matrixSize:0);
    if(rows<2)
        throw std::invalid_argument(
            "Hankel-DMD settings leave fewer than two shifted rows");
    RMatrix x(rows,options.matrixSize),y(rows,options.matrixSize);
    for(std::size_t row=0;row<rows;++row)
        for(std::size_t column=0;column<options.matrixSize;++column) {
            x(row,column)=sampled[row+column];
            y(row,column)=sampled[row+column+1];
        }
    DMDOptions dmdOptions=options;
    dmdOptions.exactModes=false;
    const auto base=dmdPair(
        x,y,deltaTime*static_cast<double>(options.sampleStride),dmdOptions,false);
    HankelDMDResult result;
    static_cast<DMDResult&>(result)=base;
    result.matrixRows=rows;
    result.matrixColumns=options.matrixSize;
    result.sourceSamples=sampled.size();
    CMatrix modes{result.modes,rows,result.rank};
    CMatrix terminal(rows,1);
    for(std::size_t row=0;row<rows;++row)
        terminal(row,0)=x(row,x.columns-1);
    result.modalAmplitudes=solve(modes,terminal).values;
    result.sourceModeIndices.resize(result.rank);
    std::iota(result.sourceModeIndices.begin(),result.sourceModeIndices.end(),0);
    if(options.terminalAmplitudeEnergySort) {
        std::stable_sort(result.sourceModeIndices.begin(),
            result.sourceModeIndices.end(),[&](std::size_t left,std::size_t right) {
                return std::abs(result.modalAmplitudes[left])
                    >std::abs(result.modalAmplitudes[right]);
            });
        const auto order=result.sourceModeIndices;
        reorder(result,order);
    }
    modes={result.modes,rows,result.rank};
    RMatrix all(rows,x.columns+1);
    for(std::size_t row=0;row<rows;++row) {
        for(std::size_t column=0;column<x.columns;++column)
            all(row,column)=x(row,column);
        all(row,x.columns)=y(row,y.columns-1);
    }
    result.coordinates=transpose(solve(modes,complexValues(all))).values;
    result.samples=all.columns;
    result.reconstructionMatrix.resize(result.rank);
    for(std::size_t mode=0;mode<result.rank;++mode)
        result.reconstructionMatrix[mode]=modes(0,mode);
    return result;
}

}  // namespace

SVDResult singularValueDecomposition(
    const std::vector<double>& values,std::size_t rows,std::size_t columns,
    std::size_t rank) {
    const auto decomposition=numerics::singularValueDecomposition(
        values,rows,columns,rank);
    return {decomposition.left,decomposition.singularValues,
        decomposition.rightTranspose,decomposition.rows,decomposition.columns,
        decomposition.rank};
}

DMDResult dynamicModeDecomposition(
    const std::vector<double>& states,std::size_t samples,
    std::size_t stateDimension,double deltaTime,const DMDOptions& options) {
    if(samples<3||!stateDimension||!(deltaTime>0.0)
        ||states.size()!=samples*stateDimension)
        throw std::invalid_argument(
            "DMD requires at least three samples and positive delta time");
    RMatrix x(stateDimension,samples-1),y(stateDimension,samples-1);
    for(std::size_t sample=0;sample+1<samples;++sample)
        for(std::size_t state=0;state<stateDimension;++state) {
            x(state,sample)=states[sample*stateDimension+state];
            y(state,sample)=states[(sample+1)*stateDimension+state];
        }
    auto result=dmdPair(x,y,deltaTime,options,false);
    const CMatrix modes{result.modes,stateDimension,result.rank};
    CMatrix stateColumns(stateDimension,samples);
    for(std::size_t sample=0;sample<samples;++sample)
        for(std::size_t state=0;state<stateDimension;++state)
            stateColumns(state,sample)=states[sample*stateDimension+state];
    result.coordinates=transpose(solve(modes,stateColumns)).values;
    result.samples=samples;
    return result;
}

HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& states,std::size_t samples,
    std::size_t stateDimension,double deltaTime,const HankelDMDOptions& options) {
    if(options.observableIndex>=stateDimension
        ||states.size()!=samples*stateDimension)
        throw std::invalid_argument("Hankel-DMD observable index is outside the state");
    std::vector<double> signal(samples);
    for(std::size_t sample=0;sample<samples;++sample)
        signal[sample]=states[sample*stateDimension+options.observableIndex];
    auto result=hankelSignal(signal,deltaTime,options);
    result.observableName="state_"+std::to_string(options.observableIndex);
    return result;
}

HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& states,std::size_t samples,
    std::size_t stateDimension,double deltaTime,const Observable& observable,
    const HankelDMDOptions& options) {
    if(observable.stateDimension()!=stateDimension
        ||observable.outputDimension()!=1)
        throw std::invalid_argument(
            "Hankel-DMD requires one scalar state-aligned observable");
    auto result=hankelSignal(
        observable.evaluateRealBatch(states,samples),deltaTime,options);
    result.observableName=observable.descriptor().name;
    return result;
}

}  // namespace dynsyslab::core

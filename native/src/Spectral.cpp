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
#include <set>
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

std::size_t optimalHardThresholdRank(std::vector<double> singular,std::size_t rows,std::size_t columns) {
    if(singular.empty())return 0;
    const double beta=double(std::min(rows,columns))/double(std::max(rows,columns));
    const double omega=0.56*beta*beta*beta-0.95*beta*beta+1.82*beta+1.43;
    auto sorted=singular;const auto mid=sorted.begin()+std::ptrdiff_t(sorted.size()/2);
    std::nth_element(sorted.begin(),mid,sorted.end());
    const double threshold=omega*(*mid);
    std::size_t k=0;while(k<singular.size()&&singular[k]>threshold)++k;
    return std::max<std::size_t>(k,1);
}

std::size_t retainedRank(const std::vector<double>& singular,std::size_t rows,
                         std::size_t columns,std::size_t requested,
                         double tolerance) {
    if(!std::isfinite(tolerance)||tolerance<0)
        throw std::invalid_argument("DMD singular value tolerance must be finite and nonnegative");
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
                  const DMDOptions& options,bool coordinates,
                  std::vector<double>* singularValues=nullptr) {
    if(x.rows<2||x.columns<2||x.rows!=y.rows||x.columns!=y.columns)
        throw std::invalid_argument(
            "DMD snapshot matrices require equal dimensions of at least 2x2");
    const auto svd=numerics::singularValueDecomposition(
        x.values,x.rows,x.columns);
    if(singularValues)*singularValues=svd.singularValues;
    std::size_t rank=retainedRank(
        svd.singularValues,x.rows,x.columns,options.rank,
        options.singularValueTolerance);
    if(options.optimalHardThreshold)
        rank=std::min(rank,optimalHardThresholdRank(svd.singularValues,x.rows,x.columns));
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
    if(coordinates) {
        // For projected DMD, U is orthonormal: (U W)^+ H = W^+ U^T H.
        // Solve in rank dimensions instead of factoring the tall complex modes.
        if(options.exactModes)result.coordinates=transpose(solve(modes,complexValues(x))).values;
        else {
            const auto projectedX=product(transpose(basis),x);
            CMatrix all(rank,x.columns+1);
            for(std::size_t row=0;row<rank;++row) {
                for(std::size_t col=0;col<x.columns;++col)all(row,col)=projectedX(row,col);
                all(row,x.columns)=projected(row,y.columns-1);
            }
            result.coordinates=transpose(solve(eigenvectors,all)).values;
        }
    }
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
    const auto coordinates=result.coordinates;
    for(std::size_t row=0;row<result.samples;++row)
        for(std::size_t col=0;col<result.rank;++col)
            result.coordinates[row*result.rank+col]=coordinates[row*result.rank+order[col]];
    result.sourceModeIndices=order;
}

HankelDMDResult hankelSignals(const std::vector<double>& signal,std::size_t samples,
    const std::vector<std::string>& names,double deltaTime,const HankelDMDOptions& options) {
    const std::size_t channels=names.size();
    if(!channels||signal.size()!=samples*channels||options.matrixSize<2
        ||options.samplingWindow<2||!options.sampleStride||!options.snapshotStride
        ||!(deltaTime>0.0)||!std::isfinite(deltaTime))
        throw std::invalid_argument("invalid multichannel Hankel-DMD settings");
    std::set<std::string> unique;
    for(const auto& name:names)
        if(name.empty()||!unique.insert(name).second)
            throw std::invalid_argument("Hankel channel names must be nonempty and unique");
    auto scales=options.channelScales;
    if(scales.empty())scales.assign(channels,1.0);
    if(scales.size()!=channels)
        throw std::invalid_argument("Hankel channel scales must match channels");
    for(double scale:scales)if(!(scale>0.0)||!std::isfinite(scale))
        throw std::invalid_argument("Hankel channel scales must be positive and finite");
    for(double value:signal)if(!std::isfinite(value))
        throw std::invalid_argument("Hankel measurements must be finite");
    const auto window=fixedPairHankelWindow(samples,0,options.samplingWindow,
        options.matrixSize,options.sampleStride,options.snapshotStride);
    std::vector<HankelChannel> selected;
    for (std::size_t c=0;c<channels;++c) selected.push_back({c,0,scales[c]});
    auto h=buildHankelMatrix(signal,samples,channels,selected,window,true);
    return fitHankelDMD(std::move(h),names,deltaTime,options);
}

} // namespace

HankelDMDResult fitHankelDMD(HankelMatrix h,const std::vector<std::string>& names,
    double deltaTime,const HankelDMDOptions& options) {
    const auto channels=h.channels.size(), delays=h.window.delays, columns=h.columns;
    if (!channels || names.size()!=channels || std::set<std::string>(names.begin(),names.end()).size()!=channels
        || std::any_of(names.begin(),names.end(),[](const auto& n){return n.empty();})
        || !(deltaTime>0.) || !std::isfinite(deltaTime) || columns<2 || !delays
        || h.rows/channels!=delays || h.rows%channels || !h.window.sampleStride || !h.sourceSamples
        || h.starts.size()!=columns || h.ends.size()!=columns || h.x.size()/columns!=h.rows
        || h.x.size()%columns || h.y.size()!=h.x.size())
        throw std::invalid_argument("invalid Hankel DMD pair or channel names");
    RMatrix x(std::move(h.x),h.rows,columns), y(std::move(h.y),h.rows,columns);
    DMDOptions dmdOptions=options;
    dmdOptions.exactModes=false;
    HankelDMDResult result;
    static_cast<DMDResult&>(result)=dmdPair(x,y,
        deltaTime*static_cast<double>(h.window.sampleStride),dmdOptions,true,
        &result.singularValues);
    result.matrixRows=h.rows;result.matrixColumns=columns;
    result.delayCount=delays;result.referenceDelay=options.alignToWindowEnd?delays-1:0;
    result.sourceSamples=1+(h.sourceSamples-1)/h.window.sampleStride;result.observableNames=names;
    result.observableName=channels==1?names.front():"joint-observables";
    result.modalAmplitudes.assign(result.coordinates.begin()+(columns-1)*result.rank,
        result.coordinates.begin()+columns*result.rank);
    result.sourceModeIndices.resize(result.rank);
    std::iota(result.sourceModeIndices.begin(),result.sourceModeIndices.end(),0);
    if(options.terminalAmplitudeEnergySort) {
        std::stable_sort(result.sourceModeIndices.begin(),result.sourceModeIndices.end(),
            [&](std::size_t left,std::size_t right) {
                return std::abs(result.modalAmplitudes[left])>std::abs(result.modalAmplitudes[right]);
            });
        const auto order=result.sourceModeIndices;reorder(result,order);
    }
    result.samples=columns+1;
    result.coordinateSampleOffsets=options.alignToWindowEnd?std::move(h.ends):std::move(h.starts);
    result.coordinateSampleOffsets.push_back(result.coordinateSampleOffsets.back()+h.window.sampleStride);
    result.reconstructionMatrix.resize(channels*result.rank);
    for(std::size_t channel=0;channel<channels;++channel)
        for(std::size_t mode=0;mode<result.rank;++mode)
            result.reconstructionMatrix[channel*result.rank+mode]=h.channels[channel].scale
                *result.modes[(result.referenceDelay*channels+channel)*result.rank+mode];
    return result;
}


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
    return hankelSignals(signal,samples,{"state_"+std::to_string(options.observableIndex)},deltaTime,options);
}

HankelDMDResult hankelDynamicModeDecomposition(
    const std::vector<double>& states,std::size_t samples,
    std::size_t stateDimension,double deltaTime,const Observable& observable,
    const HankelDMDOptions& options) {
    if(observable.stateDimension()!=stateDimension
        ||observable.outputDimension()!=1)
        throw std::invalid_argument(
            "Hankel-DMD requires one scalar state-aligned observable");
    return hankelSignals(observable.evaluateRealBatch(states,samples),samples,
        {observable.descriptor().name},deltaTime,options);
}

HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>& values,
    std::size_t samples,const std::vector<std::string>& names,double dt,
    const HankelDMDOptions& options) {
    return hankelSignals(values,samples,names,dt,options);
}

HankelDMDResult hankelDynamicModeDecomposition(const std::vector<double>& states,
    std::size_t samples,std::size_t dimension,double dt,
    const std::vector<Observable>& observables,const HankelDMDOptions& options) {
    std::vector<std::string> names;
    std::vector<double> values(samples*observables.size());
    for(std::size_t channel=0;channel<observables.size();++channel) {
        const auto& observable=observables[channel];
        if(observable.stateDimension()!=dimension||observable.outputDimension()!=1)
            throw std::invalid_argument("Joint Hankel-DMD requires scalar state-aligned observables");
        names.push_back(observable.descriptor().name);
        const auto signal=observable.evaluateRealBatch(states,samples);
        for(std::size_t sample=0;sample<samples;++sample)
            values[sample*observables.size()+channel]=signal[sample];
    }
    return hankelSignals(values,samples,names,dt,options);
}

}  // namespace dynsyslab::core

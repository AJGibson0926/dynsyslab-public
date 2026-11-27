#include <dynsyslab/core/EdmdControl.hpp>
#include <dynsyslab/core/SINDy.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
void validate(const std::vector<double>& values,std::size_t rows,std::size_t columns) {
    if(!rows||!columns||rows>std::numeric_limits<std::size_t>::max()/columns
        ||values.size()!=rows*columns
        ||!std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("EDMDc requires finite, correctly shaped snapshots");
}
void checkSize(std::size_t dimension,std::size_t degree) {
    // Bound combinatorial growth before asking the shared SINDy enumerator to
    // allocate. This is a small dense spectral model, not a raw-window lift.
    if(!dimension||!degree||degree>32)
        throw std::invalid_argument("EDMDc requires positive dimension and polynomial degree 1..32");
    long double features=1;
    for(std::size_t j=1;j<=degree;++j) {
        features*= (static_cast<long double>(dimension)+j)/j;
        if(features>4096)
            throw std::invalid_argument("EDMDc polynomial dictionary exceeds 4096 features; reduce coordinates or degree");
    }
}
void validateDictionary(const EdmdControlDictionary& d) {
    checkSize(d.dimension,d.degree);
    validate(d.center,1,d.dimension); validate(d.scale,1,d.dimension);
    if(std::any_of(d.scale.begin(),d.scale.end(),[](double s){return s<=0;})
        ||d.exponents!=polynomialExponentTable(d.dimension,d.degree,d.includeConstant))
        throw std::invalid_argument("EDMDc requires positive scales and the canonical SINDy monomial table");
}
} // namespace

std::vector<double> liftEdmdControlSnapshots(const EdmdControlDictionary& d,
    const std::vector<double>& snapshots,std::size_t samples) {
    validateDictionary(d); validate(snapshots,samples,d.dimension);
    const auto features=d.exponents.size()/d.dimension;
    if(samples>std::numeric_limits<std::size_t>::max()/features)
        throw std::invalid_argument("EDMDc lifted shape overflow");
    std::vector<double> result(samples*features),q(d.dimension);
    for(std::size_t k=0;k<samples;++k) {
        for(std::size_t j=0;j<d.dimension;++j) {
            q[j]=(snapshots[k*d.dimension+j]-d.center[j])/d.scale[j];
            if(!std::isfinite(q[j])) throw std::invalid_argument("EDMDc normalization overflow");
        }
        for(std::size_t f=0;f<features;++f) {
            double value=1;
            for(std::size_t j=0;j<d.dimension;++j)
                for(std::size_t p=0;p<d.exponents[f*d.dimension+j];++p) value*=q[j];
            if(!std::isfinite(value)) throw std::invalid_argument("EDMDc polynomial feature overflow");
            result[k*features+f]=value;
        }
    }
    return result;
}

EdmdControlResult extendedDynamicModeDecompositionWithControl(
    const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<double>& inputs,std::size_t pairs,std::size_t dimension,
    std::size_t inputDimension,double deltaTime,const EdmdControlOptions& options) {
    validate(x,pairs,dimension); validate(y,pairs,dimension);
    checkSize(dimension,options.polynomialDegree);
    EdmdControlResult out;
    auto& d=out.dictionary;
    d.dimension=dimension; d.degree=options.polynomialDegree;
    d.includeConstant=options.includeConstant;
    d.exponents=polynomialExponentTable(dimension,d.degree,d.includeConstant);
    d.center.assign(dimension,0); d.scale.assign(dimension,1);
    if(options.normalize) for(std::size_t j=0;j<dimension;++j) {
        long double sum=0,squared=0;
        for(std::size_t k=0;k<pairs;++k) sum+=x[k*dimension+j];
        d.center[j]=static_cast<double>(sum/pairs);
        for(std::size_t k=0;k<pairs;++k) {
            const long double value=static_cast<long double>(x[k*dimension+j])-d.center[j];
            squared+=value*value;
        }
        d.scale[j]=static_cast<double>(std::sqrt(squared/pairs));
        if(d.scale[j]<=1e-15) d.scale[j]=1;
    }
    const auto liftedX=liftEdmdControlSnapshots(d,x,pairs);
    const auto liftedY=liftEdmdControlSnapshots(d,y,pairs);
    out.lifted=dynamicModeDecompositionWithControl(liftedX,liftedY,inputs,pairs,
        d.exponents.size()/dimension,inputDimension,deltaTime,options.regression);
    const auto predicted=predictEdmdControl(out,x,inputs,pairs);
    long double error=0,reference=0;
    for(std::size_t i=0;i<y.size();++i) {
        const long double residual=static_cast<long double>(predicted[i])-y[i];
        error+=residual*residual; reference+=static_cast<long double>(y[i])*y[i];
    }
    out.snapshotOneStepRelativeError=static_cast<double>(std::sqrt(error/
        std::max(reference,static_cast<long double>(std::numeric_limits<double>::min()))));
    return out;
}

std::vector<std::complex<double>> evaluateEdmdControlEigenfunctions(
    const EdmdControlResult& model,const std::vector<double>& snapshots,std::size_t samples) {
    const auto theta=liftEdmdControlSnapshots(model.dictionary,snapshots,samples);
    if(model.lifted.dimension!=model.dictionary.exponents.size()/model.dictionary.dimension)
        throw std::invalid_argument("EDMDc model/dictionary dimension mismatch");
    return evaluateDmdControlEigenfunctions(model.lifted,theta,samples);
}

EigenfunctionFitResult EdmdControlResult::eigenfunctionSet() const {
    validateDictionary(dictionary);
    if(!lifted.hasEigenfunctions||lifted.spectrum.size()!=lifted.rank)
        throw std::invalid_argument("EDMDc has no exportable polynomial eigenbasis");
    FunctionDomain q;q.dimension=dictionary.dimension;
    return {FunctionExpansion({dictionary.dimension,dictionary.exponents,dictionary.center,dictionary.scale},
        lifted.eigenfunctionCoefficients,q),
        lifted.spectrum,{}};
}

std::vector<double> predictEdmdControl(const EdmdControlResult& model,
    const std::vector<double>& snapshots,const std::vector<double>& inputs,std::size_t samples) {
    const auto& d=model.dictionary;
    const auto theta=liftEdmdControlSnapshots(d,snapshots,samples);
    const auto features=d.exponents.size()/d.dimension;
    if(model.lifted.dimension!=features)
        throw std::invalid_argument("EDMDc model/dictionary dimension mismatch");
    const auto prediction=predictDmdControl(model.lifted,theta,inputs,samples);
    std::vector<double> result(samples*d.dimension);
    const std::size_t firstLinear=d.includeConstant?1:0;
    for(std::size_t k=0;k<samples;++k)
        for(std::size_t j=0;j<d.dimension;++j)
            result[k*d.dimension+j]=d.center[j]+d.scale[j]*prediction[k*features+firstLinear+j];
    return result;
}
} // namespace dynsyslab::core

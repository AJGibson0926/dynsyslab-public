#include <dynsyslab/core/Decomposition.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <dynsyslab/numerics/NormalEquations.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <dynsyslab/numerics/Svd.hpp>

#include "FlatMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using RMatrix=detail::FlatMatrix<double>;
using CMatrix=detail::FlatMatrix<std::complex<double>>;

std::vector<std::vector<std::size_t>> monomials(
    std::size_t dimension,std::size_t degree,bool constant) {
    std::vector<std::vector<std::size_t>> result;
    if(constant)result.push_back(std::vector<std::size_t>(dimension));
    std::vector<std::size_t> exponent(dimension);
    std::function<void(std::size_t,std::size_t)> enumerate=
        [&](std::size_t axis,std::size_t left) {
            if(axis+1==dimension) {
                exponent[axis]=left;result.push_back(exponent);return;
            }
            for(std::size_t value=0;value<=left;++value) {
                exponent[axis]=value;enumerate(axis+1,left-value);
            }
        };
    for(std::size_t order=1;order<=degree;++order)enumerate(0,order);
    return result;
}

RMatrix dictionary(const RMatrix& states,const DictionaryOptions& options) {
    if(options.kind==DictionaryKind::polynomial) {
        if(!options.polynomialDegree)
            throw std::invalid_argument("EDMD polynomial degree must be positive");
        const auto exponents=monomials(
            states.columns,options.polynomialDegree,options.includeConstant);
        RMatrix result(states.rows,exponents.size());
        for(std::size_t sample=0;sample<states.rows;++sample)
            for(std::size_t feature=0;feature<exponents.size();++feature) {
                double value=1.0;
                for(std::size_t axis=0;axis<states.columns;++axis)
                    value*=std::pow(states(sample,axis),
                                    static_cast<int>(exponents[feature][axis]));
                result(sample,feature)=value;
            }
        return result;
    }
    if(!options.centerCount||options.radialBasisCenters.size()
        !=options.centerCount*states.columns||!(options.radialBasisWidth>0.0)
        ||!std::isfinite(options.radialBasisWidth))
        throw std::invalid_argument("EDMD radial-basis dictionary is incomplete");
    RMatrix result(states.rows,options.centerCount);
    const double denominator=2.0*options.radialBasisWidth
        *options.radialBasisWidth;
    for(std::size_t sample=0;sample<states.rows;++sample)
        for(std::size_t center=0;center<options.centerCount;++center) {
            double distance=0.0;
            for(std::size_t axis=0;axis<states.columns;++axis) {
                const double difference=states(sample,axis)
                    -options.radialBasisCenters[center*states.columns+axis];
                distance+=difference*difference;
            }
            result(sample,center)=std::exp(-distance/denominator);
        }
    return result;
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

CMatrix product(const CMatrix& left,const CMatrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

RMatrix transpose(const RMatrix& source) {
    RMatrix result(source.columns,source.rows);
    for(std::size_t row=0;row<source.rows;++row)
        for(std::size_t column=0;column<source.columns;++column)
            result(column,row)=source(row,column);
    return result;
}

void requireUniform(const std::vector<double>& time) {
    if(time.size()<4)
        throw std::invalid_argument("HAVOK requires at least four samples");
    const double step=time[1]-time[0];
    if(!(step>0.0)||!std::isfinite(step))
        throw std::invalid_argument("HAVOK time must increase");
    for(std::size_t index=2;index<time.size();++index)
        if(std::abs((time[index]-time[index-1])-step)
            >1.0e-10*std::max(1.0,std::abs(step)))
            throw std::invalid_argument("HAVOK requires uniform time steps");
}

}  // namespace

EDMDResult extendedDynamicModeDecomposition(
    const std::vector<double>& stateValues,std::size_t samples,
    std::size_t dimension,double step,const EDMDOptions& options) {
    if(samples<3||!dimension||stateValues.size()!=samples*dimension
        ||!(step>0.0)||options.regularization<0.0
        ||!std::isfinite(options.regularization)
        ||!std::all_of(stateValues.begin(),stateValues.end(),
            [](double value){return std::isfinite(value);}))
        throw std::invalid_argument("invalid EDMD inputs");
    const RMatrix states{stateValues,samples,dimension};
    const RMatrix all=dictionary(states,options.dictionary);
    RMatrix x(samples-1,all.columns),y(samples-1,all.columns);
    for(std::size_t sample=0;sample+1<samples;++sample)
        for(std::size_t feature=0;feature<all.columns;++feature) {
            x(sample,feature)=all(sample,feature);
            y(sample,feature)=all(sample+1,feature);
        }
    const RMatrix operation{numerics::solveNormalEquations(
        x.values,x.rows,x.columns,y.values,y.columns,options.regularization),
        x.columns,y.columns};
    const auto eigen=numerics::eigensystem(operation.values,operation.rows);
    const CMatrix coefficients{eigen.rightEigenvectors,operation.rows,operation.columns};
    const CMatrix functions=product(complexValues(all),coefficients);
    numerics::LeastSquaresOptions complexOptions;
    complexOptions.method=numerics::LeastSquaresMethod::divideAndConquerSvd;
    const CMatrix stateComplex=complexValues(states);
    const CMatrix modeCoefficients{numerics::solveLeastSquares(
        functions.values,functions.rows,functions.columns,stateComplex.values,
        stateComplex.columns,complexOptions),functions.columns,stateComplex.columns};
    const CMatrix modes=transpose(modeCoefficients);
    EDMDResult result;
    result.samples=samples;
    result.stateDimension=dimension;
    result.featureCount=all.columns;
    result.deltaTime=step;
    result.spectrum=Spectrum::fromDiscreteMultipliers(eigen.eigenvalues,step);
    result.operatorMatrix=operation.values;
    result.eigenfunctionCoefficients=coefficients.values;
    if(options.dictionary.kind==DictionaryKind::polynomial)
        for(const auto& row:monomials(dimension,options.dictionary.polynomialDegree,options.dictionary.includeConstant))
            result.dictionaryExponents.insert(result.dictionaryExponents.end(),row.begin(),row.end());
    result.eigenfunctionValues=functions.values;
    result.modes=modes.values;
    return result;
}

HAVOKResult havokDecomposition(
    const std::vector<double>& stateValues,const std::vector<double>& time,
    std::size_t samples,std::size_t dimension,const HAVOKOptions& options) {
    if(stateValues.size()!=samples*dimension||time.size()!=samples
        ||options.observableIndex>=dimension||options.delays<2
        ||options.delays>=samples-1||options.rank<2)
        throw std::invalid_argument("invalid HAVOK inputs");
    requireUniform(time);
    const std::size_t columns=samples-options.delays+1;
    RMatrix hankel(options.delays,columns);
    std::vector<double> signal(samples);
    for(std::size_t sample=0;sample<samples;++sample)
        signal[sample]=stateValues[sample*dimension+options.observableIndex];
    for(std::size_t delay=0;delay<options.delays;++delay)
        for(std::size_t column=0;column<columns;++column)
            hankel(delay,column)=signal[column+delay];
    const auto svd=numerics::singularValueDecomposition(
        hankel.values,hankel.rows,hankel.columns,options.rank);
    const std::size_t rank=svd.rank;
    if(rank<2)throw std::runtime_error("HAVOK numerical rank is below two");
    const RMatrix left{svd.left,options.delays,rank};
    const RMatrix rightTranspose{svd.rightTranspose,rank,columns};
    const RMatrix coordinates=transpose(rightTranspose);
    const double step=time[1]-time[0];
    const auto order=options.derivativeOrder;
    if(order!=1&&order!=2&&order!=4)throw std::invalid_argument("HAVOK derivative order must be 1, 2 or 4");
    if(order==4&&columns<5)throw std::invalid_argument("fourth-order HAVOK derivatives need five delay columns");
    RMatrix derivative(columns,rank);
    for(std::size_t mode=0;mode<rank;++mode) {
        const auto x=[&](std::size_t k){return coordinates(k,mode);};
        if(order==1) {
            derivative(0,mode)=(x(1)-x(0))/step;
            derivative(columns-1,mode)=(x(columns-1)-x(columns-2))/step;
        } else {
            derivative(0,mode)=(-3.0*x(0)+4.0*x(1)-x(2))/(2.0*step);
            derivative(columns-1,mode)=(3.0*x(columns-1)-4.0*x(columns-2)+x(columns-3))/(2.0*step);
        }
        for(std::size_t sample=1;sample+1<columns;++sample)
            derivative(sample,mode)=order==4&&sample>=2&&sample+2<columns
                ?(-x(sample+2)+8.0*x(sample+1)-8.0*x(sample-1)+x(sample-2))/(12.0*step)
                :(x(sample+1)-x(sample-1))/(2.0*step);
    }
    RMatrix design(columns-1,rank),target(columns-1,rank-1);
    for(std::size_t sample=0;sample+1<columns;++sample) {
        for(std::size_t mode=0;mode<rank;++mode)
            design(sample,mode)=coordinates(sample,mode);
        for(std::size_t mode=0;mode+1<rank;++mode)
            target(sample,mode)=derivative(sample,mode);
    }
    const RMatrix coefficients{numerics::solveLeastSquares(
        design.values,design.rows,design.columns,target.values,target.columns),
        design.columns,target.columns};
    HAVOKResult result;
    result.delays=options.delays;
    result.rank=rank;
    result.embeddedSamples=columns;
    result.deltaTime=step;
    result.generator.resize((rank-1)*(rank-1));
    result.input.resize(rank-1);
    for(std::size_t output=0;output+1<rank;++output) {
        for(std::size_t predictor=0;predictor+1<rank;++predictor)
            result.generator[output*(rank-1)+predictor]=
                coefficients(predictor,output);
        result.input[output]=coefficients(rank-1,output);
    }
    result.forcing.resize(columns);
    for(std::size_t sample=0;sample<columns;++sample)
        result.forcing[sample]=coordinates(sample,rank-1);
    result.coordinates=coordinates.values;
    result.coordinateTime.assign(time.begin(),time.begin()+columns);
    result.leftSingularVectors=left.values;
    result.singularValues=svd.singularValues;
    result.rightSingularVectors=rightTranspose.values;
    result.observable=std::move(signal);
    return result;
}

HAVOKResult havokDecomposition(
    const std::vector<double>& states,const std::vector<double>& time,
    std::size_t samples,std::size_t dimension,const Observable& observable,
    const HAVOKOptions& options) {
    if(observable.stateDimension()!=dimension||observable.outputDimension()!=1)
        throw std::invalid_argument(
            "HAVOK requires one scalar state-aligned observable");
    auto scalar=observable.evaluateRealBatch(states,samples);
    HAVOKOptions scalarOptions=options;
    scalarOptions.observableIndex=0;
    return havokDecomposition(scalar,time,samples,1,scalarOptions);
}

}  // namespace dynsyslab::core

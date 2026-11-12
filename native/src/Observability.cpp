#include <dynsyslab/core/Observability.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Svd.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include "FlatMatrix.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
namespace dynsyslab::core {
namespace {
using Matrix=detail::FlatMatrix<double>;

Matrix product(const Matrix& left,const Matrix& right) {
    return {numerics::multiply(left.values,left.rows,left.columns,
        right.values,right.columns),left.rows,right.columns};
}

Matrix mapped(const std::vector<double>& values,std::size_t rows,
              std::size_t columns,const char* name) {
    if(!rows||!columns||values.size()!=rows*columns)
        throw std::invalid_argument(std::string(name)+" dimensions are inconsistent");
    if(!std::all_of(values.begin(),values.end(),
        [](double value){return std::isfinite(value);}))
        throw std::invalid_argument(std::string(name)+" must be finite");
    return {values,rows,columns};
}

Matrix transpose(const Matrix& source) {
    Matrix result(source.columns,source.rows);
    for(std::size_t row=0;row<source.rows;++row)
        for(std::size_t column=0;column<source.columns;++column)
            result(column,row)=source(row,column);
    return result;
}

Matrix identity(std::size_t dimension) {
    Matrix result(dimension,dimension);
    for(std::size_t index=0;index<dimension;++index)result(index,index)=1.0;
    return result;
}

Matrix symmetric(Matrix value,const char* name) {
    for(std::size_t row=0;row<value.rows;++row)
        for(std::size_t column=row;column<value.columns;++column) {
            const double mean=0.5*(value(row,column)+value(column,row));
            if(!std::isfinite(mean))
                throw std::runtime_error(std::string(name)+" is not finite");
            value(row,column)=mean;value(column,row)=mean;
        }
    return value;
}

double maximumAbsolute(const Matrix& matrix) {
    double result=0.0;
    for(const double value:matrix.values)result=std::max(result,std::abs(value));
    return result;
}

bool allFinite(const Matrix& matrix) {
    return std::all_of(matrix.values.begin(),matrix.values.end(),
        [](double value){return std::isfinite(value);});
}

void spectralSummary(ObservabilityResult& result,const Matrix& gramian,
                     bool retainExistingSingularValues) {
    const auto eigen=numerics::symmetricEigensystem(
        gramian.values,gramian.rows);
    Matrix directions(gramian.rows,gramian.columns);
    std::vector<double> gramianStrengths(gramian.rows);
    for(std::size_t index=0;index<gramian.rows;++index) {
        const std::size_t source=gramian.rows-index-1;
        gramianStrengths[index]=std::sqrt(std::max(0.0,eigen.eigenvalues[source]));
        for(std::size_t row=0;row<gramian.rows;++row)
            directions(row,index)=eigen.eigenvectors[row*gramian.rows+source];
    }
    if(!retainExistingSingularValues)result.singularValues=gramianStrengths;
    for(double strength:gramianStrengths)result.gramianEigenvalues.push_back(strength*strength);
    result.gramian=gramian.values;
    result.rightDirections=std::move(directions.values);
}

void summarizeStrengths(ObservabilityResult& result,double relativeTolerance,
                        std::size_t toleranceScale) {
    const double largest=result.singularValues.empty()
        ?0.0:result.singularValues.front();
    result.rankTolerance=relativeTolerance
        *static_cast<double>(toleranceScale)*largest;
    double energy=0.0;
    for(const double strength:result.singularValues)energy+=strength*strength;
    double entropy=0.0;
    for(const double strength:result.singularValues) {
        const double probability=energy>std::numeric_limits<double>::min()
            ?strength*strength/energy:0.0;
        if(probability>0.0)entropy-=probability*std::log(probability);
        if(strength>result.rankTolerance)++result.rank;
    }
    result.effectiveRank=energy>std::numeric_limits<double>::min()
        ?std::exp(entropy):0.0;
    result.resolvedConditionNumber=std::numeric_limits<double>::infinity();
    for(auto iterator=result.singularValues.rbegin();
        iterator!=result.singularValues.rend();++iterator)
        if(*iterator>result.rankTolerance) {
            result.resolvedConditionNumber=largest/(*iterator);break;
        }
    result.conditionNumber=result.rank==result.stateDimension
        ?result.resolvedConditionNumber:std::numeric_limits<double>::infinity();
}

}  // namespace

ObservabilityResult analyzeObservability(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs,std::size_t blocks,
    double relativeTolerance) {
    ObservabilityOptions options;
    options.blockCount=blocks;
    options.relativeTolerance=relativeTolerance;
    return analyzeObservability(generator,output,states,outputs,options);
}

const char* observabilityMethodName(ObservabilityMethod method) noexcept {
    switch(method) {
    case ObservabilityMethod::classicalMatrix:return "matrix";
    case ObservabilityMethod::discreteInfiniteGramian:return "discrete-infinite-gramian";
    case ObservabilityMethod::discreteGramian:return "discrete-gramian";
    case ObservabilityMethod::continuousFiniteGramian:return "continuous-gramian";
    case ObservabilityMethod::continuousInfiniteGramian:
        return "continuous-infinite-gramian";
    }
    return "unknown";
}

ObservabilityResult analyzeObservability(
    const std::vector<double>& generator,const std::vector<double>& output,
    std::size_t states,std::size_t outputs,const ObservabilityOptions& options) {
    if(!std::isfinite(options.relativeTolerance)||options.relativeTolerance<=0.0)
        throw std::invalid_argument("observability tolerance must be positive");
    const Matrix A=mapped(generator,states,states,"observability generator");
    const Matrix C=mapped(output,outputs,states,"observability output");
    ObservabilityResult result;
    result.method=options.method;
    result.stateDimension=states;
    result.outputDimension=outputs;
    if(options.method==ObservabilityMethod::classicalMatrix) {
        const std::size_t blocks=options.blockCount?options.blockCount:states;
        if(blocks>(std::numeric_limits<std::size_t>::max)()/outputs)
            throw std::invalid_argument("observability matrix dimensions overflow");
        if(!std::isfinite(options.maximumAbsoluteEntry)
            ||options.maximumAbsoluteEntry<=0.0)
            throw std::invalid_argument("observability safety bound must be positive");
        Matrix observability(outputs*blocks,states),power=identity(states);
        for(std::size_t block=0;block<blocks;++block) {
            const Matrix row=product(C,power);
            if(!allFinite(row)||maximumAbsolute(row)>options.maximumAbsoluteEntry)
                throw std::runtime_error(
                    "observability matrix exceeded its numerical safety bound");
            std::copy(row.values.begin(),row.values.end(),
                observability.values.begin()
                    +static_cast<std::ptrdiff_t>(block*outputs*states));
            if(block+1<blocks)power=product(power,A);
        }
        const auto svd=numerics::singularValueDecomposition(
            observability.values,observability.rows,observability.columns);
        result.matrix=observability.values;
        result.singularValues=svd.singularValues;
        result.leftDirections=svd.left;
        result.blockCount=blocks;
        result.horizon=static_cast<double>(blocks);
        spectralSummary(result,symmetric(product(
            transpose(observability),observability),"observability Gramian"),true);
        // Use the matrix SVD itself: squaring a poorly conditioned matrix loses
        // weak directions long before the SVD does.
        for(std::size_t row=0;row<states;++row)
            for(std::size_t mode=0;mode<svd.rank;++mode)
                result.rightDirections[row*states+mode]=svd.rightTranspose[mode*states+row];
    } else {
        Matrix gramian;
        if(options.method==ObservabilityMethod::discreteGramian) {
            const std::size_t blocks=options.blockCount?options.blockCount:states;
            if(options.sampleTime<0.0||!std::isfinite(options.sampleTime))
                throw std::invalid_argument("Gramian sample time cannot be negative");
            Matrix transition=A;
            if(options.sampleTime>0.0)
                transition={numerics::matrixExponential(
                    A.values,states,options.sampleTime),states,states};
            gramian=Matrix(states,states);
            Matrix block=C;
            for(std::size_t index=0;index<blocks;++index) {
                const Matrix contribution=product(transpose(block),block);
                for(std::size_t entry=0;entry<gramian.values.size();++entry)
                    gramian.values[entry]+=contribution.values[entry];
                block=product(block,transition);
                if(!allFinite(gramian)||!allFinite(block))
                    throw std::runtime_error(
                        "discrete observability Gramian overflowed");
            }
            gramian=symmetric(std::move(gramian),"discrete observability Gramian");
            result.blockCount=blocks;
            result.horizon=static_cast<double>(blocks);
        } else if(options.method==ObservabilityMethod::discreteInfiniteGramian) {
            if(options.sampleTime<0.0||!std::isfinite(options.sampleTime))
                throw std::invalid_argument("Gramian sample time cannot be negative");
            const auto transition=options.sampleTime>0.0
                ?numerics::matrixExponential(generator,states,options.sampleTime):generator;
            gramian={numerics::discreteInfiniteObservabilityGramian(
                transition,output,states,outputs),states,states};
            result.horizon=std::numeric_limits<double>::infinity();
        } else if(options.method==ObservabilityMethod::continuousFiniteGramian) {
            gramian={numerics::finiteObservabilityGramian(
                generator,output,states,outputs,options.horizon),states,states};
            result.horizon=options.horizon;
        } else {
            gramian={numerics::infiniteObservabilityGramian(
                generator,output,states,outputs),states,states};
            result.horizon=std::numeric_limits<double>::infinity();
        }
        spectralSummary(result,gramian,false);
    }
    const std::size_t toleranceScale=
        options.method==ObservabilityMethod::classicalMatrix
        ?std::max(states,outputs*result.blockCount):std::max(states,outputs);
    summarizeStrengths(result,options.relativeTolerance,toleranceScale);
    return result;
}


} // namespace dynsyslab::core

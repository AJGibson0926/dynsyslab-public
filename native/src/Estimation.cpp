#include <dynsyslab/core/Estimation.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <dynsyslab/numerics/Control.hpp>
#include <dynsyslab/numerics/Eigensystem.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Product.hpp>
#include <dynsyslab/numerics/Svd.hpp>

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

bool realMode(std::complex<double> value) {
    return std::abs(value.imag())<=1.0e-12;
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
    if(!retainExistingSingularValues)result.singularValues=std::move(gramianStrengths);
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
    result.conditionNumber=std::numeric_limits<double>::infinity();
    for(auto iterator=result.singularValues.rbegin();
        iterator!=result.singularValues.rend();++iterator)
        if(*iterator>result.rankTolerance) {
            result.conditionNumber=largest/(*iterator);break;
        }
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
            if(maximumAbsolute(row)>options.maximumAbsoluteEntry)
                throw std::runtime_error(
                    "observability matrix exceeded its numerical safety bound");
            std::copy(row.values.begin(),row.values.end(),
                observability.values.begin()
                    +static_cast<std::ptrdiff_t>(block*outputs*states));
            power=product(power,A);
        }
        const auto svd=numerics::singularValueDecomposition(
            observability.values,observability.rows,observability.columns);
        result.matrix=observability.values;
        result.singularValues=svd.singularValues;
        result.blockCount=blocks;
        result.horizon=static_cast<double>(blocks);
        spectralSummary(result,symmetric(product(
            transpose(observability),observability),"observability Gramian"),true);
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

std::vector<double> reconstructionSplitMatrix(
    const KoopmanModel& model,const ReconstructionResult& reconstruction) {
    if(reconstruction.modes!=model.liftedDimension()
        ||reconstruction.matrix.size()!=reconstruction.outputs*reconstruction.modes)
        throw std::invalid_argument("reconstruction does not match KROM");
    std::vector<double> result(reconstruction.outputs*model.splitDimension());
    for(std::size_t output=0;output<reconstruction.outputs;++output) {
        std::size_t split=0;
        for(std::size_t mode=0;mode<model.liftedDimension();++mode) {
            const auto coefficient=
                reconstruction.matrix[output*reconstruction.modes+mode];
            result[output*model.splitDimension()+split++]=coefficient.real();
            if(!realMode(model.continuousEigenvalues()[mode]))
                result[output*model.splitDimension()+split++]=-coefficient.imag();
        }
    }
    return result;
}

ContinuousLQESolution solveContinuousLQE(
    const std::vector<double>& generator,const std::vector<double>& measurement,
    const std::vector<double>& process,const std::vector<double>& noise,
    std::size_t states,std::size_t measurements,bool requireDetectable) {
    const Matrix A=mapped(generator,states,states,"LQE generator");
    const Matrix C=mapped(measurement,measurements,states,"LQE measurement");
    mapped(process,states,states,"LQE process covariance");
    mapped(noise,measurements,measurements,"LQE measurement covariance");
    const bool detectable=numerics::isDetectable(
        generator,measurement,states,measurements);
    if(requireDetectable&&!detectable)
        throw std::invalid_argument("continuous LQE pair (A,C) is not detectable");
    const auto dual=numerics::solveContinuousLqr(
        transpose(A).values,transpose(C).values,process,noise,states,measurements);
    const Matrix gainDual{dual.gain,measurements,states};
    const Matrix gain=transpose(gainDual);
    Matrix closed=A;
    const Matrix feedback=product(gain,C);
    for(std::size_t index=0;index<closed.values.size();++index)
        closed.values[index]-=feedback.values[index];
    const auto poles=numerics::eigensystem(closed.values,states,false);
    ContinuousLQESolution result;
    result.stateDimension=states;
    result.measurementDimension=measurements;
    result.detectable=detectable;
    result.relativeResidual=dual.relativeResidual;
    result.covariance=dual.riccati;
    result.gain=gain.values;
    result.observerPoles=poles.eigenvalues;
    return result;
}

KoopmanLQEDesign designKoopmanLQE(
    const KoopmanModel& model,const std::string& measurementName,
    const std::string& stateName,const std::vector<double>& process,
    const std::vector<double>& noise,bool requireDetectable) {
    if(!model.hasReconstruction(measurementName)
        ||!model.hasReconstruction(stateName))
        throw std::invalid_argument(
            "KLQE requires named measurement and state reconstructions");
    const auto& state=model.reconstruction(stateName);
    if(state.outputs!=model.stateDimension())
        throw std::invalid_argument(
            "KLQE state reconstruction must output the full dynamical state");
    auto result=designKoopmanLQEFromMeasurements(model,measurementName,process,noise,stateName,requireDetectable);
    result.physicalStateDimension=model.stateDimension();
    return result;
}

KoopmanLQEDesign designKoopmanLQEFromMeasurements(
    const KoopmanModel& model,const std::string& measurementName,
    const std::vector<double>& process,const std::vector<double>& noise,
    const std::string& outputName,bool requireDetectable) {
    if(!model.hasReconstruction(measurementName))
        throw std::invalid_argument("KLQE requires a named measurement reconstruction");
    const auto& measurement=model.reconstruction(measurementName);
    KoopmanLQEDesign result;
    result.modelIdentity=koopmanModelIdentity(model);
    result.measurementObservable=measurementName;
    result.stateObservable=outputName;
    result.generator=model.generatorSplit();
    result.measurementMatrix=reconstructionSplitMatrix(model,measurement);
    if(!outputName.empty()) {
        const auto& decoder=model.reconstruction(outputName);
        result.stateReconstructionMatrix=reconstructionSplitMatrix(model,decoder);
        result.outputDimension=decoder.outputs;
    }
    result.solution=solveContinuousLQE(
        result.generator,result.measurementMatrix,process,noise,
        model.splitDimension(),measurement.outputs,requireDetectable);
    return result;
}

KoopmanLQEResult runKoopmanLQE(
    const KoopmanModel& model,const KoopmanLQEDesign& design,
    const std::vector<double>& time,const std::vector<double>& observations,
    const std::vector<double>& initialPhysicalState,const std::vector<double>& truth) {
    if(initialPhysicalState.size()!=model.stateDimension())
        throw std::invalid_argument("KLQE initial physical state dimension mismatch");
    return runKoopmanLQEFromLifted(model,design,time,observations,
        model.evaluateSplit(initialPhysicalState.data()),truth);
}

KoopmanLQEResult runKoopmanLQEFromLifted(
    const KoopmanModel& model,const KoopmanLQEDesign& design,
    const std::vector<double>& time,const std::vector<double>& observations,
    const std::vector<double>& initialLifted,const std::vector<double>& truth) {
    const std::size_t samples=time.size(),states=model.splitDimension();
    const std::size_t measurements=design.solution.measurementDimension;
    const std::size_t physical=design.outputDimension?design.outputDimension:design.physicalStateDimension;
    if(!design.modelIdentity.empty()&&design.modelIdentity!=koopmanModelIdentity(model))
        throw std::invalid_argument("KLQE design belongs to a different KROM");
    if(!samples||observations.size()!=samples*measurements
        ||initialLifted.size()!=states
        ||(!truth.empty()&&truth.size()!=samples*physical))
        throw std::invalid_argument("KLQE run arrays have inconsistent dimensions");
    for(const auto* values:{&observations,&initialLifted,&truth})
        for(double v:*values)if(!std::isfinite(v))throw std::invalid_argument("KLQE values must be finite");
    for(std::size_t index=0;index<samples;++index)
        if(!std::isfinite(time[index])
            ||(index&&time[index]<=time[index-1]))
            throw std::invalid_argument(
                "KLQE time must be finite and strictly increasing");
    Matrix A=mapped(design.generator,states,states,"KLQE generator");
    const Matrix C=mapped(
        design.measurementMatrix,measurements,states,"KLQE measurement");
    const Matrix decoder=physical?mapped(
        design.stateReconstructionMatrix,physical,states,"KLQE output decoder"):Matrix{};
    if(!physical&&!design.stateReconstructionMatrix.empty())
        throw std::invalid_argument("KLQE decoder has no output dimension");
    const Matrix gain=mapped(
        design.solution.gain,states,measurements,"KLQE gain");
    const Matrix feedback=product(gain,C);
    for(std::size_t index=0;index<A.values.size();++index)
        A.values[index]-=feedback.values[index];
    std::vector<double> lifted=initialLifted;
    KoopmanLQEResult result;
    result.time=time;
    result.samples=samples;
    result.liftedDimension=states;
    result.stateDimension=physical;
    result.measurementDimension=measurements;
    result.liftedEstimate.resize(samples*states);
    result.stateEstimate.resize(samples*physical);
    result.measurementEstimate.resize(samples*measurements);
    result.innovation.resize(samples*measurements);
    if(!truth.empty())result.stateError.resize(samples*physical);
    result.stateRmse.assign(physical,0.0);
    result.innovationRmse.assign(measurements,0.0);
    for(std::size_t sample=0;sample<samples;++sample) {
        if(sample) {
            const auto held=numerics::zeroOrderHold(
                A.values,gain.values,states,measurements,
                time[sample]-time[sample-1]);
            std::vector<double> updated(states);
            for(std::size_t row=0;row<states;++row) {
                for(std::size_t column=0;column<states;++column)
                    updated[row]+=held.transition[row*states+column]*lifted[column];
                for(std::size_t measurement=0;measurement<measurements;++measurement)
                    updated[row]+=held.response[row*measurements+measurement]
                        *observations[(sample-1)*measurements+measurement];
            }
            lifted=std::move(updated);
        }
        std::vector<double> estimate(physical),predicted(measurements);
        for(std::size_t row=0;row<physical;++row)
            for(std::size_t column=0;column<states;++column)
                estimate[row]+=decoder(row,column)*lifted[column];
        for(std::size_t row=0;row<measurements;++row)
            for(std::size_t column=0;column<states;++column)
                predicted[row]+=C(row,column)*lifted[column];
        for(std::size_t index=0;index<states;++index)
            result.liftedEstimate[sample*states+index]=lifted[index];
        for(std::size_t index=0;index<physical;++index) {
            result.stateEstimate[sample*physical+index]=estimate[index];
            if(!truth.empty()) {
                const double error=truth[sample*physical+index]-estimate[index];
                result.stateError[sample*physical+index]=error;
                result.stateRmse[index]+=error*error;
            }
        }
        for(std::size_t index=0;index<measurements;++index) {
            const double innovation=
                observations[sample*measurements+index]-predicted[index];
            result.measurementEstimate[sample*measurements+index]=predicted[index];
            result.innovation[sample*measurements+index]=innovation;
            result.innovationRmse[index]+=innovation*innovation;
        }
    }
    for(double& value:result.stateRmse)
        value=std::sqrt(value/static_cast<double>(samples));
    for(double& value:result.innovationRmse)
        value=std::sqrt(value/static_cast<double>(samples));
    return result;
}

}  // namespace dynsyslab::core
